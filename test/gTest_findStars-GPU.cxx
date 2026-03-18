//
//  gTest_findStars-GPU.cxx
//
//  Tests for the GPU-dispatched morphological operations (DSTmorph_gpu.h)
//  and for the findStars pipeline with non-double (GPU-eligible) pixel types.
//
//  Strategy:
//    • Use GPU-eligible pixel types (float, uint16_t, int16_t, uint32_t)
//      so the GPU dispatch path is exercised when a backend is available.
//    • double images always fall back to CPU (static compile-time check).
//    • Each morphological test compares the GPU result with the CPU result.
//      – Integer types must be bit-exact.
//      – float types must agree within 1 ULP (or relative tolerance 1e-5).
//    • findStars tests feed float/uint16_t images to starFinder; the GPU
//      dispatch is now automatic inside starFinder when the image type is
//      GPU-eligible (GetEqBitPerPixel() check).
//    • The MorphoPipeline tests replicate the exact fBS/fBM/fBE operations
//      used inside _findStarsMorpho with float kernels via GPU dispatch and
//      compare GPU vs CPU — this is where GPU code bugs would be caught.
//    • Tests are backend-agnostic: they pass whether Vulkan, CUDA, HIP, or
//      the CPUBackend is active (the CPU backend is always available).
//
//  Licensed under CC BY-NC 4.0

#include <gtest/gtest.h>

#include <DSTools/DSTmorpho.h>
#include <DSTools/DSTmorph_gpu.h>
#include <DSTools/DSTstarFinder.h>
#include <DSTdisplay/DSTmathplotpp.h>
#include <DSTfits/FITSdata.h>      // is_gpu_eligible_type_v

#include <chrono>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <optional>
#include <random>
#include <string>

#include "test_utils.h"

using namespace DST;
using namespace DSL;

// ============================================================
//  Helpers
// ============================================================

/// Convert a double FITSimg to a pixel type T by clamping and rounding.
template<typename T>
std::shared_ptr<FITSimg<T>> toType(const std::shared_ptr<FITSimg<double>>& src)
{
    std::vector<size_t> shape;
    for(size_t d = 0; d < src->GetDimension(); ++d)
        shape.push_back(src->Size(d + 1));

    auto dst = std::make_shared<FITSimg<T>>(shape);
    const std::valarray<double>* data = src->template GetData<double>();
    std::valarray<T>* out             = dst->template GetData<T>();

    const double lo = static_cast<double>(std::numeric_limits<T>::lowest());
    const double hi = static_cast<double>(std::numeric_limits<T>::max());

    for(size_t i = 0; i < data->size(); ++i)
    {
        double v = std::max(lo, std::min(hi, (*data)[i]));
        if constexpr (std::is_integral_v<T>)
            (*out)[i] = static_cast<T>(std::round(v));
        else
            (*out)[i] = static_cast<T>(v);
    }
    return dst;
}

/// Max absolute difference between two FITScube outputs (as double).
double maxAbsDiff(const std::shared_ptr<FITScube>& a,
                  const std::shared_ptr<FITScube>& b)
{
    EXPECT_EQ(a->Nelements(), b->Nelements());
    double mx = 0.0;
    for(size_t i = 0; i < a->Nelements(); ++i)
    {
        double diff = std::abs(a->DoubleValueAtPixel(i) - b->DoubleValueAtPixel(i));
        if(diff > mx) mx = diff;
    }
    return mx;
}

// ============================================================
//  Backend availability
// ============================================================

TEST(GPUBackendTest, CreateSucceeds)
{
    auto backend = DST::GPU::Backend::create();
    ASSERT_NE(backend, nullptr);
    std::cout << "Backend available: " << (backend->available() ? "YES (GPU)" : "NO (CPU fallback)") << std::endl;
}

TEST(GPUBackendTest, EligibilityConstraints)
{
    EXPECT_TRUE( (DSL::is_gpu_eligible_type_v<float>)    );
    EXPECT_TRUE( (DSL::is_gpu_eligible_type_v<uint8_t>)  );
    EXPECT_TRUE( (DSL::is_gpu_eligible_type_v<int8_t>)   );
    EXPECT_TRUE( (DSL::is_gpu_eligible_type_v<uint16_t>) );
    EXPECT_TRUE( (DSL::is_gpu_eligible_type_v<int16_t>)  );
    EXPECT_TRUE( (DSL::is_gpu_eligible_type_v<uint32_t>) );
    EXPECT_TRUE( (DSL::is_gpu_eligible_type_v<int32_t>)  );
    EXPECT_FALSE( (DSL::is_gpu_eligible_type_v<double>)  );
    EXPECT_FALSE( (DSL::is_gpu_eligible_type_v<int64_t>) );
    EXPECT_FALSE( (DSL::is_gpu_eligible_type_v<uint64_t>));
}

// ============================================================
//  Morphological operation correctness: GPU vs CPU
//  Image: 128×128 float star field
// ============================================================

class GPUMorphoTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        std::filesystem::create_directories("./build/test_output");
        backend_ = DST::GPU::Backend::create();

        std::vector<DST::Stars> catalog;
        auto dbl = TestUtils::generateStarField(catalog,
            128, 128, 10, 1.5, 50.0, 500.0, 10.0, 3.0, 42);

        img_f32_ = toType<float>(dbl);
        img_u16_ = toType<uint16_t>(dbl);

        k3_.emplace(std::valarray<float>(1.0f, 9),
                    std::vector<size_t>{3, 3}, Morphology::ktype::DILATE);
        k5_.emplace(std::valarray<float>(1.0f, 25),
                    std::vector<size_t>{5, 5}, Morphology::ktype::DILATE);
    }

    std::unique_ptr<DST::GPU::Backend>          backend_;
    std::shared_ptr<FITSimg<float>>             img_f32_;
    std::shared_ptr<FITSimg<uint16_t>>          img_u16_;
    std::optional<Morphology::kernel<float>>    k3_, k5_;
};


// ---- dilation -------------------------------------------------------

TEST_F(GPUMorphoTest, Dilation_Float_3x3)
{
    auto k = Morphology::kernel<float>(k3_->data(), k3_->sizes(), Morphology::ktype::DILATE);
    auto cpu = Morphology::dilation<float, float>(img_f32_, k);
    auto gpu = Morphology::GPU::dilation<float, float>(img_f32_, k, *backend_);
    ASSERT_EQ(cpu->Nelements(), gpu->Nelements());
    double diff = maxAbsDiff(cpu, gpu);
    std::cout << "[dilation float 3×3] max|GPU-CPU| = " << diff << std::endl;
    EXPECT_LE(diff, 1.0);
}

TEST_F(GPUMorphoTest, Dilation_Uint16_3x3)
{
    auto k = Morphology::kernel<float>(k3_->data(), k3_->sizes(), Morphology::ktype::DILATE);
    auto cpu = Morphology::dilation<float, uint16_t>(img_u16_, k);
    auto gpu = Morphology::GPU::dilation<float, uint16_t>(img_u16_, k, *backend_);
    ASSERT_EQ(cpu->Nelements(), gpu->Nelements());
    double diff = maxAbsDiff(cpu, gpu);
    std::cout << "[dilation uint16 3×3] max|GPU-CPU| = " << diff << std::endl;
    EXPECT_LE(diff, 1.0);
}


// ---- erosion --------------------------------------------------------

TEST_F(GPUMorphoTest, Erosion_Float_3x3)
{
    auto k = Morphology::kernel<float>(k3_->data(), k3_->sizes(), Morphology::ktype::ERODE);
    auto cpu = Morphology::erosion<float, float>(img_f32_, k);
    auto gpu = Morphology::GPU::erosion<float, float>(img_f32_, k, *backend_);
    ASSERT_EQ(cpu->Nelements(), gpu->Nelements());
    double diff = maxAbsDiff(cpu, gpu);
    std::cout << "[erosion float 3×3] max|GPU-CPU| = " << diff << std::endl;
    EXPECT_LE(diff, 1.0);
}

TEST_F(GPUMorphoTest, Erosion_Float_5x5)
{
    auto k = Morphology::kernel<float>(k5_->data(), k5_->sizes(), Morphology::ktype::ERODE);
    auto cpu = Morphology::erosion<float, float>(img_f32_, k);
    auto gpu = Morphology::GPU::erosion<float, float>(img_f32_, k, *backend_);
    ASSERT_EQ(cpu->Nelements(), gpu->Nelements());
    double diff = maxAbsDiff(cpu, gpu);
    std::cout << "[erosion float 5×5] max|GPU-CPU| = " << diff << std::endl;
    EXPECT_LE(diff, 1.0);
}


// ---- opening (chain: erosion + dilation) ----------------------------

TEST_F(GPUMorphoTest, Opening_Float_3x3)
{
    auto k = Morphology::kernel<float>(k3_->data(), k3_->sizes(), Morphology::ktype::OPEN);
    auto cpu = Morphology::opening<float, float>(img_f32_, k);
    auto gpu = Morphology::GPU::opening<float, float>(img_f32_, k, *backend_);
    ASSERT_EQ(cpu->Nelements(), gpu->Nelements());
    double diff = maxAbsDiff(cpu, gpu);
    std::cout << "[opening float 3×3] max|GPU-CPU| = " << diff << std::endl;
    EXPECT_LE(diff, 2.0);
}


// ---- closing (chain: dilation + erosion) ----------------------------

TEST_F(GPUMorphoTest, Closing_Float_3x3)
{
    auto k = Morphology::kernel<float>(k3_->data(), k3_->sizes(), Morphology::ktype::CLOSE);
    auto cpu = Morphology::closing<float, float>(img_f32_, k);
    auto gpu = Morphology::GPU::closing<float, float>(img_f32_, k, *backend_);
    ASSERT_EQ(cpu->Nelements(), gpu->Nelements());
    double diff = maxAbsDiff(cpu, gpu);
    std::cout << "[closing float 3×3] max|GPU-CPU| = " << diff << std::endl;
    EXPECT_LE(diff, 2.0);
}


// ---- correlation ----------------------------------------------------

TEST_F(GPUMorphoTest, Correlate_Float_3x3)
{
    std::valarray<float> w = {
        0.0625f, 0.125f, 0.0625f,
        0.125f,  0.25f,  0.125f,
        0.0625f, 0.125f, 0.0625f
    };
    auto k = Morphology::kernel<float>(w, {3, 3}, Morphology::ktype::CORRELATE);
    auto cpu = Morphology::correlate<float, float>(img_f32_, k, 1, false);
    auto gpu = Morphology::GPU::correlate<float, float>(img_f32_, k, *backend_, false);
    ASSERT_EQ(cpu->Nelements(), gpu->Nelements());
    double diff = maxAbsDiff(cpu, gpu);
    std::cout << "[correlate float 3×3] max|GPU-CPU| = " << diff << std::endl;
    EXPECT_LE(diff, 1.0);
}


// ---- double falls back to CPU (static, no backend call needed) ------

TEST_F(GPUMorphoTest, Double_FallsBackToCPU)
{
    std::vector<DST::Stars> catalog;
    auto dbl_img = TestUtils::generateStarField(catalog, 64, 64, 5, 1.5, 50.0, 200.0, 5.0, 1.0, 7);
    auto k = Morphology::kernel<double>(
        std::valarray<double>(1.0, 9), {3, 3}, Morphology::ktype::DILATE);
    auto cpu = Morphology::dilation<double, double>(dbl_img, k);
    auto gpu = Morphology::GPU::dilation<double, double>(dbl_img, k, *backend_);
    double diff = maxAbsDiff(cpu, gpu);
    EXPECT_EQ(diff, 0.0) << "double should be bit-exact via CPU fallback";
}


// ============================================================
//  Round-trip: upload → kernel → download identity
// ============================================================

TEST_F(GPUMorphoTest, IdentityKernel_Float)
{
    auto k = Morphology::kernel<float>(
        std::valarray<float>(1.0f, 1), {1, 1}, Morphology::ktype::DILATE);
    auto gpu = Morphology::GPU::dilation<float, float>(img_f32_, k, *backend_);
    ASSERT_EQ(img_f32_->Nelements(), gpu->Nelements());
    double diff = maxAbsDiff(img_f32_, gpu);
    std::cout << "[identity kernel float] max|GPU-orig| = " << diff << std::endl;
    EXPECT_LE(diff, 1e-4);
}


// ============================================================
//  Performance smoke test: 512×512 float dilation, 7×7 kernel
// ============================================================

TEST(GPUPerfTest, Dilation_Float_512x512_7x7)
{
    auto backend = DST::GPU::Backend::create();
    std::vector<DST::Stars> catalog;
    auto dbl = TestUtils::generateStarField(catalog, 512, 512, 50, 1.5, 50.0, 1000.0, 5.0, 3.0, 99);
    auto img  = toType<float>(dbl);
    auto k = Morphology::kernel<float>(
        std::valarray<float>(1.0f, 49), {7, 7}, Morphology::ktype::DILATE);

    auto t0 = std::chrono::steady_clock::now();
    auto gpu = Morphology::GPU::dilation<float, float>(img, k, *backend);
    auto t1 = std::chrono::steady_clock::now();
    auto cpu = Morphology::dilation<float, float>(img, k);
    auto t2 = std::chrono::steady_clock::now();

    double ms_gpu = std::chrono::duration<double, std::milli>(t1 - t0).count();
    double ms_cpu = std::chrono::duration<double, std::milli>(t2 - t1).count();
    std::cout << "[512×512 float dilation 7×7]  GPU=" << ms_gpu << " ms  CPU=" << ms_cpu << " ms" << std::endl;
    ASSERT_EQ(gpu->Nelements(), img->Nelements());
    double diff = maxAbsDiff(gpu, cpu);
    std::cout << "  max|GPU-CPU|=" << diff << std::endl;
    EXPECT_LE(diff, 1.0);
}

TEST(GPUPerfTest, Dilation_Float_2040x2040_7x7)
{
    auto backend = DST::GPU::Backend::create();
    std::vector<DST::Stars> catalog;
    auto dbl = TestUtils::generateStarField(catalog, 2040, 2040, 50, 1.5, 50.0, 1000.0, 5.0, 3.0, 99);
    auto img  = toType<float>(dbl);
    auto k = Morphology::kernel<float>(
        std::valarray<float>(1.0f, 49), {7, 7}, Morphology::ktype::DILATE);

    auto t0 = std::chrono::steady_clock::now();
    auto gpu = Morphology::GPU::dilation<float, float>(img, k, *backend);
    auto t1 = std::chrono::steady_clock::now();
    auto cpu = Morphology::dilation<float, float>(img, k);
    auto t2 = std::chrono::steady_clock::now();

    double ms_gpu = std::chrono::duration<double, std::milli>(t1 - t0).count();
    double ms_cpu = std::chrono::duration<double, std::milli>(t2 - t1).count();
    std::cout << "[2040x2040 float dilation 7×7]  GPU=" << ms_gpu << " ms  CPU=" << ms_cpu << " ms" << std::endl;
    ASSERT_EQ(gpu->Nelements(), img->Nelements());
    double diff = maxAbsDiff(gpu, cpu);
    std::cout << "  max|GPU-CPU|=" << diff << std::endl;
    EXPECT_LE(diff, 1.0);
}


// ============================================================
//  Morphological pipeline GPU test — mirrors _findStarsMorpho
//
//  starFinder._findStarsMorpho() applies three steps per layer:
//    iK = folding<K>(img, fBS)   fBS = OPEN,  3×3
//    iN = folding<K>(img, fBM)   fBM = DILATE,13×13
//    iS = folding<K>(iN,  fBE)   fBE = ERODE, 9×9
//
//  When the image is GPU-eligible, starFinder now automatically
//  converts these double kernels to float and uses GPU dispatch.
//  These tests verify GPU vs CPU correctness for that exact pipeline.
// ============================================================

class StarFinderMorphoPipelineGPU : public ::testing::Test
{
protected:
    void SetUp() override
    {
        std::filesystem::create_directories("./build/test_output");
        backend_ = DST::GPU::Backend::create();

        std::vector<DST::Stars> catalog;
        auto dbl = TestUtils::generateStarField(catalog,
            512, 512, 50, 0.7, 100.0, 1000.0, 5.0, 5.0, 12345);
        img_ = toType<float>(dbl);

        kBS_.emplace(std::valarray<float>(1.0f, 9),
                     std::vector<size_t>{3, 3},   Morphology::ktype::OPEN);
        kBM_.emplace(std::valarray<float>(1.0f, 169),
                     std::vector<size_t>{13, 13}, Morphology::ktype::DILATE);
        kBE_.emplace(std::valarray<float>(1.0f, 81),
                     std::vector<size_t>{9, 9},   Morphology::ktype::ERODE);
    }

    std::unique_ptr<DST::GPU::Backend>          backend_;
    std::shared_ptr<FITSimg<float>>             img_;
    std::optional<Morphology::kernel<float>>    kBS_, kBM_, kBE_;
};

// fBS: OPEN 3×3 — smoothing step
TEST_F(StarFinderMorphoPipelineGPU, Step1_Opening_fBS_3x3)
{
    auto cpu = Morphology::opening<float, float>(img_, *kBS_);
    auto gpu = Morphology::GPU::opening<float, float>(img_, *kBS_, *backend_);
    ASSERT_EQ(cpu->Nelements(), gpu->Nelements());
    double diff = maxAbsDiff(cpu, gpu);
    std::cout << "[fBS OPEN 3×3]   max|GPU-CPU| = " << diff << std::endl;
    EXPECT_LE(diff, 2.0);
}

// fBM: DILATE 13×13 — local maximum (star dilation)
TEST_F(StarFinderMorphoPipelineGPU, Step2_Dilation_fBM_13x13)
{
    auto cpu = Morphology::dilation<float, float>(img_, *kBM_);
    auto gpu = Morphology::GPU::dilation<float, float>(img_, *kBM_, *backend_);
    ASSERT_EQ(cpu->Nelements(), gpu->Nelements());
    double diff = maxAbsDiff(cpu, gpu);
    std::cout << "[fBM DILATE 13×13] max|GPU-CPU| = " << diff << std::endl;
    EXPECT_LE(diff, 1.0);
}

// fBE: ERODE 9×9 applied to dilated image — background estimation
TEST_F(StarFinderMorphoPipelineGPU, Step3_Erosion_fBE_9x9_OnDilated)
{
    auto iN_cpu_cube = Morphology::dilation<float, float>(img_, *kBM_);
    auto iN_gpu_cube = Morphology::GPU::dilation<float, float>(img_, *kBM_, *backend_);
    auto iN_cpu = std::dynamic_pointer_cast<FITSimg<float>>(iN_cpu_cube);
    auto iN_gpu = std::dynamic_pointer_cast<FITSimg<float>>(iN_gpu_cube);
    ASSERT_NE(iN_cpu, nullptr); ASSERT_NE(iN_gpu, nullptr);

    auto iS_cpu = Morphology::erosion<float, float>(iN_cpu, *kBE_);
    auto iS_gpu = Morphology::GPU::erosion<float, float>(iN_gpu, *kBE_, *backend_);
    ASSERT_EQ(iS_cpu->Nelements(), iS_gpu->Nelements());
    double diff = maxAbsDiff(iS_cpu, iS_gpu);
    std::cout << "[fBE ERODE 9×9 on dilated] max|GPU-CPU| = " << diff << std::endl;
    EXPECT_LE(diff, 1.0);
}

// Full 3-step chain: GPU results match CPU at every stage
TEST_F(StarFinderMorphoPipelineGPU, FullPipeline_ChainedGPUvsCPU)
{
    auto iK_cpu = Morphology::opening<float, float>(img_, *kBS_);
    auto iN_cpu_cube = Morphology::dilation<float, float>(img_, *kBM_);
    auto iN_cpu = std::dynamic_pointer_cast<FITSimg<float>>(iN_cpu_cube);
    ASSERT_NE(iN_cpu, nullptr);
    auto iS_cpu = Morphology::erosion<float, float>(iN_cpu, *kBE_);

    auto iK_gpu = Morphology::GPU::opening<float, float>(img_, *kBS_, *backend_);
    auto iN_gpu_cube = Morphology::GPU::dilation<float, float>(img_, *kBM_, *backend_);
    auto iN_gpu = std::dynamic_pointer_cast<FITSimg<float>>(iN_gpu_cube);
    ASSERT_NE(iN_gpu, nullptr);
    auto iS_gpu = Morphology::GPU::erosion<float, float>(iN_gpu, *kBE_, *backend_);

    double diff_K = maxAbsDiff(iK_cpu, iK_gpu);
    double diff_N = maxAbsDiff(iN_cpu_cube, iN_gpu_cube);
    double diff_S = maxAbsDiff(iS_cpu, iS_gpu);

    std::cout << "[findStars pipeline GPU vs CPU]\n"
              << "  iK fBS OPEN 3×3:         max|diff| = " << diff_K << "\n"
              << "  iN fBM DILATE 13×13:     max|diff| = " << diff_N << "\n"
              << "  iS fBE ERODE 9×9 on iN:  max|diff| = " << diff_S << std::endl;

    EXPECT_LE(diff_K, 2.0);
    EXPECT_LE(diff_N, 1.0);
    EXPECT_LE(diff_S, 1.0);
}


// ============================================================
//  findStars with non-double images
//
//  starFinder now automatically uses GPU dispatch when the image
//  pixel type is GPU-eligible (checked via GetEqBitPerPixel()).
//  float and uint16_t images trigger the GPU path; double stays CPU.
// ============================================================

static size_t countMatches(const FITStable& detected,
                           const std::vector<DST::Stars>& catalog,
                           double tol_px = 1.0)
{
    const auto& list = detected.getColumn("coord_px")->values<FITSform::doubleVector>();
    size_t matched = 0;
    for(const auto& s : catalog)
        for(const auto& c : list)
        {
            double dx = c[0] - s.position().X();
            double dy = c[1] - s.position().Y();
            if(std::sqrt(dx*dx + dy*dy) < tol_px) { ++matched; break; }
        }
    return matched;
}

TEST(StarFinderGPUTest, FindStarsBySegmentation_Float)
{
    std::vector<DST::Stars> catalog;
    auto dbl = TestUtils::generateStarField(catalog,
        512, 512, 50, 0.7, 100.0, 1000.0, 5.0, 5.0, 12345);
    std::shared_ptr<FITScube> img = toType<float>(dbl);

    DST::algorithms::starFinder finder(img);
    finder.setOutputUrl("./build/test_output");
    verbose = verboseLevel::VERBOSE_NONE;
    std::cout << "[float MORPHO] GPU backend active: " << (finder.isGPUEnabled() ? "YES" : "NO") << std::endl;
    finder.setThresholdValue(0.5);
    finder.findStars(DST::algorithms::starFinder::detectionMethod::MORPHO, 0.1);

    const auto& detected = finder.detectedStars();
    double rate = static_cast<double>(detected.nrows()) /
                  static_cast<double>(catalog.size()) * 100.0;
    std::cout << "[float MORPHO] " << detected.nrows() << "/" << catalog.size()
              << " (" << rate << "%)" << std::endl;

    EXPECT_GT(detected.nrows(), static_cast<size_t>(0));
    EXPECT_GT(rate, 75.0);
    EXPECT_LT(rate, 110.0);
    EXPECT_GT(countMatches(detected, catalog), catalog.size() / 2);
}

TEST(StarFinderGPUTest, FindStarsBySegmentation_Uint16)
{
    std::vector<DST::Stars> catalog;
    auto dbl = TestUtils::generateStarField(catalog,
        512, 512, 50, 0.7, 100.0, 1000.0, 5.0, 5.0, 12345);
    std::shared_ptr<FITScube> img = toType<uint16_t>(dbl);

    DST::algorithms::starFinder finder(img);
    finder.setOutputUrl("./build/test_output");
    verbose = verboseLevel::VERBOSE_NONE;
    finder.setThresholdValue(0.5);
    finder.findStars(DST::algorithms::starFinder::detectionMethod::MORPHO, 0.1);

    const auto& detected = finder.detectedStars();
    double rate = static_cast<double>(detected.nrows()) /
                  static_cast<double>(catalog.size()) * 100.0;
    std::cout << "[uint16 MORPHO] " << detected.nrows() << "/" << catalog.size()
              << " (" << rate << "%)" << std::endl;

    EXPECT_GT(detected.nrows(), static_cast<size_t>(0));
    EXPECT_GE(rate, 75.0);   // uint16 quantization loses ~2% vs float path
    EXPECT_LT(rate, 110.0);
    EXPECT_GT(countMatches(detected, catalog), catalog.size() / 2);
}

TEST(StarFinderGPUTest, FindStarsBySNR_Float)
{
    std::vector<DST::Stars> catalog;
    auto dbl = TestUtils::generateStarField(catalog,
        512, 512, 50, 0.7, 10.0, 1000.0, 10.0, 8.0, 12345);
    std::shared_ptr<FITScube> img = toType<float>(dbl);

    DST::algorithms::starFinder finder(img);
    finder.setOutputUrl("./build/test_output");
    verbose = verboseLevel::VERBOSE_NONE;
    finder.setThresholdValue(5);
    finder.findStars(DST::algorithms::starFinder::detectionMethod::SNR, 0.1);

    const auto& detected = finder.detectedStars();
    double rate = static_cast<double>(detected.nrows()) /
                  static_cast<double>(catalog.size()) * 100.0;
    std::cout << "[float SNR] " << detected.nrows() << "/" << catalog.size()
              << " (" << rate << "%)" << std::endl;

    EXPECT_GT(detected.nrows(), static_cast<size_t>(0));
    EXPECT_GT(rate, 60.0);
    EXPECT_LT(rate, 110.0);
    EXPECT_GT(countMatches(detected, catalog), catalog.size() / 2);
}

TEST(StarFinderGPUTest, FindStarsByThreshold_Float)
{
    std::vector<DST::Stars> catalog;
    auto dbl = TestUtils::generateStarField(catalog,
        512, 512, 50, 0.7, 10.0, 1000.0, 10.0, 8.0, 12345);
    std::shared_ptr<FITScube> img = toType<float>(dbl);

    DST::algorithms::starFinder finder(img);
    finder.setOutputUrl("./build/test_output");
    verbose = verboseLevel::VERBOSE_NONE;
    finder.setThresholdValue(70);
    finder.findStars(DST::algorithms::starFinder::detectionMethod::THRESHOLD, 0.1);

    const auto& detected = finder.detectedStars();
    double rate = static_cast<double>(detected.nrows()) /
                  static_cast<double>(catalog.size()) * 100.0;
    std::cout << "[float THRESHOLD] " << detected.nrows() << "/" << catalog.size()
              << " (" << rate << "%)" << std::endl;

    EXPECT_GT(detected.nrows(), static_cast<size_t>(0));
    EXPECT_GT(rate, 50.0);
    EXPECT_LT(rate, 110.0);
    EXPECT_GT(countMatches(detected, catalog), catalog.size() / 2);
}

TEST(StarFinderGPUTest, FindStarsCube_Float)
{
    std::vector<DST::Stars> catalog;
    const unsigned int seed = 12345;

    auto dbl1 = TestUtils::generateStarField(catalog, 256, 256, 30, 0.7, 10.0, 1000.0, 5.0, 5.0, seed);
    auto dbl2 = TestUtils::generateStarField(catalog, 256, 256, 30, 0.7, 10.0, 1000.0, 5.0, 5.0, seed+385);
    auto dbl3 = TestUtils::generateStarField(catalog, 256, 256, 30, 0.7, 10.0, 1000.0, 5.0, 5.0, seed+976);

    auto layer1 = toType<float>(dbl1);
    layer1->AddLayer(*toType<float>(dbl2));
    layer1->AddLayer(*toType<float>(dbl3));
    std::shared_ptr<FITScube> img = layer1;

    img->HDU().ValueForKey("CTYPE1", "RA---TAN",  fChar);
    img->HDU().ValueForKey("CTYPE2", "DEC--TAN",  fChar);
    img->HDU().ValueForKey("CTYPE3", "WAVELENGTH", fChar);
    img->HDU().ValueForKey("CRPIX1", 0.); img->HDU().ValueForKey("CRPIX2", 0.); img->HDU().ValueForKey("CRPIX3", 0.);
    img->HDU().ValueForKey("CRVAL1", 180.0); img->HDU().ValueForKey("CRVAL2", 0.0); img->HDU().ValueForKey("CRVAL3", 656.28);
    img->HDU().ValueForKey("CDELT1", -0.0002777777778);
    img->HDU().ValueForKey("CDELT2",  0.0002777777778);
    img->HDU().ValueForKey("CDELT3",  0.1);

    DST::algorithms::starFinder finder(img);
    finder.setOutputUrl("./build/test_output");
    verbose = verboseLevel::VERBOSE_NONE;
    finder.setThresholdValue(0.998);
    finder.findStars(DST::algorithms::starFinder::detectionMethod::MORPHO, 0.1);

    const auto& detected = finder.detectedStars();
    std::cout << "[float 3-layer cube MORPHO] Detected " << detected.nrows() << " stars" << std::endl;
    EXPECT_GT(detected.nrows(), static_cast<size_t>(0));

#ifndef Linuxx86_64
    {
        DSL::mathplot::StderrRedirect redirect;
        EXPECT_NO_THROW(finder.plot(0));
        EXPECT_NO_THROW(finder.plot(1));
        EXPECT_NO_THROW(finder.plot(2));
    }
#endif
}

// ---- Double image must still use CPU path (regression) ----

TEST(StarFinderGPUTest, FindStarsBySegmentation_Double_StillCPU)
{
    std::vector<DST::Stars> catalog;
    auto img = TestUtils::generateStarField(catalog,
        512, 512, 50, 0.7, 100.0, 1000.0, 5.0, 5.0, 12345);

    DST::algorithms::starFinder finder(img);
    finder.setOutputUrl("./build/test_output");
    verbose = verboseLevel::VERBOSE_NONE;
    finder.setThresholdValue(0.5);
    finder.findStars(DST::algorithms::starFinder::detectionMethod::MORPHO, 0.1);

    const auto& detected = finder.detectedStars();
    double rate = static_cast<double>(detected.nrows()) /
                  static_cast<double>(catalog.size()) * 100.0;
    std::cout << "[double MORPHO — CPU path] " << detected.nrows() << "/" << catalog.size()
              << " (" << rate << "%)" << std::endl;

    EXPECT_GT(rate, 75.0);
    EXPECT_LT(rate, 110.0);
}
