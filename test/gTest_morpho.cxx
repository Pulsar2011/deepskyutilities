#include <gtest/gtest.h>
#include <DSTools/DSTmorpho.h>
#include <valarray>
#include <vector>

using namespace DST;
using namespace Morphology;

// Define all types to test
typedef ::testing::Types<uint8_t, uint16_t, uint32_t, uint64_t, 
                         int8_t, int16_t, int32_t, int64_t, 
                         float, double> NumericTypes;

// Test fixture for constructor tests
template <typename T>
class KernelConstructorTest : public ::testing::Test
{
    protected:
        std::valarray<T> test_data;
        std::vector<size_t> test_size;

        void SetUp() override
        {
            test_data = {T(1), T(2), T(3), T(4), T(5), T(6)};
            test_size = {2, 3};
        }
};

TYPED_TEST_SUITE(KernelConstructorTest, NumericTypes);

TYPED_TEST(KernelConstructorTest, VectorConstructor)
{
    kernel<TypeParam> k(this->test_data, this->test_size);
    
    EXPECT_EQ(k.data().size(), 6);
    EXPECT_EQ(k.sizes().size(), 2);
    EXPECT_EQ(k.size(0), 2);
    EXPECT_EQ(k.size(1), 3);
    EXPECT_EQ(k.data()[0], TypeParam(1));
    EXPECT_EQ(k.data()[5], TypeParam(6));
}

TYPED_TEST(KernelConstructorTest, InitializerListConstructor)
{
    kernel<TypeParam> k(this->test_data, {2, 3});
    
    EXPECT_EQ(k.data().size(), 6);
    EXPECT_EQ(k.sizes().size(), 2);
    EXPECT_EQ(k.size(0), 2);
    EXPECT_EQ(k.size(1), 3);
}

TYPED_TEST(KernelConstructorTest, CopyConstructor)
{
    kernel<TypeParam> k1(this->test_data, this->test_size);
    kernel<TypeParam> k2(k1);
    
    EXPECT_EQ(k2.data().size(), k1.data().size());
    EXPECT_EQ(k2.sizes().size(), k1.sizes().size());
    EXPECT_EQ(k2.size(0), k1.size(0));
    EXPECT_EQ(k2.size(1), k1.size(1));
    
    for(size_t i = 0; i < k1.data().size(); ++i) {
        EXPECT_EQ(k2.data()[i], k1.data()[i]);
    }
}

TYPED_TEST(KernelConstructorTest, InvalidSizeThrows)
{
    std::valarray<TypeParam> bad_data = {TypeParam(1), TypeParam(2), TypeParam(3)};
    std::vector<size_t> bad_size = {2, 3}; // Size mismatch: 3 elements vs 2*3=6
    
    EXPECT_THROW(kernel<TypeParam>(bad_data, bad_size), std::invalid_argument);
}

TYPED_TEST(KernelConstructorTest, OutOfRangeSizeAccess)
{
    kernel<TypeParam> k(this->test_data, this->test_size);
    
    EXPECT_THROW(k.size(5), std::out_of_range);
}

// Test fixture for operator tests
template <typename T>
class KernelOperatorTest : public ::testing::Test
{
    protected:
        std::valarray<T> data1;
        std::valarray<T> data2;
        std::vector<size_t> size_vec;

        void SetUp() override
        {
            data1 = {T(2), T(4), T(6), T(8)};
            data2 = {T(1), T(2), T(3), T(4)};
            size_vec = {2, 2};
        }
};

TYPED_TEST_SUITE(KernelOperatorTest, NumericTypes);

TYPED_TEST(KernelOperatorTest, AssignmentOperator)
{
    kernel<TypeParam> k1(this->data1, this->size_vec);
    kernel<TypeParam> k2(this->data2, this->size_vec);
    
    k2 = k1;
    
    EXPECT_EQ(k2.data().size(), k1.data().size());
    EXPECT_EQ(k2.type(), k1.type());
    for(size_t i = 0; i < k1.data().size(); ++i) {
        EXPECT_EQ(k2.data()[i], k1.data()[i]);
    }
}

TYPED_TEST(KernelOperatorTest, KernelAdditionOperator)
{
    kernel<TypeParam> k1(this->data1, this->size_vec);
    kernel<TypeParam> k2(this->data2, this->size_vec);
    
    k1 += k2;
    
    EXPECT_EQ(k1.data()[0], TypeParam(3));
    EXPECT_EQ(k1.data()[1], TypeParam(6));
    EXPECT_EQ(k1.data()[2], TypeParam(9));
    EXPECT_EQ(k1.data()[3], TypeParam(12));
}

TYPED_TEST(KernelOperatorTest, KernelSubtractionOperator)
{
    kernel<TypeParam> k1(this->data1, this->size_vec);
    kernel<TypeParam> k2(this->data2, this->size_vec);
    
    k1 -= k2;
    
    EXPECT_EQ(k1.data()[0], TypeParam(1));
    EXPECT_EQ(k1.data()[1], TypeParam(2));
    EXPECT_EQ(k1.data()[2], TypeParam(3));
    EXPECT_EQ(k1.data()[3], TypeParam(4));
}

TYPED_TEST(KernelOperatorTest, KernelMultiplicationOperator)
{
    kernel<TypeParam> k1(this->data1, this->size_vec);
    kernel<TypeParam> k2(this->data2, this->size_vec);
    
    k1 *= k2;
    
    EXPECT_EQ(k1.data()[0], TypeParam(2));
    EXPECT_EQ(k1.data()[1], TypeParam(8));
    EXPECT_EQ(k1.data()[2], TypeParam(18));
    EXPECT_EQ(k1.data()[3], TypeParam(32));
}

TYPED_TEST(KernelOperatorTest, KernelDivisionOperator)
{
    kernel<TypeParam> k1(this->data1, this->size_vec);
    kernel<TypeParam> k2(this->data2, this->size_vec);
    
    k1 /= k2;
    
    EXPECT_EQ(k1.data()[0], TypeParam(2));
    EXPECT_EQ(k1.data()[1], TypeParam(2));
    EXPECT_EQ(k1.data()[2], TypeParam(2));
    EXPECT_EQ(k1.data()[3], TypeParam(2));
}

TYPED_TEST(KernelOperatorTest, KernelOperatorSizeMismatchThrows)
{
    kernel<TypeParam> k1(this->data1, this->size_vec);
    std::valarray<TypeParam> small_data = {TypeParam(1), TypeParam(2)};
    kernel<TypeParam> k2(small_data, {1, 2});
    
    EXPECT_THROW(k1 += k2, std::invalid_argument);
    EXPECT_THROW(k1 -= k2, std::invalid_argument);
    EXPECT_THROW(k1 *= k2, std::invalid_argument);
    EXPECT_THROW(k1 /= k2, std::invalid_argument);
}

TYPED_TEST(KernelOperatorTest, ArrayAdditionOperator)
{
    kernel<TypeParam> k(this->data1, this->size_vec);
    
    k += this->data2;
    
    EXPECT_EQ(k.data()[0], TypeParam(3));
    EXPECT_EQ(k.data()[3], TypeParam(12));
}

TYPED_TEST(KernelOperatorTest, ArraySubtractionOperator)
{
    kernel<TypeParam> k(this->data1, this->size_vec);
    
    k -= this->data2;
    
    EXPECT_EQ(k.data()[0], TypeParam(1));
    EXPECT_EQ(k.data()[3], TypeParam(4));
}

TYPED_TEST(KernelOperatorTest, ArrayMultiplicationOperator)
{
    kernel<TypeParam> k(this->data1, this->size_vec);
    
    k *= this->data2;
    
    EXPECT_EQ(k.data()[0], TypeParam(2));
    EXPECT_EQ(k.data()[3], TypeParam(32));
}

TYPED_TEST(KernelOperatorTest, ArrayDivisionOperator)
{
    kernel<TypeParam> k(this->data1, this->size_vec);
    
    k /= this->data2;
    
    EXPECT_EQ(k.data()[0], TypeParam(2));
    EXPECT_EQ(k.data()[3], TypeParam(2));
}

TYPED_TEST(KernelOperatorTest, ArrayOperatorSizeMismatchThrows)
{
    kernel<TypeParam> k(this->data1, this->size_vec);
    std::valarray<TypeParam> small_array = {TypeParam(1), TypeParam(2)};
    
    EXPECT_THROW(k += small_array, std::invalid_argument);
    EXPECT_THROW(k -= small_array, std::invalid_argument);
    EXPECT_THROW(k *= small_array, std::invalid_argument);
    EXPECT_THROW(k /= small_array, std::invalid_argument);
}

TYPED_TEST(KernelOperatorTest, ScalarAdditionOperator)
{
    kernel<TypeParam> k(this->data1, this->size_vec);
    
    k += TypeParam(2);
    
    EXPECT_EQ(k.data()[0], TypeParam(4));
    EXPECT_EQ(k.data()[3], TypeParam(10));
}

TYPED_TEST(KernelOperatorTest, ScalarSubtractionOperator)
{
    kernel<TypeParam> k(this->data1, this->size_vec);
    
    k -= TypeParam(1);
    
    EXPECT_EQ(k.data()[0], TypeParam(1));
    EXPECT_EQ(k.data()[3], TypeParam(7));
}

TYPED_TEST(KernelOperatorTest, ScalarMultiplicationOperator)
{
    kernel<TypeParam> k(this->data1, this->size_vec);
    
    k *= TypeParam(2);
    
    EXPECT_EQ(k.data()[0], TypeParam(4));
    EXPECT_EQ(k.data()[3], TypeParam(16));
}

TYPED_TEST(KernelOperatorTest, ScalarDivisionOperator)
{
    kernel<TypeParam> k(this->data1, this->size_vec);
    
    k /= TypeParam(2);
    
    EXPECT_EQ(k.data()[0], TypeParam(1));
    EXPECT_EQ(k.data()[3], TypeParam(4));
}

// Test fixture for static methods
template <typename T>
class KernelStaticMethodTest : public ::testing::Test {};

TYPED_TEST_SUITE(KernelStaticMethodTest, NumericTypes);

TYPED_TEST(KernelStaticMethodTest, OneKernelVector)
{
    std::vector<size_t> size_vec = {3, 4};
    auto k = kernel<TypeParam>::one(size_vec);
    
    EXPECT_EQ(k.data().size(), 12);
    EXPECT_EQ(k.size(0), 3);
    EXPECT_EQ(k.size(1), 4);
    
    for(size_t i = 0; i < k.data().size(); ++i) {
        EXPECT_EQ(k.data()[i], TypeParam(1));
    }
}

TYPED_TEST(KernelStaticMethodTest, OneKernelInitializerList)
{
    auto k = kernel<TypeParam>::one({2, 5});
    
    EXPECT_EQ(k.data().size(), 10);
    
    for(size_t i = 0; i < k.data().size(); ++i) {
        EXPECT_EQ(k.data()[i], TypeParam(1));
    }
}

TYPED_TEST(KernelStaticMethodTest, ZeroKernelVector)
{
    std::vector<size_t> size_vec = {3, 4};
    auto k = kernel<TypeParam>::zero(size_vec);
    
    EXPECT_EQ(k.data().size(), 12);
    EXPECT_EQ(k.size(0), 3);
    EXPECT_EQ(k.size(1), 4);
    
    for(size_t i = 0; i < k.data().size(); ++i) {
        EXPECT_EQ(k.data()[i], TypeParam(0));
    }
}

TYPED_TEST(KernelStaticMethodTest, ZeroKernelInitializerList)
{
    auto k = kernel<TypeParam>::zero({2, 5});
    
    EXPECT_EQ(k.data().size(), 10);
    
    for(size_t i = 0; i < k.data().size(); ++i) {
        EXPECT_EQ(k.data()[i], TypeParam(0));
    }
}

// New tests for folding
template <typename T>
class KernelFoldingTest : public ::testing::Test {};

TYPED_TEST_SUITE(KernelFoldingTest, NumericTypes);

TYPED_TEST(KernelFoldingTest, DimensionMatching)
{
    using T = TypeParam;

    // Create 5x5 image initialized to zero
    size_t maxSize = 11;
    size_t depth = 5;
    auto img = std::make_shared< FITSimg<T> >(std::vector<size_t>{maxSize,maxSize,depth});
    auto data_ptr = img->template GetData<T>();
    (*data_ptr) *= 0;

    for(size_t i = 0; i < maxSize; i+=5)
    for(size_t j = 0; j < maxSize; j+=5)
    for(size_t k = 0; k < depth; k++)
        (*data_ptr)[img->PixelIndex({i,j,k})] = T(1);

    // 3x3 cross kernel: [0,1,0, 1,1,1, 0,1,0]
    std::valarray<T> kdata = { T(0), T(1), T(0),
                               T(1), T(1), T(1),
                               T(0), T(1), T(0) };
    kernel<T> k(kdata, std::vector<size_t>{3,3});

    // Call folding wrapper
    EXPECT_ANY_THROW(folding<T>(std::dynamic_pointer_cast< FITScube >(img), k, 1));
}

TYPED_TEST(KernelFoldingTest, FoldDeltaWithCrossKernel_OneThreads)
{
    using T = TypeParam;

    // Create 5x5 image initialized to zero
    size_t maxSize = 101;
    auto img = std::make_shared< FITSimg<T> >(std::vector<size_t>{maxSize,maxSize});
    auto data_ptr = img->template GetData<T>();
    (*data_ptr) *= 0;

    for(size_t i = 5; i <= maxSize-5; i+=5)
    for(size_t j = 5; j <= maxSize-5; j+=5)
    {
        (*data_ptr)[img->PixelIndex({i,j})] = T(1);
        if(i==(maxSize-1)/2 && j==(maxSize-1)/2)
            img->MaskPixels( {img->PixelIndex({i,j})} ); // Mask center pixel
    }

    // 3x3 cross kernel: [0,1,0, 1,1,1, 0,1,0]
    std::valarray<T> kdata = { T(0), T(1), T(0),
                               T(1), T(1), T(1),
                               T(0), T(1), T(0) };
    kernel<T> k(kdata, std::vector<size_t>{3,3});

    // Call folding wrapper
    auto out_ptr = folding<T>(std::dynamic_pointer_cast< FITScube >(img), k, 1);

    // Cast result back to image of same pixel type
    auto out_img = std::dynamic_pointer_cast< FITSimg<T> >(out_ptr);
    ASSERT_NE(out_img, nullptr);

    auto out_data = out_img->template GetData<T>();
    EXPECT_EQ(out_data->size(), static_cast<size_t>(maxSize*maxSize));

    // Expect ones at center and its 4-neighbors (linear indices)
    std::vector<size_t> expected_ones;
    for(size_t i = 5; i <= maxSize-5; i+=5)
    for(size_t j = 5; j <= maxSize-5; j+=5)
    {
        if(i==(maxSize-1)/2 && j==(maxSize-1)/2)
            continue; // Skip center pixel which is masked

        expected_ones.push_back(img->PixelIndex({i  ,j  }));
        expected_ones.push_back(img->PixelIndex({i-1,j  }));
        expected_ones.push_back(img->PixelIndex({i+1,j  }));
        expected_ones.push_back(img->PixelIndex({i  ,j-1}));
        expected_ones.push_back(img->PixelIndex({i  ,j+1}));
    }

    for(size_t i = 0; i < out_data->size(); ++i)
    {
        bool should_be_one = (std::find(expected_ones.begin(), expected_ones.end(), i) != expected_ones.end());
        if(should_be_one)
            EXPECT_EQ((*out_data)[i], T(1));
        else
            EXPECT_EQ((*out_data)[i], T(0));
    }

    EXPECT_TRUE(out_img->HDU().Exists("DST_PROC"));
    EXPECT_TRUE(out_img->HDU().Exists("DST_VER"));
    EXPECT_TRUE(out_img->HDU().Exists("DST_SHA"));
    EXPECT_EQ( out_img->HDU().GetValueForKey("DST_PROC"), "FOLDING");
    EXPECT_EQ( out_img->HDU().GetValueForKey("DST_VER"), DST::gGIT::version_short());
    EXPECT_EQ( out_img->HDU().GetValueForKey("DST_SHA"), DST::gGIT::sha1());

    out_img->Write("./build/KernelFoldingOne_"+demangle(typeid(T).name())+".fits",true);
}

TYPED_TEST(KernelFoldingTest, FoldDeltaWithCrossKernel_MultiThreads)
{
    using T = TypeParam;

    // Create 5x5 image initialized to zero
    size_t maxSize = 101;
    auto img = std::make_shared< FITSimg<T> >(std::vector<size_t>{maxSize,maxSize});
    auto data_ptr = img->template GetData<T>();
    (*data_ptr) *= 0;

    for(size_t i = 0; i < maxSize; i+=5)
    for(size_t j = 0; j < maxSize; j+=5)
    {
        (*data_ptr)[img->PixelIndex({i,j})] = T(1);
        if(i==(maxSize-1)/2 && j==(maxSize-1)/2)
            img->MaskPixels( {img->PixelIndex({i,j})} ); // Mask center pixel
    }

    // 3x3 cross kernel: [0,1,0, 1,1,1, 0,1,0]
    std::valarray<T> kdata = { T(0), T(1), T(0),
                               T(1), T(1), T(1),
                               T(0), T(1), T(0) };
    kernel<T> k(kdata, std::vector<size_t>{3,3});

    // Call folding wrapper
    auto out_ptr = folding<T>(std::dynamic_pointer_cast< FITScube >(img), k, 8);

    // Cast result back to image of same pixel type
    auto out_img = std::dynamic_pointer_cast< FITSimg<T> >(out_ptr);
    ASSERT_NE(out_img, nullptr);

    auto out_data = out_img->template GetData<T>();
    EXPECT_EQ(out_data->size(), static_cast<size_t>(maxSize*maxSize));

    // Expect ones at center and its 4-neighbors (linear indices)
    std::vector<size_t> expected_ones;
    for(size_t i = 0; i < maxSize; i+=5)
    for(size_t j = 0; j < maxSize; j+=5)
    {
        if(i==(maxSize-1)/2 && j==(maxSize-1)/2)
            continue; // Skip center pixel which is masked

        expected_ones.push_back(img->PixelIndex({i  ,j  }));
        expected_ones.push_back(img->PixelIndex({i-1,j  }));
        expected_ones.push_back(img->PixelIndex({i+1,j  }));
        expected_ones.push_back(img->PixelIndex({i  ,j-1}));
        expected_ones.push_back(img->PixelIndex({i  ,j+1}));
    }

    for(size_t i = 0; i < out_data->size(); ++i)
    {
        bool should_be_one = (std::find(expected_ones.begin(), expected_ones.end(), i) != expected_ones.end());
        if(should_be_one)
            EXPECT_EQ((*out_data)[i], T(1));
        else
            EXPECT_EQ((*out_data)[i], T(0));
    }

    EXPECT_TRUE(out_img->HDU().Exists("DST_PROC"));
    EXPECT_TRUE(out_img->HDU().Exists("DST_VER"));
    EXPECT_TRUE(out_img->HDU().Exists("DST_SHA"));
    EXPECT_EQ( out_img->HDU().GetValueForKey("DST_PROC"), "FOLDING");
    EXPECT_EQ( out_img->HDU().GetValueForKey("DST_VER"), DST::gGIT::version_short());
    EXPECT_EQ( out_img->HDU().GetValueForKey("DST_SHA"), DST::gGIT::sha1());
    
    out_img->Write("./build/KernelFoldingMultiple_"+demangle(typeid(T).name())+".fits",true);
}

TYPED_TEST(KernelFoldingTest, FoldDeltaWithCrossKernel_3D)
{
    using T = TypeParam;

    // Create 5x5 image initialized to zero
    size_t maxSize = 11;
    size_t depth = 5;
    auto img = std::make_shared< FITSimg<T> >(std::vector<size_t>{maxSize,maxSize,depth});
    auto data_ptr = img->template GetData<T>();
    (*data_ptr) *= 0;

    for(size_t i = 0; i < maxSize; i+=5)
    for(size_t j = 0; j < maxSize; j+=5)
    for(size_t k = 0; k < depth; k++)
        (*data_ptr)[img->PixelIndex({i,j,k})] = T(1);

    // 3x3 cross kernel: [0,1,0, 1,1,1, 0,1,0]
    std::valarray<T> kdata = { T(0), T(1), T(0),
                               T(1), T(1), T(1),
                               T(0), T(1), T(0) };
    kernel<T> k(kdata, std::vector<size_t>{3,3,1});

    // Call folding wrapper
    auto out_ptr = folding<T>(std::dynamic_pointer_cast< FITScube >(img), k, 8);

    // Cast result back to image of same pixel type
    auto out_img = std::dynamic_pointer_cast< FITSimg<T> >(out_ptr);
    ASSERT_NE(out_img, nullptr);

    auto out_data = out_img->template GetData<T>();
    EXPECT_EQ(out_data->size(), static_cast<size_t>(maxSize*maxSize*depth));

    // Expect ones at center and its 4-neighbors (linear indices)
    std::vector<size_t> expected_ones;
    for(size_t i = 0; i < maxSize; i+=5)
    for(size_t j = 0; j < maxSize; j+=5)
    for(size_t k = 0; k < depth; k++)
    {
        expected_ones.push_back(img->PixelIndex({i  ,j  ,k}));
        expected_ones.push_back(img->PixelIndex({i-1,j  ,k}));
        expected_ones.push_back(img->PixelIndex({i+1,j  ,k}));
        expected_ones.push_back(img->PixelIndex({i  ,j-1,k}));
        expected_ones.push_back(img->PixelIndex({i  ,j+1,k}));
    }

    for(size_t i = 0; i < out_data->size(); ++i)
    {
        bool should_be_one = (std::find(expected_ones.begin(), expected_ones.end(), i) != expected_ones.end());
        if(should_be_one)
            EXPECT_EQ((*out_data)[i], T(1));
        else
            EXPECT_EQ((*out_data)[i], T(0));
    }

    EXPECT_TRUE(out_img->HDU().Exists("DST_PROC"));
    EXPECT_TRUE(out_img->HDU().Exists("DST_VER"));
    EXPECT_TRUE(out_img->HDU().Exists("DST_SHA"));
    EXPECT_EQ( out_img->HDU().GetValueForKey("DST_PROC"), "FOLDING");
    EXPECT_EQ( out_img->HDU().GetValueForKey("DST_VER"), DST::gGIT::version_short());
    EXPECT_EQ( out_img->HDU().GetValueForKey("DST_SHA"), DST::gGIT::sha1());

    out_img->Write("./build/KernelFolding3D_"+demangle(typeid(T).name())+".fits",true);
}

TEST(MorphologyCorrelateTest, OnesKernel_NoNormalize_And_Normalize)
{
    // Build a simple deterministic 5x5 image: value = x + 10*y (y=0 bottom)
    const size_t w = 5;
    const size_t h = 5;
    auto img = std::make_shared< FITSimg<double> >(std::vector<size_t>{w,h});
    auto data_ptr = img->template GetData<double>();
    (*data_ptr) *= 0.0;

    for(size_t y = 0; y < h; ++y)
    for(size_t x = 0; x < w; ++x)
        (*data_ptr)[img->PixelIndex({x,y})] = static_cast<double>(x) + 10.0*static_cast<double>(y);

    // 3x3 ones kernel
    std::valarray<double> kdata = {
        1,1,1,
        1,1,1,
        1,1,1
    };
    kernel<double> k(kdata, std::vector<size_t>{3,3});

    // Correlate without normalization: central pixel should be 3x3 sum
    auto out_sum_ptr = correlate<double>(std::dynamic_pointer_cast<FITScube>(img), k, 1, false);
    auto out_sum_img = std::dynamic_pointer_cast< FITSimg<double> >(out_sum_ptr);
    ASSERT_NE(out_sum_img, nullptr);

    const size_t cx = 2;
    const size_t cy = 2;
    double expected_sum = 0.0;
    for(int dy = -1; dy <= 1; ++dy)
    for(int dx = -1; dx <= 1; ++dx)
        expected_sum += (*data_ptr)[img->PixelIndex({static_cast<size_t>(cx+dx), static_cast<size_t>(cy+dy)})];

    EXPECT_DOUBLE_EQ((*out_sum_img->GetData<double>())[out_sum_img->PixelIndex({cx,cy})], expected_sum);
    EXPECT_TRUE(out_sum_img->HDU().Exists("DST_PROC"));
    EXPECT_EQ(out_sum_img->HDU().GetValueForKey("DST_PROC"), "CORRELATE");
    EXPECT_TRUE(out_sum_img->HDU().Exists("DST_NORM"));
    EXPECT_EQ(out_sum_img->HDU().GetValueForKey("DST_NORM"), "F");

    // Correlate with normalization: central pixel should be 3x3 mean (divide by 9)
    auto out_mean_ptr = correlate<double>(std::dynamic_pointer_cast<FITScube>(img), k, 1, true);
    auto out_mean_img = std::dynamic_pointer_cast< FITSimg<double> >(out_mean_ptr);
    ASSERT_NE(out_mean_img, nullptr);
    EXPECT_DOUBLE_EQ((*out_mean_img->GetData<double>())[out_mean_img->PixelIndex({cx,cy})], expected_sum / 9.0);
    EXPECT_EQ(out_mean_img->HDU().GetValueForKey("DST_NORM"), "T");

    // Border behavior with normalization: only overlapping pixels count
    const size_t bx = 0;
    const size_t by = 0;
    double border_sum = 0.0;
    double border_count = 0.0;
    for(int dy = -1; dy <= 1; ++dy)
    for(int dx = -1; dx <= 1; ++dx)
    {
        const int tx = static_cast<int>(bx) + dx;
        const int ty = static_cast<int>(by) + dy;
        if(tx < 0 || ty < 0 || tx >= static_cast<int>(w) || ty >= static_cast<int>(h))
            continue;
        border_sum += (*data_ptr)[img->PixelIndex({static_cast<size_t>(tx), static_cast<size_t>(ty)})];
        border_count += 1.0;
    }
    EXPECT_DOUBLE_EQ((*out_mean_img->GetData<double>())[out_mean_img->PixelIndex({bx,by})], border_sum / border_count);
}

TEST(MorphologyCorrelateTest, NormalizeSkipsMaskedPixels)
{
    const size_t w = 5;
    const size_t h = 5;
    auto img = std::make_shared< FITSimg<double> >(std::vector<size_t>{w,h});
    auto data_ptr = img->template GetData<double>();
    (*data_ptr) *= 0.0;

    for(size_t y = 0; y < h; ++y)
    for(size_t x = 0; x < w; ++x)
        (*data_ptr)[img->PixelIndex({x,y})] = 1.0;

    // Mask one neighbor in the 3x3 around the center
    const size_t cx = 2;
    const size_t cy = 2;
    img->MaskPixels({img->PixelIndex({cx+1, cy})});

    std::valarray<double> kdata = {
        1,1,1,
        1,1,1,
        1,1,1
    };
    kernel<double> k(kdata, std::vector<size_t>{3,3});

    auto out_ptr = correlate<double>(std::dynamic_pointer_cast<FITScube>(img), k, 1, true);
    auto out_img = std::dynamic_pointer_cast< FITSimg<double> >(out_ptr);
    ASSERT_NE(out_img, nullptr);

    // With all ones, normalized correlation is the mean of valid samples.
    // 3x3 neighborhood would have 9 samples, but one is masked => 8 valid => expected mean = 1.
    EXPECT_DOUBLE_EQ((*out_img->GetData<double>())[out_img->PixelIndex({cx,cy})], 1.0);
}

TYPED_TEST(KernelFoldingTest, Dilation)
{
    using T = TypeParam;

    // Use explicit 11x11 input/expected arrays (rows provided top-to-bottom).
    const size_t width  = 11;
    const size_t height = 11;

    auto img = std::make_shared< FITSimg<T> >(std::vector<size_t>{width,height});
    auto data_ptr = img->template GetData<T>();
    (*data_ptr) *= 0;

    // Input rows (top-to-bottom)
    const int inputRows[11][11] = {
        {0,0,0,0,0,0,0,0,0,0,0},
        {0,1,1,1,1,0,0,1,1,1,0},
        {0,1,1,1,1,0,0,1,1,1,0},
        {0,1,1,1,1,1,1,1,1,1,0},
        {0,1,1,1,1,1,1,1,1,1,0},
        {0,1,1,0,0,0,1,1,1,1,0},
        {0,1,1,0,0,0,1,1,1,1,0},
        {0,1,1,0,0,0,1,1,1,1,0},
        {0,1,1,1,1,1,1,1,0,0,0},
        {0,1,1,1,1,1,1,1,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0}
    };

    // Expected output rows (top-to-bottom)
    const int expectedRows[11][11] = {
        {1,1,1,1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1,1,1,1},
        {1,1,1,1,0,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1,1,0,0},
        {1,1,1,1,1,1,1,1,1,0,0}
    };

    // Fill image: inputRows are top-to-bottom; PixelIndex expects y=0 as bottom => map accordingly.
    for(size_t row = 0; row < height; ++row)
    {
        for(size_t col = 0; col < width; ++col)
        {
            const int v = inputRows[row][col];
            size_t x = col;
            size_t y = height - 1 - row;
            (*data_ptr)[ img->PixelIndex({x,y}) ] = static_cast<T>(v);
        }
    }

    std::valarray<T> kdata = { T(1), T(1), T(1),
                               T(1), T(1), T(1),
                               T(1), T(1), T(1) };
    kernel<T> k(kdata, std::vector<size_t>{3,3});

    // Call dilation
    auto out_ptr = dilation<T>(std::dynamic_pointer_cast< FITScube >(img), k, 1);

    auto out_img = std::dynamic_pointer_cast< FITSimg<T> >(out_ptr);
    ASSERT_NE(out_img, nullptr);

    auto out_data = out_img->template GetData<T>();
    EXPECT_EQ(out_data->size(), static_cast<size_t>(width*height));

    // Validate output against expectedRows (top-to-bottom)
    for(size_t row = 0; row < height; ++row)
    {
        for(size_t col = 0; col < width; ++col)
        {
            size_t x = col;
            size_t y = height - 1 - row;
            size_t idx = out_img->PixelIndex({x,y});
            T expected = static_cast<T>( expectedRows[row][col] );
            EXPECT_EQ( (*out_data)[idx], expected ) << "Mismatch at (col=" << col << ", row(top)=" << row << ")";
        }
    }

    EXPECT_TRUE(out_img->HDU().Exists("DST_PROC"));
    EXPECT_TRUE(out_img->HDU().Exists("DST_VER"));
    EXPECT_TRUE(out_img->HDU().Exists("DST_SHA"));
    EXPECT_EQ( out_img->HDU().GetValueForKey("DST_PROC"), "DILATION");
    EXPECT_EQ( out_img->HDU().GetValueForKey("DST_VER"), DST::gGIT::version_short());
    EXPECT_EQ( out_img->HDU().GetValueForKey("DST_SHA"), DST::gGIT::sha1());

    out_img->Write("./build/KernelDilationOne_"+demangle(typeid(T).name())+".fits",true);
}

TYPED_TEST(KernelFoldingTest, Erosion)
{
    using T = TypeParam;

    // Use explicit 13x13 input/expected arrays (rows provided top-to-bottom).
    const size_t width  = 13;
    const size_t height = 13;

    auto img = std::make_shared< FITSimg<T> >(std::vector<size_t>{width,height});
    auto data_ptr = img->template GetData<T>();
    (*data_ptr) *= 0;

    // Input rows (top-to-bottom)
    const int inputRows[13][13] = {
        {1,1,1,1,1,1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,0,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1,1,1,1,1,1},
        {0,1,1,1,1,1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1,1,1,1,1,1}
    };

    // Expected output rows (top-to-bottom)
    const int expectedRows[13][13] = {
        {1,1,1,1,1,0,0,0,1,1,1,1,1},
        {1,1,1,1,1,0,0,0,1,1,1,1,1},
        {1,1,1,1,1,0,0,0,1,1,1,1,1},
        {1,1,1,1,1,1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1,1,1,1,1,1},
        {0,0,1,1,1,1,1,1,1,1,1,1,1},
        {0,0,1,1,1,1,1,1,1,1,1,1,1},
        {0,0,1,1,1,1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1,1,1,1,1,1}
    };

    // Fill image: inputRows are top-to-bottom; PixelIndex expects y=0 bottom => map accordingly.
    for(size_t row = 0; row < height; ++row)
    {
        for(size_t col = 0; col < width; ++col)
        {
            const int v = inputRows[row][col];
            size_t x = col;
            size_t y = height - 1 - row;
            (*data_ptr)[ img->PixelIndex({x,y}) ] = static_cast<T>(v);
        }
    }

    std::valarray<T> kdata = { T(1), T(1), T(1),
                               T(1), T(1), T(1),
                               T(1), T(1), T(1) };
    kernel<T> k(kdata, std::vector<size_t>{3,3});

    // Call erosion
    auto out_ptr = erosion<T>(std::dynamic_pointer_cast< FITScube >(img), k, 1);

    auto out_img = std::dynamic_pointer_cast< FITSimg<T> >(out_ptr);
    ASSERT_NE(out_img, nullptr);

    auto out_data = out_img->template GetData<T>();
    EXPECT_EQ(out_data->size(), static_cast<size_t>(width*height));

    // Validate output against expectedRows (top-to-bottom)
    for(size_t row = 0; row < height; ++row)
    {
        for(size_t col = 0; col < width; ++col)
        {
            size_t x = col;
            size_t y = height - 1 - row;
            size_t idx = out_img->PixelIndex({x,y});
            T expected = static_cast<T>( expectedRows[row][col] );
            EXPECT_EQ( (*out_data)[idx], expected ) << "Mismatch at (col=" << col << ", row(top)=" << row << ")";
        }
    }

    EXPECT_TRUE(out_img->HDU().Exists("DST_PROC"));
    EXPECT_TRUE(out_img->HDU().Exists("DST_VER"));
    EXPECT_TRUE(out_img->HDU().Exists("DST_SHA"));
    EXPECT_EQ( out_img->HDU().GetValueForKey("DST_PROC"), "EROSION");
    EXPECT_EQ( out_img->HDU().GetValueForKey("DST_VER"), DST::gGIT::version_short());
    EXPECT_EQ( out_img->HDU().GetValueForKey("DST_SHA"), DST::gGIT::sha1());

    out_img->Write("./build/KernelErosionOne_"+demangle(typeid(T).name())+".fits",true);
}

TYPED_TEST(KernelFoldingTest, Opening)
{
    using T = TypeParam;

    // Replace previous random setup with the explicit input/expected arrays (provided top-to-bottom).
    const size_t width  = 16;
    const size_t height = 15;
    size_t maxSizeX = width;
    size_t maxSizeY = height;

    auto img = std::make_shared< FITSimg<T> >(std::vector<size_t>{maxSizeX,maxSizeY});
    auto data_ptr = img->template GetData<T>();
    (*data_ptr) *= 0;

    // Input rows given top-to-bottom (row 0 = top)
    const int inputRows[15][16] = {
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0},
        {0,0,1,1,1,1,1,0,0,0,0,1,1,1,0,0},
        {0,0,1,1,1,1,1,0,0,0,1,1,1,1,0,0},
        {0,0,1,1,1,1,0,0,0,1,1,1,1,1,0,0},
        {0,0,0,1,1,0,0,0,1,1,1,1,1,0,0,0},
        {0,0,0,0,0,0,0,1,1,1,1,1,0,0,0,0},
        {0,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0},
        {0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0},
        {0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0},
        {0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0},
        {0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0},
        {0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
    };

    // Expected rows given top-to-bottom
    const int expectedRows[15][16] = {
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0},
        {0,0,1,1,1,1,0,0,0,0,0,1,1,1,0,0},
        {0,0,1,1,1,1,0,0,0,0,1,1,1,1,0,0},
        {0,0,1,1,1,1,0,0,0,1,1,1,1,1,0,0},
        {0,0,0,0,0,0,0,0,1,1,1,1,1,0,0,0},
        {0,0,0,0,0,0,0,1,1,1,1,1,0,0,0,0},
        {0,0,0,0,0,0,1,1,1,1,1,0,0,0,0,0},
        {0,0,0,0,0,1,1,1,1,1,0,0,0,0,0,0},
        {0,0,0,0,1,1,1,1,1,0,0,0,0,0,0,0},
        {0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0},
        {0,0,0,0,1,1,1,1,1,1,1,0,0,0,0,0},
        {0,0,0,0,0,1,1,1,1,1,1,0,0,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0}
    };

    // Fill image: provided rows are top-to-bottom, PixelIndex expects y=0 bottom row -> map accordingly.
    for(size_t row = 0; row < height; ++row)
    {
        for(size_t col = 0; col < width; ++col)
        {
            const int v = inputRows[row][col];
            size_t x = col;
            size_t y = height - 1 - row; // map top-row index to image y
            (*data_ptr)[ img->PixelIndex({x,y}) ] = static_cast<T>(v);
        }
    }

    std::valarray<T> kdata = { T(1), T(1), T(1),
                               T(1), T(1), T(1),
                               T(1), T(1), T(1) };
    kernel<T> k(kdata, std::vector<size_t>{3,3});

    // Call opening wrapper
    auto out_ptr = opening<T>(std::dynamic_pointer_cast< FITScube >(img), k, 1);

    auto out_img = std::dynamic_pointer_cast< FITSimg<T> >(out_ptr);
    ASSERT_NE(out_img, nullptr);

    auto out_data = out_img->template GetData<T>();
    EXPECT_EQ(out_data->size(), static_cast<size_t>(maxSizeX*maxSizeY));

    // Validate output against expectedRows (provided top-to-bottom)
    for(size_t row = 0; row < height; ++row)
    {
        for(size_t col = 0; col < width; ++col)
        {
            size_t x = col;
            size_t y = height - 1 - row;
            size_t idx = out_img->PixelIndex({x,y});
            T expected = static_cast<T>( expectedRows[row][col] );
            EXPECT_EQ( (*out_data)[idx], expected ) << "Mismatch at (col=" << col << ", row(top)=" << row << ")";
        }
    }

    EXPECT_TRUE(out_img->HDU().Exists("DST_PROC"));
    EXPECT_TRUE(out_img->HDU().Exists("DST_VER"));
    EXPECT_TRUE(out_img->HDU().Exists("DST_SHA"));
    EXPECT_EQ( out_img->HDU().GetValueForKey("DST_PROC"), "OPENING");
    EXPECT_EQ( out_img->HDU().GetValueForKey("DST_VER"), DST::gGIT::version_short());
    EXPECT_EQ( out_img->HDU().GetValueForKey("DST_SHA"), DST::gGIT::sha1());

    out_img->Write("./build/KernelOpeningOne_"+demangle(typeid(T).name())+".fits",true);
}

TYPED_TEST(KernelFoldingTest, Closing)
{
    using T = TypeParam;

    // Replace previous random setup with explicit input/expected arrays (provided top-to-bottom).
    const size_t width  = 11;
    const size_t height = 11;

    auto img = std::make_shared< FITSimg<T> >(std::vector<size_t>{width,height});
    auto data_ptr = img->template GetData<T>();
    (*data_ptr) *= 0;

    // Input rows given top-to-bottom
    const int inputRows[11][11] = {
        {0,0,0,0,0,0,0,0,0,0,0},
        {0,1,1,1,1,0,0,1,1,1,0},
        {0,1,1,1,1,0,0,1,1,1,0},
        {0,1,1,1,1,1,1,1,1,1,0},
        {0,1,1,1,1,0,0,0,1,1,0},
        {0,1,1,1,1,0,0,0,1,1,0},
        {0,1,0,0,1,0,0,0,1,1,0},
        {0,1,0,0,1,1,1,1,1,1,0},
        {0,1,1,1,1,1,1,1,0,0,0},
        {0,1,1,1,1,1,1,1,0,0,0},
        {0,0,0,0,0,0,0,0,0,0,0}
    };

    // Expected rows given top-to-bottom
    const int expectedRows[11][11] = {
        {1,1,1,1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1,1,1,1},
        {1,1,1,1,1,0,0,0,1,1,1},
        {1,1,1,1,1,0,0,0,1,1,1},
        {1,1,1,1,1,0,0,0,1,1,1},
        {1,1,1,1,1,1,1,1,1,1,1},
        {1,1,1,1,1,1,1,1,0,0,0},
        {1,1,1,1,1,1,1,1,0,0,0},
        {1,1,1,1,1,1,1,1,0,0,0}
    };

    // Fill image: provided rows are top-to-bottom, PixelIndex expects y=0 bottom row -> map accordingly.
    for(size_t row = 0; row < height; ++row)
    {
        for(size_t col = 0; col < width; ++col)
        {
            const int v = inputRows[row][col];
            size_t x = col;
            size_t y = height - 1 - row; // map top-row index to image y
            (*data_ptr)[ img->PixelIndex({x,y}) ] = static_cast<T>(v);
        }
    }

    std::valarray<T> kdata = { T(1), T(1), T(1),
                               T(1), T(1), T(1),
                               T(1), T(1), T(1) };
    kernel<T> k(kdata, std::vector<size_t>{3,3});

    // Call closing wrapper
    auto out_ptr = closing<T>(std::dynamic_pointer_cast< FITScube >(img), k, 1);

    auto out_img = std::dynamic_pointer_cast< FITSimg<T> >(out_ptr);
    ASSERT_NE(out_img, nullptr);

    auto out_data = out_img->template GetData<T>();
    EXPECT_EQ(out_data->size(), static_cast<size_t>(width*height));

    // Validate output against expectedRows (top-to-bottom)
    for(size_t row = 0; row < height; ++row)
    {
        for(size_t col = 0; col < width; ++col)
        {
            size_t x = col;
            size_t y = height - 1 - row;
            size_t idx = out_img->PixelIndex({x,y});
            T expected = static_cast<T>( expectedRows[row][col] );
            EXPECT_EQ( (*out_data)[idx], expected ) << "Mismatch at (col=" << col << ", row(top)=" << row << ")";
        }
    }

    EXPECT_TRUE(out_img->HDU().Exists("DST_PROC"));
    EXPECT_TRUE(out_img->HDU().Exists("DST_VER"));
    EXPECT_TRUE(out_img->HDU().Exists("DST_SHA"));
    EXPECT_EQ( out_img->HDU().GetValueForKey("DST_PROC"), "CLOSING");
    EXPECT_EQ( out_img->HDU().GetValueForKey("DST_VER"), DST::gGIT::version_short());
    EXPECT_EQ( out_img->HDU().GetValueForKey("DST_SHA"), DST::gGIT::sha1());

    out_img->Write("./build/KernelClosingOne_"+demangle(typeid(T).name())+".fits",true);
}

TYPED_TEST(KernelFoldingTest, bTopHat)
{
    using T = TypeParam;

    // Use explicit 5x5 input (rows provided top-to-bottom).
    const size_t width  = 5;
    const size_t height = 5;

    auto img = std::make_shared< FITSimg<T> >(std::vector<size_t>{width,height});
    auto data_ptr = img->template GetData<T>();
    (*data_ptr) *= 0;

    const int inputRows[5][5] = {
        {7,6,6,6,7},
        {6,5,4,5,6},
        {6,4,0,4,6},
        {6,5,4,5,6},
        {7,6,6,6,7}
    };

    // Fill image (map top-to-bottom rows to y=height-1-row)
    for(size_t row = 0; row < height; ++row)
    for(size_t col = 0; col < width;  ++col)
    {
        size_t x = col;
        size_t y = height - 1 - row;
        (*data_ptr)[ img->PixelIndex({x,y}) ] = static_cast<T>(inputRows[row][col]);
    }

    // 3x3 all-ones kernel
    std::valarray<T> kdata = { T(1), T(1), T(1),
                               T(1), T(1), T(1),
                               T(1), T(1), T(1) };
    kernel<T> k(kdata, std::vector<size_t>{3,3});

    // bTopHat(out) = closing(in) - in
    auto out_ptr   = bTopHat<T>(std::dynamic_pointer_cast< FITScube >(img), k, 1);
    auto out_img   = std::dynamic_pointer_cast< FITSimg<T> >(out_ptr);
    ASSERT_NE(out_img, nullptr);
    auto out_data  = out_img->template GetData<T>();

    auto close_ptr = closing<T>(std::dynamic_pointer_cast< FITScube >(img), k, 1);
    auto close_img = std::dynamic_pointer_cast< FITSimg<T> >(close_ptr);
    ASSERT_NE(close_img, nullptr);
    auto close_data = close_img->template GetData<T>();

    ASSERT_EQ(out_data->size(), data_ptr->size());
    for(size_t i = 0; i < out_data->size(); ++i)
        EXPECT_EQ((*out_data)[i], static_cast<T>((*close_data)[i] - (*data_ptr)[i])) << "Mismatch at linear index " << i;

    EXPECT_TRUE(out_img->HDU().Exists("DST_PROC"));
    EXPECT_TRUE(out_img->HDU().Exists("DST_VER"));
    EXPECT_TRUE(out_img->HDU().Exists("DST_SHA"));
    EXPECT_EQ( out_img->HDU().GetValueForKey("DST_PROC"), "black Top-Hat");
    EXPECT_EQ( out_img->HDU().GetValueForKey("DST_VER"), DST::gGIT::version_short());
    EXPECT_EQ( out_img->HDU().GetValueForKey("DST_SHA"), DST::gGIT::sha1());

    out_img->Write("./build/KernelbTopHatOne_"+demangle(typeid(T).name())+".fits",true);
}

TYPED_TEST(KernelFoldingTest, wTopHat)
{
    using T = TypeParam;

    // Use explicit 5x5 input (rows provided top-to-bottom).
    const size_t width  = 5;
    const size_t height = 5;

    auto img = std::make_shared< FITSimg<T> >(std::vector<size_t>{width,height});
    auto data_ptr = img->template GetData<T>();
    (*data_ptr) *= 0;

    const int inputRows[5][5] = {
        {7,6,6,6,7},
        {6,5,4,5,6},
        {6,4,0,4,6},
        {6,5,4,5,6},
        {7,6,6,6,7}
    };

    for(size_t row = 0; row < height; ++row)
    for(size_t col = 0; col < width;  ++col)
    {
        size_t x = col;
        size_t y = height - 1 - row;
        (*data_ptr)[ img->PixelIndex({x,y}) ] = static_cast<T>(inputRows[row][col]);
    }

    std::valarray<T> kdata = { T(1), T(1), T(1),
                               T(1), T(1), T(1),
                               T(1), T(1), T(1) };
    kernel<T> k(kdata, std::vector<size_t>{3,3});

    // wTopHat(out) = in - opening(in)
    auto out_ptr   = wTopHat<T>(std::dynamic_pointer_cast< FITScube >(img), k, 1);
    auto out_img   = std::dynamic_pointer_cast< FITSimg<T> >(out_ptr);
    ASSERT_NE(out_img, nullptr);
    auto out_data  = out_img->template GetData<T>();

    auto open_ptr  = opening<T>(std::dynamic_pointer_cast< FITScube >(img), k, 1);
    auto open_img  = std::dynamic_pointer_cast< FITSimg<T> >(open_ptr);
    ASSERT_NE(open_img, nullptr);
    auto open_data = open_img->template GetData<T>();

    ASSERT_EQ(out_data->size(), data_ptr->size());
    for(size_t i = 0; i < out_data->size(); ++i)
        EXPECT_EQ((*out_data)[i], static_cast<T>((*data_ptr)[i] - (*open_data)[i])) << "Mismatch at linear index " << i;

    EXPECT_TRUE(out_img->HDU().Exists("DST_PROC"));
    EXPECT_TRUE(out_img->HDU().Exists("DST_VER"));
    EXPECT_TRUE(out_img->HDU().Exists("DST_SHA"));
    EXPECT_EQ( out_img->HDU().GetValueForKey("DST_PROC"), "White Top-Hat");
    EXPECT_EQ( out_img->HDU().GetValueForKey("DST_VER"), DST::gGIT::version_short());
    EXPECT_EQ( out_img->HDU().GetValueForKey("DST_SHA"), DST::gGIT::sha1());

    out_img->Write("./build/KernelwTopHatOne_"+demangle(typeid(T).name())+".fits",true);
}

// Simple, unambiguous sanity checks:

TYPED_TEST(KernelFoldingTest, wTopHat_SinglePixel)
{
    using T = TypeParam;
    const size_t W = 7, H = 7;

    auto img = std::make_shared< FITSimg<T> >(std::vector<size_t>{W,H});
    auto data = img->template GetData<T>();
    (*data) *= 0;

    // single bright pixel at center
    size_t cx = W/2, cy = H/2;
    (*data)[ img->PixelIndex({cx,cy}) ] = T(1);

    std::valarray<T> kdata = { T(1), T(1), T(1),
                               T(1), T(1), T(1),
                               T(1), T(1), T(1) };
    kernel<T> k(kdata, std::vector<size_t>{3,3});

    auto out_ptr = wTopHat<T>(std::dynamic_pointer_cast< FITScube >(img), k, 1);
    auto out_img = std::dynamic_pointer_cast< FITSimg<T> >(out_ptr);
    ASSERT_NE(out_img, nullptr);
    auto out = out_img->template GetData<T>();

    ASSERT_EQ(out->size(), data->size());
    for(size_t i = 0; i < out->size(); ++i)
        EXPECT_EQ((*out)[i], (*data)[i]); // opening removes isolated bright pixel -> wTopHat == original
}

// ============================================================
// ktype::wMEDIAN tests
//
// Motivation: MEDIAN on a sparse kernel (e.g. a diagonal line in a large
// bounding box) returns 0 because zero-weight cells are multiplied into the
// sorted buffer (ttmp = pixel * 0 = 0).  For an N-element kernel whose support
// is only M << N/2 pixels, the majority of sorted values are 0, so the median
// is 0.  wMEDIAN fixes this by only collecting pixels where kernel weight > 0.
// ============================================================

// 1. Dense all-ones kernel: wMEDIAN must equal MEDIAN at interior pixels
//    (no OOB zero-padding difference at borders).
TYPED_TEST(KernelFoldingTest, wMEDIAN_DenseKernelMatchesMedian)
{
    using T = TypeParam;
    const size_t W = 11, H = 11;

    auto img = std::make_shared<FITSimg<T>>(std::vector<size_t>{W, H});
    auto data = img->template GetData<T>();
    for(size_t y = 0; y < H; ++y)
    for(size_t x = 0; x < W; ++x)
        (*data)[img->PixelIndex({x,y})] = static_cast<T>((x*3 + y*7) % 13 + 1);

    std::valarray<T> kdata = {T(1),T(1),T(1),
                               T(1),T(1),T(1),
                               T(1),T(1),T(1)};
    kernel<T> km (kdata, {3,3}); km .setType(ktype::MEDIAN);
    kernel<T> kwm(kdata, {3,3}); kwm.setType(ktype::wMEDIAN);

    auto med_ptr  = folding<T>(std::dynamic_pointer_cast<FITScube>(img), km,  1);
    auto wmed_ptr = folding<T>(std::dynamic_pointer_cast<FITScube>(img), kwm, 1);
    auto med_img  = std::dynamic_pointer_cast<FITSimg<T>>(med_ptr);
    auto wmed_img = std::dynamic_pointer_cast<FITSimg<T>>(wmed_ptr);
    ASSERT_NE(med_img,  nullptr);
    ASSERT_NE(wmed_img, nullptr);

    auto med_data  = med_img->template GetData<T>();
    auto wmed_data = wmed_img->template GetData<T>();

    // Skip the 1-pixel border: OOB pixels are zero-padded by MEDIAN but
    // excluded by wMEDIAN, so border results intentionally differ.
    for(size_t y = 1; y < H-1; ++y)
    for(size_t x = 1; x < W-1; ++x)
    {
        const size_t idx = med_img->PixelIndex({x,y});
        EXPECT_EQ((*med_data)[idx], (*wmed_data)[idx])
            << "Mismatch at (" << x << "," << y << ")";
    }
}

// 2. Multi-thread and single-thread must produce identical results.
TYPED_TEST(KernelFoldingTest, wMEDIAN_MultiThread_MatchesSingleThread)
{
    using T = TypeParam;
    const size_t W = 51, H = 51;

    auto img = std::make_shared<FITSimg<T>>(std::vector<size_t>{W, H});
    auto data = img->template GetData<T>();
    for(size_t y = 0; y < H; ++y)
    for(size_t x = 0; x < W; ++x)
        (*data)[img->PixelIndex({x,y})] = static_cast<T>((x*3 + y*7 + 1) % 20 + 1);

    // 5×5 diagonal kernel: 5 non-zero cells out of 25 (flat index = x + y*5).
    std::valarray<T> kdata(T(0), 25);
    kdata[0]  = T(1);   // kIdx=(0,0) → offset (-2,-2)
    kdata[6]  = T(1);   // kIdx=(1,1) → offset (-1,-1)
    kdata[12] = T(1);   // kIdx=(2,2) → offset ( 0, 0) — center
    kdata[18] = T(1);   // kIdx=(3,3) → offset (+1,+1)
    kdata[24] = T(1);   // kIdx=(4,4) → offset (+2,+2)
    kernel<T> kwm(kdata, {5,5});
    kwm.setType(ktype::wMEDIAN);

    auto st_ptr = folding<T>(std::dynamic_pointer_cast<FITScube>(img), kwm, 1);
    auto mt_ptr = folding<T>(std::dynamic_pointer_cast<FITScube>(img), kwm, 4);

    auto st_img = std::dynamic_pointer_cast<FITSimg<T>>(st_ptr);
    auto mt_img = std::dynamic_pointer_cast<FITSimg<T>>(mt_ptr);
    ASSERT_NE(st_img, nullptr);
    ASSERT_NE(mt_img, nullptr);

    auto st_data = st_img->template GetData<T>();
    auto mt_data = mt_img->template GetData<T>();
    for(size_t i = 0; i < st_data->size(); ++i)
        EXPECT_EQ((*st_data)[i], (*mt_data)[i]) << "Mismatch at index " << i;
}

// ---- Non-typed wMEDIAN tests (double, specific values) ----

// 3. Core regression: MEDIAN returns 0 on a sparse kernel (>50 % zero cells),
//    wMEDIAN returns the true median of the support.
//
//    9×1 kernel: only positions 0, 4, 8 are non-zero (3 out of 9 = 67 % zeros).
//    Uniform image = 10.
//    MEDIAN: sorted buffer = {0,0,0,0,0,0,10,10,10} → median = 0.
//    wMEDIAN: buffer = {10,10,10} → median = 10.
TEST(MorphologywMEDIANTest, SparseKernel_MedianGivesZero_wMEDIANGivesCorrect)
{
    const size_t W = 20, H = 1;
    auto img = std::make_shared<FITSimg<double>>(std::vector<size_t>{W,H});
    (*img->template GetData<double>()) = 10.0;

    std::valarray<double> kdata(0.0, 9);
    kdata[0] = 1.0;
    kdata[4] = 1.0;
    kdata[8] = 1.0;

    kernel<double> km (kdata, {9,1}); km .setType(ktype::MEDIAN);
    kernel<double> kwm(kdata, {9,1}); kwm.setType(ktype::wMEDIAN);

    auto med_ptr  = folding<double>(std::dynamic_pointer_cast<FITScube>(img), km,  1);
    auto wmed_ptr = folding<double>(std::dynamic_pointer_cast<FITScube>(img), kwm, 1);
    auto med_img  = std::dynamic_pointer_cast<FITSimg<double>>(med_ptr);
    auto wmed_img = std::dynamic_pointer_cast<FITSimg<double>>(wmed_ptr);
    ASSERT_NE(med_img,  nullptr);
    ASSERT_NE(wmed_img, nullptr);

    // x=10: kCenter=4, support at x={6,10,14} — all in-bounds.
    const size_t cx = 10;
    EXPECT_DOUBLE_EQ((*med_img ->GetData<double>())[med_img ->PixelIndex(std::vector<size_t>({cx,0}))],  0.0);
    EXPECT_DOUBLE_EQ((*wmed_img->GetData<double>())[wmed_img->PixelIndex(std::vector<size_t>({cx,0}))], 10.0);
}

// 4. Diagonal line kernel (realistic NISP dispersion kernel scenario):
//    5×5 kernel, only the main diagonal set (5 ones out of 25 = 80 % zeros).
//    Uniform image = 42.
//    MEDIAN: 20 zeros dominate → median = 0.
//    wMEDIAN: 5 support pixels, all 42 → median = 42.
TEST(MorphologywMEDIANTest, DiagonalKernel_MedianZero_wMEDIANNonZero)
{
    const size_t KW = 5, KH = 5;
    std::valarray<double> kdata(0.0, KW * KH);
    kdata[0 + 0*KW] = 1.0;   // kIdx=(0,0) → offset (-2,-2)
    kdata[1 + 1*KW] = 1.0;   // kIdx=(1,1) → offset (-1,-1)
    kdata[2 + 2*KW] = 1.0;   // kIdx=(2,2) → center
    kdata[3 + 3*KW] = 1.0;   // kIdx=(3,3) → offset (+1,+1)
    kdata[4 + 4*KW] = 1.0;   // kIdx=(4,4) → offset (+2,+2)

    kernel<double> km (kdata, {KW,KH}); km .setType(ktype::MEDIAN);
    kernel<double> kwm(kdata, {KW,KH}); kwm.setType(ktype::wMEDIAN);

    const size_t W = 21, H = 21;
    auto img = std::make_shared<FITSimg<double>>(std::vector<size_t>{W,H});
    (*img->template GetData<double>()) = 42.0;

    auto med_ptr  = folding<double>(std::dynamic_pointer_cast<FITScube>(img), km,  1);
    auto wmed_ptr = folding<double>(std::dynamic_pointer_cast<FITScube>(img), kwm, 1);
    auto med_img  = std::dynamic_pointer_cast<FITSimg<double>>(med_ptr);
    auto wmed_img = std::dynamic_pointer_cast<FITSimg<double>>(wmed_ptr);
    ASSERT_NE(med_img,  nullptr);
    ASSERT_NE(wmed_img, nullptr);

    const size_t cx = W/2, cy = H/2;
    EXPECT_DOUBLE_EQ((*med_img ->GetData<double>())[med_img ->PixelIndex({cx,cy})],  0.0);
    EXPECT_DOUBLE_EQ((*wmed_img->GetData<double>())[wmed_img->PixelIndex({cx,cy})], 42.0);
}

// 5. Odd support count: wMEDIAN returns the true middle value.
//
//    1D image, values 1..9.
//    5×1 sparse kernel: weight > 0 at positions 0, 2, 4 (kCenter = 2).
//    For pixel x=4: support at x={2,4,6} with values {3,5,7} → median = 5.
TEST(MorphologywMEDIANTest, OddSupportCount_CorrectMedian)
{
    const size_t W = 9, H = 1;
    auto img = std::make_shared<FITSimg<double>>(std::vector<size_t>{W,H});
    auto data = img->template GetData<double>();
    for(size_t x = 0; x < W; ++x)
        (*data)[img->PixelIndex(std::vector<size_t>({x,0}))] = static_cast<double>(x + 1);

    std::valarray<double> kdata = {1.0, 0.0, 1.0, 0.0, 1.0};
    kernel<double> kwm(kdata, {5,1});
    kwm.setType(ktype::wMEDIAN);

    auto out_ptr = folding<double>(std::dynamic_pointer_cast<FITScube>(img), kwm, 1);
    auto out_img = std::dynamic_pointer_cast<FITSimg<double>>(out_ptr);
    ASSERT_NE(out_img, nullptr);

    EXPECT_DOUBLE_EQ((*out_img->GetData<double>())[out_img->PixelIndex({4,0})], 5.0);
}

// 6. Even support count: wMEDIAN returns the average of the two middle values.
//
//    1D image, values 1..9.
//    7×1 sparse kernel: weight > 0 at positions 0, 2, 4, 6 (kCenter = 3).
//    For pixel x=4: support at x={1,3,5,7} with values {2,4,6,8}.
//    Sorted: {2,4,6,8} → median = (4+6)/2 = 5.0.
TEST(MorphologywMEDIANTest, EvenSupportCount_MiddleAverage)
{
    const size_t W = 9, H = 1;
    auto img = std::make_shared<FITSimg<double>>(std::vector<size_t>{W,H});
    auto data = img->template GetData<double>();
    for(size_t x = 0; x < W; ++x)
        (*data)[img->PixelIndex(std::vector<size_t>({x,0}))] = static_cast<double>(x + 1);

    std::valarray<double> kdata = {1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0};
    kernel<double> kwm(kdata, {7,1});
    kwm.setType(ktype::wMEDIAN);

    auto out_ptr = folding<double>(std::dynamic_pointer_cast<FITScube>(img), kwm, 1);
    auto out_img = std::dynamic_pointer_cast<FITSimg<double>>(out_ptr);
    ASSERT_NE(out_img, nullptr);

    EXPECT_DOUBLE_EQ((*out_img->GetData<double>())[out_img->PixelIndex({4,0})], 5.0);
}

// 7. Masked pixels must be excluded from the wMEDIAN support.
//
//    Same 5×1 sparse kernel as test 5; mask x=6 (value 7).
//    Support reduces to {x=2(val=3), x=4(val=5)} → even count → median = (3+5)/2 = 4.0.
TEST(MorphologywMEDIANTest, MaskedPixelsExcluded)
{
    const size_t W = 9, H = 1;
    auto img = std::make_shared<FITSimg<double>>(std::vector<size_t>{W,H});
    auto data = img->template GetData<double>();
    for(size_t x = 0; x < W; ++x)
        (*data)[img->PixelIndex(std::vector<size_t>({x,0}))] = static_cast<double>(x + 1);

    img->MaskPixels({img->PixelIndex({6,0})});

    std::valarray<double> kdata = {1.0, 0.0, 1.0, 0.0, 1.0};
    kernel<double> kwm(kdata, {5,1});
    kwm.setType(ktype::wMEDIAN);

    auto out_ptr = folding<double>(std::dynamic_pointer_cast<FITScube>(img), kwm, 1);
    auto out_img = std::dynamic_pointer_cast<FITSimg<double>>(out_ptr);
    ASSERT_NE(out_img, nullptr);

    EXPECT_DOUBLE_EQ((*out_img->GetData<double>())[out_img->PixelIndex({4,0})], 4.0);
}

// 8. When all support pixels are masked, wMEDIAN must return 0.
TEST(MorphologywMEDIANTest, AllSupportMasked_ReturnsZero)
{
    const size_t W = 9, H = 1;
    auto img = std::make_shared<FITSimg<double>>(std::vector<size_t>{W,H});
    (*img->template GetData<double>()) = 10.0;

    img->MaskPixels({img->PixelIndex({2,0}),
                     img->PixelIndex({4,0}),
                     img->PixelIndex({6,0})});

    std::valarray<double> kdata = {1.0, 0.0, 1.0, 0.0, 1.0};
    kernel<double> kwm(kdata, {5,1});
    kwm.setType(ktype::wMEDIAN);

    auto out_ptr = folding<double>(std::dynamic_pointer_cast<FITScube>(img), kwm, 1);
    auto out_img = std::dynamic_pointer_cast<FITSimg<double>>(out_ptr);
    ASSERT_NE(out_img, nullptr);

    EXPECT_DOUBLE_EQ((*out_img->GetData<double>())[out_img->PixelIndex({4,0})], 0.0);
}

TYPED_TEST(KernelFoldingTest, bTopHat_SingleHole)
{
    using T = TypeParam;
    const size_t W = 7, H = 7;

    auto img = std::make_shared< FITSimg<T> >(std::vector<size_t>{W,H});
    auto data = img->template GetData<T>();
    (*data) = T(1);

    // single dark hole at center
    size_t cx = W/2, cy = H/2;
    (*data)[ img->PixelIndex({cx,cy}) ] = T(0);

    std::valarray<T> kdata = { T(1), T(1), T(1),
                               T(1), T(1), T(1),
                               T(1), T(1), T(1) };
    kernel<T> k(kdata, std::vector<size_t>{3,3});

    auto closing_ptr = closing<T>(std::dynamic_pointer_cast< FITScube >(img), k, 1);
    auto closing_img = std::dynamic_pointer_cast< FITSimg<T> >(closing_ptr);
    ASSERT_NE(closing_img, nullptr);
    auto closing_data = closing_img->template GetData<T>();
    for(size_t y = 0; y < H; ++y)
    for(size_t x = 0; x < W; ++x)
    {
        size_t idx = closing_img->PixelIndex({x,y});
        T expected = (x==0 || x+1 == W || y==0 || y+1==H) ? T(1) : T(1);
        EXPECT_EQ((*closing_data)[idx], expected) << "Mismatch at ("<<x<<","<<y<<")";
    }

    auto out_ptr = bTopHat<T>(std::dynamic_pointer_cast< FITScube >(img), k, 1);
    auto out_img = std::dynamic_pointer_cast< FITSimg<T> >(out_ptr);
    ASSERT_NE(out_img, nullptr);
    auto out = out_img->template GetData<T>();

    // Expect a single 1 at the hole, zeros elsewhere
    for(size_t y = 0; y < H; ++y)
    for(size_t x = 0; x < W; ++x)
    {
        size_t idx = out_img->PixelIndex({x,y});
        T expected = ((x==cx && y==cy)) ? T(1) :((x==0 || x+1 == W || y==0 || y+1==H)? T(0) : T(0));
        EXPECT_EQ((*out)[idx], expected) << "Mismatch at ("<<x<<","<<y<<")";
    }
}

