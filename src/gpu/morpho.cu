//
//  morpho.cu
//
//  CUDA __global__ kernels for morphological operations.
//
//  KEY DESIGN: This file #includes DSTmorph_ops.h directly.
//  The same per-pixel functions (pixel_dilation, pixel_erosion,
//  pixel_correlate, pixel_fold) used by the CPU path are compiled
//  into GPU kernels here by nvcc — zero algorithmic duplication.
//  DST_HOST_DEVICE expands to __host__ __device__ under nvcc, making
//  every helper in DSTmorph_ops.h available in both host and device code.
//
//  Computation precision:
//    All intermediate arithmetic is done in float32, using the
//    pixel_xxx<float, PixT> template instantiation where:
//      • float = computation type (kernel weights are always float32 on GPU)
//      • PixT  = input image pixel type
//    The result is cast back to PixT before writing to the output buffer.
//    Precision is exact for all pixel types with ≤ 24-bit value range
//    (int8/uint8/int16/uint16).  For int32/uint32 with values > 2^24, the
//    float32 intermediate may differ by ≤1 ULP in min/max comparisons;
//    use the CPU path if exact int32 comparisons are required.
//
//  MEDIAN handling:
//    pixel_fold<float, PixT> for MEDIAN allocates a float buf[MAX_KERNEL_SIZE]
//    stack array (≤2 KB per thread).  Kernels with kSize > MAX_KERNEL_SIZE
//    are rejected by cuda_launch_kernel() — the caller (GPUContext::run) must
//    fall back to CPU for oversized MEDIAN kernels.
//
//  Licensed under PolyForm Noncommercial 1.0.0

#include "../include/DSTools/DSTmorph_ops.h"
#include "CUDABackend.h"

#include <cuda_runtime.h>
#include <cstdint>
#include <stdexcept>
#include <string>

// ------------------------------------------------------------
// Tuning constants
// ------------------------------------------------------------
static constexpr unsigned CUDA_BLOCK_SIZE = 256;

// ------------------------------------------------------------
// __global__ kernel — templated on the image pixel type PixT.
// One thread per output pixel.  The operation type (op) selects
// the appropriate per-pixel function from DSTmorph_ops.h.
// ------------------------------------------------------------
template<typename PixT>
__global__
void morph_dispatch(const PixT*      in,
                    const bool*      mask,      // may be nullptr
                    PixT*            out,
                    const float*     kw,        // kernel weights (float32)
                    const size_t*    kSizes,    // kernel dimension sizes
                    const ptrdiff_t* kCenter,   // kernel center offsets
                    size_t           kDim,
                    size_t           kSize,
                    const size_t*    imgSizes,  // image dimension sizes
                    size_t           imgDim,
                    size_t           nPixels,
                    int              op,        // DST::Morphology::ktype as int
                    bool             normalize)
{
    const size_t i = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if(i >= nPixels) return;

    // All per-pixel functions are instantiated with:
    //   U = float  (computation and kernel-weight type)
    //   K = PixT   (input image pixel type)
    // The result (float) is cast back to PixT before storage.

    using namespace DST::Morphology;

    float result = 0.0f;

    if(op == static_cast<int>(ktype::DILATE))
    {
        result = ops::pixel_dilation<float, PixT>(
            in, mask, kw, kSizes, kDim, kSize, kCenter, imgSizes, imgDim, i);
    }
    else if(op == static_cast<int>(ktype::ERODE))
    {
        result = ops::pixel_erosion<float, PixT>(
            in, mask, kw, kSizes, kDim, kSize, kCenter, imgSizes, imgDim, i);
    }
    else if(op == static_cast<int>(ktype::CORRELATE) ||
            op == static_cast<int>(ktype::nCORRELATE))
    {
        const bool norm = (op == static_cast<int>(ktype::nCORRELATE)) || normalize;
        result = ops::pixel_correlate<float, PixT>(
            in, mask, kw, kSizes, kDim, kSize, kCenter, imgSizes, imgDim, i, norm);
    }
    else
    {
        // Fold operations: MIN, MAX, MEAN, MEDIAN, SUM, RMS, RMSE
        // ktype enum values 0–6 map directly to ops::fold_op (enforced by static_cast
        // compatibility designed in DSTmorph_ops.h).
        result = ops::pixel_fold<float, PixT>(
            in, mask, kw, kSizes, kDim, kSize, kCenter, imgSizes, imgDim, i,
            static_cast<ops::fold_op>(op));
    }

    out[i] = static_cast<PixT>(result);
}

// ------------------------------------------------------------
// Helper: launch morph_dispatch<PixT> and synchronise.
// ------------------------------------------------------------
template<typename PixT>
static void launch(const void* in, const bool* mask, void* out,
                   const float* kw,
                   const size_t* kSizes, const ptrdiff_t* kCenter,
                   size_t kDim, size_t kSize,
                   const size_t* imgSizes, size_t imgDim, size_t nPixels,
                   int op, bool normalize)
{
    const unsigned nBlocks = static_cast<unsigned>(
        (nPixels + CUDA_BLOCK_SIZE - 1) / CUDA_BLOCK_SIZE);

    morph_dispatch<PixT><<<nBlocks, CUDA_BLOCK_SIZE>>>(
        static_cast<const PixT*>(in),
        mask,
        static_cast<PixT*>(out),
        kw, kSizes, kCenter,
        kDim, kSize,
        imgSizes, imgDim, nPixels,
        op, normalize);

    cudaError_t err = cudaDeviceSynchronize();
    if(err != cudaSuccess)
        throw std::runtime_error(
            std::string("morpho.cu: CUDA kernel failed: ")
            + cudaGetErrorString(err));
}

// ------------------------------------------------------------
// cuda_launch_kernel — public entry point called by CUDABackend.
// Dispatches to the correct PixT template instantiation based on
// the pixelType tag.
// ------------------------------------------------------------
namespace DST { namespace GPU {

void cuda_launch_kernel(Morphology::ktype  op,
                        const void*        in,
                        const bool*        mask,
                        void*              out,
                        const float*       kw,
                        const size_t*      kSizes,
                        const ptrdiff_t*   kCenter,
                        size_t             kDim,
                        size_t             kSize,
                        const size_t*      imgSizes,
                        size_t             imgDim,
                        size_t             nPixels,
                        PixelType          pixelType,
                        bool               normalize)
{
    // MEDIAN with oversized kernel cannot use the stack buffer — reject early.
    if(op == Morphology::ktype::MEDIAN &&
       kSize > Morphology::ops::MAX_KERNEL_SIZE)
    {
        throw std::runtime_error(
            "cuda_launch_kernel: MEDIAN kernel exceeds MAX_KERNEL_SIZE; "
            "use the CPU path for this kernel size.");
    }

    const int op_i = static_cast<int>(op);

    switch(pixelType)
    {
        case PixelType::UInt8:
            launch<uint8_t> (in, mask, out, kw, kSizes, kCenter, kDim, kSize, imgSizes, imgDim, nPixels, op_i, normalize);
            break;
        case PixelType::Int8:
            launch<int8_t>  (in, mask, out, kw, kSizes, kCenter, kDim, kSize, imgSizes, imgDim, nPixels, op_i, normalize);
            break;
        case PixelType::UInt16:
            launch<uint16_t>(in, mask, out, kw, kSizes, kCenter, kDim, kSize, imgSizes, imgDim, nPixels, op_i, normalize);
            break;
        case PixelType::Int16:
            launch<int16_t> (in, mask, out, kw, kSizes, kCenter, kDim, kSize, imgSizes, imgDim, nPixels, op_i, normalize);
            break;
        case PixelType::UInt32:
            launch<uint32_t>(in, mask, out, kw, kSizes, kCenter, kDim, kSize, imgSizes, imgDim, nPixels, op_i, normalize);
            break;
        case PixelType::Int32:
            launch<int32_t> (in, mask, out, kw, kSizes, kCenter, kDim, kSize, imgSizes, imgDim, nPixels, op_i, normalize);
            break;
        case PixelType::Float32:
            launch<float>   (in, mask, out, kw, kSizes, kCenter, kDim, kSize, imgSizes, imgDim, nPixels, op_i, normalize);
            break;
        default:
            throw std::invalid_argument("cuda_launch_kernel: unsupported PixelType");
    }
}

}} // namespace DST::GPU
