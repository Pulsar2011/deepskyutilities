//
//  morpho.hip
//
//  HIP __global__ kernels for morphological operations (AMD ROCm).
//
//  This file is the HIP analogue of morpho.cu and follows the same design:
//    • #includes DSTmorph_ops.h directly — same per-pixel algorithm,
//      compiled for AMD GPU threads by hipcc.
//    • DST_HOST_DEVICE expands to __host__ __device__ under hipcc.
//    • All computation is done in float32 (pixel_xxx<float, PixT>).
//
//  The only differences from morpho.cu are the API prefix (hip* instead
//  of cuda*), the include of <hip/hip_runtime.h>, and the function name
//  hip_launch_kernel instead of cuda_launch_kernel.
//
//  Licensed under PolyForm Noncommercial 1.0.0

#include "../include/DSTools/DSTmorph_ops.h"
#include "HIPBackend.h"

#include <hip/hip_runtime.h>
#include <cstdint>
#include <stdexcept>
#include <string>

// ------------------------------------------------------------
// Tuning constants
// ------------------------------------------------------------
static constexpr unsigned HIP_BLOCK_SIZE = 256;

// ------------------------------------------------------------
// __global__ kernel — templated on the image pixel type PixT.
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
    const size_t i = static_cast<size_t>(hipBlockIdx_x) * hipBlockDim_x + hipThreadIdx_x;
    if(i >= nPixels) return;

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
        (nPixels + HIP_BLOCK_SIZE - 1) / HIP_BLOCK_SIZE);

    hipLaunchKernelGGL(morph_dispatch<PixT>,
        dim3(nBlocks), dim3(HIP_BLOCK_SIZE), 0, 0,
        static_cast<const PixT*>(in),
        mask,
        static_cast<PixT*>(out),
        kw, kSizes, kCenter,
        kDim, kSize,
        imgSizes, imgDim, nPixels,
        op, normalize);

    hipError_t err = hipDeviceSynchronize();
    if(err != hipSuccess)
        throw std::runtime_error(
            std::string("morpho.hip: HIP kernel failed: ")
            + hipGetErrorString(err));
}

// ------------------------------------------------------------
// hip_launch_kernel — public entry point called by HIPBackend.
// ------------------------------------------------------------
namespace DST { namespace GPU {

void hip_launch_kernel(Morphology::ktype  op,
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
    if(op == Morphology::ktype::MEDIAN &&
       kSize > Morphology::ops::MAX_KERNEL_SIZE)
    {
        throw std::runtime_error(
            "hip_launch_kernel: MEDIAN kernel exceeds MAX_KERNEL_SIZE; "
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
            throw std::invalid_argument("hip_launch_kernel: unsupported PixelType");
    }
}

}} // namespace DST::GPU
