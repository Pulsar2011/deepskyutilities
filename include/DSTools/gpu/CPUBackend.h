//
//  CPUBackend.h
//
//  CPU implementation of the Backend interface.
//  "Device" memory is ordinary heap memory; upload/download are memcpy.
//  runKernel() calls the portable per-pixel functions from DSTmorph_ops.h
//  in a plain for-loop, giving a CPU path through the full Backend abstraction.
//
//  This backend is always available and serves as:
//    1. The permanent fallback for platforms without a GPU backend.
//    2. The test vehicle for GPUContext before real GPU backends exist.
//    3. The mandatory CPU path for float64 / int64 / uint64 / size_t images.
//
//  Licensed under CC BY-NC 4.0

#ifndef __DST_CPUBackend__
#define __DST_CPUBackend__

#include "Backend.h"
#include <DSTools/DSTmorph_ops.h>

#include <cstring>    // memcpy
#include <stdexcept>

namespace DST { namespace GPU {

class CPUBackend : public Backend
{
public:
    bool available() const override { return true; }

    // ------------------------------------------------------------------
    // Memory management — host heap acts as "device"
    // ------------------------------------------------------------------
    DeviceBuffer allocate(size_t bytes) override
    {
        DeviceBuffer b;
        b.ptr   = ::operator new(bytes);
        b.bytes = bytes;
        return b;
    }

    void upload(DeviceBuffer dst, const void* src, size_t bytes) override
    {
        std::memcpy(dst.ptr, src, bytes);
    }

    void download(void* dst, DeviceBuffer src, size_t bytes) override
    {
        std::memcpy(dst, src.ptr, bytes);
    }

    void free(DeviceBuffer buf) override
    {
        ::operator delete(buf.ptr);
    }

    // ------------------------------------------------------------------
    // runKernel — dispatches to per-pixel ops functions
    // ------------------------------------------------------------------
    void runKernel(Morphology::ktype  op,
                   DeviceBuffer       in,
                   DeviceBuffer       mask,
                   DeviceBuffer       out,
                   DeviceBuffer       kernelW,
                   DeviceBuffer       kSizes,
                   DeviceBuffer       kCenter,
                   size_t             kDim,
                   size_t             kSize,
                   DeviceBuffer       imgSizes,
                   size_t             imgDim,
                   size_t             nPixels,
                   PixelType          pixelType,
                   bool               normalize) override
    {
        switch(pixelType)
        {
            case PixelType::UInt8:
                runTyped<uint8_t>(op, in, mask, out, kernelW, kSizes, kCenter,
                                  kDim, kSize, imgSizes, imgDim, nPixels, normalize);
                break;
            case PixelType::Int8:
                runTyped<int8_t>(op, in, mask, out, kernelW, kSizes, kCenter,
                                 kDim, kSize, imgSizes, imgDim, nPixels, normalize);
                break;
            case PixelType::UInt16:
                runTyped<uint16_t>(op, in, mask, out, kernelW, kSizes, kCenter,
                                   kDim, kSize, imgSizes, imgDim, nPixels, normalize);
                break;
            case PixelType::Int16:
                runTyped<int16_t>(op, in, mask, out, kernelW, kSizes, kCenter,
                                  kDim, kSize, imgSizes, imgDim, nPixels, normalize);
                break;
            case PixelType::UInt32:
                runTyped<uint32_t>(op, in, mask, out, kernelW, kSizes, kCenter,
                                   kDim, kSize, imgSizes, imgDim, nPixels, normalize);
                break;
            case PixelType::Int32:
                runTyped<int32_t>(op, in, mask, out, kernelW, kSizes, kCenter,
                                  kDim, kSize, imgSizes, imgDim, nPixels, normalize);
                break;
            case PixelType::Float32:
                runTyped<float>(op, in, mask, out, kernelW, kSizes, kCenter,
                                kDim, kSize, imgSizes, imgDim, nPixels, normalize);
                break;
            default:
                throw std::invalid_argument("CPUBackend::runKernel: unsupported PixelType");
        }
    }

private:
    // The kernel weights are always stored as float32 in DeviceBuffer (GPU convention).
    // For integer pixel types, U = pixel type (the per-pixel ops promote internally).
    // Here we use the pixel type as both K and U for simplicity on the CPU backend
    // (the ops functions cast kernel weights via static_cast<U> from float internally
    // by reading kw as float* and interpreting as U in the template).
    //
    // Actually, per-pixel ops take `const U* kernelWeights`. On the CPU backend we
    // re-cast the kernel weights from float back to U. For morphological footprint ops
    // (dilation, erosion), only the sign of the weight matters (>0 test). For fold/
    // correlate, the weight value is used in arithmetic. We therefore store a U-typed
    // weight copy locally for correctness.
    template<typename T>
    void runTyped(Morphology::ktype op,
                  DeviceBuffer in, DeviceBuffer mask, DeviceBuffer out,
                  DeviceBuffer kernelW, DeviceBuffer kSizes, DeviceBuffer kCenter,
                  size_t kDim, size_t kSize,
                  DeviceBuffer imgSizes, size_t imgDim,
                  size_t nPixels, bool normalize)
    {
        const T*         in_ptr   = static_cast<const T*>(in.ptr);
        T*               out_ptr  = static_cast<T*>(out.ptr);
        const bool*      mask_ptr = (mask.ptr && mask.bytes > 0)
                                        ? static_cast<const bool*>(mask.ptr)
                                        : nullptr;
        const float*     kw_f    = static_cast<const float*>(kernelW.ptr);
        const size_t*    ks_ptr  = static_cast<const size_t*>(kSizes.ptr);
        const ptrdiff_t* kc_ptr  = static_cast<const ptrdiff_t*>(kCenter.ptr);
        const size_t*    is_ptr  = static_cast<const size_t*>(imgSizes.ptr);

        // Build T-typed kernel weights (promoted from float storage)
        // using a local heap buffer to avoid VLA and stack size issues.
        std::vector<T> kw_typed(kSize);
        for(size_t i = 0; i < kSize; ++i)
            kw_typed[i] = static_cast<T>(kw_f[i]);

        const T* kw_ptr = kw_typed.data();

        if(op == Morphology::ktype::DILATE)
        {
            for(size_t i = 0; i < nPixels; ++i)
                out_ptr[i] = Morphology::ops::pixel_dilation<T,T>(
                    in_ptr, mask_ptr, kw_ptr, ks_ptr, kDim, kSize, kc_ptr, is_ptr, imgDim, i);
        }
        else if(op == Morphology::ktype::ERODE)
        {
            for(size_t i = 0; i < nPixels; ++i)
                out_ptr[i] = Morphology::ops::pixel_erosion<T,T>(
                    in_ptr, mask_ptr, kw_ptr, ks_ptr, kDim, kSize, kc_ptr, is_ptr, imgDim, i);
        }
        else if(op == Morphology::ktype::CORRELATE || op == Morphology::ktype::nCORRELATE)
        {
            const bool norm = (op == Morphology::ktype::nCORRELATE) || normalize;
            for(size_t i = 0; i < nPixels; ++i)
                out_ptr[i] = Morphology::ops::pixel_correlate<T,T>(
                    in_ptr, mask_ptr, kw_ptr, ks_ptr, kDim, kSize, kc_ptr, is_ptr, imgDim, i, norm);
        }
        else
        {
            const Morphology::ops::fold_op fop = static_cast<Morphology::ops::fold_op>(op);
            for(size_t i = 0; i < nPixels; ++i)
                out_ptr[i] = Morphology::ops::pixel_fold<T,T>(
                    in_ptr, mask_ptr, kw_ptr, ks_ptr, kDim, kSize, kc_ptr, is_ptr, imgDim, i, fop);
        }
    }
};

}} // namespace DST::GPU

#endif // __DST_CPUBackend__
