//
//  HIPBackend.h
//
//  HIP (AMD ROCm) implementation of the Backend interface.
//  Include this header only from HIP translation units (.hip files),
//  or from host code guarded with #ifdef USE_HIP.
//
//  Memory management delegates directly to the HIP runtime:
//    allocate → hipMalloc
//    upload   → hipMemcpy (Host → Device)
//    download → hipMemcpy (Device → Host)
//    free     → hipFree
//
//  runKernel() forwards to hip_launch_kernel() defined in morpho.hip,
//  which contains the actual __global__ kernel dispatch.
//
//  The design mirrors CUDABackend.h exactly; only the API prefix changes
//  from cuda* to hip*.
//
//  Licensed under PolyForm Noncommercial 1.0.0

#ifndef __DST_HIPBackend__
#define __DST_HIPBackend__

#include "Backend.h"

#include <hip/hip_runtime.h>
#include <stdexcept>
#include <string>

namespace DST { namespace GPU {

// ------------------------------------------------------------
// hip_launch_kernel — declared here, defined in morpho.hip.
// Selects the __global__ kernel specialisation for pixelType,
// launches it, and synchronises before returning.
// ------------------------------------------------------------
void hip_launch_kernel(Morphology::ktype  op,
                       const void*        in,        // device ptr, typed by pixelType
                       const bool*        mask,      // device ptr, may be nullptr
                       void*              out,       // device ptr, typed by pixelType
                       const float*       kw,        // device ptr, float32 weights
                       const size_t*      kSizes,    // device ptr, size_t[kDim]
                       const ptrdiff_t*   kCenter,   // device ptr, ptrdiff_t[kDim]
                       size_t             kDim,
                       size_t             kSize,
                       const size_t*      imgSizes,  // device ptr, size_t[imgDim]
                       size_t             imgDim,
                       size_t             nPixels,
                       PixelType          pixelType,
                       bool               normalize);

// ------------------------------------------------------------
// HIPBackend
// ------------------------------------------------------------
class HIPBackend : public Backend
{
public:
    bool available() const override
    {
        int nDevices = 0;
        hipError_t err = hipGetDeviceCount(&nDevices);
        return (err == hipSuccess && nDevices > 0);
    }

    // ------------------------------------------------------------------
    // Memory management
    // ------------------------------------------------------------------
    DeviceBuffer allocate(size_t bytes) override
    {
        DeviceBuffer b;
        b.bytes = bytes;
        hipError_t err = hipMalloc(&b.ptr, bytes);
        if(err != hipSuccess)
            throw std::runtime_error(
                std::string("HIPBackend::allocate: hipMalloc failed: ")
                + hipGetErrorString(err));
        return b;
    }

    void upload(DeviceBuffer dst, const void* src, size_t bytes) override
    {
        hipError_t err = hipMemcpy(dst.ptr, src, bytes, hipMemcpyHostToDevice);
        if(err != hipSuccess)
            throw std::runtime_error(
                std::string("HIPBackend::upload: hipMemcpy H→D failed: ")
                + hipGetErrorString(err));
    }

    void download(void* dst, DeviceBuffer src, size_t bytes) override
    {
        hipError_t err = hipMemcpy(dst, src.ptr, bytes, hipMemcpyDeviceToHost);
        if(err != hipSuccess)
            throw std::runtime_error(
                std::string("HIPBackend::download: hipMemcpy D→H failed: ")
                + hipGetErrorString(err));
    }

    void free(DeviceBuffer buf) override
    {
        hipFree(buf.ptr);  // hipFree ignores nullptr — safe
    }

    // ------------------------------------------------------------------
    // Kernel execution — delegates to hip_launch_kernel (morpho.hip)
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
                   bool               normalize = false) override
    {
        hip_launch_kernel(
            op,
            in.ptr,
            static_cast<const bool*>(mask.ptr),
            out.ptr,
            static_cast<const float*>(kernelW.ptr),
            static_cast<const size_t*>(kSizes.ptr),
            static_cast<const ptrdiff_t*>(kCenter.ptr),
            kDim, kSize,
            static_cast<const size_t*>(imgSizes.ptr),
            imgDim, nPixels,
            pixelType, normalize);
    }
};

}} // namespace DST::GPU

#endif // __DST_HIPBackend__
