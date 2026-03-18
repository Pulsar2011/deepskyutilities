//
//  CUDABackend.h
//
//  CUDA implementation of the Backend interface.
//  Include this header only from CUDA translation units (.cu files),
//  or from host code guarded with #ifdef USE_CUDA.
//
//  Memory management delegates directly to the CUDA runtime:
//    allocate → cudaMalloc
//    upload   → cudaMemcpy (Host → Device)
//    download → cudaMemcpy (Device → Host)
//    free     → cudaFree
//
//  runKernel() forwards to cuda_launch_kernel() defined in morpho.cu,
//  which contains the actual __global__ kernel dispatch.  This keeps
//  the __global__ kernels in a separate translation unit compiled by nvcc,
//  while this header file can be #included by the Backend::create() factory.
//
//  Licensed under CC BY-NC 4.0

#ifndef __DST_CUDABackend__
#define __DST_CUDABackend__

#include "Backend.h"

#include <cuda_runtime.h>
#include <stdexcept>
#include <string>

namespace DST { namespace GPU {

// ------------------------------------------------------------
// cuda_launch_kernel — declared here, defined in morpho.cu.
// Selects the __global__ kernel specialisation for pixelType,
// launches it, and synchronises before returning.
// ------------------------------------------------------------
void cuda_launch_kernel(Morphology::ktype  op,
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
// CUDABackend
// ------------------------------------------------------------
class CUDABackend : public Backend
{
public:
    bool available() const override
    {
        int nDevices = 0;
        cudaError_t err = cudaGetDeviceCount(&nDevices);
        return (err == cudaSuccess && nDevices > 0);
    }

    // ------------------------------------------------------------------
    // Memory management
    // ------------------------------------------------------------------
    DeviceBuffer allocate(size_t bytes) override
    {
        DeviceBuffer b;
        b.bytes = bytes;
        cudaError_t err = cudaMalloc(&b.ptr, bytes);
        if(err != cudaSuccess)
            throw std::runtime_error(
                std::string("CUDABackend::allocate: cudaMalloc failed: ")
                + cudaGetErrorString(err));
        return b;
    }

    void upload(DeviceBuffer dst, const void* src, size_t bytes) override
    {
        cudaError_t err = cudaMemcpy(dst.ptr, src, bytes, cudaMemcpyHostToDevice);
        if(err != cudaSuccess)
            throw std::runtime_error(
                std::string("CUDABackend::upload: cudaMemcpy H→D failed: ")
                + cudaGetErrorString(err));
    }

    void download(void* dst, DeviceBuffer src, size_t bytes) override
    {
        cudaError_t err = cudaMemcpy(dst, src.ptr, bytes, cudaMemcpyDeviceToHost);
        if(err != cudaSuccess)
            throw std::runtime_error(
                std::string("CUDABackend::download: cudaMemcpy D→H failed: ")
                + cudaGetErrorString(err));
    }

    void free(DeviceBuffer buf) override
    {
        cudaFree(buf.ptr);  // cudaFree ignores nullptr — safe
    }

    // ------------------------------------------------------------------
    // Kernel execution — delegates to cuda_launch_kernel (morpho.cu)
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
        cuda_launch_kernel(
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

#endif // __DST_CUDABackend__
