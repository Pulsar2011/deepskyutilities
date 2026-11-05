//
//  Backend.cpp
//
//  Backend::create() factory — selects the best available backend.
//  Compile-time guards control which backends are compiled in.
//  Licensed under PolyForm Noncommercial 1.0.0

#include "Backend.h"
#include "CPUBackend.h"

#ifdef USE_CUDA
#  include "CUDABackend.h"
#endif
#ifdef USE_HIP
#  include "HIPBackend.h"
#endif
#ifdef USE_VULKAN
#  include "VulkanBackend.h"
#endif
#ifdef USE_MPI
#  include "MPIBackend.h"
#endif

namespace DST { namespace GPU {

std::unique_ptr<Backend> Backend::create()
{
#ifdef USE_MPI
    {
        // MPI must already be initialised by the caller.
        // MPIBackend is only a useful accelerator when nRanks > 1.
        auto b = std::make_unique<MPIBackend>();
        if(b->available()) return b;
    }
#endif

#ifdef USE_CUDA
    {
        auto b = std::make_unique<CUDABackend>();
        if(b->available()) return b;
    }
#endif

#ifdef USE_HIP
    {
        auto b = std::make_unique<HIPBackend>();
        if(b->available()) return b;
    }
#endif

#ifdef USE_VULKAN
    {
        auto b = std::make_unique<VulkanBackend>();
        if(b->available()) return b;
    }
#endif

    // CPU is always available
    return std::make_unique<CPUBackend>();
}

}} // namespace DST::GPU
