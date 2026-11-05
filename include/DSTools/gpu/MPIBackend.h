//
//  MPIBackend.h
//
//  MPI domain-decomposition backend for morphological operations.
//
//  Architecture:
//    • Splits the image along its slowest-varying axis (highest axis index)
//      into N_rank tiles of approximately equal size.
//    • Adds a halo of haloWidth = max(kCenter) rows on each tile boundary
//      so that every output pixel within the tile can be computed correctly
//      from local data alone.
//    • Each MPI rank runs its tile through its own per-node backend
//      (selected by Backend::create() — whichever GPU or CPU is fastest
//      on that node).
//    • After the kernel, halo rows in the output tile are discarded;
//      the valid interior rows are gathered to rank 0 via MPI_Gatherv.
//
//  Memory management:
//    "Device" memory for the MPI backend is ordinary host heap (same as
//    CPUBackend).  The per-node backend owns device-side allocations for
//    the actual computation.
//
//  Usage:
//    Only rank 0 calls allocate/upload/download as expected; other ranks
//    participate only during runKernel() via MPI collective calls.
//    runKernel() is a collective operation — all ranks must call it.
//
//  Compile guard: #ifdef USE_MPI.  MPI must be initialised by the caller
//  (MPI_Init) before constructing this object.
//
//  Licensed under PolyForm Noncommercial 1.0.0

#ifndef __DST_MPIBackend__
#define __DST_MPIBackend__

#include "Backend.h"

#include <mpi.h>
#include <memory>
#include <vector>
#include <stdexcept>
#include <cstring>
#include <cstdint>

namespace DST { namespace GPU {

class MPIBackend : public Backend
{
public:
    /// Construct an MPIBackend.
    /// \param comm  MPI communicator (default: MPI_COMM_WORLD).
    /// \param perNodeBackend  Override the per-node backend (null = auto-create
    ///                        via Backend::create() on each rank).
    explicit MPIBackend(MPI_Comm comm              = MPI_COMM_WORLD,
                        Backend* perNodeBackend    = nullptr);
    ~MPIBackend() override = default;

    bool available() const override;

    // ------------------------------------------------------------------
    // Memory management — host heap (rank 0 only for image data)
    // ------------------------------------------------------------------
    DeviceBuffer allocate(size_t bytes) override;
    void upload  (DeviceBuffer dst, const void* src, size_t bytes) override;
    void download(void* dst, DeviceBuffer src, size_t bytes) override;
    void free    (DeviceBuffer buf) override;

    // ------------------------------------------------------------------
    // runKernel — COLLECTIVE: all ranks must call simultaneously.
    // Performs scatter → per-node kernel → gather.
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
                   bool               normalize = false) override;

private:
    MPI_Comm comm_;
    int      rank_ = 0;
    int      nRanks_ = 1;

    // Per-node backend (owns device-side memory for the tile computation)
    std::unique_ptr<Backend> perNode_;

    // Utility
    static size_t pixelSize(PixelType pt);
    static MPI_Datatype mpiType(PixelType pt);
};

}} // namespace DST::GPU

#endif // __DST_MPIBackend__
