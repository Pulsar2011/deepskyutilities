//
//  MPIBackend.cpp
//
//  MPI domain-decomposition backend implementation.
//  Licensed under PolyForm Noncommercial 1.0.0

#include "MPIBackend.h"
#include "CPUBackend.h"

#include <algorithm>
#include <cstring>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

namespace DST { namespace GPU {

// ============================================================
//  Pixel type helpers
// ============================================================
size_t MPIBackend::pixelSize(PixelType pt)
{
    switch(pt)
    {
        case PixelType::UInt8:  case PixelType::Int8:   return 1;
        case PixelType::UInt16: case PixelType::Int16:  return 2;
        case PixelType::UInt32: case PixelType::Int32:
        case PixelType::Float32:                        return 4;
        default:
            throw std::invalid_argument("MPIBackend::pixelSize: unknown PixelType");
    }
}

MPI_Datatype MPIBackend::mpiType(PixelType pt)
{
    switch(pt)
    {
        case PixelType::UInt8:  return MPI_UINT8_T;
        case PixelType::Int8:   return MPI_INT8_T;
        case PixelType::UInt16: return MPI_UINT16_T;
        case PixelType::Int16:  return MPI_INT16_T;
        case PixelType::UInt32: return MPI_UINT32_T;
        case PixelType::Int32:  return MPI_INT32_T;
        case PixelType::Float32:return MPI_FLOAT;
        default:
            throw std::invalid_argument("MPIBackend::mpiType: unknown PixelType");
    }
}

// ============================================================
//  Constructor
// ============================================================
MPIBackend::MPIBackend(MPI_Comm comm, Backend* perNodeBackend)
    : comm_(comm)
{
    MPI_Comm_rank(comm_, &rank_);
    MPI_Comm_size(comm_, &nRanks_);

    if(perNodeBackend)
        perNode_.reset(perNodeBackend);   // takes ownership
    else
        perNode_ = Backend::create();     // auto-select on each node
}

// ============================================================
//  available()
// ============================================================
bool MPIBackend::available() const
{
    return nRanks_ > 1 && perNode_ && perNode_->available();
}

// ============================================================
//  Memory management — host heap, rank 0 only for image buffers.
//  All ranks call allocate() because the interface is collective in
//  runKernel(); non-zero ranks may receive nullptr safely.
// ============================================================
DeviceBuffer MPIBackend::allocate(size_t bytes)
{
    DeviceBuffer b;
    b.bytes = bytes;
    b.ptr   = (bytes > 0) ? ::operator new(bytes) : nullptr;
    return b;
}

void MPIBackend::upload(DeviceBuffer dst, const void* src, size_t bytes)
{
    if(dst.ptr && src && bytes > 0)
        std::memcpy(dst.ptr, src, bytes);
}

void MPIBackend::download(void* dst, DeviceBuffer src, size_t bytes)
{
    if(dst && src.ptr && bytes > 0)
        std::memcpy(dst, src.ptr, bytes);
}

void MPIBackend::free(DeviceBuffer buf)
{
    ::operator delete(buf.ptr);
}

// ============================================================
//  runKernel — collective across all ranks
//
//  Domain decomposition strategy:
//    • The image is N-dimensional; decomposition is along the LAST
//      (slowest-varying) axis (axis imgDim-1), i.e. the "rows" axis
//      in a 2-D image, or "planes" in a 3-D image.
//    • haloWidth = max(kCenter[kDim-1])  (the half-width of the kernel
//      along the decomposition axis; other axes are always local).
//    • Each rank receives a contiguous block of slices of the input +
//      halo slices from neighbouring ranks.
//    • After local computation, the halo rows of the OUTPUT are discarded;
//      only the valid interior rows are gathered to rank 0.
//
//  Kernel geometry assumption:
//    The kernel's last dimension covers the decomposition axis.
//    If kDim < imgDim, axes 0..(imgDim-kDim-1) are not covered by the
//    kernel, so no halo is needed on those axes — only the last kernel
//    axis matters for the halo.
// ============================================================
void MPIBackend::runKernel(Morphology::ktype  op,
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
                            bool               normalize)
{
    if(nRanks_ == 1)
    {
        // Single-rank fallback: run directly via per-node backend
        DeviceBuffer d_in     = perNode_->allocate(in.bytes);
        DeviceBuffer d_mask   = (mask.ptr && mask.bytes > 0)
                                    ? perNode_->allocate(mask.bytes)
                                    : DeviceBuffer{nullptr, 0};
        DeviceBuffer d_out    = perNode_->allocate(out.bytes);
        DeviceBuffer d_kw     = perNode_->allocate(kernelW.bytes);
        DeviceBuffer d_ks     = perNode_->allocate(kSizes.bytes);
        DeviceBuffer d_kc     = perNode_->allocate(kCenter.bytes);
        DeviceBuffer d_is     = perNode_->allocate(imgSizes.bytes);

        perNode_->upload(d_in, in.ptr, in.bytes);
        if(d_mask.ptr) perNode_->upload(d_mask, mask.ptr, mask.bytes);
        perNode_->upload(d_kw, kernelW.ptr, kernelW.bytes);
        perNode_->upload(d_ks, kSizes.ptr, kSizes.bytes);
        perNode_->upload(d_kc, kCenter.ptr, kCenter.bytes);
        perNode_->upload(d_is, imgSizes.ptr, imgSizes.bytes);

        perNode_->runKernel(op, d_in, d_mask, d_out, d_kw, d_ks, d_kc,
                             kDim, kSize, d_is, imgDim, nPixels, pixelType, normalize);

        perNode_->download(out.ptr, d_out, out.bytes);

        perNode_->free(d_in);
        if(d_mask.ptr) perNode_->free(d_mask);
        perNode_->free(d_out);
        perNode_->free(d_kw);
        perNode_->free(d_ks);
        perNode_->free(d_kc);
        perNode_->free(d_is);
        return;
    }

    // ------------------------------------------------------------------
    // Read geometry from host-side DeviceBuffers
    // ------------------------------------------------------------------
    const size_t* imgSzHost = static_cast<const size_t*>(imgSizes.ptr);
    const ptrdiff_t* kcHost = static_cast<const ptrdiff_t*>(kCenter.ptr);

    // Halo width = half-size of the kernel along the decomposition axis.
    // kCenter[kDim-1] is the center offset of the last kernel axis.
    const size_t decompAxis  = imgDim - 1;          // slowest-varying image axis
    const size_t kDecompAxis = kDim - 1;             // corresponding kernel axis
    const size_t haloWidth   = (kDim > 0)
        ? static_cast<size_t>(kcHost[kDecompAxis])
        : 0;

    // Number of "slices" along the decomposition axis
    const size_t totalSlices = imgSzHost[decompAxis];

    // Pixels per slice (all image axes except the decomposition axis)
    size_t slicePixels = 1;
    for(size_t d = 0; d < decompAxis; ++d)
        slicePixels *= imgSzHost[d];

    const size_t pxBytes = pixelSize(pixelType);
    const size_t sliceBytes = slicePixels * pxBytes;

    // ------------------------------------------------------------------
    // Compute per-rank slice assignment
    // ------------------------------------------------------------------
    std::vector<size_t> rankSliceCount(nRanks_);
    std::vector<size_t> rankSliceStart(nRanks_);
    {
        size_t base = totalSlices / nRanks_;
        size_t rem  = totalSlices % nRanks_;
        size_t offset = 0;
        for(int r = 0; r < nRanks_; ++r)
        {
            rankSliceCount[r] = base + (static_cast<size_t>(r) < rem ? 1 : 0);
            rankSliceStart[r] = offset;
            offset += rankSliceCount[r];
        }
    }

    const size_t mySliceStart  = rankSliceStart[rank_];
    const size_t mySliceCount  = rankSliceCount[rank_];

    // With halo
    const size_t haloStart = (mySliceStart >= haloWidth)
                                ? mySliceStart - haloWidth : 0;
    const size_t haloEnd   = std::min(mySliceStart + mySliceCount + haloWidth,
                                       totalSlices);
    const size_t haloSlices = haloEnd - haloStart;

    // Interior slice offset within the halo tile
    const size_t interiorOffset = mySliceStart - haloStart;

    // ------------------------------------------------------------------
    // Broadcast kernel metadata to all ranks
    // ------------------------------------------------------------------
    std::vector<size_t>    ks_bcast(kDim);
    std::vector<ptrdiff_t> kc_bcast(kDim);
    std::vector<float>     kw_bcast(kSize);
    std::vector<size_t>    is_bcast(imgDim);

    if(rank_ == 0)
    {
        std::memcpy(ks_bcast.data(), kSizes.ptr,   kDim  * sizeof(size_t));
        std::memcpy(kc_bcast.data(), kCenter.ptr,  kDim  * sizeof(ptrdiff_t));
        std::memcpy(kw_bcast.data(), kernelW.ptr,  kSize * sizeof(float));
        std::memcpy(is_bcast.data(), imgSizes.ptr, imgDim* sizeof(size_t));
    }
    MPI_Bcast(ks_bcast.data(), static_cast<int>(kDim  * sizeof(size_t)),
              MPI_BYTE, 0, comm_);
    MPI_Bcast(kc_bcast.data(), static_cast<int>(kDim  * sizeof(ptrdiff_t)),
              MPI_BYTE, 0, comm_);
    MPI_Bcast(kw_bcast.data(), static_cast<int>(kSize * sizeof(float)),
              MPI_BYTE, 0, comm_);
    MPI_Bcast(is_bcast.data(), static_cast<int>(imgDim* sizeof(size_t)),
              MPI_BYTE, 0, comm_);

    // ------------------------------------------------------------------
    // Scatter image data + mask with halo via MPI_Scatterv
    // ------------------------------------------------------------------
    std::vector<uint8_t> tileData(haloSlices * sliceBytes);
    std::vector<uint8_t> tileMask;
    const bool doMask = (mask.ptr && mask.bytes > 0);
    if(doMask) tileMask.resize(haloSlices * slicePixels);

    // Build scatter counts and displacements (in bytes) for each rank
    std::vector<int> scatterCounts(nRanks_), scatterDispls(nRanks_);
    for(int r = 0; r < nRanks_; ++r)
    {
        size_t rHaloStart = (rankSliceStart[r] >= haloWidth)
                                ? rankSliceStart[r] - haloWidth : 0;
        size_t rHaloEnd   = std::min(rankSliceStart[r] + rankSliceCount[r] + haloWidth,
                                      totalSlices);
        size_t rHalo      = rHaloEnd - rHaloStart;
        scatterCounts[r]  = static_cast<int>(rHalo * sliceBytes);
        scatterDispls[r]  = static_cast<int>(rHaloStart * sliceBytes);
    }

    MPI_Scatterv(rank_ == 0 ? static_cast<const uint8_t*>(in.ptr) : nullptr,
                 scatterCounts.data(), scatterDispls.data(), MPI_BYTE,
                 tileData.data(), static_cast<int>(haloSlices * sliceBytes),
                 MPI_BYTE, 0, comm_);

    if(doMask)
    {
        std::vector<int> mScatterCounts(nRanks_), mScatterDispls(nRanks_);
        for(int r = 0; r < nRanks_; ++r)
        {
            size_t rHaloStart = (rankSliceStart[r] >= haloWidth)
                                    ? rankSliceStart[r] - haloWidth : 0;
            size_t rHaloEnd   = std::min(rankSliceStart[r] + rankSliceCount[r] + haloWidth,
                                          totalSlices);
            size_t rHalo      = rHaloEnd - rHaloStart;
            mScatterCounts[r] = static_cast<int>(rHalo * slicePixels);
            mScatterDispls[r] = static_cast<int>(rHaloStart * slicePixels);
        }
        MPI_Scatterv(rank_ == 0 ? static_cast<const uint8_t*>(mask.ptr) : nullptr,
                     mScatterCounts.data(), mScatterDispls.data(), MPI_BYTE,
                     tileMask.data(), static_cast<int>(haloSlices * slicePixels),
                     MPI_BYTE, 0, comm_);
    }

    // ------------------------------------------------------------------
    // Build tile image sizes for the per-node kernel
    // ------------------------------------------------------------------
    std::vector<size_t> tileImgSizes(is_bcast.begin(), is_bcast.end());
    tileImgSizes[decompAxis] = haloSlices;
    const size_t tilePixels  = haloSlices * slicePixels;
    const size_t tileBytes   = tilePixels * pxBytes;

    // ------------------------------------------------------------------
    // Upload tile to per-node backend and run
    // ------------------------------------------------------------------
    DeviceBuffer d_in   = perNode_->allocate(tileBytes);
    DeviceBuffer d_mask = doMask
                            ? perNode_->allocate(haloSlices * slicePixels * sizeof(bool))
                            : DeviceBuffer{nullptr, 0};
    DeviceBuffer d_out  = perNode_->allocate(tileBytes);
    DeviceBuffer d_kw   = perNode_->allocate(kSize * sizeof(float));
    DeviceBuffer d_ks   = perNode_->allocate(kDim  * sizeof(size_t));
    DeviceBuffer d_kc   = perNode_->allocate(kDim  * sizeof(ptrdiff_t));
    DeviceBuffer d_is   = perNode_->allocate(imgDim * sizeof(size_t));

    perNode_->upload(d_in, tileData.data(), tileBytes);
    if(doMask)
        perNode_->upload(d_mask, tileMask.data(), haloSlices * slicePixels * sizeof(bool));
    perNode_->upload(d_kw, kw_bcast.data(), kSize * sizeof(float));
    perNode_->upload(d_ks, ks_bcast.data(), kDim  * sizeof(size_t));
    perNode_->upload(d_kc, kc_bcast.data(), kDim  * sizeof(ptrdiff_t));
    perNode_->upload(d_is, tileImgSizes.data(), imgDim * sizeof(size_t));

    perNode_->runKernel(op, d_in, d_mask, d_out, d_kw, d_ks, d_kc,
                         kDim, kSize, d_is, imgDim, tilePixels,
                         pixelType, normalize);

    // Download local tile result
    std::vector<uint8_t> tileOut(tileBytes);
    perNode_->download(tileOut.data(), d_out, tileBytes);

    perNode_->free(d_in);
    if(d_mask.ptr) perNode_->free(d_mask);
    perNode_->free(d_out);
    perNode_->free(d_kw);
    perNode_->free(d_ks);
    perNode_->free(d_kc);
    perNode_->free(d_is);

    // ------------------------------------------------------------------
    // Extract interior rows (discard halo rows from output)
    // ------------------------------------------------------------------
    const size_t interiorSlices = mySliceCount;
    const size_t interiorBytes  = interiorSlices * sliceBytes;
    const uint8_t* interiorPtr  = tileOut.data() + interiorOffset * sliceBytes;

    // ------------------------------------------------------------------
    // Gather interior output to rank 0 via MPI_Gatherv
    // ------------------------------------------------------------------
    std::vector<int> gatherCounts(nRanks_), gatherDispls(nRanks_);
    for(int r = 0; r < nRanks_; ++r)
    {
        gatherCounts[r] = static_cast<int>(rankSliceCount[r] * sliceBytes);
        gatherDispls[r] = static_cast<int>(rankSliceStart[r] * sliceBytes);
    }

    MPI_Gatherv(interiorPtr, static_cast<int>(interiorBytes), MPI_BYTE,
                rank_ == 0 ? static_cast<uint8_t*>(out.ptr) : nullptr,
                gatherCounts.data(), gatherDispls.data(), MPI_BYTE,
                0, comm_);
}

}} // namespace DST::GPU
