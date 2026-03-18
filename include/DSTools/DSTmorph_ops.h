//
//  DSTmorph_ops.h
//
//  File created by GILLARD William
//  Centre de Physic des Particules de Marseille
//  Licensed under CC BY-NC 4.0
//  You may share and adapt this code with attribution,
//  but not for commercial purposes.
//  Licence text: https://creativecommons.org/licenses/by-nc/4.0/
//
//  Per-pixel morphological operations — portable across CPU and GPU.
//
//  Design rules:
//    • No STL types in per-pixel function scope (no std::valarray, std::vector).
//    • Only <cstddef>, <limits>, <cmath> in scope.
//    • All functions marked DST_HOST_DEVICE so nvcc/hipcc compile them for device.
//    • C++17 compatible (no C++20/23 features).
//    • Index helpers replicate FITScube::PixelCoordinates / PixelIndex exactly.
//    • Kernel decode direction matches the original DSTmorpho.h per operation:
//        dilation / folding / correlate → ascending (axis 0 fastest, forward loop)
//        erosion                        → descending (axis kDim-1 fastest, reverse loop)

#ifndef __DST_morph_ops__
#define __DST_morph_ops__

#include <cstddef>    // size_t, ptrdiff_t
#include <limits>     // std::numeric_limits
#include <cmath>      // std::sqrt

// ---------------------------------------------------------------------------
// GPU portability macro
// ---------------------------------------------------------------------------
#if defined(__CUDACC__) || defined(__HIPCC__)
#  define DST_HOST_DEVICE __host__ __device__
#else
#  define DST_HOST_DEVICE
#endif

namespace DST { namespace Morphology { namespace ops {

// ---------------------------------------------------------------------------
// Compile-time limits
// ---------------------------------------------------------------------------

/// Maximum image / kernel dimensionality supported via stack arrays.
static constexpr size_t MAX_KDIM = 8;

/// Maximum total kernel elements for MEDIAN stack buffer.
/// Covers up to a ~22×22 kernel (512 elements, 2 KB per GPU thread).
/// Kernels exceeding this limit with MEDIAN fall back to the CPU path.
static constexpr size_t MAX_KERNEL_SIZE = 512;

// ---------------------------------------------------------------------------
// fold_op: integer encoding of folding operations.
// Values match ktype enum in DSTmorpho.h for direct static_cast compatibility.
// ---------------------------------------------------------------------------
enum class fold_op : int
{
    MIN    = 0,
    MAX    = 1,
    MEAN   = 2,
    MEDIAN = 3,
    SUM    = 4,
    RMS    = 5,
    RMSE   = 6
};

// ---------------------------------------------------------------------------
// Index helpers — replicate FITScube::PixelCoordinates and PixelIndex exactly
// ---------------------------------------------------------------------------

/// Flat pixel index → per-axis coordinates.
/// Matches FITScube::PixelCoordinates(k): coords[i] = (k / stride_i) % imgSizes[i]
/// where stride_i = product(imgSizes[0..i-1]).  Axis 0 is fastest varying.
DST_HOST_DEVICE inline
void decode_pixel_coords(size_t k,
                         const size_t* imgSizes, size_t imgDim,
                         size_t* coords)
{
    for(size_t i = 0; i < imgDim; ++i)
    {
        size_t stride = 1;
        for(size_t l = 0; l < i; ++l) stride *= imgSizes[l];
        coords[i] = (k / stride) % imgSizes[i];
    }
}

/// Per-axis coordinates → flat pixel index.
/// Matches FITScube::PixelIndex(iPx): index = sum_i( coords[i] * stride_i )
/// where stride_i = product(imgSizes[0..i-1]).
DST_HOST_DEVICE inline
size_t compute_flat_index(const size_t* coords,
                          const size_t* imgSizes, size_t imgDim)
{
    size_t index = 0;
    for(size_t naxe = 0; naxe < imgDim; ++naxe)
    {
        size_t arg = coords[naxe];
        for(size_t iX = 0; iX < naxe; ++iX) arg *= imgSizes[iX];
        index += arg;
    }
    return index;
}

/// Kernel flat index → per-axis kernel index, ASCENDING order (axis 0 fastest).
/// Used by dilation, folding, correlate — matches those functions in DSTmorpho.h.
DST_HOST_DEVICE inline
void decode_kernel_idx_asc(size_t kj,
                            const size_t* kSizes, size_t kDim,
                            size_t* kIdx)
{
    size_t rem = kj;
    for(size_t ka = 0; ka < kDim; ++ka)
    {
        kIdx[ka] = rem % kSizes[ka];
        rem /= kSizes[ka];
    }
}

/// Kernel flat index → per-axis kernel index, DESCENDING order (axis kDim-1 fastest).
/// Used by erosion — preserves the original erosion traversal order in DSTmorpho.h.
DST_HOST_DEVICE inline
void decode_kernel_idx_desc(size_t kj,
                             const size_t* kSizes, size_t kDim,
                             size_t* kIdx)
{
    size_t rem = kj;
    for(ptrdiff_t ka = static_cast<ptrdiff_t>(kDim) - 1; ka >= 0; --ka)
    {
        kIdx[static_cast<size_t>(ka)] = rem % kSizes[static_cast<size_t>(ka)];
        rem /= kSizes[static_cast<size_t>(ka)];
    }
}

// ---------------------------------------------------------------------------
// Neighbor lookup helper
// ---------------------------------------------------------------------------

/// Given the coordinates of the current pixel and the decoded kernel multi-index,
/// compute the flat index of the neighbor addressed by this kernel element.
///
/// The kernel maps to the LAST kDim axes of the image
/// (kernel axis ka → image axis imgDim - kDim + ka), matching DSTmorpho.h.
///
/// Returns true and sets *tgtIndex_out if the neighbor is in bounds.
/// Returns false (out of bounds) otherwise.
DST_HOST_DEVICE inline
bool neighbor_index(const size_t* pxCoords,
                    const size_t* kIdx,
                    const ptrdiff_t* kCenter,
                    const size_t* imgSizes, size_t imgDim, size_t kDim,
                    size_t* tgtIndex_out)
{
    const size_t kAxis0 = imgDim - kDim;

    size_t tgtCoords[MAX_KDIM];
    for(size_t a = 0; a < imgDim; ++a) tgtCoords[a] = pxCoords[a];

    for(size_t ka = 0; ka < kDim; ++ka)
    {
        const size_t imgAxis = kAxis0 + ka;
        if(imgAxis >= imgDim) return false;

        const ptrdiff_t offset = static_cast<ptrdiff_t>(kIdx[ka]) - kCenter[ka];
        const ptrdiff_t tgt   = static_cast<ptrdiff_t>(pxCoords[imgAxis]) + offset;

        if(tgt < 0 || static_cast<size_t>(tgt) >= imgSizes[imgAxis]) return false;

        tgtCoords[imgAxis] = static_cast<size_t>(tgt);
    }

    *tgtIndex_out = compute_flat_index(tgtCoords, imgSizes, imgDim);
    return true;
}

// ---------------------------------------------------------------------------
// Insertion sort for small arrays (MEDIAN on GPU)
// ---------------------------------------------------------------------------
template<typename T>
DST_HOST_DEVICE inline
void sort_small(T* arr, size_t n)
{
    for(size_t i = 1; i < n; ++i)
    {
        T key = arr[i];
        ptrdiff_t j = static_cast<ptrdiff_t>(i) - 1;
        while(j >= 0 && arr[j] > key)
        {
            arr[j + 1] = arr[j];
            --j;
        }
        arr[j + 1] = key;
    }
}

// ---------------------------------------------------------------------------
// DILATION — per-pixel
// ---------------------------------------------------------------------------
/// Computes the morphological dilation at pixel pixelIdx.
/// Running max over the kernel footprint (kernelWeights[kj] > 0).
/// OOB and masked pixels are excluded. Returns U(0) if footprint is empty.
/// Kernel decode: ascending (matches dilation in DSTmorpho.h).
///
/// \tparam U  Promoted computation type (kernel<T>::promote_with<K>)
/// \tparam K  Input image pixel type
template<typename U, typename K>
DST_HOST_DEVICE
U pixel_dilation(const K*         in_data,
                 const bool*      mask,
                 const U*         kernelWeights,
                 const size_t*    kSizes,
                 size_t           kDim,
                 size_t           kSize,
                 const ptrdiff_t* kCenter,
                 const size_t*    imgSizes,
                 size_t           imgDim,
                 size_t           pixelIdx)
{
    size_t pxCoords[MAX_KDIM];
    decode_pixel_coords(pixelIdx, imgSizes, imgDim, pxCoords);

    U    result = std::numeric_limits<U>::lowest();
    bool found  = false;

    for(size_t kj = 0; kj < kSize; ++kj)
    {
        if(!(kernelWeights[kj] > U(0))) continue;

        size_t kIdx[MAX_KDIM];
        decode_kernel_idx_asc(kj, kSizes, kDim, kIdx);

        size_t tgtIndex;
        if(!neighbor_index(pxCoords, kIdx, kCenter, imgSizes, imgDim, kDim, &tgtIndex))
            continue;

        if(mask && mask[tgtIndex]) continue;

        const U v = static_cast<U>(in_data[tgtIndex]);
        if(!found || v > result) { result = v; found = true; }
    }

    return found ? result : U(0);
}

// ---------------------------------------------------------------------------
// EROSION — per-pixel
// ---------------------------------------------------------------------------
/// Computes the morphological erosion at pixel pixelIdx.
/// Running min over the kernel footprint (kernelWeights[kj] > 0).
/// OOB and masked pixels are excluded. Returns U(0) if footprint is empty.
/// Kernel decode: DESCENDING (matches erosion in DSTmorpho.h exactly).
///
/// \tparam U  Promoted computation type
/// \tparam K  Input image pixel type
template<typename U, typename K>
DST_HOST_DEVICE
U pixel_erosion(const K*         in_data,
                const bool*      mask,
                const U*         kernelWeights,
                const size_t*    kSizes,
                size_t           kDim,
                size_t           kSize,
                const ptrdiff_t* kCenter,
                const size_t*    imgSizes,
                size_t           imgDim,
                size_t           pixelIdx)
{
    size_t pxCoords[MAX_KDIM];
    decode_pixel_coords(pixelIdx, imgSizes, imgDim, pxCoords);

    U    result = std::numeric_limits<U>::max();
    bool found  = false;

    for(size_t kj = 0; kj < kSize; ++kj)
    {
        if(!(kernelWeights[kj] > U(0))) continue;

        size_t kIdx[MAX_KDIM];
        decode_kernel_idx_desc(kj, kSizes, kDim, kIdx);   // descending — matches original

        size_t tgtIndex;
        if(!neighbor_index(pxCoords, kIdx, kCenter, imgSizes, imgDim, kDim, &tgtIndex))
            continue;

        if(mask && mask[tgtIndex]) continue;

        const U v = static_cast<U>(in_data[tgtIndex]);
        if(!found || v < result) { result = v; found = true; }
    }

    return found ? result : U(0);
}

// ---------------------------------------------------------------------------
// CORRELATE — per-pixel
// ---------------------------------------------------------------------------
/// Weighted sum: out(p) = sum_x in(p+x) * k(x).
/// If normalize=true, divides by sum(|w|) over valid samples.
/// OOB and masked pixels are skipped.
/// Kernel decode: ascending (matches correlate in DSTmorpho.h).
///
/// \tparam U  Promoted computation type
/// \tparam K  Input image pixel type
template<typename U, typename K>
DST_HOST_DEVICE
U pixel_correlate(const K*         in_data,
                  const bool*      mask,
                  const U*         kernelWeights,
                  const size_t*    kSizes,
                  size_t           kDim,
                  size_t           kSize,
                  const ptrdiff_t* kCenter,
                  const size_t*    imgSizes,
                  size_t           imgDim,
                  size_t           pixelIdx,
                  bool             normalize)
{
    size_t pxCoords[MAX_KDIM];
    decode_pixel_coords(pixelIdx, imgSizes, imgDim, pxCoords);

    U acc  = U(0);
    U wsum = U(0);

    for(size_t kj = 0; kj < kSize; ++kj)
    {
        size_t kIdx[MAX_KDIM];
        decode_kernel_idx_asc(kj, kSizes, kDim, kIdx);

        size_t tgtIndex;
        if(!neighbor_index(pxCoords, kIdx, kCenter, imgSizes, imgDim, kDim, &tgtIndex))
            continue;

        if(mask && mask[tgtIndex]) continue;

        const U w = kernelWeights[kj];
        const U v = static_cast<U>(in_data[tgtIndex]);
        acc += v * w;

        if(normalize)
        {
            const U aw = (w < U(0)) ? U(-w) : w;
            if(aw == aw) wsum += aw;   // NaN guard on kernel weight
        }
    }

    if(normalize && wsum != U(0)) acc /= wsum;

    return acc;
}

// ---------------------------------------------------------------------------
// FOLD — per-pixel
// ---------------------------------------------------------------------------
/// General folding operation dispatched by fold_op.
/// Replicates folding() in DSTmorpho.h:
///   - OOB pixels:    tmp[kj] = U(0),  ttmp[kj] = 0 * w  (zero-pad boundary)
///   - Masked pixels: tmp[kj] = quiet_NaN(), ttmp[kj] = NaN*w → NaN → filtered
///   - Valid pixels:  tmp[kj] = pixel value,  ttmp[kj] = value * w
///   Filter: ttmp_val == ttmp_val (IEEE-754 NaN check, filters masked floats)
///
/// MEDIAN only works when kSize <= MAX_KERNEL_SIZE; caller must ensure this
/// or dispatch to a CPU-only path.
/// Kernel decode: ascending (matches folding in DSTmorpho.h).
///
/// \tparam U  Promoted computation type
/// \tparam K  Input image pixel type
template<typename U, typename K>
DST_HOST_DEVICE
U pixel_fold(const K*         in_data,
             const bool*      mask,
             const U*         kernelWeights,
             const size_t*    kSizes,
             size_t           kDim,
             size_t           kSize,
             const ptrdiff_t* kCenter,
             const size_t*    imgSizes,
             size_t           imgDim,
             size_t           pixelIdx,
             fold_op          op)
{
    size_t pxCoords[MAX_KDIM];
    decode_pixel_coords(pixelIdx, imgSizes, imgDim, pxCoords);

    // Running accumulators
    U      vsum   = U(0);
    U      vsumsq = U(0);
    U      vmin   = std::numeric_limits<U>::max();
    U      vmax   = std::numeric_limits<U>::lowest();
    size_t count  = 0;

    // MEDIAN stack buffer (only populated when op == MEDIAN)
    U      medBuf[MAX_KERNEL_SIZE];
    size_t medCount = 0;

    for(size_t kj = 0; kj < kSize; ++kj)
    {
        const U w = kernelWeights[kj];

        size_t kIdx[MAX_KDIM];
        decode_kernel_idx_asc(kj, kSizes, kDim, kIdx);

        // Compute weighted value (ttmp[kj] = tmp[kj] * kernelWeight[kj])
        U ttmp_val;

        size_t tgtIndex;
        if(!neighbor_index(pxCoords, kIdx, kCenter, imgSizes, imgDim, kDim, &tgtIndex))
        {
            // OOB: tmp[kj] = U(0)  →  ttmp[kj] = U(0) * w
            // For finite w this is 0; NaN*0 is NaN (filtered below).
            ttmp_val = U(0) * w;
        }
        else if(mask && mask[tgtIndex])
        {
            // Masked: tmp[kj] = quiet_NaN()  →  ttmp[kj] = NaN*w = NaN (floats)
            // For integer U quiet_NaN() == 0, so masked pixels contribute 0
            // (same behaviour as the original code).
            ttmp_val = std::numeric_limits<U>::quiet_NaN() * w;
        }
        else
        {
            ttmp_val = static_cast<U>(in_data[tgtIndex]) * w;
        }

        // IEEE-754 NaN filter: ttmp[ttmp==ttmp] — matches the original valarray idiom
        if(ttmp_val != ttmp_val) continue;

        ++count;
        vsum   += ttmp_val;
        if(ttmp_val < vmin) vmin = ttmp_val;
        if(ttmp_val > vmax) vmax = ttmp_val;
        vsumsq += ttmp_val * ttmp_val;

        if(op == fold_op::MEDIAN && medCount < MAX_KERNEL_SIZE)
            medBuf[medCount++] = ttmp_val;
    }

    if(count == 0) return U(0);

    switch(op)
    {
        case fold_op::MIN:  return vmin;
        case fold_op::MAX:  return vmax;
        case fold_op::SUM:  return vsum;

        case fold_op::MEAN:
            return vsum / static_cast<U>(count);

        case fold_op::RMS:
            return static_cast<U>(
                std::sqrt(static_cast<double>(vsumsq) /
                          static_cast<double>(count)));

        case fold_op::RMSE:
        {
            // Original: sumsq = sum(v_i^2 - mean^2), sqrt(sumsq)
            const U mean = vsum / static_cast<U>(count);
            const U sumsq = vsumsq - static_cast<U>(count) * mean * mean;
            return static_cast<U>(std::sqrt(static_cast<double>(sumsq)));
        }

        case fold_op::MEDIAN:
        {
            if(medCount == 0) return U(0);
            sort_small(medBuf, medCount);
            if(medCount % 2 == 0)
                return (medBuf[medCount / 2 - 1] + medBuf[medCount / 2]) / U(2);
            else
                return medBuf[medCount / 2];
        }

        default:
            return vsum;
    }
}

}}} // namespace DST::Morphology::ops

#endif // __DST_morph_ops__
