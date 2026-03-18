//
//  DSTmorph_gpu.h
//
//  Host-only GPU dispatch wrappers for morphological operations.
//
//  Each function mirrors the DST::Morphology:: CPU API but adds a
//  DST::GPU::Backend& parameter.  At compile time, if constexpr routes:
//    • Non-GPU-eligible pixel types (double, int64, uint64, size_t)
//      → CPU path (DST::Morphology:: functions), nThreads is used.
//    • GPU-eligible pixel types (all ≤ 32-bit integers + float)
//      → GPUContext path if backend.available(), else CPU fallback.
//
//  Output type on GPU path is always K (the input pixel type).
//  Output type on CPU path is the promoted type U = kernel<T>::promote_with<K>,
//  which may differ from K when the kernel type T promotes the result.
//
//  Primitive operations  (dilation, erosion, correlate, folding):
//    One GPUContext, one run(), one download.
//
//  Composition operations (opening, closing):
//    One GPUContext, two chained run() calls — data stays device-resident
//    between the two morphological passes.
//
//  Arithmetic compositions (wTopHat, bTopHat, gradiant):
//    GPU for the morphological passes; pixel arithmetic done on the host
//    after a single download.
//
//  Licensed under CC BY-NC 4.0

#ifndef __DST_morph_gpu__
#define __DST_morph_gpu__

#include "DSTmorpho.h"            // CPU path + kernel<T>, ktype, FITScube …

#include <DSTfits/FITSdata.h>     // is_gpu_eligible_type_v

#include <DSTools/gpu/Backend.h>    // DST::GPU::Backend
#include <DSTools/gpu/GPUContext.h> // DST::GPU::GPUContext<T>

#include <memory>
#include <vector>
#include <stdexcept>

// ============================================================
// Convenience macro: build a typed kernel with a different ktype
// (same weights and sizes, only the operation tag changes).
// ============================================================
#define DST_GPU_RETYPE_KERNEL(k, new_ktype) \
    DST::Morphology::kernel<T>((k).data(), (k).sizes(), (new_ktype))

namespace DST { namespace Morphology { namespace GPU {

// ------------------------------------------------------------
// Internal helpers
// ------------------------------------------------------------
namespace detail {

/// Allocate an output FITSimg<K> with the same shape and HDU as \a in.
template<typename K>
std::shared_ptr<FITScube> make_output(const FITScube& in,
                                       const char*     proc_name)
{
    std::vector<size_t> shape;
    shape.reserve(in.GetDimension());
    for(size_t d = 0; d < in.GetDimension(); ++d)
        shape.push_back(in.Size(d + 1));

    FITSimg<K>* out = new FITSimg<K>(shape);

    for(FITSDictionary::const_iterator it = in.HDU().begin();
        it != in.HDU().end(); ++it)
    {
        if(!out->HDU().Exists(it->first))
            out->HDU().ValueForKey(it->first,
                                   it->second.value(),
                                   it->second.type(),
                                   it->second.comment());
    }
    out->HDU().ValueForKey("DST_PROC", proc_name,
                            key_type::fChar,
                            "Processing history keyword");
    out->HDU().ValueForKey("DST_VER",
                            DST::gGIT::version_short(),
                            key_type::fChar,
                            "DST framwork version library");
    out->HDU().ValueForKey("DST_SHA",
                            DST::gGIT::sha1(),
                            key_type::fChar,
                            "DST framwork GIT identifier");

    return std::shared_ptr<FITScube>(out);
}

/// Run a single morphological kernel on the GPU and return the result.
/// K = image pixel type (GPU-eligible, enforced by caller).
/// T = kernel weight type.
template<typename K, typename T>
std::shared_ptr<FITScube> gpu_single(DST::GPU::Backend&                         backend,
                                      const FITScube&                             in,
                                      const DST::Morphology::kernel<T>&           k,
                                      const char*                                 proc_name)
{
    DST::GPU::GPUContext<K> ctx(backend, in);
    ctx.run(k);
    std::shared_ptr<FITScube> out = make_output<K>(in, proc_name);
    ctx.download(*out);
    return out;
}

/// Run two chained morphological kernels on the GPU and return the result.
/// Data stays device-resident between k1 and k2 — no PCIe round-trip.
template<typename K, typename T>
std::shared_ptr<FITScube> gpu_chain2(DST::GPU::Backend&               backend,
                                      const FITScube&                  in,
                                      const DST::Morphology::kernel<T>& k1,
                                      const DST::Morphology::kernel<T>& k2,
                                      const char*                       proc_name)
{
    DST::GPU::GPUContext<K> ctx(backend, in);
    ctx.run(k1);
    ctx.run(k2);
    std::shared_ptr<FITScube> out = make_output<K>(in, proc_name);
    ctx.download(*out);
    return out;
}

} // namespace detail


// ============================================================
//  dilation
// ============================================================

/// GPU-dispatched dilation (typed image variant).
template<typename T, typename K>
std::shared_ptr<FITScube> dilation(const std::shared_ptr< FITSimg<K> >& in,
                                    const DST::Morphology::kernel<T>&    k,
                                    DST::GPU::Backend&                   backend,
                                    size_t                               nThreads = 1)
{
    if constexpr (!DSL::is_gpu_eligible_type_v<K>)
        return DST::Morphology::dilation<T, K>(in, k, nThreads);
    else
    {
        if(backend.available())
            return detail::gpu_single<K, T>(backend, *in, k, "DILATION");
        else
            return DST::Morphology::dilation<T, K>(in, k, nThreads);
    }
}

/// GPU-dispatched dilation (runtime-typed FITScube variant).
template<typename T>
std::shared_ptr<FITScube> dilation(const std::shared_ptr<FITScube>& in,
                                    const DST::Morphology::kernel<T>& k,
                                    DST::GPU::Backend&                backend,
                                    size_t                            nThreads = 1)
{
    if(!in)
        throw FITSexception(BAD_DIMEN,"DST::GPU","dilation","Input pointer is null.");
    switch(in->GetEqBitPerPixel())
    {
        case SBYTE_IMG:    return dilation<T,int8_t>  (std::dynamic_pointer_cast<FITSimg<int8_t>  >(in), k, backend, nThreads);
        case BYTE_IMG:     return dilation<T,uint8_t> (std::dynamic_pointer_cast<FITSimg<uint8_t> >(in), k, backend, nThreads);
        case SHORT_IMG:    return dilation<T,int16_t> (std::dynamic_pointer_cast<FITSimg<int16_t> >(in), k, backend, nThreads);
        case USHORT_IMG:   return dilation<T,uint16_t>(std::dynamic_pointer_cast<FITSimg<uint16_t>>(in), k, backend, nThreads);
        case LONG_IMG:     return dilation<T,int32_t> (std::dynamic_pointer_cast<FITSimg<int32_t> >(in), k, backend, nThreads);
        case ULONG_IMG:    return dilation<T,uint32_t>(std::dynamic_pointer_cast<FITSimg<uint32_t>>(in), k, backend, nThreads);
        case LONGLONG_IMG: return dilation<T,int64_t> (std::dynamic_pointer_cast<FITSimg<int64_t> >(in), k, backend, nThreads);
        case ULONGLONG_IMG:return dilation<T,uint64_t>(std::dynamic_pointer_cast<FITSimg<uint64_t>>(in), k, backend, nThreads);
        case FLOAT_IMG:    return dilation<T,float>   (std::dynamic_pointer_cast<FITSimg<float>   >(in), k, backend, nThreads);
        case DOUBLE_IMG:   return dilation<T,double>  (std::dynamic_pointer_cast<FITSimg<double>  >(in), k, backend, nThreads);
        default: break;
    }
    throw FITSexception(BAD_DIMEN,"DST::GPU","dilation","Unsupported pixel type.");
}


// ============================================================
//  erosion
// ============================================================

template<typename T, typename K>
std::shared_ptr<FITScube> erosion(const std::shared_ptr< FITSimg<K> >& in,
                                   const DST::Morphology::kernel<T>&    k,
                                   DST::GPU::Backend&                   backend,
                                   size_t                               nThreads = 1)
{
    if constexpr (!DSL::is_gpu_eligible_type_v<K>)
        return DST::Morphology::erosion<T, K>(in, k, nThreads);
    else
    {
        if(backend.available())
            return detail::gpu_single<K, T>(backend, *in, k, "EROSION");
        else
            return DST::Morphology::erosion<T, K>(in, k, nThreads);
    }
}

template<typename T>
std::shared_ptr<FITScube> erosion(const std::shared_ptr<FITScube>& in,
                                   const DST::Morphology::kernel<T>& k,
                                   DST::GPU::Backend&                backend,
                                   size_t                            nThreads = 1)
{
    if(!in)
        throw FITSexception(BAD_DIMEN,"DST::GPU","erosion","Input pointer is null.");
    switch(in->GetEqBitPerPixel())
    {
        case SBYTE_IMG:    return erosion<T,int8_t>  (std::dynamic_pointer_cast<FITSimg<int8_t>  >(in), k, backend, nThreads);
        case BYTE_IMG:     return erosion<T,uint8_t> (std::dynamic_pointer_cast<FITSimg<uint8_t> >(in), k, backend, nThreads);
        case SHORT_IMG:    return erosion<T,int16_t> (std::dynamic_pointer_cast<FITSimg<int16_t> >(in), k, backend, nThreads);
        case USHORT_IMG:   return erosion<T,uint16_t>(std::dynamic_pointer_cast<FITSimg<uint16_t>>(in), k, backend, nThreads);
        case LONG_IMG:     return erosion<T,int32_t> (std::dynamic_pointer_cast<FITSimg<int32_t> >(in), k, backend, nThreads);
        case ULONG_IMG:    return erosion<T,uint32_t>(std::dynamic_pointer_cast<FITSimg<uint32_t>>(in), k, backend, nThreads);
        case LONGLONG_IMG: return erosion<T,int64_t> (std::dynamic_pointer_cast<FITSimg<int64_t> >(in), k, backend, nThreads);
        case ULONGLONG_IMG:return erosion<T,uint64_t>(std::dynamic_pointer_cast<FITSimg<uint64_t>>(in), k, backend, nThreads);
        case FLOAT_IMG:    return erosion<T,float>   (std::dynamic_pointer_cast<FITSimg<float>   >(in), k, backend, nThreads);
        case DOUBLE_IMG:   return erosion<T,double>  (std::dynamic_pointer_cast<FITSimg<double>  >(in), k, backend, nThreads);
        default: break;
    }
    throw FITSexception(BAD_DIMEN,"DST::GPU","erosion","Unsupported pixel type.");
}


// ============================================================
//  correlate
// ============================================================

template<typename T, typename K>
std::shared_ptr<FITScube> correlate(const std::shared_ptr< FITSimg<K> >& in,
                                     const DST::Morphology::kernel<T>&    k,
                                     DST::GPU::Backend&                   backend,
                                     size_t                               nThreads  = 1,
                                     bool                                 normalize = false)
{
    if constexpr (!DSL::is_gpu_eligible_type_v<K>)
        return DST::Morphology::correlate<T, K>(in, k, nThreads, normalize);
    else
    {
        if(backend.available())
        {
            // Build a kernel tagged with the correct correlate ktype
            const auto ktype_tag = normalize ? DST::Morphology::ktype::nCORRELATE
                                             : DST::Morphology::ktype::CORRELATE;
            const DST::Morphology::kernel<T> k_corr(k.data(), k.sizes(), ktype_tag);
            return detail::gpu_single<K, T>(backend, *in, k_corr,
                                             normalize ? "nCORRELATE" : "CORRELATE");
        }
        else
            return DST::Morphology::correlate<T, K>(in, k, nThreads, normalize);
    }
}

template<typename T>
std::shared_ptr<FITScube> correlate(const std::shared_ptr<FITScube>& in,
                                     const DST::Morphology::kernel<T>& k,
                                     DST::GPU::Backend&                backend,
                                     size_t                            nThreads  = 1,
                                     bool                              normalize = false)
{
    if(!in)
        throw FITSexception(BAD_DIMEN,"DST::GPU","correlate","Input pointer is null.");
    switch(in->GetEqBitPerPixel())
    {
        case SBYTE_IMG:    return correlate<T,int8_t>  (std::dynamic_pointer_cast<FITSimg<int8_t>  >(in), k, backend, nThreads, normalize);
        case BYTE_IMG:     return correlate<T,uint8_t> (std::dynamic_pointer_cast<FITSimg<uint8_t> >(in), k, backend, nThreads, normalize);
        case SHORT_IMG:    return correlate<T,int16_t> (std::dynamic_pointer_cast<FITSimg<int16_t> >(in), k, backend, nThreads, normalize);
        case USHORT_IMG:   return correlate<T,uint16_t>(std::dynamic_pointer_cast<FITSimg<uint16_t>>(in), k, backend, nThreads, normalize);
        case LONG_IMG:     return correlate<T,int32_t> (std::dynamic_pointer_cast<FITSimg<int32_t> >(in), k, backend, nThreads, normalize);
        case ULONG_IMG:    return correlate<T,uint32_t>(std::dynamic_pointer_cast<FITSimg<uint32_t>>(in), k, backend, nThreads, normalize);
        case LONGLONG_IMG: return correlate<T,int64_t> (std::dynamic_pointer_cast<FITSimg<int64_t> >(in), k, backend, nThreads, normalize);
        case ULONGLONG_IMG:return correlate<T,uint64_t>(std::dynamic_pointer_cast<FITSimg<uint64_t>>(in), k, backend, nThreads, normalize);
        case FLOAT_IMG:    return correlate<T,float>   (std::dynamic_pointer_cast<FITSimg<float>   >(in), k, backend, nThreads, normalize);
        case DOUBLE_IMG:   return correlate<T,double>  (std::dynamic_pointer_cast<FITSimg<double>  >(in), k, backend, nThreads, normalize);
        default: break;
    }
    throw FITSexception(BAD_DIMEN,"DST::GPU","correlate","Unsupported pixel type.");
}


// ============================================================
//  opening  (erode → dilate, device-resident between passes)
// ============================================================

template<typename T, typename K>
std::shared_ptr<FITScube> opening(const std::shared_ptr< FITSimg<K> >& in,
                                   const DST::Morphology::kernel<T>&    k,
                                   DST::GPU::Backend&                   backend,
                                   size_t                               nThreads = 1)
{
    if constexpr (!DSL::is_gpu_eligible_type_v<K>)
        return DST::Morphology::opening<T, K>(in, k, nThreads);
    else
    {
        if(backend.available())
        {
            const DST::Morphology::kernel<T> k_ero(k.data(), k.sizes(),
                                                    DST::Morphology::ktype::ERODE);
            const DST::Morphology::kernel<T> k_dil(k.data(), k.sizes(),
                                                    DST::Morphology::ktype::DILATE);
            return detail::gpu_chain2<K, T>(backend, *in, k_ero, k_dil, "OPENING");
        }
        else
            return DST::Morphology::opening<T, K>(in, k, nThreads);
    }
}

template<typename T>
std::shared_ptr<FITScube> opening(const std::shared_ptr<FITScube>& in,
                                   const DST::Morphology::kernel<T>& k,
                                   DST::GPU::Backend&                backend,
                                   size_t                            nThreads = 1)
{
    if(!in)
        throw FITSexception(BAD_DIMEN,"DST::GPU","opening","Input pointer is null.");
    switch(in->GetEqBitPerPixel())
    {
        case SBYTE_IMG:    return opening<T,int8_t>  (std::dynamic_pointer_cast<FITSimg<int8_t>  >(in), k, backend, nThreads);
        case BYTE_IMG:     return opening<T,uint8_t> (std::dynamic_pointer_cast<FITSimg<uint8_t> >(in), k, backend, nThreads);
        case SHORT_IMG:    return opening<T,int16_t> (std::dynamic_pointer_cast<FITSimg<int16_t> >(in), k, backend, nThreads);
        case USHORT_IMG:   return opening<T,uint16_t>(std::dynamic_pointer_cast<FITSimg<uint16_t>>(in), k, backend, nThreads);
        case LONG_IMG:     return opening<T,int32_t> (std::dynamic_pointer_cast<FITSimg<int32_t> >(in), k, backend, nThreads);
        case ULONG_IMG:    return opening<T,uint32_t>(std::dynamic_pointer_cast<FITSimg<uint32_t>>(in), k, backend, nThreads);
        case LONGLONG_IMG: return opening<T,int64_t> (std::dynamic_pointer_cast<FITSimg<int64_t> >(in), k, backend, nThreads);
        case ULONGLONG_IMG:return opening<T,uint64_t>(std::dynamic_pointer_cast<FITSimg<uint64_t>>(in), k, backend, nThreads);
        case FLOAT_IMG:    return opening<T,float>   (std::dynamic_pointer_cast<FITSimg<float>   >(in), k, backend, nThreads);
        case DOUBLE_IMG:   return opening<T,double>  (std::dynamic_pointer_cast<FITSimg<double>  >(in), k, backend, nThreads);
        default: break;
    }
    throw FITSexception(BAD_DIMEN,"DST::GPU","opening","Unsupported pixel type.");
}


// ============================================================
//  closing  (dilate → erode, device-resident between passes)
// ============================================================

template<typename T, typename K>
std::shared_ptr<FITScube> closing(const std::shared_ptr< FITSimg<K> >& in,
                                   const DST::Morphology::kernel<T>&    k,
                                   DST::GPU::Backend&                   backend,
                                   size_t                               nThreads = 1)
{
    if constexpr (!DSL::is_gpu_eligible_type_v<K>)
        return DST::Morphology::closing<T, K>(in, k, nThreads);
    else
    {
        if(backend.available())
        {
            const DST::Morphology::kernel<T> k_dil(k.data(), k.sizes(),
                                                    DST::Morphology::ktype::DILATE);
            const DST::Morphology::kernel<T> k_ero(k.data(), k.sizes(),
                                                    DST::Morphology::ktype::ERODE);
            return detail::gpu_chain2<K, T>(backend, *in, k_dil, k_ero, "CLOSING");
        }
        else
            return DST::Morphology::closing<T, K>(in, k, nThreads);
    }
}

template<typename T>
std::shared_ptr<FITScube> closing(const std::shared_ptr<FITScube>& in,
                                   const DST::Morphology::kernel<T>& k,
                                   DST::GPU::Backend&                backend,
                                   size_t                            nThreads = 1)
{
    if(!in)
        throw FITSexception(BAD_DIMEN,"DST::GPU","closing","Input pointer is null.");
    switch(in->GetEqBitPerPixel())
    {
        case SBYTE_IMG:    return closing<T,int8_t>  (std::dynamic_pointer_cast<FITSimg<int8_t>  >(in), k, backend, nThreads);
        case BYTE_IMG:     return closing<T,uint8_t> (std::dynamic_pointer_cast<FITSimg<uint8_t> >(in), k, backend, nThreads);
        case SHORT_IMG:    return closing<T,int16_t> (std::dynamic_pointer_cast<FITSimg<int16_t> >(in), k, backend, nThreads);
        case USHORT_IMG:   return closing<T,uint16_t>(std::dynamic_pointer_cast<FITSimg<uint16_t>>(in), k, backend, nThreads);
        case LONG_IMG:     return closing<T,int32_t> (std::dynamic_pointer_cast<FITSimg<int32_t> >(in), k, backend, nThreads);
        case ULONG_IMG:    return closing<T,uint32_t>(std::dynamic_pointer_cast<FITSimg<uint32_t>>(in), k, backend, nThreads);
        case LONGLONG_IMG: return closing<T,int64_t> (std::dynamic_pointer_cast<FITSimg<int64_t> >(in), k, backend, nThreads);
        case ULONGLONG_IMG:return closing<T,uint64_t>(std::dynamic_pointer_cast<FITSimg<uint64_t>>(in), k, backend, nThreads);
        case FLOAT_IMG:    return closing<T,float>   (std::dynamic_pointer_cast<FITSimg<float>   >(in), k, backend, nThreads);
        case DOUBLE_IMG:   return closing<T,double>  (std::dynamic_pointer_cast<FITSimg<double>  >(in), k, backend, nThreads);
        default: break;
    }
    throw FITSexception(BAD_DIMEN,"DST::GPU","closing","Unsupported pixel type.");
}


// ============================================================
//  wTopHat  (input − opening)
//  GPU path: opening is computed device-resident (chain2),
//  then downloaded; pixel subtraction is done on the host.
// ============================================================

template<typename T, typename K>
std::shared_ptr<FITScube> wTopHat(const std::shared_ptr< FITSimg<K> >& in,
                                   const DST::Morphology::kernel<T>&    k,
                                   DST::GPU::Backend&                   backend,
                                   size_t                               nThreads = 1)
{
    if constexpr (!DSL::is_gpu_eligible_type_v<K>)
        return DST::Morphology::wTopHat<T, K>(in, k, nThreads);
    else
    {
        if(!backend.available())
            return DST::Morphology::wTopHat<T, K>(in, k, nThreads);

        // Compute opening on GPU (erode → dilate, device-resident)
        const DST::Morphology::kernel<T> k_ero(k.data(), k.sizes(),
                                                DST::Morphology::ktype::ERODE);
        const DST::Morphology::kernel<T> k_dil(k.data(), k.sizes(),
                                                DST::Morphology::ktype::DILATE);
        std::shared_ptr<FITScube> opened =
            detail::gpu_chain2<K, T>(backend, *in, k_ero, k_dil, "temp");

        // result = in − opened  (host arithmetic)
        std::shared_ptr<FITScube> result = detail::make_output<K>(*in, "White Top-Hat");
        K*       res_raw    = result->template raw_data<K>();
        const K* in_raw     = in->template raw_data<K>();
        const K* opened_raw = opened->template raw_data<K>();
        const size_t N      = in->Nelements();
        for(size_t i = 0; i < N; ++i)
            res_raw[i] = in_raw[i] - opened_raw[i];

        return result;
    }
}

template<typename T>
std::shared_ptr<FITScube> wTopHat(const std::shared_ptr<FITScube>& in,
                                   const DST::Morphology::kernel<T>& k,
                                   DST::GPU::Backend&                backend,
                                   size_t                            nThreads = 1)
{
    if(!in)
        throw FITSexception(BAD_DIMEN,"DST::GPU","wTopHat","Input pointer is null.");
    switch(in->GetEqBitPerPixel())
    {
        case SBYTE_IMG:    return wTopHat<T,int8_t>  (std::dynamic_pointer_cast<FITSimg<int8_t>  >(in), k, backend, nThreads);
        case BYTE_IMG:     return wTopHat<T,uint8_t> (std::dynamic_pointer_cast<FITSimg<uint8_t> >(in), k, backend, nThreads);
        case SHORT_IMG:    return wTopHat<T,int16_t> (std::dynamic_pointer_cast<FITSimg<int16_t> >(in), k, backend, nThreads);
        case USHORT_IMG:   return wTopHat<T,uint16_t>(std::dynamic_pointer_cast<FITSimg<uint16_t>>(in), k, backend, nThreads);
        case LONG_IMG:     return wTopHat<T,int32_t> (std::dynamic_pointer_cast<FITSimg<int32_t> >(in), k, backend, nThreads);
        case ULONG_IMG:    return wTopHat<T,uint32_t>(std::dynamic_pointer_cast<FITSimg<uint32_t>>(in), k, backend, nThreads);
        case LONGLONG_IMG: return wTopHat<T,int64_t> (std::dynamic_pointer_cast<FITSimg<int64_t> >(in), k, backend, nThreads);
        case ULONGLONG_IMG:return wTopHat<T,uint64_t>(std::dynamic_pointer_cast<FITSimg<uint64_t>>(in), k, backend, nThreads);
        case FLOAT_IMG:    return wTopHat<T,float>   (std::dynamic_pointer_cast<FITSimg<float>   >(in), k, backend, nThreads);
        case DOUBLE_IMG:   return wTopHat<T,double>  (std::dynamic_pointer_cast<FITSimg<double>  >(in), k, backend, nThreads);
        default: break;
    }
    throw FITSexception(BAD_DIMEN,"DST::GPU","wTopHat","Unsupported pixel type.");
}


// ============================================================
//  bTopHat  (closing − input)
//  GPU path: closing computed device-resident (chain2),
//  then pixel subtraction on the host.
// ============================================================

template<typename T, typename K>
std::shared_ptr<FITScube> bTopHat(const std::shared_ptr< FITSimg<K> >& in,
                                   const DST::Morphology::kernel<T>&    k,
                                   DST::GPU::Backend&                   backend,
                                   size_t                               nThreads = 1)
{
    if constexpr (!DSL::is_gpu_eligible_type_v<K>)
        return DST::Morphology::bTopHat<T, K>(in, k, nThreads);
    else
    {
        if(!backend.available())
            return DST::Morphology::bTopHat<T, K>(in, k, nThreads);

        // Compute closing on GPU (dilate → erode, device-resident)
        const DST::Morphology::kernel<T> k_dil(k.data(), k.sizes(),
                                                DST::Morphology::ktype::DILATE);
        const DST::Morphology::kernel<T> k_ero(k.data(), k.sizes(),
                                                DST::Morphology::ktype::ERODE);
        std::shared_ptr<FITScube> closed =
            detail::gpu_chain2<K, T>(backend, *in, k_dil, k_ero, "temp");

        // result = closed − in  (host arithmetic)
        std::shared_ptr<FITScube> result = detail::make_output<K>(*in, "Black Top-Hat");
        K*       res_raw    = result->template raw_data<K>();
        const K* in_raw     = in->template raw_data<K>();
        const K* closed_raw = closed->template raw_data<K>();
        const size_t N      = in->Nelements();
        for(size_t i = 0; i < N; ++i)
            res_raw[i] = closed_raw[i] - in_raw[i];

        return result;
    }
}

template<typename T>
std::shared_ptr<FITScube> bTopHat(const std::shared_ptr<FITScube>& in,
                                   const DST::Morphology::kernel<T>& k,
                                   DST::GPU::Backend&                backend,
                                   size_t                            nThreads = 1)
{
    if(!in)
        throw FITSexception(BAD_DIMEN,"DST::GPU","bTopHat","Input pointer is null.");
    switch(in->GetEqBitPerPixel())
    {
        case SBYTE_IMG:    return bTopHat<T,int8_t>  (std::dynamic_pointer_cast<FITSimg<int8_t>  >(in), k, backend, nThreads);
        case BYTE_IMG:     return bTopHat<T,uint8_t> (std::dynamic_pointer_cast<FITSimg<uint8_t> >(in), k, backend, nThreads);
        case SHORT_IMG:    return bTopHat<T,int16_t> (std::dynamic_pointer_cast<FITSimg<int16_t> >(in), k, backend, nThreads);
        case USHORT_IMG:   return bTopHat<T,uint16_t>(std::dynamic_pointer_cast<FITSimg<uint16_t>>(in), k, backend, nThreads);
        case LONG_IMG:     return bTopHat<T,int32_t> (std::dynamic_pointer_cast<FITSimg<int32_t> >(in), k, backend, nThreads);
        case ULONG_IMG:    return bTopHat<T,uint32_t>(std::dynamic_pointer_cast<FITSimg<uint32_t>>(in), k, backend, nThreads);
        case LONGLONG_IMG: return bTopHat<T,int64_t> (std::dynamic_pointer_cast<FITSimg<int64_t> >(in), k, backend, nThreads);
        case ULONGLONG_IMG:return bTopHat<T,uint64_t>(std::dynamic_pointer_cast<FITSimg<uint64_t>>(in), k, backend, nThreads);
        case FLOAT_IMG:    return bTopHat<T,float>   (std::dynamic_pointer_cast<FITSimg<float>   >(in), k, backend, nThreads);
        case DOUBLE_IMG:   return bTopHat<T,double>  (std::dynamic_pointer_cast<FITSimg<double>  >(in), k, backend, nThreads);
        default: break;
    }
    throw FITSexception(BAD_DIMEN,"DST::GPU","bTopHat","Unsupported pixel type.");
}


// ============================================================
//  gradiant  (dilation − erosion)
//  GPU path: two separate GPUContexts (one per primitive);
//  pixel subtraction on the host.  Unlike the CPU path, the
//  two GPU operations are not concurrent (they share the same
//  backend/device queue), but both avoid PCIe round-trips
//  for intermediate results.
// ============================================================

template<typename T, typename K>
std::shared_ptr<FITScube> gradiant(const std::shared_ptr< FITSimg<K> >& in,
                                    const DST::Morphology::kernel<T>&    k,
                                    DST::GPU::Backend&                   backend,
                                    size_t                               nThreads = 1)
{
    if constexpr (!DSL::is_gpu_eligible_type_v<K>)
        return DST::Morphology::gradiant<T, K>(in, k, nThreads);
    else
    {
        if(!backend.available())
            return DST::Morphology::gradiant<T, K>(in, k, nThreads);

        const DST::Morphology::kernel<T> k_dil(k.data(), k.sizes(),
                                                DST::Morphology::ktype::DILATE);
        const DST::Morphology::kernel<T> k_ero(k.data(), k.sizes(),
                                                DST::Morphology::ktype::ERODE);

        std::shared_ptr<FITScube> dilated =
            detail::gpu_single<K, T>(backend, *in, k_dil, "temp");
        std::shared_ptr<FITScube> eroded  =
            detail::gpu_single<K, T>(backend, *in, k_ero, "temp");

        // result = dilated − eroded  (host arithmetic)
        std::shared_ptr<FITScube> result = detail::make_output<K>(*in, "GRADIANT");
        K*       res_raw = result->template raw_data<K>();
        const K* dil_raw = dilated->template raw_data<K>();
        const K* ero_raw = eroded->template raw_data<K>();
        const size_t N   = in->Nelements();
        for(size_t i = 0; i < N; ++i)
            res_raw[i] = dil_raw[i] - ero_raw[i];

        return result;
    }
}

template<typename T>
std::shared_ptr<FITScube> gradiant(const std::shared_ptr<FITScube>& in,
                                    const DST::Morphology::kernel<T>& k,
                                    DST::GPU::Backend&                backend,
                                    size_t                            nThreads = 1)
{
    if(!in)
        throw FITSexception(BAD_DIMEN,"DST::GPU","gradiant","Input pointer is null.");
    switch(in->GetEqBitPerPixel())
    {
        case SBYTE_IMG:    return gradiant<T,int8_t>  (std::dynamic_pointer_cast<FITSimg<int8_t>  >(in), k, backend, nThreads);
        case BYTE_IMG:     return gradiant<T,uint8_t> (std::dynamic_pointer_cast<FITSimg<uint8_t> >(in), k, backend, nThreads);
        case SHORT_IMG:    return gradiant<T,int16_t> (std::dynamic_pointer_cast<FITSimg<int16_t> >(in), k, backend, nThreads);
        case USHORT_IMG:   return gradiant<T,uint16_t>(std::dynamic_pointer_cast<FITSimg<uint16_t>>(in), k, backend, nThreads);
        case LONG_IMG:     return gradiant<T,int32_t> (std::dynamic_pointer_cast<FITSimg<int32_t> >(in), k, backend, nThreads);
        case ULONG_IMG:    return gradiant<T,uint32_t>(std::dynamic_pointer_cast<FITSimg<uint32_t>>(in), k, backend, nThreads);
        case LONGLONG_IMG: return gradiant<T,int64_t> (std::dynamic_pointer_cast<FITSimg<int64_t> >(in), k, backend, nThreads);
        case ULONGLONG_IMG:return gradiant<T,uint64_t>(std::dynamic_pointer_cast<FITSimg<uint64_t>>(in), k, backend, nThreads);
        case FLOAT_IMG:    return gradiant<T,float>   (std::dynamic_pointer_cast<FITSimg<float>   >(in), k, backend, nThreads);
        case DOUBLE_IMG:   return gradiant<T,double>  (std::dynamic_pointer_cast<FITSimg<double>  >(in), k, backend, nThreads);
        default: break;
    }
    throw FITSexception(BAD_DIMEN,"DST::GPU","gradiant","Unsupported pixel type.");
}


// ============================================================
//  folding  — dispatches to the appropriate primitive
// ============================================================

template<typename T, typename K>
std::shared_ptr<FITScube> folding(const std::shared_ptr< FITSimg<K> >& in,
                                   const DST::Morphology::kernel<T>&    k,
                                   DST::GPU::Backend&                   backend,
                                   size_t                               nThreads = 1)
{
    using namespace DST::Morphology;
    switch(k.type())
    {
        case ktype::DILATE:       return dilation<T, K>(in, k, backend, nThreads);
        case ktype::ERODE:        return erosion <T, K>(in, k, backend, nThreads);
        case ktype::OPEN:         return opening <T, K>(in, k, backend, nThreads);
        case ktype::CLOSE:        return closing <T, K>(in, k, backend, nThreads);
        case ktype::WHITE_TOPHAT: return wTopHat <T, K>(in, k, backend, nThreads);
        case ktype::BLACK_TOPHAT: return bTopHat <T, K>(in, k, backend, nThreads);
        case ktype::CORRELATE:    return correlate<T, K>(in, k, backend, nThreads, false);
        case ktype::nCORRELATE:   return correlate<T, K>(in, k, backend, nThreads, true);
        default: break;
    }

    // Remaining fold operations (MIN/MAX/MEAN/MEDIAN/SUM/RMS/RMSE)
    if constexpr (!DSL::is_gpu_eligible_type_v<K>)
        return DST::Morphology::folding<T, K>(in, k, nThreads);
    else
    {
        if(backend.available())
            return detail::gpu_single<K, T>(backend, *in, k, "FOLD");
        else
            return DST::Morphology::folding<T, K>(in, k, nThreads);
    }
}

template<typename T>
std::shared_ptr<FITScube> folding(const std::shared_ptr<FITScube>& in,
                                   const DST::Morphology::kernel<T>& k,
                                   DST::GPU::Backend&                backend,
                                   size_t                            nThreads = 1)
{
    if(!in)
        throw FITSexception(BAD_DIMEN,"DST::GPU","folding","Input pointer is null.");
    switch(in->GetEqBitPerPixel())
    {
        case SBYTE_IMG:    return folding<T,int8_t>  (std::dynamic_pointer_cast<FITSimg<int8_t>  >(in), k, backend, nThreads);
        case BYTE_IMG:     return folding<T,uint8_t> (std::dynamic_pointer_cast<FITSimg<uint8_t> >(in), k, backend, nThreads);
        case SHORT_IMG:    return folding<T,int16_t> (std::dynamic_pointer_cast<FITSimg<int16_t> >(in), k, backend, nThreads);
        case USHORT_IMG:   return folding<T,uint16_t>(std::dynamic_pointer_cast<FITSimg<uint16_t>>(in), k, backend, nThreads);
        case LONG_IMG:     return folding<T,int32_t> (std::dynamic_pointer_cast<FITSimg<int32_t> >(in), k, backend, nThreads);
        case ULONG_IMG:    return folding<T,uint32_t>(std::dynamic_pointer_cast<FITSimg<uint32_t>>(in), k, backend, nThreads);
        case LONGLONG_IMG: return folding<T,int64_t> (std::dynamic_pointer_cast<FITSimg<int64_t> >(in), k, backend, nThreads);
        case ULONGLONG_IMG:return folding<T,uint64_t>(std::dynamic_pointer_cast<FITSimg<uint64_t>>(in), k, backend, nThreads);
        case FLOAT_IMG:    return folding<T,float>   (std::dynamic_pointer_cast<FITSimg<float>   >(in), k, backend, nThreads);
        case DOUBLE_IMG:   return folding<T,double>  (std::dynamic_pointer_cast<FITSimg<double>  >(in), k, backend, nThreads);
        default: break;
    }
    throw FITSexception(BAD_DIMEN,"DST::GPU","folding","Unsupported pixel type.");
}

}}} // namespace DST::Morphology::GPU

#undef DST_GPU_RETYPE_KERNEL

#endif // __DST_morph_gpu__
