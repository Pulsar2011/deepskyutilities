//
//  GPUContext.h
//
//  Device-residency context for morphological operations.
//
//  Keeps image data on the device across multiple chained operations,
//  avoiding expensive PCIe round-trips for every intermediate result.
//
//  Usage:
//    auto backend = Backend::create();
//    GPUContext<float> ctx(*backend, *fits_image);  // upload once
//    ctx.run(k_erode);                              // stays on device
//    ctx.run(k_dilate);                             // stays on device
//    ctx.download(*output_fits);                    // download once
//
//  Internals:
//    • Two device buffers (d_a_, d_b_) are double-buffered: each run()
//      reads from the "input" buffer and writes to the "output" buffer,
//      then swaps them. No host memory is touched between construction
//      and download().
//    • The boolean mask is uploaded once in the constructor and reused
//      across all run() calls.
//    • The kernel weights are uploaded as float32 before each run().
//
//  Template parameter T: image pixel type. Must satisfy is_gpu_eligible_type_v<T>.
//  Licensed under CC BY-NC 4.0

#ifndef __DST_GPUContext__
#define __DST_GPUContext__

#include "Backend.h"

#include <DSTfits/FITSdata.h>   // is_gpu_eligible_type_v
#include <DSTfits/FITSimg.h>    // FITScube

#include <vector>
#include <cstddef>
#include <stdexcept>

namespace DST { namespace GPU {

template<typename T>
class GPUContext
{
    static_assert(DSL::is_gpu_eligible_type_v<T>,
        "GPUContext: pixel type T is not GPU-eligible. "
        "Use the CPU path for double/int64/uint64/size_t.");

public:
    /// Upload the image (and its mask) to the device.
    /// \param backend  Backend to use for all subsequent operations.
    /// \param img      Source image cube (must contain type T data).
    explicit GPUContext(Backend& backend, const FITScube& img)
        : backend_(backend)
    {
        const T* host_data = img.raw_data<T>();
        if(!host_data)
            throw std::invalid_argument("GPUContext: image has no data of the requested type");

        nPixels_ = img.Nelements();
        imgDim_  = img.GetDimension();

        imgSizes_.resize(imgDim_);
        for(size_t d = 0; d < imgDim_; ++d)
            imgSizes_[d] = img.Size(d + 1);

        const size_t dataBytes = nPixels_ * sizeof(T);
        const size_t maskBytes = nPixels_ * sizeof(bool);

        // Allocate double buffers
        d_a_ = backend_.allocate(dataBytes);
        d_b_ = backend_.allocate(dataBytes);

        // Upload image into d_a_ (current "input" buffer)
        backend_.upload(d_a_, host_data, dataBytes);

        // Upload mask (may be empty)
        const bool* host_mask = img.raw_mask();
        if(host_mask && img.GetMask().size() == nPixels_)
        {
            d_mask_ = backend_.allocate(maskBytes);
            backend_.upload(d_mask_, host_mask, maskBytes);
        }
        else
        {
            d_mask_ = {nullptr, 0};  // no mask
        }

        a_is_input_ = true;
    }

    ~GPUContext()
    {
        backend_.free(d_a_);
        backend_.free(d_b_);
        if(d_mask_.ptr) backend_.free(d_mask_);
    }

    // Non-copyable, non-movable (owns device resources)
    GPUContext(const GPUContext&)            = delete;
    GPUContext& operator=(const GPUContext&) = delete;

    /// Run a morphological operation on the device-resident data.
    /// The kernel weights are uploaded as float32 for the GPU (they are
    /// always converted to float32 before crossing the host/device boundary).
    template<typename K>
    void run(const Morphology::kernel<K>& k)
    {
        const size_t kSize = k.size();
        const size_t kDim  = k.sizes().size();

        // Build kernel center offsets
        std::vector<ptrdiff_t> kCenter(kDim);
        for(size_t a = 0; a < kDim; ++a)
            kCenter[a] = static_cast<ptrdiff_t>(k.sizes()[a] / 2);

        // Upload kernel weights as float32
        std::vector<float> kw_f(kSize);
        for(size_t i = 0; i < kSize; ++i)
            kw_f[i] = static_cast<float>(k.data()[i]);

        const size_t kSizeBytes   = kDim * sizeof(size_t);
        const size_t kCenterBytes = kDim * sizeof(ptrdiff_t);
        const size_t kWBytes      = kSize * sizeof(float);
        const size_t iSizeBytes   = imgDim_ * sizeof(size_t);

        DeviceBuffer d_kw      = backend_.allocate(kWBytes);
        DeviceBuffer d_ks      = backend_.allocate(kSizeBytes);
        DeviceBuffer d_kc      = backend_.allocate(kCenterBytes);
        DeviceBuffer d_is      = backend_.allocate(iSizeBytes);

        backend_.upload(d_kw, kw_f.data(),        kWBytes);
        backend_.upload(d_ks, k.sizes().data(),   kSizeBytes);
        backend_.upload(d_kc, kCenter.data(),     kCenterBytes);
        backend_.upload(d_is, imgSizes_.data(),   iSizeBytes);

        DeviceBuffer& d_in  = a_is_input_ ? d_a_ : d_b_;
        DeviceBuffer& d_out = a_is_input_ ? d_b_ : d_a_;

        backend_.runKernel(
            k.type(),
            d_in, d_mask_, d_out,
            d_kw, d_ks, d_kc,
            kDim, kSize,
            d_is, imgDim_,
            nPixels_,
            PixelTypeOf<T>::value,
            /*normalize=*/false);

        // Swap input/output buffers
        a_is_input_ = !a_is_input_;

        backend_.free(d_kw);
        backend_.free(d_ks);
        backend_.free(d_kc);
        backend_.free(d_is);
    }

    /// Download the current device-resident result into a host FITScube.
    /// The output cube must already be allocated with the correct shape and type T.
    void download(FITScube& out) const
    {
        T* host_out = out.raw_data<T>();
        if(!host_out)
            throw std::invalid_argument("GPUContext::download: output image has wrong type or no data");

        const DeviceBuffer& d_in = a_is_input_ ? d_a_ : d_b_;
        backend_.download(host_out, d_in, nPixels_ * sizeof(T));
    }

private:
    Backend&             backend_;
    DeviceBuffer         d_a_;          // double-buffer A
    DeviceBuffer         d_b_;          // double-buffer B
    DeviceBuffer         d_mask_;       // mask (upload-once)
    std::vector<size_t>  imgSizes_;     // image dimension sizes (host copy)
    size_t               imgDim_    = 0;
    size_t               nPixels_   = 0;
    bool                 a_is_input_ = true;
};

}} // namespace DST::GPU

#endif // __DST_GPUContext__
