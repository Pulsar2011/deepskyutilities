//
//  Backend.h
//
//  Abstract GPU/CPU backend interface for morphological operations.
//  Licensed under CC BY-NC 4.0
//
//  Usage:
//    auto b = Backend::create();      // factory — returns best available backend
//    DeviceBuffer d = b->allocate(n); // allocate device memory
//    b->upload(d, src, n);            // copy host → device
//    b->runKernel(...);               // execute morphological op
//    b->download(dst, d, n);         // copy device → host
//    b->free(d);                      // release device memory

#ifndef __DST_Backend__
#define __DST_Backend__

#include <cstddef>
#include <cstdint>
#include <memory>

#include <DSTools/DSTmorpho.h>   // ktype, DST::Morphology::ops::fold_op

namespace DST { namespace GPU {

/// Opaque handle to a device-side memory allocation.
/// `ptr` is the device address (type depends on backend: void* for CPU/CUDA,
/// id<MTLBuffer> for Metal, etc.).  `bytes` is the allocation size.
struct DeviceBuffer
{
    void*  ptr   = nullptr;
    size_t bytes = 0;
};

/// Pixel type tag passed to runKernel so it can select the right template instantiation.
enum class PixelType : int
{
    UInt8   = 0,
    Int8    = 1,
    UInt16  = 2,
    Int16   = 3,
    UInt32  = 4,
    Int32   = 5,
    Float32 = 6
};

/// Map a C++ type to the PixelType enum.
template<typename T> struct PixelTypeOf;
template<> struct PixelTypeOf<uint8_t>  { static constexpr PixelType value = PixelType::UInt8;   };
template<> struct PixelTypeOf<int8_t>   { static constexpr PixelType value = PixelType::Int8;    };
template<> struct PixelTypeOf<uint16_t> { static constexpr PixelType value = PixelType::UInt16;  };
template<> struct PixelTypeOf<int16_t>  { static constexpr PixelType value = PixelType::Int16;   };
template<> struct PixelTypeOf<uint32_t> { static constexpr PixelType value = PixelType::UInt32;  };
template<> struct PixelTypeOf<int32_t>  { static constexpr PixelType value = PixelType::Int32;   };
template<> struct PixelTypeOf<float>    { static constexpr PixelType value = PixelType::Float32; };

/// Abstract backend interface.
///
/// Concrete implementations: CPUBackend (always available),
/// CUDABackend (#ifdef USE_CUDA), HIPBackend (#ifdef USE_HIP),
/// VulkanBackend (#ifdef USE_VULKAN), MPIBackend (#ifdef USE_MPI).
class Backend
{
public:
    virtual ~Backend() = default;

    /// Returns true if the backend is actually available at runtime
    /// (e.g. CUDA device present, Vulkan driver loaded).
    virtual bool available() const = 0;

    // ------------------------------------------------------------------
    // Memory management
    // ------------------------------------------------------------------
    virtual DeviceBuffer allocate(size_t bytes) = 0;
    virtual void upload  (DeviceBuffer dst, const void* src, size_t bytes) = 0;
    virtual void download(void* dst, DeviceBuffer src, size_t bytes) = 0;
    virtual void free    (DeviceBuffer buf) = 0;

    // ------------------------------------------------------------------
    // Kernel execution
    // ------------------------------------------------------------------
    /// Run a morphological operation on device memory.
    ///
    /// \param op           Operation type (DILATE/ERODE handled via
    ///                     pixel_dilation/pixel_erosion; fold ops via pixel_fold;
    ///                     CORRELATE via pixel_correlate).
    /// \param in           Device buffer — input image (flat array of pixelType)
    /// \param mask         Device buffer — boolean mask (bool[nPixels]), may be {nullptr,0}
    /// \param out          Device buffer — output image (flat array of pixelType)
    /// \param kernelW      Device buffer — kernel weights (float, always float32 on GPU)
    /// \param kSizes       Device buffer — kernel dimension sizes (size_t[kDim])
    /// \param kCenter      Device buffer — kernel center offsets (ptrdiff_t[kDim])
    /// \param kDim         Number of kernel dimensions
    /// \param kSize        Total kernel elements
    /// \param imgSizes     Device buffer — image dimension sizes (size_t[imgDim])
    /// \param imgDim       Number of image dimensions
    /// \param nPixels      Total pixels in image
    /// \param pixelType    Pixel storage type (selects template instantiation)
    /// \param normalize    For CORRELATE: divide by sum(|w|)
    virtual void runKernel(Morphology::ktype  op,
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
                           bool               normalize = false) = 0;

    // ------------------------------------------------------------------
    // Factory
    // ------------------------------------------------------------------
    /// Create the best available backend for the current platform.
    /// Selection order: CUDA → HIP → Vulkan → CPU.
    static std::unique_ptr<Backend> create();
};

}} // namespace DST::GPU

#endif // __DST_Backend__
