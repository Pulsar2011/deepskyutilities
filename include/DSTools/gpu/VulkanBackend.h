//
//  VulkanBackend.h
//
//  Vulkan Compute implementation of the Backend interface.
//  Include only from translation units that also link against Vulkan.
//
//  Architecture:
//    • One compute pipeline loaded from morpho.spv (compiled from
//      shaders/morpho.comp by glslangValidator at build time).
//    • Push constants encode the operation type, pixel type, and
//      all kernel/image geometry needed by the GLSL shader.
//    • 8 descriptor bindings (in, mask, out, kw, kSizes, kCenter,
//      imgSizes, and a hasMask flag encoded in PC) per dispatch.
//    • Output buffer is zero-initialised before each kernel launch
//      (required for the atomic-OR sub-word writes in the shader).
//    • Compute queue, command pool, and pipeline are created lazily
//      on the first call to runKernel().
//
//  Pixel data storage on device:
//    Device buffers store raw bytes identical to the host layout
//    (same as CUDA/HIP).  The shader unpacks bytes from uint32 words.
//    Sub-32-bit writes use atomicOr; output must be zero-initialised.
//
//  SPIR-V path:
//    At build time set DST_VULKAN_SPIRV_PATH via CMake or compile with
//        -DDST_VULKAN_SPIRV_PATH=\"/path/to/morpho.spv\"
//    Default: looks for morpho.spv in the same directory as the executable.
//
//  Licensed under PolyForm Noncommercial 1.0.0

#ifndef __DST_VulkanBackend__
#define __DST_VulkanBackend__

#include "Backend.h"

#include <vulkan/vulkan.h>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace DST { namespace GPU {

class VulkanBackend : public Backend
{
public:
    VulkanBackend();
    ~VulkanBackend() override;

    // Non-copyable
    VulkanBackend(const VulkanBackend&)            = delete;
    VulkanBackend& operator=(const VulkanBackend&) = delete;

    bool available() const override;

    // ------------------------------------------------------------------
    // Memory management — device-local buffers via VkBuffer/VkDeviceMemory
    // ------------------------------------------------------------------
    DeviceBuffer allocate(size_t bytes) override;
    void upload  (DeviceBuffer dst, const void* src, size_t bytes) override;
    void download(void* dst, DeviceBuffer src, size_t bytes) override;
    void free    (DeviceBuffer buf) override;

    // ------------------------------------------------------------------
    // Kernel execution
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
    // Vulkan object handles — all owned by this instance
    VkInstance               instance_       = VK_NULL_HANDLE;
    VkPhysicalDevice         physDev_        = VK_NULL_HANDLE;
    VkDevice                 device_         = VK_NULL_HANDLE;
    VkQueue                  queue_          = VK_NULL_HANDLE;
    uint32_t                 queueFamily_    = 0;

    VkCommandPool            cmdPool_        = VK_NULL_HANDLE;
    VkDescriptorSetLayout    dsLayout_       = VK_NULL_HANDLE;
    VkPipelineLayout         pipeLayout_     = VK_NULL_HANDLE;
    VkPipeline               pipeline_       = VK_NULL_HANDLE;

    bool                     initialised_    = false;
    bool                     available_      = false;

    // Staging buffer (host-visible, used for upload/download)
    struct StagingBuf {
        VkBuffer       buf  = VK_NULL_HANDLE;
        VkDeviceMemory mem  = VK_NULL_HANDLE;
        size_t         size = 0;
        void*          mapped = nullptr;
    } staging_;

    // Private helpers
    void            init();
    void            createInstance();
    void            selectPhysicalDevice();
    void            createDevice();
    void            createCommandPool();
    void            loadPipeline();
    void            ensureStaging(size_t bytes);
    void            destroyStaging();

    VkBuffer        createBuffer(size_t bytes, VkBufferUsageFlags usage,
                                 VkMemoryPropertyFlags props,
                                 VkDeviceMemory& outMem);
    uint32_t        findMemType(uint32_t typeBits, VkMemoryPropertyFlags flags);

    VkCommandBuffer beginOneShot();
    void            endOneShot(VkCommandBuffer cmd);
    void            fillBuffer(VkBuffer buf, VkDeviceSize offset,
                               VkDeviceSize size, uint32_t data);

    // Push constant layout must match the GLSL push_constant block exactly
    struct PushConstants
    {
        uint32_t nPixels;
        uint32_t kDim;
        uint32_t kSize;
        uint32_t imgDim;
        int32_t  op;
        int32_t  pixelType;
        int32_t  normalize;
        int32_t  hasMask;
    };
};

}} // namespace DST::GPU

#endif // __DST_VulkanBackend__
