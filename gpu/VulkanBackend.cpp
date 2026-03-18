//
//  VulkanBackend.cpp
//
//  Vulkan Compute backend implementation.
//  Licensed under CC BY-NC 4.0

#include "VulkanBackend.h"

#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdlib>   // getenv

// ============================================================
//  SPIR-V path resolution
// ============================================================
#ifndef DST_VULKAN_SPIRV_PATH
#  define DST_VULKAN_SPIRV_PATH "morpho.spv"
#endif

static std::vector<uint32_t> loadSPIRV(const std::string& path)
{
    std::ifstream f(path, std::ios::ate | std::ios::binary);
    if(!f.is_open())
        throw std::runtime_error("VulkanBackend: cannot open SPIR-V: " + path);
    const size_t bytes = static_cast<size_t>(f.tellg());
    if(bytes % 4 != 0)
        throw std::runtime_error("VulkanBackend: SPIR-V size not 4-byte aligned");
    std::vector<uint32_t> code(bytes / 4);
    f.seekg(0);
    f.read(reinterpret_cast<char*>(code.data()), static_cast<std::streamsize>(bytes));
    return code;
}

// ============================================================
//  Helper macros
// ============================================================
#define VK_CHECK(expr) \
    do { \
        VkResult _r = (expr); \
        if(_r != VK_SUCCESS) \
            throw std::runtime_error(std::string("Vulkan error at " #expr ": ") \
                                     + std::to_string(static_cast<int>(_r))); \
    } while(0)

namespace DST { namespace GPU {

// ============================================================
//  Constructor / Destructor
// ============================================================
VulkanBackend::VulkanBackend()
{
    try {
        init();
    } catch(...) {
        // Swallow — available_ stays false; Backend::create() will fall through
    }
}

VulkanBackend::~VulkanBackend()
{
    if(!device_) return;
    vkDeviceWaitIdle(device_);

    destroyStaging();

    if(pipeline_)   vkDestroyPipeline      (device_, pipeline_,   nullptr);
    if(pipeLayout_) vkDestroyPipelineLayout(device_, pipeLayout_, nullptr);
    if(dsLayout_)   vkDestroyDescriptorSetLayout(device_, dsLayout_, nullptr);
    if(cmdPool_)    vkDestroyCommandPool   (device_, cmdPool_,    nullptr);
    if(device_)     vkDestroyDevice        (device_,              nullptr);
    if(instance_)   vkDestroyInstance      (instance_,            nullptr);
}

// ============================================================
//  available()
// ============================================================
bool VulkanBackend::available() const { return available_; }

// ============================================================
//  Initialisation
// ============================================================
void VulkanBackend::init()
{
    createInstance();
    selectPhysicalDevice();
    createDevice();
    createCommandPool();
    loadPipeline();
    available_   = true;
    initialised_ = true;
}

void VulkanBackend::createInstance()
{
    VkApplicationInfo appInfo{};
    appInfo.sType            = VK_STRUCTURE_TYPE_APPLICATION_INFO;
    appInfo.pApplicationName = "DSTmorph";
    appInfo.apiVersion       = VK_API_VERSION_1_1;

    VkInstanceCreateInfo ci{};
    ci.sType            = VK_STRUCTURE_TYPE_INSTANCE_CREATE_INFO;
    ci.pApplicationInfo = &appInfo;
    VK_CHECK(vkCreateInstance(&ci, nullptr, &instance_));
}

void VulkanBackend::selectPhysicalDevice()
{
    uint32_t count = 0;
    vkEnumeratePhysicalDevices(instance_, &count, nullptr);
    if(count == 0)
        throw std::runtime_error("VulkanBackend: no Vulkan physical devices found");

    std::vector<VkPhysicalDevice> devs(count);
    vkEnumeratePhysicalDevices(instance_, &count, devs.data());

    // Prefer discrete GPU; fall back to any device with a compute queue
    VkPhysicalDevice fallback = VK_NULL_HANDLE;
    for(auto& d : devs)
    {
        VkPhysicalDeviceProperties props;
        vkGetPhysicalDeviceProperties(d, &props);

        uint32_t qCount = 0;
        vkGetPhysicalDeviceQueueFamilyProperties(d, &qCount, nullptr);
        std::vector<VkQueueFamilyProperties> qProps(qCount);
        vkGetPhysicalDeviceQueueFamilyProperties(d, &qCount, qProps.data());

        for(uint32_t qi = 0; qi < qCount; ++qi)
        {
            if(qProps[qi].queueFlags & VK_QUEUE_COMPUTE_BIT)
            {
                if(!fallback) { fallback = d; queueFamily_ = qi; }
                if(props.deviceType == VK_PHYSICAL_DEVICE_TYPE_DISCRETE_GPU)
                {
                    physDev_      = d;
                    queueFamily_  = qi;
                }
                break;
            }
        }
    }

    if(!physDev_) physDev_ = fallback;
    if(!physDev_)
        throw std::runtime_error("VulkanBackend: no compute-capable Vulkan device found");
}

void VulkanBackend::createDevice()
{
    const float prio = 1.0f;
    VkDeviceQueueCreateInfo qci{};
    qci.sType            = VK_STRUCTURE_TYPE_DEVICE_QUEUE_CREATE_INFO;
    qci.queueFamilyIndex = queueFamily_;
    qci.queueCount       = 1;
    qci.pQueuePriorities = &prio;

    VkDeviceCreateInfo dci{};
    dci.sType            = VK_STRUCTURE_TYPE_DEVICE_CREATE_INFO;
    dci.queueCreateInfoCount = 1;
    dci.pQueueCreateInfos    = &qci;
    VK_CHECK(vkCreateDevice(physDev_, &dci, nullptr, &device_));
    vkGetDeviceQueue(device_, queueFamily_, 0, &queue_);
}

void VulkanBackend::createCommandPool()
{
    VkCommandPoolCreateInfo ci{};
    ci.sType            = VK_STRUCTURE_TYPE_COMMAND_POOL_CREATE_INFO;
    ci.queueFamilyIndex = queueFamily_;
    ci.flags            = VK_COMMAND_POOL_CREATE_RESET_COMMAND_BUFFER_BIT;
    VK_CHECK(vkCreateCommandPool(device_, &ci, nullptr, &cmdPool_));
}

void VulkanBackend::loadPipeline()
{
    // ------------------------------------------------------------------
    // Descriptor set layout: 7 storage buffer bindings
    // ------------------------------------------------------------------
    const uint32_t N_BINDINGS = 7;
    VkDescriptorSetLayoutBinding bindings[N_BINDINGS] = {};
    for(uint32_t i = 0; i < N_BINDINGS; ++i)
    {
        bindings[i].binding         = i;
        bindings[i].descriptorType  = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
        bindings[i].descriptorCount = 1;
        bindings[i].stageFlags      = VK_SHADER_STAGE_COMPUTE_BIT;
    }
    VkDescriptorSetLayoutCreateInfo dsci{};
    dsci.sType        = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_LAYOUT_CREATE_INFO;
    dsci.bindingCount = N_BINDINGS;
    dsci.pBindings    = bindings;
    VK_CHECK(vkCreateDescriptorSetLayout(device_, &dsci, nullptr, &dsLayout_));

    // ------------------------------------------------------------------
    // Pipeline layout (descriptor set + push constants)
    // ------------------------------------------------------------------
    VkPushConstantRange pcRange{};
    pcRange.stageFlags = VK_SHADER_STAGE_COMPUTE_BIT;
    pcRange.offset     = 0;
    pcRange.size       = sizeof(PushConstants);

    VkPipelineLayoutCreateInfo plci{};
    plci.sType                  = VK_STRUCTURE_TYPE_PIPELINE_LAYOUT_CREATE_INFO;
    plci.setLayoutCount         = 1;
    plci.pSetLayouts            = &dsLayout_;
    plci.pushConstantRangeCount = 1;
    plci.pPushConstantRanges    = &pcRange;
    VK_CHECK(vkCreatePipelineLayout(device_, &plci, nullptr, &pipeLayout_));

    // ------------------------------------------------------------------
    // Load SPIR-V
    // ------------------------------------------------------------------
    const char* envPath = std::getenv("DST_VULKAN_SPIRV_PATH");
    std::string spirvPath = envPath ? envPath : DST_VULKAN_SPIRV_PATH;
    std::vector<uint32_t> code = loadSPIRV(spirvPath);

    VkShaderModuleCreateInfo smci{};
    smci.sType    = VK_STRUCTURE_TYPE_SHADER_MODULE_CREATE_INFO;
    smci.codeSize = code.size() * sizeof(uint32_t);
    smci.pCode    = code.data();
    VkShaderModule shaderModule;
    VK_CHECK(vkCreateShaderModule(device_, &smci, nullptr, &shaderModule));

    // ------------------------------------------------------------------
    // Compute pipeline
    // ------------------------------------------------------------------
    VkComputePipelineCreateInfo cpci{};
    cpci.sType        = VK_STRUCTURE_TYPE_COMPUTE_PIPELINE_CREATE_INFO;
    cpci.layout       = pipeLayout_;
    cpci.stage.sType  = VK_STRUCTURE_TYPE_PIPELINE_SHADER_STAGE_CREATE_INFO;
    cpci.stage.stage  = VK_SHADER_STAGE_COMPUTE_BIT;
    cpci.stage.module = shaderModule;
    cpci.stage.pName  = "main";
    VK_CHECK(vkCreateComputePipelines(device_, VK_NULL_HANDLE, 1, &cpci, nullptr, &pipeline_));
    vkDestroyShaderModule(device_, shaderModule, nullptr);
}

// ============================================================
//  Memory helpers
// ============================================================
uint32_t VulkanBackend::findMemType(uint32_t typeBits, VkMemoryPropertyFlags flags)
{
    VkPhysicalDeviceMemoryProperties props;
    vkGetPhysicalDeviceMemoryProperties(physDev_, &props);
    for(uint32_t i = 0; i < props.memoryTypeCount; ++i)
        if((typeBits & (1u << i)) &&
           (props.memoryTypes[i].propertyFlags & flags) == flags)
            return i;
    throw std::runtime_error("VulkanBackend: no suitable memory type");
}

VkBuffer VulkanBackend::createBuffer(size_t bytes, VkBufferUsageFlags usage,
                                      VkMemoryPropertyFlags props,
                                      VkDeviceMemory& outMem)
{
    VkBufferCreateInfo bci{};
    bci.sType       = VK_STRUCTURE_TYPE_BUFFER_CREATE_INFO;
    bci.size        = bytes;
    bci.usage       = usage;
    bci.sharingMode = VK_SHARING_MODE_EXCLUSIVE;
    VkBuffer buf;
    VK_CHECK(vkCreateBuffer(device_, &bci, nullptr, &buf));

    VkMemoryRequirements mr;
    vkGetBufferMemoryRequirements(device_, buf, &mr);

    VkMemoryAllocateInfo ai{};
    ai.sType           = VK_STRUCTURE_TYPE_MEMORY_ALLOCATE_INFO;
    ai.allocationSize  = mr.size;
    ai.memoryTypeIndex = findMemType(mr.memoryTypeBits, props);
    VK_CHECK(vkAllocateMemory(device_, &ai, nullptr, &outMem));
    vkBindBufferMemory(device_, buf, outMem, 0);
    return buf;
}

void VulkanBackend::ensureStaging(size_t bytes)
{
    if(staging_.size >= bytes) return;
    destroyStaging();

    VkDeviceMemory mem;
    staging_.buf  = createBuffer(bytes,
                                  VK_BUFFER_USAGE_TRANSFER_SRC_BIT |
                                  VK_BUFFER_USAGE_TRANSFER_DST_BIT,
                                  VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT |
                                  VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
                                  mem);
    staging_.mem  = mem;
    staging_.size = bytes;
    vkMapMemory(device_, staging_.mem, 0, bytes, 0, &staging_.mapped);
}

void VulkanBackend::destroyStaging()
{
    if(!staging_.buf) return;
    vkUnmapMemory(device_, staging_.mem);
    vkDestroyBuffer(device_, staging_.buf, nullptr);
    vkFreeMemory   (device_, staging_.mem, nullptr);
    staging_ = {};
}

// ============================================================
//  One-shot command buffer helpers
// ============================================================
VkCommandBuffer VulkanBackend::beginOneShot()
{
    VkCommandBufferAllocateInfo ai{};
    ai.sType              = VK_STRUCTURE_TYPE_COMMAND_BUFFER_ALLOCATE_INFO;
    ai.commandPool        = cmdPool_;
    ai.level              = VK_COMMAND_BUFFER_LEVEL_PRIMARY;
    ai.commandBufferCount = 1;
    VkCommandBuffer cmd;
    vkAllocateCommandBuffers(device_, &ai, &cmd);

    VkCommandBufferBeginInfo bi{};
    bi.sType = VK_STRUCTURE_TYPE_COMMAND_BUFFER_BEGIN_INFO;
    bi.flags = VK_COMMAND_BUFFER_USAGE_ONE_TIME_SUBMIT_BIT;
    vkBeginCommandBuffer(cmd, &bi);
    return cmd;
}

void VulkanBackend::endOneShot(VkCommandBuffer cmd)
{
    vkEndCommandBuffer(cmd);

    VkSubmitInfo si{};
    si.sType              = VK_STRUCTURE_TYPE_SUBMIT_INFO;
    si.commandBufferCount = 1;
    si.pCommandBuffers    = &cmd;

    VkFenceCreateInfo fci{};
    fci.sType = VK_STRUCTURE_TYPE_FENCE_CREATE_INFO;
    VkFence fence;
    vkCreateFence(device_, &fci, nullptr, &fence);

    vkQueueSubmit(queue_, 1, &si, fence);
    vkWaitForFences(device_, 1, &fence, VK_TRUE, UINT64_MAX);

    vkDestroyFence(device_, fence, nullptr);
    vkFreeCommandBuffers(device_, cmdPool_, 1, &cmd);
}

void VulkanBackend::fillBuffer(VkBuffer buf, VkDeviceSize offset,
                                VkDeviceSize size, uint32_t data)
{
    VkCommandBuffer cmd = beginOneShot();
    vkCmdFillBuffer(cmd, buf, offset, size, data);
    endOneShot(cmd);
}

// ============================================================
//  Memory management (public)
//  DeviceBuffer::ptr holds a heap-allocated VkBuffer/VkDeviceMemory pair.
// ============================================================
struct VkBufPair { VkBuffer buf; VkDeviceMemory mem; };

DeviceBuffer VulkanBackend::allocate(size_t bytes)
{
    if(!available_)
        throw std::runtime_error("VulkanBackend::allocate: backend not available");

    auto* pair = new VkBufPair();
    pair->buf  = createBuffer(bytes,
                               VK_BUFFER_USAGE_STORAGE_BUFFER_BIT |
                               VK_BUFFER_USAGE_TRANSFER_SRC_BIT   |
                               VK_BUFFER_USAGE_TRANSFER_DST_BIT,
                               VK_MEMORY_PROPERTY_DEVICE_LOCAL_BIT,
                               pair->mem);
    return DeviceBuffer{ pair, bytes };
}

void VulkanBackend::upload(DeviceBuffer dst, const void* src, size_t bytes)
{
    ensureStaging(bytes);
    std::memcpy(staging_.mapped, src, bytes);

    auto* pair  = static_cast<VkBufPair*>(dst.ptr);
    VkCommandBuffer cmd = beginOneShot();
    VkBufferCopy region{ 0, 0, bytes };
    vkCmdCopyBuffer(cmd, staging_.buf, pair->buf, 1, &region);
    endOneShot(cmd);
}

void VulkanBackend::download(void* dst, DeviceBuffer src, size_t bytes)
{
    ensureStaging(bytes);

    auto* pair = static_cast<VkBufPair*>(src.ptr);
    VkCommandBuffer cmd = beginOneShot();
    VkBufferCopy region{ 0, 0, bytes };
    vkCmdCopyBuffer(cmd, pair->buf, staging_.buf, 1, &region);
    endOneShot(cmd);

    std::memcpy(dst, staging_.mapped, bytes);
}

void VulkanBackend::free(DeviceBuffer buf)
{
    if(!buf.ptr) return;
    auto* pair = static_cast<VkBufPair*>(buf.ptr);
    vkDestroyBuffer(device_, pair->buf, nullptr);
    vkFreeMemory   (device_, pair->mem, nullptr);
    delete pair;
}

// ============================================================
//  runKernel
// ============================================================
void VulkanBackend::runKernel(Morphology::ktype  op,
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
    if(!available_)
        throw std::runtime_error("VulkanBackend::runKernel: backend not available");

    // ------------------------------------------------------------------
    // 1. Zero-initialise output buffer (required for atomicOr pixel writes)
    // ------------------------------------------------------------------
    auto* outPair = static_cast<VkBufPair*>(out.ptr);
    fillBuffer(outPair->buf, 0, out.bytes, 0u);

    // ------------------------------------------------------------------
    // 2. Create descriptor pool for this dispatch
    // ------------------------------------------------------------------
    const uint32_t N_BINDINGS = 7;
    VkDescriptorPoolSize poolSize{};
    poolSize.type            = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
    poolSize.descriptorCount = N_BINDINGS;

    VkDescriptorPoolCreateInfo dpci{};
    dpci.sType         = VK_STRUCTURE_TYPE_DESCRIPTOR_POOL_CREATE_INFO;
    dpci.maxSets       = 1;
    dpci.poolSizeCount = 1;
    dpci.pPoolSizes    = &poolSize;
    VkDescriptorPool dpool;
    VK_CHECK(vkCreateDescriptorPool(device_, &dpci, nullptr, &dpool));

    // ------------------------------------------------------------------
    // 3. Allocate and update descriptor set
    // ------------------------------------------------------------------
    VkDescriptorSetAllocateInfo dsai{};
    dsai.sType              = VK_STRUCTURE_TYPE_DESCRIPTOR_SET_ALLOCATE_INFO;
    dsai.descriptorPool     = dpool;
    dsai.descriptorSetCount = 1;
    dsai.pSetLayouts        = &dsLayout_;
    VkDescriptorSet ds;
    VK_CHECK(vkAllocateDescriptorSets(device_, &dsai, &ds));

    // Provide a dummy 4-byte buffer for the mask when there is none
    DeviceBuffer dummyMask{ nullptr, 0 };
    VkBufPair* dummyPair = nullptr;
    const bool hasMask = (mask.ptr != nullptr && mask.bytes > 0);
    if(!hasMask)
    {
        dummyMask = allocate(4);
        dummyPair = static_cast<VkBufPair*>(dummyMask.ptr);
    }

    auto buf_of = [&](DeviceBuffer db) -> VkBuffer {
        return static_cast<VkBufPair*>(db.ptr)->buf;
    };

    VkDescriptorBufferInfo bi[N_BINDINGS];
    bi[0] = { buf_of(in),                      0, VK_WHOLE_SIZE };
    bi[1] = { hasMask ? buf_of(mask)
                       : dummyPair->buf,         0, VK_WHOLE_SIZE };
    bi[2] = { buf_of(out),                     0, VK_WHOLE_SIZE };
    bi[3] = { buf_of(kernelW),                  0, VK_WHOLE_SIZE };
    bi[4] = { buf_of(kSizes),                   0, VK_WHOLE_SIZE };
    bi[5] = { buf_of(kCenter),                  0, VK_WHOLE_SIZE };
    bi[6] = { buf_of(imgSizes),                 0, VK_WHOLE_SIZE };

    VkWriteDescriptorSet writes[N_BINDINGS];
    for(uint32_t i = 0; i < N_BINDINGS; ++i)
    {
        writes[i]                  = {};
        writes[i].sType            = VK_STRUCTURE_TYPE_WRITE_DESCRIPTOR_SET;
        writes[i].dstSet           = ds;
        writes[i].dstBinding       = i;
        writes[i].descriptorCount  = 1;
        writes[i].descriptorType   = VK_DESCRIPTOR_TYPE_STORAGE_BUFFER;
        writes[i].pBufferInfo      = &bi[i];
    }
    vkUpdateDescriptorSets(device_, N_BINDINGS, writes, 0, nullptr);

    // ------------------------------------------------------------------
    // 4. Record and submit command buffer
    // ------------------------------------------------------------------
    PushConstants pc{};
    pc.nPixels   = static_cast<uint32_t>(nPixels);
    pc.kDim      = static_cast<uint32_t>(kDim);
    pc.kSize     = static_cast<uint32_t>(kSize);
    pc.imgDim    = static_cast<uint32_t>(imgDim);
    pc.op        = static_cast<int32_t>(op);
    pc.pixelType = static_cast<int32_t>(pixelType);
    pc.normalize = normalize ? 1 : 0;
    pc.hasMask   = hasMask   ? 1 : 0;

    const uint32_t groups = static_cast<uint32_t>((nPixels + 255u) / 256u);

    VkCommandBuffer cmd = beginOneShot();
    vkCmdBindPipeline(cmd, VK_PIPELINE_BIND_POINT_COMPUTE, pipeline_);
    vkCmdBindDescriptorSets(cmd, VK_PIPELINE_BIND_POINT_COMPUTE,
                             pipeLayout_, 0, 1, &ds, 0, nullptr);
    vkCmdPushConstants(cmd, pipeLayout_,
                        VK_SHADER_STAGE_COMPUTE_BIT,
                        0, sizeof(PushConstants), &pc);
    vkCmdDispatch(cmd, groups, 1, 1);
    endOneShot(cmd);

    // ------------------------------------------------------------------
    // 5. Cleanup
    // ------------------------------------------------------------------
    if(!hasMask && dummyPair)
        free(dummyMask);

    vkDestroyDescriptorPool(device_, dpool, nullptr);
}

}} // namespace DST::GPU
