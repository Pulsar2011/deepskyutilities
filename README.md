# DeepSkyUtilities

A collection of astronomical image processing utilities for the DeepSky Project. Provides star detection, morphological operators, FITS I/O, and optional GPU-accelerated processing.

---

## Building

### Requirements

- CMake >= 3.16
- C++17 or later (C++23 preferred)
- [DSTFits](../deepskyfits), [DSTMath](../deepskymath), [DSTDisplay](../deepskydisplay) libraries
- cfitsio, Minuit2, matplot

### Basic build (CPU only)

```bash
mkdir build && cd build
cmake ..
make -j$(nproc)
```

---

## GPU Support

The library ships a multi-backend GPU dispatch layer for morphological operators used by `DST::algorithms::starFinder`. The backend is selected at CMake configure time. At runtime, dispatch is **automatic**: float and integer images are processed on GPU; double-precision images always fall back to CPU.

### Backend status

| Backend | CMake option | Validated | Notes |
|---|---|---|---|
| **CPU** (reference) | always on | **Yes** — full test suite | Multi-threaded C++ fallback |
| **Vulkan Compute** | `USE_VULKAN=ON` | **Yes** — validated against CPU | Tested on macOS (MoltenVK) and Linux |
| CUDA | *(not yet wired)* | No — lacking test facilities | Source in `gpu/morpho.cu` |
| HIP/ROCm | *(not yet wired)* | No — lacking test facilities | Source in `gpu/morpho.hip` |
| MPI | `USE_MPI=ON` | No — lacking test facilities | Domain-decomposition for distributed clusters |

### Enabling Vulkan Compute

**Prerequisites:** Install the [Vulkan SDK](https://vulkan.lunarg.com) (provides `glslc` or `glslangValidator` for shader compilation).

```bash
cmake .. -DUSE_VULKAN=ON
make -j$(nproc)
```

CMake will automatically compile the GLSL compute shader (`gpu/shaders/morpho.comp`) to SPIR-V at build time. If neither `glslc` nor `glslangValidator` is found, CMake will abort with an informative error.

**macOS:** Vulkan is provided via MoltenVK, included in the Vulkan SDK. No additional steps needed.

**Linux:** Install the Vulkan SDK or `vulkan-tools` + `glslang-tools` from your package manager:
```bash
# Debian/Ubuntu
sudo apt install vulkan-tools glslang-tools libvulkan-dev
```

### Enabling MPI (experimental, not yet validated)

```bash
cmake .. -DUSE_MPI=ON
make -j$(nproc)
```

Requires an MPI implementation (`OpenMPI`, `MPICH`, etc.) installed and findable by `find_package(MPI)`.

### Combining backends

```bash
cmake .. -DUSE_VULKAN=ON -DUSE_MPI=ON
```

### Runtime behaviour

No code changes are needed to benefit from GPU acceleration. `starFinder` automatically uses the GPU when:
- A real GPU backend is active (`isGPUEnabled()` returns `true`)
- The image pixel type is GPU-eligible (float, uint8, uint16, int16, int32, uint32)

Double-precision images always use the CPU path regardless of backend.

```cpp
// GPU dispatch is fully automatic — no API change required
DST::algorithms::starFinder finder(skyImage);
finder.findStars();   // uses GPU if available and image type is eligible

// Inspect which backend is active
if (finder.isGPUEnabled())
    std::cout << "Running on GPU\n";
```

---

## Running the tests

### CPU regression tests

```bash
./build/exe/gTest_findStars
```

### GPU tests (requires `USE_VULKAN=ON`)

```bash
./build/exe/gTest_findStars-GPU
```

The GPU test suite covers:
- Backend creation and eligibility checks (`GPUBackendTest`)
- Morphological correctness vs. CPU reference for float and uint16 images (`GPUMorphoTest`)
- Full `_findStarsMorpho` pipeline GPU vs. CPU (`StarFinderMorphoPipelineGPU`)
- End-to-end star detection — MORPHO, SNR, THRESHOLD methods (`StarFinderGPUTest`)
- Performance benchmarks on 512×512 and 2040×2040 images (`GPUPerfTest`)

---

## License

Licensed under [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/).
Centre de Physique des Particules de Marseille — GILLARD William.
