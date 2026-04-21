# DeepSkyUtilities

A collection of astronomical image processing utilities for the DeepSky Project. Provides star detection, morphological operators, FITS I/O, and optional GPU-accelerated processing.

---

## Building

### Requirements

- CMake >= 3.16
- C++17 or later (C++23 preferred)
- [DSTFits](../deepskyfits), [DSTMath](../deepskymath), [DSTDisplay](../deepskydisplay) libraries
- cfitsio, Minuit2

The following system packages must be installed before running CMake:

**Debian/Ubuntu:**
```bash
sudo apt install \
  cmake make git curl wget \
  pkg-config build-essential \
  libgsl-dev \
  libblas-dev liblapack-dev libopenblas-dev liblapacke-dev \
  libfftw3-dev \
  autoconf automake libtool \
  astrometry.net \
  qt6-base-dev \
  doxygen graphviz \
  gnuplot-nox libcairo2 libpango-1.0-0 libpangocairo-1.0-0 \
  libboost-filesystem-dev libboost-chrono-dev libboost-thread-dev libboost-regex-dev
```

**macOS (Homebrew):**
```bash
brew install cmake git fftw openblas automake astrometry-net qt doxygen boost gnuplot
```

> **Note:** `nlohmann/json` and `scamp` are included as Git submodules and built automatically by CMake — no separate installation is needed.  
> Run `git submodule update --init --recursive` after cloning.

### Basic build (CPU only)

```bash
git submodule update --init --recursive
mkdir build && cd build
cmake ..
make -j$(nproc)
```

---

## Astrometry

The library provides `DST::astrometry::solver` — a C++ interface to two external plate-solving engines:

| Engine | What it does |
|---|---|
| **Astrometry.net** (`solve-field`) | Blind plate-solve: identifies a star field from the image itself against pre-built index files and writes WCS keywords back into the FITS header. |
| **SCAMP** | Precise WCS refinement: cross-matches a SExtractor-format source catalogue against an online reference catalogue (default: GAIA-DR3 via VizieR) and can model higher-order optical distortions. |

SCAMP is built automatically from the `scamp` submodule during `make`. `astrometry.net` must be installed system-wide (see Requirements above).

The solver is configured via a JSON file:

```json
{
  "astrometry.net": {
    "path": "/usr/bin",
    "parameters": {
      "scale-units": "arcsecperpix",
      "scale-low": 0.5,
      "scale-high": 1.2,
      "skip-solved": false
    }
  },
  "scamp": {
    "path": "",
    "parameters": {
      "astref_catalog": "GAIA-DR3",
      "distort_degrees": 3
    }
  }
}
```

**Astrometry.net (blind solve):**
```cpp
DST::astrometry::solver s("/path/to/config.json");
s.setAstrometryParam("scale-units", std::string("arcsecperpix"));
s.setAstrometryParam("scale-low",  0.5f);
s.setAstrometryParam("scale-high", 1.2f);
s.solve("image.fits");   // WCS keywords written back into image.fits
```

**SCAMP refinement from a SExtractor catalogue:**
```cpp
DST::astrometry::solver s("/path/to/config.json");
s.setScampParam("astref_catalog", std::string("GAIA-DR3"));
auto img = std::make_shared<FITScube>(...);
s.solve(img, "sources.ldac");   // runs SCAMP on sources.ldac
```

`starFinder::toFITS_LDAC()` exports a detected star catalogue directly to the FITS-LDAC format expected by SCAMP:
```cpp
DST::algorithms::starFinder finder(skyImage);
finder.findStars();
finder.toFITS_LDAC("sources.ldac");
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

### Astrometry tests

```bash
./build/exe/gTest_astrometry
```

Tests are organised in five suites:

| Suite | What is tested | Requires |
|---|---|---|
| 1 — Constructor / config | JSON parsing, invalid config rejection | nothing |
| 2 — Parameter API | `setAstrometryParam` / `setScampParam` type safety | nothing |
| 3 — FITS-LDAC export | `toFITS_LDAC()` output format | nothing |
| 4 — Live Astrometry.net solve | Full blind solve + WCS keywords | `solve-field` binary + index file |
| 5 — WCS accuracy | RA/Dec error < 1 arcsec vs. reference WCS | `solve-field` + M11 DSS2 reference FITS |

Suites 4 and 5 are **automatically skipped** when `solve-field` is not present. The M11 reference FITS and the astrometry.net index file (`index-5206-30.fits`) are downloaded at build time by CMake.

---

## License

Licensed under [CC BY-NC 4.0](https://creativecommons.org/licenses/by-nc/4.0/).
Centre de Physique des Particules de Marseille — GILLARD William.
