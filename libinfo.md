My Main Page                         {#mainpage}
============

# <img src="DeepSkyUtilitiesLogo.png" style="width:10%"> DeepSkyUtilities

A C++ library of astronomical image-processing utilities for the DeepSky Project, providing morphological operators, star detection algorithms, and an optional multi-backend GPU dispatch layer, all built on top of DeepSkyFits, DeepSkyMath, and DeepSkyDisplay.

- Repo: https://gitlab.in2p3.fr/deepskyproject/deepskyutilities.git
- C++ Standard: prefers C++23/20, falls back to C++17 (as configured by CMake)
- Platforms: Linux and macOS

--------------------------------------------------------------------------------

## Build and install

Full prerequisites and build options (CPU-only, Vulkan GPU, MPI) are described in the [README](README.md).

Quick start — CPU-only build:

```bash
mkdir -p build
cmake -S . -B build
cmake --build build -j
```

Enable Vulkan GPU acceleration:

```bash
cmake -S . -B build -DUSE_VULKAN=ON
cmake --build build -j
```

--------------------------------------------------------------------------------

## Library

DeepSkyUtilities exposes three tightly related layers: the **morphological kernel and operators** (`DST::Morphology`), the **GPU dispatch layer** that transparently accelerates those operators, and the **star-finder algorithm** (`DST::algorithms::starFinder`) that orchestrates both. A small **version-management singleton** (`DST::gGIT`) rounds out the public API. This section introduces each layer, its design choices, and typical usage patterns.

### DST::Morphology::kernel\<T\>: the convolution kernel

`kernel<T>` is the central building block of the morphological pipeline. It is a valarray-backed N-dimensional array of arithmetic type `T` paired with an operation type (`ktype`) that determines how it participates in a fold or correlation.

**Supported operation types (`ktype` enum):**

| Enum value      | Meaning |
|-----------------|---------|
| `MIN`           | Minimum reduction over kernel footprint |
| `MAX`           | Maximum reduction (= dilation kernel) |
| `MEAN`          | Mean reduction |
| `MEDIAN`        | Median reduction (sorted, stack-allocated ≤ 512 elements) |
| `SUM`           | Sum reduction |
| `RMS`           | Root-mean-square reduction |
| `RMSE`          | Root-mean-square error (needs reference) |
| `DILATE`        | Morphological dilation (alias for MAX) |
| `ERODE`         | Morphological erosion (alias for MIN) |
| `OPEN`          | Erosion followed by dilation (removes dark spots) |
| `CLOSE`         | Dilation followed by erosion (fills bright holes) |
| `WHITE_TOPHAT`  | Original − Opening (enhances bright features) |
| `BLACK_TOPHAT`  | Closing − Original (enhances dark features) |
| `CORRELATE`     | Weighted sum (linear filter) |
| `nCORRELATE`    | Normalised weighted sum |

`kernel<T>` supports compound assignment operators (`+=`, `-=`, `*=`, `/=`) with scalars, arrays, and other kernels, as well as a **type-promotion** mechanism (`promote_with<Other>`) so that output images may carry a higher-precision type than the input.

```c++
// 3×3 flat dilation kernel (all ones, DILATE type)
DST::Morphology::kernel<double> k({1,1,1, 1,1,1, 1,1,1}, {3,3}, DST::Morphology::DILATE);

// 5×5 Gaussian-shaped correlation kernel
DST::Morphology::kernel<float> gk(gValues, {5,5}, DST::Morphology::CORRELATE);

// Scale the kernel in-place
gk *= 2.0f;
```

--------------------------------------------------------------------------------

### Morphological operators

All operators share a common signature pattern:

- Accept a typed `FITSimg<K>` **or** a runtime-typed `FITScube`.
- Return an image whose pixel type may be promoted by the kernel's `promote_with<K>` rule.
- Accept an optional `nThreads` argument for multi-threaded CPU execution.
- Preserve and update FITS header metadata (HISTORY, BSCALE/BZERO, NAXISn).
- Throw `FITSexception` on CFITSIO or dimension errors.

#### dilation — morphological maximum

Replaces each pixel with the maximum of its kernel neighbourhood. Bright structures grow; small dark gaps fill. The result type is `kernel<K>::promote_with<K>`.

```c++
auto sky = fm.GetPrimary();                        // shared_ptr<FITScube>
DST::Morphology::kernel<double> k3x3({...}, {3,3}, DST::Morphology::DILATE);

auto dilated = DST::Morphology::dilation(*sky, k3x3, /*nThreads=*/4);
dilated->Write("dilated.fits", /*replace=*/true);
```

#### erosion — morphological minimum

Replaces each pixel with the minimum of its kernel neighbourhood. Bright structures shrink; thin bright filaments disappear.

```c++
auto eroded = DST::Morphology::erosion(*sky, k3x3, /*nThreads=*/4);
```

#### opening — erosion then dilation

Removes bright protrusions smaller than the kernel while preserving larger structures. Useful for background estimation.

```c++
auto opened = DST::Morphology::opening(*sky, k3x3, /*nThreads=*/4);
```

#### closing — dilation then erosion

Fills dark holes smaller than the kernel while preserving larger voids. Useful for filling gaps in detected structures.

```c++
auto closed = DST::Morphology::closing(*sky, k3x3, /*nThreads=*/4);
```

#### wTopHat — white top-hat transform

`wTopHat = original − opening`. Isolates bright features (stars, compact sources) that are smaller than the kernel against a smoothly varying background.

```c++
auto tophat = DST::Morphology::wTopHat(*sky, k3x3, /*nThreads=*/4);
```

#### bTopHat — black top-hat transform

`bTopHat = closing − original`. Isolates dark features (absorption features, dust lanes) smaller than the kernel.

```c++
auto btophat = DST::Morphology::bTopHat(*sky, k3x3, /*nThreads=*/4);
```

#### gradiant — morphological gradient

`gradiant = dilation − erosion`. Highlights edges and boundaries of structures. (Named `gradiant` in the API to preserve the original spelling.)

```c++
auto grad = DST::Morphology::gradiant(*sky, k3x3, /*nThreads=*/4);
```

#### folding — generic reduce

The primitive that underlies dilation and erosion: applies any `ktype` reduction over the kernel footprint. Useful when you need MIN, MAX, MEAN, MEDIAN, SUM, RMS, or RMSE with full control.

```c++
DST::Morphology::kernel<double> medKernel({...}, {5,5}, DST::Morphology::MEDIAN);
auto localMedian = DST::Morphology::folding(*sky, medKernel, /*nThreads=*/4);
```

#### correlate — linear filter

Computes a weighted sum of the neighbourhood, optionally normalised (`nCORRELATE`). Equivalent to discrete convolution when the kernel is not separable.

```c++
DST::Morphology::kernel<float> gaussian(gValues, {7,7}, DST::Morphology::CORRELATE);
auto smoothed = DST::Morphology::correlate(*sky, gaussian, /*nThreads=*/4);
```

--------------------------------------------------------------------------------

### GPU dispatch layer

The GPU layer sits transparently between the morphological API and the hardware. No API change is required: the same `dilation()`, `erosion()`, etc. calls are dispatched to GPU or CPU depending on the active backend and the image pixel type.

**Backend selection at a glance:**

| Backend | CMake flag | Validated | Notes |
|---|---|---|---|
| **CPU** (reference) | always on | Yes — full test suite | Multi-threaded C++; handles all types including `double` |
| **Vulkan Compute** | `USE_VULKAN=ON` | Yes — validated vs. CPU | Tested on macOS (MoltenVK) and Linux |
| CUDA | *(not yet wired)* | No | Source in `gpu/morpho.cu` |
| HIP/ROCm | *(not yet wired)* | No | Source in `gpu/morpho.hip` |
| MPI | `USE_MPI=ON` | No — experimental | Domain-decomposition for distributed clusters |

**Runtime dispatch rules:**

- GPU-eligible pixel types: `float`, `uint8_t`, `int8_t`, `uint16_t`, `int16_t`, `uint32_t`, `int32_t`.
- Non-eligible types (`double`, `int64_t`, `uint64_t`, `size_t`) always route to CPU.
- GPU path: single upload → one or more chained kernels → single download (minimises PCIe traffic).
- Composition operations (opening, closing, top-hat, gradient): intermediate results stay device-resident via `GPUContext<T>`.

The `GPUContext<T>` RAII object manages double-buffering and device memory lifetime; users never interact with it directly.

--------------------------------------------------------------------------------

### DST::algorithms::starFinder: star detection

`starFinder` wraps a `FITScube` sky image and exposes a configurable, GPU-aware pipeline for detecting and measuring point sources.

**Detection methods (`detectionMethod` enum):**

| Method | Description |
|--------|-------------|
| `MORPHO` | Morphological white top-hat + watershed flooding. Default. Best for crowded fields. |
| `THRESHOLD` | Simple pixel-value threshold followed by connected-component labelling. |
| `SNR` | Threshold expressed as a multiple of the local noise (signal-to-noise ratio). |
| `CUSTOM` | User-defined pre-processing via the virtual `_customPreprocessing()` hook. |

**Key configuration:**

| Setter | Effect |
|--------|--------|
| `setDetectionMethod(m)` | Switch detection method |
| `setThresholdValue(v)` | Threshold for THRESHOLD/SNR methods |
| `switchfBSkernel(k)` | Replace the background-subtraction kernel |
| `switchfBEkernel(k)` | Replace the erosion kernel |
| `switchfBMkernel(k)` | Replace the morphological-matching kernel |
| `setOutputUrl(path)` | Write intermediate products to disk |

**Detection and measurement flow:**

```c++
// 1. Open the sky image
FITSmanager fm("sky.fits");
auto sky = fm.GetPrimary();   // shared_ptr<FITScube>

// 2. Create the finder (GPU backend auto-selected at construction)
DST::algorithms::starFinder finder(sky);

// 3. Optionally choose a detection method and threshold
finder.setDetectionMethod(DST::algorithms::starFinder::detectionMethod::SNR);
finder.setThresholdValue(5.0);   // 5-sigma threshold

// 4. Run detection
finder.findStars();              // GPU dispatched when available

// 5. Run photometric measurement on all detected stars
finder.measureStars();

// 6. Export results as a FITS binary table
FITStable results = finder.detectedStars("STARS");
results.write("stars.fits", 0, /*replace=*/true);

// 7. Inspect GPU status
if (finder.isGPUEnabled())
    std::cout << "Detection ran on GPU\n";
```

**Customising the morphological kernels:**

```c++
// Replace the default 3×3 binary erosion kernel with a 5×5 disk
DST::Morphology::kernel<double> disk5({0,1,1,1,0,
                                        1,1,1,1,1,
                                        1,1,1,1,1,
                                        1,1,1,1,1,
                                        0,1,1,1,0}, {5,5}, DST::Morphology::ERODE);
finder.switchfBEkernel(disk5);
finder.findStars();
```

**Evaluating stellar profile parameters on a single star:**

```c++
// Fit a 2D Moffat profile to star #42
finder.EvalStellarProperties(42, DST::Stars::evalType::MOFFAT2D);

FITStable tbl = finder.detectedStars();
// Table columns include: ID, X, Y, FWHM_X, FWHM_Y, ANGLE, FLUX, SNR, BKG, ...
```

--------------------------------------------------------------------------------

### DST::Stars: stellar properties and profile fitting

`Stars` stores all measured properties of a single detected source and drives the profile-fitting engine.

**Measured properties:**

| Accessor | Type | Description |
|----------|------|-------------|
| `ID()` | `size_t` | Unique star identifier |
| `position()` | `point` | Centroid (x, y) in pixel coordinates |
| `maxima()` | `double` | Peak intensity |
| `maximaPosition()` | `point` | Position of peak pixel |
| `FWHM()` | `point` | Full width at half maximum (x, y) |
| `angle()` | `double` | Orientation angle of the PSF |
| `totalFlux()` | `double` | Integrated photometric flux |
| `SNR()` | `double` | Signal-to-noise ratio |
| `background()` | `double` | Local background level |
| `evalMethod()` | `evalType` | Profile model used |

**Profile evaluation methods (`evalType` enum):**

| Method | Description |
|--------|-------------|
| `MOMENT` | Intensity-weighted centroid and second moment (fast, no fitting) |
| `GAUSS` | Circular 1D Gaussian fit |
| `GAUSS2D` | Elliptical 2D Gaussian fit (provides FWHM_X, FWHM_Y, angle) |
| `MOFFAT` | Circular Moffat profile fit (better PSF model for optical systems) |
| `MOFFAT2D` | Elliptical 2D Moffat fit (most accurate; recommended for astrometry) |

```c++
// Iterate and print a catalogue
FITStable tbl = finder.detectedStars();
// or access individual stars via EvalStellarProperties
finder.EvalStellarProperties(1, DST::Stars::evalType::GAUSS2D);
```

`dump()` prints a formatted summary of all properties to any `std::ostream`:

```c++
// Print star #3 summary to stdout
Stars s = finder.detectedStarsMap().at(3);
s.dump(std::cout);
```

--------------------------------------------------------------------------------

### DST::gGIT: version management

`gGIT` is a lightweight singleton that exposes the library's git-tagged version, short version string, and commit SHA1. It enables runtime version checks and ordered comparisons, useful when multiple DeepSky libraries are loaded together.

```c++
// Print version information
std::cout << "DeepSkyUtilities version : " << DST::gGIT::version()       << "\n";
std::cout << "Short version            : " << DST::gGIT::version_short()  << "\n";
std::cout << "Commit SHA1              : " << DST::gGIT::sha1()           << "\n";

// Guard against known-incompatible version
if (DST::gGIT::isOlderThan("2.0.0"))
    throw std::runtime_error("DeepSkyUtilities >= 2.0.0 required");

// Exact-match check
bool ok = DST::gGIT::checkVersion("1.3.2");
```

Comparison operators (`<`, `<=`, `>`, `>=`, `==`, `!=`) operate on semantic version strings (major.minor.patch), so `isNewerThan`, `isOlderThan`, `isEqual`, and `isNotEqual` static helpers cover the most common guard patterns.

--------------------------------------------------------------------------------

DeepSkyUtilities Library

Authors:
- William GILLARD, Aix-Marseille Université, CNRS, CPPM, France

[<img src="by-nc.png" style="width:10%">](licence)
