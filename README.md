# DeepSkyUtilities

A collection of astronomical image processing utilities for the DeepSky Project. Provides star detection, morphological operators, FITS I/O, optional GPU-accelerated processing, and [FitsViewer](#fitsviewer), an interactive viewer for FITS images and tables.

---

## Building

### Requirements

- CMake >= 3.24
- C++17 or later (C++23 preferred)
- The DeepSky libraries below, installed under `DeepSkyROOT` (the environment variable of that name, or by default the parent of this source directory)
- cfitsio, Minuit2

| Library | Project | Minimum version |
|---|---|---|
| DSTCore | [deepskysys](../deepskysys) | 2.0.1 |
| DSTmath | [deepskymath](../deepskymath) | 0.1.3 |
| DSTfits | [deepskyfits](../deepskyfits) | 1.5.0 |
| DSTdisplay | [deepskydisplay](../deepskydisplay) | 1.2.4 |

CMake only looks for them in `DeepSkyROOT` and stops if one is missing or too old. These versions must be installed from a source tree that provides the CMake package (`lib/cmake/<Library>`).

Each installed library carries the version of its git tag, for example `libDSTools.2.2.0.dylib` with the links `libDSTools.2.dylib` and `libDSTools.dylib`. Installing a new version removes the older ones, so programs built against a previous major version must be rebuilt.

The following system packages must be installed before running CMake:

**Debian/Ubuntu:**
```bash
sudo apt install \
  cmake make git curl wget clang \
  pkg-config build-essential \
  libgsl-dev \
  libblas-dev liblapack-dev libopenblas-dev liblapacke-dev \
  libfftw3-dev \
  autoconf automake libtool \
  astrometry.net \
  qt6-base-dev \
  doxygen graphviz \
  gnuplot-nox libcairo2 libpango-1.0-0 libpangocairo-1.0-0
```

**macOS (Homebrew):**
```bash
brew install cmake git fftw openblas automake astrometry-net qt doxygen gnuplot
```

> **Note:** `nlohmann/json` and `scamp` are included as Git submodules and built automatically by CMake — no separate installation is needed. scamp is built from a copy of its submodule in `build/external`, so that its `autogen.sh` never writes into the submodule.  
> Run `git submodule update --init --recursive` after cloning.

### Basic build (CPU only)

Everything is built out of source: configuring in the source directory (`cmake .`) is refused.

```bash
git submodule update --init --recursive
cmake --preset default        # Release, in build/      (or: cmake -B build)
cmake --build build -j
cmake --build build --target install
```

Other presets: `debug` (Debug, in `build-debug/`) and `asan` (Debug + AddressSanitizer, in `build-asan/`):
`cmake --preset asan && cmake --build --preset asan -j`.

Options: `-DCMAKE_INSTALL_PREFIX=<path>` (default: `$DeepSkyROOT`, or the parent directory), `-DWITH_ASAN=ON`,
`-DBUILD_TESTS=OFF`, and the GPU backends below. `-DWITH_DEBUG=ON` is deprecated: it means Debug + AddressSanitizer.
Executables (FitsViewer, scamp, the tests) are in `build/bin`, the library in `build/lib`.

Former in-source builds left files in the sources (`CMakeCache.txt`, `Makefile`, `lib/`, the generated
`include/DSTools/DSTastrometry_config.h`, scamp's autotools files, ...); remove them once with
`cmake -P cmake/CleanInSourceBuild.cmake`.

---

## FitsViewer

`FitsViewer` opens FITS images, data cubes and tables for interactive inspection. Its windows need Qt6, and the install places it in `<prefix>/exe`.

If [asdf2fits](../asdf2fits) is installed under `DeepSkyROOT` when CMake runs, `FitsViewer` also opens ASDF files, recognised by their content. This dependency is optional: without it, the build goes on and `FitsViewer` reads FITS files only.

```bash
FitsViewer image.fits                 # list the HDUs of the file
FitsViewer --hdu 2 image.fits         # also open HDU 2 directly
FitsViewer -z 0 1000 image.fits       # start with this colour range
FitsViewer red.fits green.fits blue.fits   # go straight to colour composition
FitsViewer --nodraw --info image.fits # print statistics, no window
```

### Command line

| Option | Effect |
|---|---|
| `--hdu N` | Open HDU `N` as soon as the list is shown. |
| `--bins N`, `-b N` | Number of bins of the image histogram (default 100). |
| `--zscale MIN MAX`, `-z MIN MAX` | Initial colour range, in pixel values. |
| `--roi X0 NX Y0 NY` | Display only this region of a 2-D image, in pixels. |
| `--nodraw`, `--draw` | Open no window. Both disable drawing, on one file and one image HDU (`--hdu`, HDU 1 by default). |
| `--info`, `-i` | With `--nodraw`: print mean, median, quartiles, RMS, min, max and sum. |
| `--print` | With `--nodraw`: print the pixel values. |

The display scale is chosen in the viewer's toolbar: the former `--opt`/`++opt` option no longer exists.

### HDU list

The main window lists each HDU with its index, extension name and type. Double-click a line to open it. The viewer depends on the HDU:

| HDU | Viewer |
|---|---|
| Header only | FITS header table |
| 1-D or 2-D image | [Image viewer](#image-viewer) |
| 3-D image | Image viewer with a layer slider |
| ASCII or binary table | [Table viewer](#table-viewer) |

Opening another HDU replaces the current viewer and keeps its window position. The **Color ▸ Compose color image…** menu opens the [colour composite](#colour-composite).

### Image viewer

The window shows the image with its colour bar, the X and Y profiles through the cursor's pixel, and the pixel histogram with its cumulative distribution, percentiles and statistics. The status bar gives the pixel, its value and, when the image has a WCS, its world coordinates.

- **Navigate:** drag to pan and use the wheel to zoom. Dragging or zooming on the colour bar changes the colour range. Profiles and the status bar pause during a drag.
- **Header**, **WCS Grid** and **Save** (PNG, PDF or JPEG) are in the toolbar. The WCS grid checkbox appears when the image has a grid to show.
- **Scale:** Linear, Sqrt, Log10, Quadratic or Sinh. Changing the scale keeps the zoom and overlaid catalogues, and resets the colour range. Pixels a scale can't represent, such as the log of a negative value, take the lowest colour. Sinh is a normalised stretch, `sinh(3x)/sinh(3)` with `x` between the image's limits, so its colour bar runs from 0 to 1. The histogram and profiles always show raw pixel values.
- **Palette:** Night (default), Grayscale, Viridis, Hot, Cold, Thermal, Ion, Candy, Geography, Polar, Spectrum, Jet or Hues.
- **Layer slider (3-D only):** scale, palette, colour range and zoom are kept from one layer to the next.

#### Catalogue overlay

**Catalogue...** overlays the positions of a FITS table on a 2-D image or a colour composite. The dialog asks for:

- the file and table HDU;
- what the coordinates are: world coordinates, converted with the image's WCS, or pixels, either FITS (first pixel = 1) or array (first pixel = 0);
- the X and Y columns. With Y on *(from X column)*, a complex column gives real = X and imaginary = Y, and a 2-element vector column gives `[0]` = X and `[1]` = Y;
- an optional [row filter](#filters). The last filter used in the window is offered again.

World coordinates are not converted to the WCS's units: a warning appears if the columns' units differ. Positions outside the image are not drawn.

Tick **Keep previous** to add a catalogue over those already displayed, up to 100 at once. Each one gets its own colour and marker, and consecutive catalogues differ by both. Unticked, a new catalogue replaces the others. **Catalogues (n)** shows which marker belongs to which catalogue, **Clear catalogues** removes them all, and clicking a marker shows its table row.

On a colour composite, world coordinates use the WCS of the first channel.

### Table viewer

The left panel lists the columns. Double-click one to open a window with its properties, statistics and values. Right-click a cell to copy it, or a column header to copy its name.

Above the plot:

- **Filter:** restricts what is plotted, with the [syntax below](#filters).
- **New column:** `NAME = expression`, computed on every row and kept in memory until the table is reopened. The file is never changed. Operators are `+ - * / ^ ( )`, functions are `sqrt exp log log10 sin cos tan asin acos atan abs pow(x, y)`, and complex columns also take `real img abs arg`. Vector columns work element by element.

The plot tabs are:

| Tab | Plot |
|---|---|
| Scatter | Y against X or the row index |
| Density | Heatmap of two columns or of a complex one, optionally weighted by a third, with log axes and colour scale |
| Histogram | Distribution of a column, with log bins and log Y axis |
| Vector row | The elements of a vector column for one row, with Previous and Next |

- **Category:** splits Scatter and Histogram into one coloured series per category of a column. On float and double columns, each category spans value ± epsilon; an epsilon of 0 gives one category per distinct value.
- **Retain plot:** adds the next plot to the current one instead of replacing it. In the Density tab it draws contours.
- **Select region:** while on, dragging on a scatter plot or heatmap draws a box and lists the rows inside it.

### Filters

Table and catalogue filters share one syntax:

```text
(OBJNAME=~WD* || SCSID>0) && RA >= 10.5
```

- A condition reads `COLUMN operator value`, with `==`, `!=`, `<`, `<=`, `>`, `>=`. A single `=` means `==`.
- `&&` and `||` combine conditions, `&&` binding tighter, and parentheses group them.
- On string columns, `=~` and `!~` match a wildcard: `*` for any characters, `?` for exactly one.
- Logical columns take `==` and `!=`, with `T`, `F`, `true`, `false`, `1` or `0`.
- Quotes keep a value containing spaces or operators in one piece.
- A value must fit its column's type: `1.5` on an integer column is an error. Numeric columns also accept hexadecimal, such as `0xFF`.
- Vector and complex columns can't be filtered.

### Colour composite

Started from the **Color** menu, or by giving several files on the command line, with one file per channel by default.

- **Model:** RGB (3 channels) or CMYK (4 channels).
- **Source:** channels from any HDUs or files, or, in RGB, three layers of one 3-D cube. All channels must have the same pixel dimensions.

Each channel has a colour bar beside the image and a group in the **Channels** panel with its minimum, maximum and Linear or Logarithmic stretch. In RGB it also has a gradient: Default hue, Viridis or Grayscale. **Reset levels** restores each channel's full range, the status bar gives every channel's value under the cursor, and the [catalogue overlay](#catalogue-overlay) is available in the toolbar.

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
cmake --preset default -DUSE_VULKAN=ON
cmake --build build -j
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
cmake --preset default -DUSE_MPI=ON
cmake --build build -j
```

Requires an MPI implementation (`OpenMPI`, `MPICH`, etc.) installed and findable by `find_package(MPI)`.

### Combining backends

```bash
cmake --preset default -DUSE_VULKAN=ON -DUSE_MPI=ON
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

All the tests: `ctest --preset default` (or `asan`, `debug`). They run from the build directory, read their
data in `build/testdata` and write their files there.

### CPU regression tests

```bash
cd build && ./bin/gTest_findStars
```

### GPU tests (requires `USE_VULKAN=ON`)

```bash
cd build && ./bin/gTest_findStars-GPU
```

The GPU test suite covers:
- Backend creation and eligibility checks (`GPUBackendTest`)
- Morphological correctness vs. CPU reference for float and uint16 images (`GPUMorphoTest`)
- Full `_findStarsMorpho` pipeline GPU vs. CPU (`StarFinderMorphoPipelineGPU`)
- End-to-end star detection — MORPHO, SNR, THRESHOLD methods (`StarFinderGPUTest`)
- Performance benchmarks on 512×512 and 2040×2040 images (`GPUPerfTest`)

### Astrometry tests

```bash
cd build && ./bin/gTest_astrometry
```

Tests are organised in five suites:

| Suite | What is tested | Requires |
|---|---|---|
| 1 — Constructor / config | JSON parsing, invalid config rejection | nothing |
| 2 — Parameter API | `setAstrometryParam` / `setScampParam` type safety | nothing |
| 3 — FITS-LDAC export | `toFITS_LDAC()` output format | nothing |
| 4 — Live Astrometry.net solve | Full blind solve + WCS keywords | `solve-field` binary + index file |
| 5 — WCS accuracy | RA/Dec error < 1 arcsec vs. reference WCS | `solve-field` + M11 DSS2 reference FITS |

Suites 4 and 5 are **automatically skipped** when `solve-field` is not present. The M11 reference FITS and the astrometry.net index file (`index-5206-30.fits`) are downloaded at build time into `build/testdata` (a failed download is not retried until its `.stamp` file is deleted).

---

## License

Licensed under the [PolyForm Noncommercial License 1.0.0](https://polyformproject.org/licenses/noncommercial/1.0.0/): free for noncommercial use. Commercial use requires a separate licence from the author. See [LICENSE.md](LICENSE.md).
Copyright (c) 2026 William Gillard.
