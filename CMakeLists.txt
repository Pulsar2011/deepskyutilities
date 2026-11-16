#-------------------- CMAKE Configuration
cmake_minimum_required(VERSION 3.24)

# -Wno-dev by default (developer and deprecation warnings, mostly from the submodules): no need to
# pass it, nor to alias cmake with it, which breaks 'cmake --build' and 'cmake --install'.
# 'cmake -Wdev ...' shows them again.
if(NOT DEFINED CACHE{CMAKE_SUPPRESS_DEVELOPER_WARNINGS})
    set(CMAKE_SUPPRESS_DEVELOPER_WARNINGS ON CACHE INTERNAL "Suppress developer warnings (-Wno-dev)")
endif()

# Everything is generated in a build directory, never in the sources: from the project directory,
# 'cmake --preset default' (or 'cmake -B build'); from a build directory, 'cmake <path to project>'.
get_filename_component(_src "${CMAKE_SOURCE_DIR}" REALPATH)
get_filename_component(_bin "${CMAKE_BINARY_DIR}" REALPATH)
if(_src STREQUAL _bin)
    message(FATAL_ERROR
        "In-source builds are not allowed. Remove what this run (and any former in-source build) created:\n"
        "    cmake -P cmake/CleanInSourceBuild.cmake\n"
        "then, from the project directory:\n"
        "    cmake --preset default        (or: cmake -B build)\n"
        "    cmake --build build -j\n"
        "    cmake --build build --target install\n"
        "Other presets: debug (build-debug/), asan (build-asan/).\n"
        "Or from a build directory: mkdir build && cd build && cmake .. && make")
endif()

project(DeepSkyUtilities LANGUAGES C CXX)

if(NOT WIN32)
    string(ASCII 27 Esc)
    set(ColourReset "${Esc}[m")
    set(ColourBold  "${Esc}[1m")
    set(Red         "${Esc}[31m")
    set(Green       "${Esc}[32m")
    set(Yellow      "${Esc}[33m")
    set(Blue        "${Esc}[34m")
    set(Magenta     "${Esc}[35m")
    set(Cyan        "${Esc}[36m")
    set(White       "${Esc}[37m")
    set(BoldRed     "${Esc}[1;31m")
    set(BoldGreen   "${Esc}[1;32m")
    set(BoldYellow  "${Esc}[1;33m")
    set(BoldBlue    "${Esc}[1;34m")
    set(BoldMagenta "${Esc}[1;35m")
    set(BoldCyan    "${Esc}[1;36m")
    set(BoldWhite   "${Esc}[1;37m")
endif()

message("${BoldCyan}----------- DeepSkyUtilities (DSTools) -----------${ColourReset}")
message("${Green}-- Project source directory :${ColourReset} ${CMAKE_SOURCE_DIR}")
message("${Green}-- Project build directory  :${ColourReset} ${CMAKE_BINARY_DIR}")

if(DEFINED ENV{DeepSkyROOT})
    set(DeepSkyROOT $ENV{DeepSkyROOT})
else()
    set(DeepSkyROOT ${CMAKE_SOURCE_DIR}/..)
endif()
get_filename_component(DeepSkyROOT "${DeepSkyROOT}" ABSOLUTE)

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${DeepSkyROOT} CACHE PATH "default install path" FORCE)
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/modules")

#-------------------- Options
option(WITH_ASAN   "Build DSTools and its tests with AddressSanitizer" OFF)
option(BUILD_TESTS "Build the GoogleTest suite" ON)
option(WITH_DEBUG  "Deprecated: same as CMAKE_BUILD_TYPE=Debug and WITH_ASAN=ON (use the asan preset)" OFF)

if(WITH_DEBUG)
    message(WARNING "WITH_DEBUG is deprecated: use 'cmake --preset asan' (Debug + AddressSanitizer, in build-asan/)")
    set(CMAKE_BUILD_TYPE Debug)
    set(WITH_ASAN ON)
endif()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Products of this build: executables in bin/, libraries in lib/ (objects stay in CMakeFiles/)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# The highest standard the compiler supports, from C++23 down
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED OFF)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
# compile_commands.json in the build directory, for clangd / IntelliSense
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

#-------------------- COMPILER
# These flags apply to every C++ target of this directory, the googletest submodule included
# (not to the C of scamp, built apart).
set(_release "$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:Release>>")
set(_debug   "$<AND:$<COMPILE_LANGUAGE:CXX>,$<CONFIG:Debug>>")

if(CMAKE_SYSTEM_NAME MATCHES "Linux")
    set(OperatingSystem "Linux")
    add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-Wall;-pipe;-pthread;-fexceptions;-Wno-unknown-pragmas;-Wmisleading-indentation;-mtune=native;-ffp-contract=on>")
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|amd64")
        add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-m64>")
        add_compile_definitions(Linuxx86_64)
    elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm64")
        add_compile_definitions(Linuxaarch64)
    endif()
    add_compile_options("$<${_release}:-flto>")
    add_link_options(-pthread "$<$<CONFIG:Release>:-flto>")
    # Interprocedural optimization (LTO) of the release build, when supported
    include(CheckIPOSupported)
    check_ipo_supported(RESULT ipo_supported OUTPUT ipo_error)
    if(ipo_supported)
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION_RELEASE ON)
    endif()
elseif(CMAKE_SYSTEM_NAME MATCHES "Darwin")
    set(OperatingSystem "Mac OS X")
    add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-Wall;-mcpu=native;-mtune=native;-ffp-contract=on>")
    add_compile_definitions(Darwin)
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "arm64|aarch64")
        add_compile_definitions(Darwinarm64)
        set(CMAKE_OSX_ARCHITECTURES "arm64")
    elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|amd64")
        add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-m64>")
        add_compile_definitions(Darwinx86_64)
        set(CMAKE_OSX_ARCHITECTURES "x86_64")
    endif()
    add_compile_options("$<${_release}:-flto=thin>")
    add_link_options("$<$<CONFIG:Release>:-flto=thin>")
endif()

add_compile_options("$<${_release}:-O3;-fomit-frame-pointer>")
# No LTO in debug: on macOS it replaces the .o files with LLVM bitcode, leaving the debug map with
# empty OSO entries that LLDB cannot follow.
add_compile_options("$<${_debug}:-O0;-g>")
add_compile_definitions("$<$<CONFIG:Debug>:__DEBUG__>")

if(WITH_ASAN)
    add_compile_options("$<$<COMPILE_LANGUAGE:CXX>:-fsanitize=address;-fno-omit-frame-pointer;-g>")
    add_link_options(-fsanitize=address)
endif()

message("${BoldYellow}-- CMAKE INSTALL PREFIX : ${ColourReset}${CMAKE_INSTALL_PREFIX}")
message("${Green}-- DeepSkyROOT :${ColourReset} ${DeepSkyROOT}")
message("${Green}-- BUILD TYPE :${ColourReset} ${CMAKE_BUILD_TYPE}")
message("${Green}-- CXX COMPILER :${ColourReset} ${CMAKE_CXX_COMPILER}")
message("${Green}-- USE CXX STANDARD :${ColourReset} ${CMAKE_CXX_STANDARD}")
if(WITH_ASAN)
    message("${Blue}-- AddressSanitizer ACTIVATED${ColourReset}")
endif()

#-------------------- VERSION FROM THE GIT TAG
find_package(Git)
include(GetGitRevisionDescription)
git_describe(VERSION --always --tags --long)

string(REGEX REPLACE "^v([0-9]+)\\..*" "\\1" VERSION_MAJOR "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.([0-9]+).*" "\\1" VERSION_MINOR "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" VERSION_PATCH "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.[0-9]+.(-*)" "\\1" VERSION_SHA1 "${VERSION}")
set(VERSION_SHORT "v${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
message("-- PKG VERSION : ${VERSION} (${VERSION_SHORT}, sha1 ${VERSION_SHA1})")

set(CMAKE_PROJECT_VERSION       "${VERSION}")
set(CMAKE_PROJECT_VERSION_MAJOR "${VERSION_MAJOR}")
set(CMAKE_PROJECT_VERSION_MINOR "${VERSION_MINOR}")
set(CMAKE_PROJECT_VERSION_PATCH "${VERSION_PATCH}")

# A clone without a reachable tag (e.g. a shallow CI clone) has no numeric version: the library
# stays unversioned and no CMake package is installed.
if(VERSION MATCHES "^v[0-9]+\\.[0-9]+\\.[0-9]+")
    set(DSTools_NUMERIC_VERSION "${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
    set(PROJECT_VERSION "${DSTools_NUMERIC_VERSION}")
else()
    message(WARNING "No version tag found (git describe gave '${VERSION}'): DSTools will not be versioned")
endif()

# Generated sources stay in the build directory
set(GENERATED_DIR ${CMAKE_BINARY_DIR}/generated)
configure_file(${CMAKE_SOURCE_DIR}/modules/DST_version.cxx.in ${GENERATED_DIR}/src/DST_version.cxx)
configure_file(${CMAKE_SOURCE_DIR}/modules/gTest_version.h.in ${GENERATED_DIR}/test/gTest_version.h)

#-------------------- DeepSky libraries (installed in DeepSkyROOT)
message("${Red}----------- SEARCH DeepSky libraries -----------${ColourReset}")
# Minimum versions of the DeepSky libraries: the latest tag of each package.
set(DSTCore_MIN_VERSION    2.0.1)
set(DSTmath_MIN_VERSION    0.1.3)
set(DSTfits_MIN_VERSION    1.5.0)
set(DSTdisplay_MIN_VERSION 1.2.4)

# Only look in DeepSkyROOT: the RPATH points there, and CMake would otherwise fall back on another
# install (e.g. found through PATH) when the one in DeepSkyROOT is too old.
foreach(_dst_pkg DSTCore DSTmath DSTfits DSTdisplay)
    # A <pkg>_DIR cached by an earlier configure is used as is, whatever PATHS and NO_DEFAULT_PATH say:
    # forget it, so that a copy installed elsewhere (e.g. in /usr) cannot come back from the cache.
    unset(${_dst_pkg}_DIR CACHE)
    find_package(${_dst_pkg} ${${_dst_pkg}_MIN_VERSION} REQUIRED PATHS ${DeepSkyROOT} NO_DEFAULT_PATH)
    # Imported targets are system targets by default, which turns ${DeepSkyROOT}/include into -isystem:
    # it then loses its -I priority over other installed copies and hides warnings from our own headers.
    # Needs CMake >= 3.25; older versions ignore the property and keep -isystem.
    set_target_properties(${_dst_pkg}::${_dst_pkg} PROPERTIES SYSTEM OFF)
    message("-- ${Green}Found${ColourReset} ${_dst_pkg} ${${_dst_pkg}_VERSION}: ${${_dst_pkg}_DIR}")
endforeach()

# The DSTfits headers include cfitsio and Minuit2 headers, installed beside them by deepskyfits
find_library(MINUIT2_LIBRARY Minuit2 PATHS ${DeepSkyROOT}/lib NO_DEFAULT_PATH REQUIRED)
foreach(_dir cfitsio Minuit2)
    if(NOT EXISTS ${DeepSkyROOT}/include/${_dir})
        message(FATAL_ERROR "${DeepSkyROOT}/include/${_dir} not found: install deepskyfits first")
    endif()
endforeach()
message("-- ${Green}Found${ColourReset} Minuit2: ${MINUIT2_LIBRARY}")
add_library(deepsky_externals INTERFACE)
target_include_directories(deepsky_externals INTERFACE ${DeepSkyROOT}/include ${DeepSkyROOT}/include/cfitsio ${DeepSkyROOT}/include/Minuit2)
target_link_libraries(deepsky_externals INTERFACE ${MINUIT2_LIBRARY})

# Optional: without asdf2fits, FitsViewer is built for FITS files only. A package file left behind
# by a removed install must not break the link, so the files themselves are checked too. The exact
# directory is given because a search from DeepSkyROOT also looks into the asdf2fits source folder,
# where an in-source build leaves a package file with a wrong prefix.
unset(asdf2fits_DIR CACHE)
# 0.2 is the first version telling where an HDU is in the tree (ASDFmanager::GetPathAtIndex).
find_package(asdf2fits 0.2 QUIET PATHS ${DeepSkyROOT}/lib/cmake/asdf2fits NO_DEFAULT_PATH)
set(ASDF2FITS_USABLE OFF)
if(asdf2fits_FOUND
   AND EXISTS ${DeepSkyROOT}/lib/libasdf2fits${CMAKE_SHARED_LIBRARY_SUFFIX}
   AND EXISTS ${DeepSkyROOT}/lib/libasdf${CMAKE_SHARED_LIBRARY_SUFFIX}
   AND EXISTS ${DeepSkyROOT}/include/asdf2fits/ASDFmanager.h)
    set(ASDF2FITS_USABLE ON)
    message("-- ${Green}Found${ColourReset} asdf2fits ${asdf2fits_VERSION}: FitsViewer will read ASDF files")
else()
    message("${Yellow}-- asdf2fits >= 0.2 not found: FitsViewer will read FITS files only${ColourReset}")
endif()

#-------------------- nlohmann/json (header-only submodule)
if(NOT EXISTS ${CMAKE_SOURCE_DIR}/json/include/nlohmann/json.hpp)
    message(FATAL_ERROR "${BoldRed}json is missing: run 'git submodule update --init'${ColourReset}")
endif()
message("${Green}-- nlohmann/json headers:${ColourReset} ${CMAKE_SOURCE_DIR}/json/include")

#-------------------- Qt6 (optional: without it, the GUI is disabled)
find_package(Qt6 COMPONENTS Widgets Core Gui PrintSupport Concurrent)
if(Qt6_FOUND)
    message("${Blue}-- Qt6 FOUND : ${ColourReset} ${Qt6_VERSION}")
    set(QT_LIBS Qt6::Widgets Qt6::Core Qt6::Gui Qt6::PrintSupport Qt6::Concurrent)
    add_compile_definitions(__WITH_QT__)
    # Qt's code generators
    set(CMAKE_AUTOMOC ON)
    set(CMAKE_AUTOUIC ON)
    set(CMAKE_AUTORCC ON)
else()
    message("${Red}-- Qt6 NOT FOUND. GUI will be disabled${ColourReset}")
    set(QT_LIBS "")
endif()

#-------------------- RPATH
# This build comes first, then DeepSkyROOT: an installed copy must not replace what the tests are testing.
set(CMAKE_BUILD_RPATH "${CMAKE_LIBRARY_OUTPUT_DIRECTORY};${DeepSkyROOT}/lib")
# Installed: <prefix>/exe and <prefix>/lib find <prefix>/lib, then DeepSkyROOT
if(APPLE)
    set(CMAKE_INSTALL_RPATH "@loader_path/../lib;${DeepSkyROOT}/lib")
else()
    set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib;${DeepSkyROOT}/lib")
endif()
# Also the directories of the external libraries (Qt from Homebrew or conda, ...)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
if(CMAKE_SYSTEM_NAME MATCHES "Linux")
    # RPATH rather than RUNPATH: also used to find the dependencies of the DeepSky libraries
    add_link_options(-Wl,--disable-new-dtags)
endif()

#-------------------- astrometry.net solve-field and scamp defaults
# The detected path is written into DSTastrometry_config.h (generated header, installed alongside
# DSTastrometry.h) so downstream consumers get the same default. DSTastrometry.h itself is never
# generated by CMake.
find_program(ASTROMETRY_NET_SOLVE_FIELD solve-field
    HINTS /opt/homebrew/bin
          /opt/homebrew/Cellar/astrometry-net/*/bin
          /usr/local/bin
          /usr/bin
)
if(ASTROMETRY_NET_SOLVE_FIELD)
    get_filename_component(ASTROMETRY_NET_DEFAULT_PATH "${ASTROMETRY_NET_SOLVE_FIELD}" DIRECTORY)
    message("${Green}-- astrometry.net solve-field:${ColourReset} ${ASTROMETRY_NET_SOLVE_FIELD}")
else()
    # Leave ASTROMETRY_NET_DEFAULT_PATH unset so #cmakedefine comments it out
    message("${Yellow}-- astrometry.net solve-field not found; ASTROMETRY_NET_DEFAULT_PATH will be undefined${ColourReset}")
endif()

# scamp install path (matches the DESTINATION of its install() below)
get_filename_component(SCAMP_DEFAULT_PATH "${CMAKE_INSTALL_PREFIX}/exe" ABSOLUTE)

configure_file(
    ${CMAKE_SOURCE_DIR}/modules/DSTastrometry_config.h.in
    ${GENERATED_DIR}/include/DSTools/DSTastrometry_config.h
)

#-------------------- scamp (Autotools, submodule)
message("${Red}----------- scamp (submodule) -----------${ColourReset}")
if(NOT EXISTS ${CMAKE_SOURCE_DIR}/scamp/configure.ac)
    message(FATAL_ERROR "${BoldRed}scamp is missing: run 'git submodule update --init'${ColourReset}")
endif()
include(ExternalProject)
set(EXTERNAL_DIR ${CMAKE_BINARY_DIR}/external)
include(ProcessorCount)
ProcessorCount(NCORES)
if(NCORES EQUAL 0)
    set(NCORES 4)
endif()
# The autotools build needs make, whatever the generator (CMAKE_MAKE_PROGRAM is ninja with Ninja)
find_program(MAKE_EXECUTABLE NAMES gmake make REQUIRED)

find_program(AUTOMAKE automake
    HINTS /opt/homebrew/opt/automake/bin   # macOS Apple Silicon (keg-only)
          /usr/local/opt/automake/bin       # macOS Intel (keg-only)
)
if(NOT AUTOMAKE)
    message(FATAL_ERROR "automake not found. Install it (e.g. brew install automake) to build scamp from the submodule.")
endif()
get_filename_component(AUTOMAKE_BIN_DIR ${AUTOMAKE} DIRECTORY)
message("${Green}-- automake found:${ColourReset} ${AUTOMAKE}")

# Locate FFTW (required by scamp) — may be keg-only on macOS
find_path(FFTW_INCLUDE_DIR fftw3.h
    HINTS /opt/homebrew/opt/fftw/include   # macOS Apple Silicon (keg-only)
          /usr/local/opt/fftw/include       # macOS Intel (keg-only)
)
find_library(FFTW_LIBRARY NAMES fftw3
    HINTS /opt/homebrew/opt/fftw/lib
          /usr/local/opt/fftw/lib
)
if(NOT FFTW_INCLUDE_DIR OR NOT FFTW_LIBRARY)
    message(FATAL_ERROR "FFTW not found. Install it (e.g. brew install fftw) to build scamp from the submodule.")
endif()
get_filename_component(FFTW_LIB_DIR ${FFTW_LIBRARY} DIRECTORY)
message("${Green}-- FFTW include:${ColourReset} ${FFTW_INCLUDE_DIR}")
message("${Green}-- FFTW lib:    ${ColourReset} ${FFTW_LIB_DIR}")

# Locate OpenBLAS (used instead of ATLAS for linear algebra) — may be keg-only on macOS.
# NO_DEFAULT_PATH prevents find_path from falling back to system SDK paths (e.g. vecLib),
# which have cblas.h but lack a working lapacke.h.
# Linux multiarch paths (CMAKE_LIBRARY_ARCHITECTURE is e.g. aarch64-linux-gnu or x86_64-linux-gnu)
# are listed explicitly so NO_DEFAULT_PATH does not block them on Debian/Ubuntu CI runners.
find_path(OPENBLAS_INCLUDE_DIR cblas.h
    HINTS $ENV{CONDA_PREFIX}/include
          $ENV{CONDA_PREFIX}/include/openblas
          /opt/homebrew/opt/openblas/include          # macOS Apple Silicon (keg-only)
          /usr/local/opt/openblas/include              # macOS Intel (keg-only)
          /usr/include/${CMAKE_LIBRARY_ARCHITECTURE}   # Debian/Ubuntu multiarch (arm64, x86_64 …)
          /usr/include/openblas
          /usr/include
    NO_DEFAULT_PATH
)
find_library(OPENBLAS_LIBRARY NAMES openblas
    HINTS $ENV{CONDA_PREFIX}/lib
          $ENV{CONDA_PREFIX}/openblas/lib
          /opt/homebrew/opt/openblas/lib
          /usr/local/opt/openblas/lib
          /usr/lib/${CMAKE_LIBRARY_ARCHITECTURE}       # Debian/Ubuntu multiarch
          /usr/lib
    NO_DEFAULT_PATH
)
if(NOT OPENBLAS_INCLUDE_DIR OR NOT OPENBLAS_LIBRARY)
    message(FATAL_ERROR "OpenBLAS not found. Install it (e.g. brew install openblas) to build scamp from the submodule.")
endif()
get_filename_component(OPENBLAS_LIB_DIR ${OPENBLAS_LIBRARY} DIRECTORY)
message("${Green}-- OpenBLAS include:${ColourReset} ${OPENBLAS_INCLUDE_DIR}")
message("${Green}-- OpenBLAS lib:    ${ColourReset} ${OPENBLAS_LIB_DIR}")

# On Debian/Ubuntu, OpenBLAS is built WITHOUT LAPACKE (NO_LAPACKE=1).
# The LAPACKE_dpotrf symbol lives in a separate liblapacke.so, and lapacke.h
# is in /usr/include/ while cblas.h is in /usr/include/<arch>/  (multiarch).
# scamp's acx_openblas.m4 assumes both headers and the library symbol live
# together under the paths given by --with-openblas-{incdir,libdir}.
#
# To satisfy scamp's configure on Debian:
#   1) Create a thin GNU ld linker script "libopenblas.so" that groups
#      the real libopenblas.so and liblapacke.so, so AC_SEARCH_LIBS finds
#      LAPACKE_dpotrf in "-lopenblas".
#   2) Create a wrapper include directory with symlinks to cblas.h (from
#      openblas) and lapacke.h (from lapacke-dev), so both headers are
#      found at the same --with-openblas-incdir path.
#
# On macOS (Homebrew) OpenBLAS already bundles LAPACKE, so this is skipped.
if(UNIX AND NOT APPLE)
    find_library(_LAPACKE_LIB NAMES lapacke
      HINTS $ENV{CONDA_PREFIX}/lib
            /usr/lib/${CMAKE_LIBRARY_ARCHITECTURE} /usr/lib
    )
    find_path(_LAPACKE_INCLUDE_DIR lapacke.h
      HINTS $ENV{CONDA_PREFIX}/include
            $ENV{CONDA_PREFIX}/include/openblas
            /usr/include/${CMAKE_LIBRARY_ARCHITECTURE}
            /usr/include
    )
    if(_LAPACKE_LIB)
        set(_OB_WRAPPER_DIR "${EXTERNAL_DIR}/openblas-wrapper")

        # -- Library wrapper (ld linker script) --
        file(MAKE_DIRECTORY "${_OB_WRAPPER_DIR}/lib")
        file(WRITE "${_OB_WRAPPER_DIR}/lib/libopenblas.so"
            "INPUT(${OPENBLAS_LIBRARY} ${_LAPACKE_LIB})\n")
        set(OPENBLAS_LIB_DIR "${_OB_WRAPPER_DIR}/lib")
        message("${Green}-- Created OpenBLAS+LAPACKE lib wrapper: ${_OB_WRAPPER_DIR}/lib/libopenblas.so${ColourReset}")

        # -- Include wrapper (symlinks) --
        # Only needed when lapacke.h and cblas.h are in different directories.
        if(_LAPACKE_INCLUDE_DIR AND NOT "${_LAPACKE_INCLUDE_DIR}" STREQUAL "${OPENBLAS_INCLUDE_DIR}")
            file(MAKE_DIRECTORY "${_OB_WRAPPER_DIR}/include")
            execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
                "${OPENBLAS_INCLUDE_DIR}/cblas.h" "${_OB_WRAPPER_DIR}/include/cblas.h")
            # lapacke_mangling.h and lapacke_config.h are pulled in by lapacke.h, when present
            foreach(_h lapacke.h lapacke_mangling.h lapacke_config.h)
                if(EXISTS "${_LAPACKE_INCLUDE_DIR}/${_h}")
                    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
                        "${_LAPACKE_INCLUDE_DIR}/${_h}" "${_OB_WRAPPER_DIR}/include/${_h}")
                endif()
            endforeach()
            set(OPENBLAS_INCLUDE_DIR "${_OB_WRAPPER_DIR}/include")
            message("${Green}-- Created OpenBLAS+LAPACKE include wrapper: ${_OB_WRAPPER_DIR}/include/${ColourReset}")
        endif()
    endif()
endif()

# autogen.sh writes into its source tree: work on a copy, so that the submodule stays clean.
# PATH gets the automake directory, so that autogen.sh finds automake/aclocal even when the package
# is keg-only and not linked into /opt/homebrew/bin. MAKEFLAGS is emptied: the jobserver of the
# calling make must not reach this independent build.
set(SCAMP_COPY    ${EXTERNAL_DIR}/scamp-src)
set(SCAMP_INSTALL ${EXTERNAL_DIR}/scamp-install)
set(SCAMP_PROGRAM ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/scamp)
ExternalProject_Add(scamp_external
    PREFIX            ${EXTERNAL_DIR}/scamp
    DOWNLOAD_COMMAND  ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/scamp ${SCAMP_COPY}
    COMMAND           ${CMAKE_COMMAND} -E rm -f ${SCAMP_COPY}/.git
    SOURCE_DIR        ${SCAMP_COPY}
    BUILD_IN_SOURCE   1
    CONFIGURE_COMMAND ${CMAKE_COMMAND} -E env PATH=${AUTOMAKE_BIN_DIR}:$ENV{PATH} sh ./autogen.sh
    COMMAND           ${CMAKE_COMMAND} -E env
                      CC=${CMAKE_C_COMPILER}
                      CFLAGS=-fPIC
                      CPPFLAGS=-I${OPENBLAS_INCLUDE_DIR}
                      ./configure
                      --prefix=${SCAMP_INSTALL}
                      --with-fftw-incdir=${FFTW_INCLUDE_DIR}
                      --with-fftw-libdir=${FFTW_LIB_DIR}
                      --enable-openblas
                      --with-openblas-incdir=${OPENBLAS_INCLUDE_DIR}
                      --with-openblas-libdir=${OPENBLAS_LIB_DIR}
    BUILD_COMMAND     ${CMAKE_COMMAND} -E env MAKEFLAGS= ${MAKE_EXECUTABLE} -j ${NCORES} -C src
    INSTALL_COMMAND   ${CMAKE_COMMAND} -E env MAKEFLAGS= ${MAKE_EXECUTABLE} -C src install
    COMMAND           ${CMAKE_COMMAND} -E copy_if_different ${SCAMP_INSTALL}/bin/scamp ${SCAMP_PROGRAM}
    BUILD_BYPRODUCTS  ${SCAMP_PROGRAM}
)
message("${Green}-- scamp will be built from a copy of the submodule:${ColourReset} ${SCAMP_PROGRAM}")

# Helper target: run 'scamp --help' once to observe its options and output format.
# Usage: cmake --build build --target scamp_probe
add_custom_target(scamp_probe
    COMMAND ${SCAMP_PROGRAM} --help || true
    DEPENDS scamp_external
    COMMENT "Running scamp --help so you can observe its output format"
    VERBATIM
)

#-------------------- GPU backends (optional)
# Backend.cpp is always compiled (contains the Backend::create() factory)
set(GPU_SRC "${CMAKE_SOURCE_DIR}/src/gpu/Backend.cpp")
set(GPU_LIBS "")
set(GPU_DEFINITIONS "")

# ---- Vulkan Compute backend ----
option(USE_VULKAN "Enable Vulkan Compute GPU backend" OFF)
if(USE_VULKAN)
    find_package(Vulkan REQUIRED)
    message("${Green}-- Vulkan FOUND:${ColourReset} ${Vulkan_LIBRARIES}")

    # Compile GLSL compute shader to SPIR-V at build time
    set(SHADER_SRC "${CMAKE_SOURCE_DIR}/src/gpu/shaders/morpho.comp")
    set(SHADER_SPV "${CMAKE_BINARY_DIR}/morpho.spv")

    find_program(GLSLC        glslc)
    find_program(GLSLANG_VAL  glslangValidator)

    if(GLSLC)
        add_custom_command(
            OUTPUT  ${SHADER_SPV}
            COMMAND ${GLSLC} -fshader-stage=compute ${SHADER_SRC} -o ${SHADER_SPV}
            DEPENDS ${SHADER_SRC}
            COMMENT "Compiling morpho.comp -> morpho.spv  [glslc]"
            VERBATIM
        )
    elseif(GLSLANG_VAL)
        add_custom_command(
            OUTPUT  ${SHADER_SPV}
            COMMAND ${GLSLANG_VAL} -V ${SHADER_SRC} -o ${SHADER_SPV}
            DEPENDS ${SHADER_SRC}
            COMMENT "Compiling morpho.comp -> morpho.spv  [glslangValidator]"
            VERBATIM
        )
    else()
        message(FATAL_ERROR
            "USE_VULKAN=ON but neither glslc nor glslangValidator was found.\n"
            "Install the Vulkan SDK (https://vulkan.lunarg.com) and re-run CMake.")
    endif()

    add_custom_target(morpho_spv DEPENDS ${SHADER_SPV})

    list(APPEND GPU_SRC "${CMAKE_SOURCE_DIR}/src/gpu/VulkanBackend.cpp")
    list(APPEND GPU_LIBS Vulkan::Vulkan)
    list(APPEND GPU_DEFINITIONS USE_VULKAN DST_VULKAN_SPIRV_PATH="${SHADER_SPV}")
    message("${Green}-- USE_VULKAN enabled.  SPIR-V -> ${SHADER_SPV}${ColourReset}")
endif()

# ---- MPI backend ----
option(USE_MPI "Enable MPI domain-decomposition backend" OFF)
if(USE_MPI)
    find_package(MPI REQUIRED)
    message("${Green}-- MPI FOUND:${ColourReset} ${MPI_CXX_LIBRARIES}")
    list(APPEND GPU_SRC "${CMAKE_SOURCE_DIR}/src/gpu/MPIBackend.cpp")
    list(APPEND GPU_LIBS MPI::MPI_CXX)
    list(APPEND GPU_DEFINITIONS USE_MPI)
    message("${Green}-- USE_MPI enabled${ColourReset}")
endif()

#-------------------- libDSTools
message("${Red}----------- libDSTools -----------${ColourReset}")
set(DSTOOLS_SOURCES
    src/DSTastrometry.cxx
    src/DSTbackground.cxx
    src/DSTimageStacking.cxx
    src/DSTstarFinder.cxx
    ${GENERATED_DIR}/src/DST_version.cxx
)
# DSTastrometry_config.h is generated: the copy an in-source build left in include/DSTools is not a source
file(GLOB DSTOOLS_PUBLIC_HEADERS ${CMAKE_SOURCE_DIR}/include/DSTools/*.h)
list(FILTER DSTOOLS_PUBLIC_HEADERS EXCLUDE REGEX "/DSTastrometry_config\\.h$")
file(GLOB DSTOOLS_GPU_HEADERS ${CMAKE_SOURCE_DIR}/include/DSTools/gpu/*.h)

add_library(DSTools SHARED ${DSTOOLS_SOURCES} ${GPU_SRC} ${DSTOOLS_PUBLIC_HEADERS})
if(DSTools_NUMERIC_VERSION)
    set_target_properties(DSTools PROPERTIES
        VERSION   "${DSTools_NUMERIC_VERSION}"
        SOVERSION "${VERSION_MAJOR}"
    )
endif()
# The generated DSTastrometry_config.h comes first: both the library and the tests find it there
target_include_directories(DSTools PUBLIC
    ${GENERATED_DIR}/include/DSTools ${GENERATED_DIR}/include
    ${CMAKE_SOURCE_DIR}/include ${CMAKE_SOURCE_DIR}/include/DSTools ${CMAKE_SOURCE_DIR}/include/DSTools/gpu
    ${CMAKE_SOURCE_DIR}/json/include
)
target_link_libraries(DSTools PUBLIC
    DSTfits::DSTfits DSTCore::DSTCore DSTdisplay::DSTdisplay DSTmath::DSTmath deepsky_externals
    ${QT_LIBS} ${GPU_LIBS}
)
target_compile_definitions(DSTools PUBLIC ${GPU_DEFINITIONS})
if(USE_VULKAN)
    add_dependencies(DSTools morpho_spv)
endif()

#-------------------- FitsViewer
add_executable(FitsViewer utilities/FitsViewer/DSTfitsViewer.cxx)
target_link_libraries(FitsViewer PRIVATE DSTools)
if(ASDF2FITS_USABLE)
    target_compile_definitions(FitsViewer PRIVATE __WITH_ASDF__)
    target_link_libraries(FitsViewer PRIVATE asdf2fits::asdf2fits)
endif()

#-------------------- Tests
if(BUILD_TESTS)
    message("${Red}----------- GoogleTest -----------${ColourReset}")
    if(NOT EXISTS ${CMAKE_SOURCE_DIR}/googletest/CMakeLists.txt)
        message(FATAL_ERROR "${BoldRed}googletest is missing: run 'git submodule update --init' (or configure with -DBUILD_TESTS=OFF)${ColourReset}")
    endif()
    set(DSTOOLS_TESTS ON)
    enable_testing()
    set(INSTALL_GTEST OFF CACHE BOOL "" FORCE)
    set(BUILD_GMOCK   OFF CACHE BOOL "" FORCE)
    set(BUILD_SHARED_LIBS OFF)
    add_subdirectory(googletest EXCLUDE_FROM_ALL)
    include(GoogleTest)

    # Test data, downloaded once into the build directory. Both downloads are optional: the tests
    # that need a missing file skip themselves. A stamp file is the output, so that a failed
    # download (no network) is not retried on every build: delete the stamp to retry.
    set(TESTDATA_DIR ${CMAKE_BINARY_DIR}/testdata)
    file(MAKE_DIRECTORY ${TESTDATA_DIR})
    set(_DL_SCRIPT "${CMAKE_SOURCE_DIR}/modules/DownloadFile.cmake")
    set(TESTDATA_STAMPS "")

    # Astrometry reference image: M11, DSS2 Red, 512x512, 0.5° FOV, from the CDS HiPS2FITS
    # service (the REST API that ESA Sky uses internally).
    set(ASTROMETRY_REF_FITS "${TESTDATA_DIR}/astrometry_reference.fits")
    set(_HIPS2FITS_URL
        "https://skyview.gsfc.nasa.gov/current/cgi/runquery.pl?Survey=DSS2+Red&Position=282.7633,-6.2717&Size=0.5&Pixels=512&Return=FITS")
    add_custom_command(
        OUTPUT  ${TESTDATA_DIR}/astrometry_reference.stamp
        COMMAND ${CMAKE_COMMAND} "-DURL=${_HIPS2FITS_URL}" "-DDEST=${ASTROMETRY_REF_FITS}"
                -DSKIP_IF_EXISTS=1 -DOPTIONAL=1 -P "${_DL_SCRIPT}"
        COMMAND ${CMAKE_COMMAND} -E touch ${TESTDATA_DIR}/astrometry_reference.stamp
        COMMENT "Downloading the M11 reference FITS (astrometry tests)"
        VERBATIM
    )
    list(APPEND TESTDATA_STAMPS ${TESTDATA_DIR}/astrometry_reference.stamp)

    # Index file for the astrometry.net integration tests, when solve-field is installed
    if(ASTROMETRY_NET_SOLVE_FIELD)
        set(ASTROMETRY_INDEX_DIR  "${TESTDATA_DIR}/astrometry_index")
        set(ASTROMETRY_INDEX_FILE "${ASTROMETRY_INDEX_DIR}/index-5206-30.fits")
        set(ASTROMETRY_INDEX_URL
            "https://portal.nersc.gov/project/cosmo/temp/dstn/index-5200/index-5206-30.fits")
        add_custom_command(
            OUTPUT  ${ASTROMETRY_INDEX_DIR}/download.stamp
            COMMAND ${CMAKE_COMMAND} -E make_directory "${ASTROMETRY_INDEX_DIR}"
            COMMAND ${CMAKE_COMMAND} "-DURL=${ASTROMETRY_INDEX_URL}" "-DDEST=${ASTROMETRY_INDEX_FILE}"
                    -DSKIP_IF_EXISTS=1 -DOPTIONAL=1 -P "${_DL_SCRIPT}"
            COMMAND ${CMAKE_COMMAND} -E touch ${ASTROMETRY_INDEX_DIR}/download.stamp
            COMMENT "Downloading astrometry.net index file index-5206-30.fits"
            VERBATIM
        )
        list(APPEND TESTDATA_STAMPS ${ASTROMETRY_INDEX_DIR}/download.stamp)
        message("${Green}-- Astrometry index will be downloaded to:${ColourReset} ${ASTROMETRY_INDEX_FILE}")
    endif()
    add_custom_target(testdata DEPENDS ${TESTDATA_STAMPS})

    # The tests run from the build directory: they read testdata/ and whatever they write stays there
    file(GLOB TEST_SOURCES CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/test/gTest_*.cxx)
    set(TEST_NAMES "")
    foreach(test_source IN LISTS TEST_SOURCES)
        get_filename_component(tname ${test_source} NAME_WE)
        list(APPEND TEST_NAMES ${tname})
        add_executable(${tname} ${test_source})
        target_link_libraries(${tname} PRIVATE DSTools gtest gtest_main)
        target_include_directories(${tname} PRIVATE ${GENERATED_DIR}/test)
        target_compile_definitions(${tname} PRIVATE SCAMP_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/scamp/tests")
        if(ASTROMETRY_NET_SOLVE_FIELD)
            target_compile_definitions(${tname} PRIVATE ASTROMETRY_NET_INDEX_FILE="${ASTROMETRY_INDEX_FILE}")
        endif()
        add_dependencies(${tname} testdata)
        gtest_discover_tests(${tname} WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
    endforeach()
    message("${Blue}-- Test executables ${ColourReset} ${TEST_NAMES} (run with ctest)")
else()
    set(DSTOOLS_TESTS OFF)
endif()

#-------------------- Doxygen documentation (target: doc, written in build/doc)
find_package(Doxygen)
if(DOXYGEN_FOUND)
    # The Doxyfile paths are relative to the sources: run there, but write in the build directory
    file(READ ${CMAKE_SOURCE_DIR}/DST.Doxyfile _doxyfile)
    string(APPEND _doxyfile "\nOUTPUT_DIRECTORY = ${CMAKE_BINARY_DIR}/doc\n")
    file(WRITE ${GENERATED_DIR}/DST.Doxyfile "${_doxyfile}")
    set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS ${CMAKE_SOURCE_DIR}/DST.Doxyfile)
    add_custom_target(doc
        COMMAND ${DOXYGEN_EXECUTABLE} ${GENERATED_DIR}/DST.Doxyfile
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        COMMENT "Generating API documentation with Doxygen in ${CMAKE_BINARY_DIR}/doc"
        VERBATIM
    )
endif()

#-------------------- Installation
install(FILES ${DSTOOLS_PUBLIC_HEADERS} ${GENERATED_DIR}/include/DSTools/DSTastrometry_config.h
        DESTINATION include/DSTools)
install(FILES ${DSTOOLS_GPU_HEADERS} DESTINATION include/DSTools/gpu)
# nlohmann/json headers (header-only, installed for downstream consumers)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/json/include/nlohmann
    DESTINATION include
    FILES_MATCHING PATTERN "*.hpp"
)
# scamp executable (built by ExternalProject, not a CMake target)
install(PROGRAMS ${SCAMP_PROGRAM} DESTINATION exe)
install(TARGETS DSTools FitsViewer
    RUNTIME DESTINATION exe
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)

if(DSTools_NUMERIC_VERSION)
    # CMake package, so that dependent projects can call find_package(DSTools <version>)
    include(CMakePackageConfigHelpers)
    set(DSTools_CMAKE_DIR lib/cmake/DSTools)
    configure_package_config_file(
        ${CMAKE_SOURCE_DIR}/modules/DSToolsConfig.cmake.in
        ${GENERATED_DIR}/DSToolsConfig.cmake
        INSTALL_DESTINATION ${DSTools_CMAKE_DIR}
    )
    write_basic_package_version_file(
        ${GENERATED_DIR}/DSToolsConfigVersion.cmake
        VERSION ${DSTools_NUMERIC_VERSION}
        COMPATIBILITY AnyNewerVersion
    )
    install(FILES
        ${GENERATED_DIR}/DSToolsConfig.cmake
        ${GENERATED_DIR}/DSToolsConfigVersion.cmake
        DESTINATION ${DSTools_CMAKE_DIR}
    )

    # CMake never removes files from the install directory: drop the versioned libraries
    # (libDSTools.2.2.0.dylib, libDSTools.so.2.2.0, ...) left behind by previous installs.
    install(CODE "set(DSTools_KEEP_VERSIONS \"${DSTools_NUMERIC_VERSION};${VERSION_MAJOR}\")")
    install(CODE [[
        file(GLOB _dst_libs "$ENV{DESTDIR}${CMAKE_INSTALL_PREFIX}/lib/libDSTools.*")
        foreach(_dst_lib ${_dst_libs})
            get_filename_component(_dst_name "${_dst_lib}" NAME)
            if(_dst_name MATCHES "^libDSTools\\.(so\\.)?([0-9.]*[0-9])(\\.dylib)?$")
                list(FIND DSTools_KEEP_VERSIONS "${CMAKE_MATCH_2}" _dst_keep)
                if(_dst_keep EQUAL -1)
                    message(STATUS "Removing old library: ${_dst_lib}")
                    file(REMOVE "${_dst_lib}")
                endif()
            endif()
        endforeach()
    ]])
endif()

#-------------------- Summary
message("${BoldGreen}----------- DeepSkyUtilities ${VERSION} configured -----------${ColourReset}")
message("${Green}-- Build directory  :${ColourReset} ${CMAKE_BINARY_DIR}")
message("${Green}-- Build type       :${ColourReset} ${CMAKE_BUILD_TYPE}")
message("${Green}-- AddressSanitizer :${ColourReset} ${WITH_ASAN}")
message("${Green}-- Tests            :${ColourReset} ${DSTOOLS_TESTS}")
message("${Green}-- DSTools          :${ColourReset} ${DSTools_NUMERIC_VERSION}")
message("${Green}-- Qt6              :${ColourReset} ${Qt6_FOUND} ${Qt6_VERSION}")
message("${Green}-- asdf2fits        :${ColourReset} ${ASDF2FITS_USABLE}")
message("${Green}-- solve-field      :${ColourReset} ${ASTROMETRY_NET_SOLVE_FIELD}")
message("${Green}-- GPU backends     :${ColourReset} Vulkan ${USE_VULKAN}, MPI ${USE_MPI}")
message("${Green}-- Install prefix   :${ColourReset} ${CMAKE_INSTALL_PREFIX}")
if(WITH_ASAN)
    message("${Yellow}--                    AddressSanitizer build: installing it replaces the release libDSTools${ColourReset}")
endif()
message("${BoldGreen}-------------------------------------------------${ColourReset}")
