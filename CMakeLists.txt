#-------------------- Configuration du projet
#-------------------- CMAKE Configuration
CMAKE_MINIMUM_REQUIRED(VERSION 3.5)

project(DeepSkyProject CXX)
message("${Green}-- Project source directory :" ${CMAKE_SOURCE_DIR} "${ColourReset}")

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_BUILD_DIRECTORY ${CMAKE_BINARY_DIR}/build)

if(DEFINED ENV{DeepSkyROOT})
    set(DeepSkyROOT $ENV{DeepSkyROOT})
else()
    set(DeepSkyROOT ${CMAKE_SOURCE_DIR}/..)
endif()

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

if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/../. CACHE PATH "default install path" FORCE)
endif()

set(EXECUTABLE_OUTPUT_PATH ./build/exe/)
set(LIBRARY_OUTPUT_PATH ./build/lib/)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/modules")

message("${BoldYellow}-- CMAKE INSTALL PREFIX : ${ColourReset}" ${CMAKE_INSTALL_PREFIX})

# Add macOS RPATH configuration
if(APPLE)
    set(CMAKE_MACOSX_RPATH ON)
    set(CMAKE_SKIP_BUILD_RPATH FALSE)
    set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)

    # RPATH when installed: DeepSkyProject/exe -> DeepSkyProject/lib
    set(CMAKE_INSTALL_RPATH "@executable_path/../lib;@loader_path/../lib")
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

    # RPATH when running from the build tree
    list(APPEND CMAKE_BUILD_RPATH
        "${CMAKE_CURRENT_BINARY_DIR}/build/lib"
        "${DeepSkyROOT}/lib"
    )
endif()

#-------------------- COMPILER
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++23" COMPILER_SUPPORTS_CXX23)
CHECK_CXX_COMPILER_FLAG("-std=c++20" COMPILER_SUPPORTS_CXX20)
CHECK_CXX_COMPILER_FLAG("-std=c++17" COMPILER_SUPPORTS_CXX17)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)

if(COMPILER_SUPPORTS_CXX23)
    message("-- USE CXX23 STANDARD")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++23")
    set(CMAKE_CXX_STANDARD "23")
elseif(COMPILER_SUPPORTS_CXX20)
    message("-- USE CXX20 STANDARD")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++20")
    set(CMAKE_CXX_STANDARD "20")
elseif(COMPILER_SUPPORTS_CXX17)
    message("-- USE CXX17 STANDARD")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17")
    set(CMAKE_CXX_STANDARD "17")
elseif(COMPILER_SUPPORTS_CXX11)
    message("-- USE CXX11 STANDARD")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
    set(CMAKE_CXX_STANDARD "11")
elseif(COMPILER_SUPPORTS_CXX0X)
    message("-- USE CX0X STANDARD")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
	message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++17 nor C++11 nor C++0x support. Please use a different C++ compiler.")
endif()

set(BUILD_SHARED_LIBS "YES")

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -v")

#-------------------- GCC specifics
IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")
    cmake_policy(SET CMP0069 NEW)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -Wall -pipe -pthread -fexceptions -fPIC -Wno-unknown-pragmas -Wmisleading-indentation -flto -mtune=native -ffp-contract=on -fomit-frame-pointer")
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|amd64")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -m64")
        ADD_DEFINITIONS(-DLinuxx86_64)
    elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm64")
        ADD_DEFINITIONS(-DLinuxaarch64)
    endif()

    include(CheckIPOSupported)
    check_ipo_supported(RESULT ipo_supported OUTPUT ipo_error)
    if(ipo_supported)
        set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)
    endif()

    set(OperatingSystem "Linux")
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    set(CMAKE_CXX_EXTENSIONS OFF)

    # ---- RPATH Linux ----
    set(CMAKE_SKIP_BUILD_RPATH FALSE)
    set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
    set(CMAKE_INSTALL_RPATH "$ORIGIN/../lib:${DeepSkyROOT}/lib")
    set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
    set(CMAKE_EXE_LINKER_FLAGS    "${CMAKE_EXE_LINKER_FLAGS} -Wl,--disable-new-dtags")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -Wl,--disable-new-dtags")
    list(APPEND CMAKE_BUILD_RPATH
        "${CMAKE_CURRENT_BINARY_DIR}/build/lib"
        "${DeepSkyROOT}/lib"
    )

ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

#-------------------- DARWIN specifics
IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    # macOS specific code with architecture distinction
    if(CMAKE_SYSTEM_PROCESSOR MATCHES "arm64|aarch64")
        # Apple Silicon (ARM64)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fPIC -mcpu=native -mtune=native -ffp-contract=on")
        set(OperatingSystem "Mac OS X")
        ADD_DEFINITIONS(-DDarwin -DDarwinarm64)
        set(CMAKE_OSX_ARCHITECTURES "arm64")
    elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|amd64")
        # Apple Intel (x86_64)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -m64 -fPIC -mcpu=native -mtune=native -ffp-contract=on")
        set(OperatingSystem "Mac OS X")
        ADD_DEFINITIONS(-DDarwin -DDarwinx86_64)
        set(CMAKE_OSX_ARCHITECTURES "x86_64")
    else()
        # Fallback
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fPIC -mcpu=native -mtune=native -ffp-contract=on")
        set(OperatingSystem "Mac OS X")
        ADD_DEFINITIONS(-DDarwin)
    endif()
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

option(WITH_DEBUG "Activate debug mode" OFF)
if(WITH_DEBUG)
	message("${Blue}-- DEBUG MODE ACTIVATED${ColourReset}")
	add_definitions(-D__DEBUG__)
	# Ensure proper CMake build type name and debug flags so lldb can see symbols.
	# Use -g3 for maximum debug info, -O0 to avoid optimization and -fno-omit-frame-pointer
	# Avoid adding any -s/--strip flags here.
	# Disable LTO in debug mode to preserve debug symbols
	set(CMAKE_INTERPROCEDURAL_OPTIMIZATION OFF)
	set(CMAKE_BUILD_TYPE "Debug" CACHE STRING "Build type" FORCE)
	# Populate the Debug-specific vars used by CMake and toolchains
	set(CMAKE_C_FLAGS_DEBUG   "${CMAKE_C_FLAGS_DEBUG} -g3 -O0 -fno-omit-frame-pointer -fsanitize=address")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g3 -O0 -fno-omit-frame-pointer -fsanitize=address")
	# Ensure linker doesn't strip in debug mode (leave linker flags untouched or add explicit)
	set(CMAKE_EXE_LINKER_FLAGS_DEBUG "${CMAKE_EXE_LINKER_FLAGS_DEBUG}")
	# Also append to global flags (keeps behaviour consistent with rest of file)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG}")
else(WITH_DEBUG)
	# For non-debug keep release-style optimizations but avoid automatic stripping here.
	# Re-enable LTO for Release builds on macOS
	if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -flto=thin -fomit-frame-pointer")
	endif()
	set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Build type" FORCE)
	set(CMAKE_C_FLAGS_RELEASE   "${CMAKE_C_FLAGS_RELEASE} -O3")
	set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE}")
endif(WITH_DEBUG)

message("${Green}-- CXX COMPILER:${ColourReset} ${CMAKE_CXX_COMPILER}")
message("${Green}-- DEFAULT CXX FLAGS:${ColourReset} ${CMAKE_CXX_FLAGS}")
message("${Green}-- CMAKE_BUILD_TYPE:${ColourReset} ${CMAKE_BUILD_TYPE}")
message("${Green}-- CMAKE_CXX_FLAGS_DEBUG:${ColourReset} ${CMAKE_CXX_FLAGS_DEBUG}")
message("${Green}-- CMAKE_EXE_LINKER_FLAGS_DEBUG:${ColourReset} ${CMAKE_EXE_LINKER_FLAGS_DEBUG}")

#-------------------- APPEND GIT VERSION TAG
# Make a version file containing the current version from git.
find_package(Git)
include(GetGitRevisionDescription)
git_describe(VERSION --always --tags --long)

#parse the version information into pieces.
string(REGEX REPLACE "^v([0-9]+)\\..*" "\\1" VERSION_MAJOR "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.([0-9]+).*" "\\1" VERSION_MINOR "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.([0-9]+).*" "\\1" VERSION_PATCH "${VERSION}")
string(REGEX REPLACE "^v[0-9]+\\.[0-9]+\\.[0-9]+.(-*)" "\\1" VERSION_SHA1 "${VERSION}")
set(VERSION_SHORT "v${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
message("-- PKG VERSION : " ${VERSION})
message("-- PKG VERSION : " ${VERSION_SHORT})
message("-- PKG MAJOR   : " ${VERSION_MAJOR})
message("-- PKG MINOR   : " ${VERSION_MINOR})
message("-- PKG PATCH   : " ${VERSION_PATCH})
message("-- PKG SHA1    : " ${VERSION_SHA1})

set(CMAKE_PROJECT_VERSION "${VERSION}")
set(CMAKE_PROJECT_VERSION_MAJOR "${VERSION_MAJOR}")
set(CMAKE_PROJECT_VERSION_MINOR "${VERSION_MINOR}")
set(CMAKE_PROJECT_VERSION_PATCH "${VERSION_PATCH}")

configure_file(${CMAKE_MODULE_PATH}/DST_version.cxx.in ${CMAKE_CURRENT_BINARY_DIR}/src/DST_version.cxx)
set(version_file "${CMAKE_CURRENT_BINARY_DIR}/src/DST_version.cxx")

configure_file(${CMAKE_MODULE_PATH}/gTest_version.h.in ${CMAKE_CURRENT_BINARY_DIR}/test/gTest_version.h)

#-------------------- Inclusion des en-têtes publics
# The generated DSTastrometry_config.h lives in the build tree; add it first
# so both the library and downstream test executables find it automatically.
INCLUDE_DIRECTORIES(
    ${CMAKE_CURRENT_BINARY_DIR}/include/DSTools
    ${CMAKE_CURRENT_BINARY_DIR}/include
    ${CMAKE_CURRENT_BINARY_DIR}/test
    ./include
    ./include/DSTools
)

#-------------------- nlohmann/json (header-only submodule)
include_directories(${CMAKE_SOURCE_DIR}/json/include)
message("${Green}-- nlohmann/json headers:${ColourReset} ${CMAKE_SOURCE_DIR}/json/include")

message("-- PROJECT DIR: " ${CMAKE_SOURCE_DIR})
message("-- PROJECT BINARY DIR: " ${PROJECT_BINARY_DIR})
message("-- CURRENT BINARY DIR: " ${CMAKE_CURRENT_BINARY_DIR})
message("-- CMAKE SOURCE DIR  : " ${CMAKE_SOURCE_DIR})


#-------------------- Find Boost dependencies
if(!${COMPILER_SUPPORTS_CXX17})
    message("${Red}----------- SEARCH BOOST Framework -----------${ColourReset}")
    find_package(Boost REQUIRED COMPONENTS filesystem chrono)
    if(${Boost_FOUND})
        message("${Blue}-- Boost INC DIR ${ColourReset} ${Boost_INCLUDE_DIRS}")
        message("${Blue}-- Boost LIBRARY DIR ${ColourReset} ${Boost_LIBRARY_DIRS}")
        message("${Blue}-- Boost LIB ${ColourReset} ${Boost_LIBRARIES}")
    endif(${Boost_FOUND})
    set(Boost_USE_MULTITHREADED ON)
endif(!${COMPILER_SUPPORTS_CXX17})

#-------------------- Find QT dependencies. NOTE THAT Qt is optional package
find_package(Qt6 COMPONENTS Widgets Core Gui PrintSupport Concurrent)
if(Qt6_FOUND)
    message("${Blue}-- Qt6 FOUND : ${ColourReset} ${Qt6_VERSION}")
    include_directories(${Qt6Widgets_INCLUDE_DIRS} ${Qt6Core_INCLUDE_DIRS} ${Qt6Gui_INCLUDE_DIRS} ${Qt6PrintSupport_INCLUDE_DIRS} ${Qt6Concurrent_INCLUDE_DIRS})
    add_definitions(${Qt6Widgets_DEFINITIONS} ${Qt6Core_DEFINITIONS} ${Qt6Gui_DEFINITIONS} ${Qt6PrintSupport_DEFINITIONS} ${Qt6Concurrent_DEFINITIONS})
    set(QT_LIBS Qt6::Widgets Qt6::Core Qt6::Gui Qt6::PrintSupport Qt6::Concurrent)
    add_definitions(-D__WITH_QT__)

    # Enable Qt's code generators
    set(CMAKE_AUTOMOC ON)
    set(CMAKE_AUTOUIC ON)
    set(CMAKE_AUTORCC ON)

else(Qt6_FOUND)
    message("${Red}-- Qt6 NOT FOUND. GUI will be disabled${ColourReset}")
    set(QT_LIBS "")
endif(Qt6_FOUND)

link_directories(${DeepSkyROOT}/lib)
include_directories(${DeepSkyROOT}/include)
find_library(DSTCore NAMES DSTCore HINTS ${DeepSkyROOT}/lib)
if(DSTCore)
    message("-- ${Green}Found${ColourReset} DSTCore: ${DSTCore}")
else()
    message(WARNING "libDSTCore not found; install the library")
endif()

find_library(DSTmath NAMES DSTmath HINTS ${DeepSkyROOT}/lib)
if(DSTmath)
    message("-- ${Green}Found${ColourReset} DSTmath: ${DSTmath}")
else()
    message(WARNING "libDSTmath not found; install the library")
endif()

link_directories(${DeepSkyROOT}/lib)
find_library(DSTfits NAMES DSTfits HINTS ${DeepSkyROOT}/lib)
if(DSTfits)
    message("-- ${Green}Found${ColourReset} DSTfits: ${DSTfits}")
else()
    message(WARNING "libDSTfits not found; install the library")
endif()

find_library(DSTdisplay NAMES DSTdisplay HINTS ${DeepSkyROOT}/lib)
if(DSTdisplay)
    message("-- ${Green}Found${ColourReset} DSTdysplay: ${DSTdisplay}")
else()
    message(WARNING "libDSTdisplay not found; install the library")
endif()

find_library(cfitsio NAMES cfitsio HINTS ${DeepSkyROOT}/lib)
if(cfitsio)
    include_directories(${DeepSkyROOT}/include/cfitsio)
    message("-- ${Green}Found${ColourReset} cfitsio: ${cfitsio}")
else()
    message(WARNING "libcfitsio not found; install the library")
endif()


find_library(Minuit2 NAMES Minuit2 HINTS ${DeepSkyROOT}/lib)
if(Minuit2)
    include_directories(${DeepSkyROOT}/include/Minuit2)
    message("-- ${Green}Found${ColourReset} Minuit2: ${Minuit2}")
else()
    message(WARNING "libMinuit2 not found; install the library")
endif()

#-------------------- Detect astrometry.net solve-field
# The detected path is written into DSTastrometry_config.h (generated header,
# installed alongside DSTastrometry.h) so downstream consumers get the same
# default.  DSTastrometry.h itself is never generated by CMake.
find_program(ASTROMETRY_NET_SOLVE_FIELD solve-field
    HINTS /opt/homebrew/bin
          /opt/homebrew/Cellar/astrometry-net/*/bin
          /usr/local/bin
          /usr/bin
)
if(ASTROMETRY_NET_SOLVE_FIELD)
    get_filename_component(ASTROMETRY_NET_DEFAULT_PATH "${ASTROMETRY_NET_SOLVE_FIELD}" DIRECTORY)
    message("${Green}-- astrometry.net solve-field:${ColourReset} ${ASTROMETRY_NET_SOLVE_FIELD}")
    message("${Green}-- ASTROMETRY_NET_DEFAULT_PATH:${ColourReset} ${ASTROMETRY_NET_DEFAULT_PATH}")

    # Download the index file needed by integration tests.
    # The file lands in build/test_output/astrometry/ so it is co-located with
    # other test artefacts and does not pollute the source tree.
    # ASTROMETRY_NET_INDEX_FILE is passed as a compile-time macro so that
    # writeMinimalConfig() in the test can inject it into the solver JSON config.
    set(ASTROMETRY_INDEX_DIR  "${CMAKE_SOURCE_DIR}/build/astrometry_index")
    set(ASTROMETRY_INDEX_FILE "${ASTROMETRY_INDEX_DIR}/index-5206-30.fits")
    set(ASTROMETRY_INDEX_URL
        "https://portal.nersc.gov/project/cosmo/temp/dstn/index-5200/index-5206-30.fits")

    # Use a stamp file as the OUTPUT so CMake does not retry the download on every
    # build invocation when the network is unavailable.  The actual index file may
    # or may not be present after the stamp is written; the astrometry test checks
    # at runtime and skips if the file is missing.
    set(ASTROMETRY_INDEX_STAMP "${ASTROMETRY_INDEX_DIR}/download.stamp")
    set(_DL_SCRIPT "${CMAKE_SOURCE_DIR}/modules/DownloadFile.cmake")
    add_custom_command(
        OUTPUT  "${ASTROMETRY_INDEX_STAMP}"
        COMMAND ${CMAKE_COMMAND} -E make_directory "${ASTROMETRY_INDEX_DIR}"
        COMMAND ${CMAKE_COMMAND} -E echo "[download_astrometry_index] cmake = ${CMAKE_COMMAND}"
        COMMAND ${CMAKE_COMMAND} -E echo "[download_astrometry_index] script = ${_DL_SCRIPT}"
        COMMAND ${CMAKE_COMMAND} -E echo "[download_astrometry_index] URL = ${ASTROMETRY_INDEX_URL}"
        COMMAND ${CMAKE_COMMAND} -E echo "[download_astrometry_index] DEST = ${ASTROMETRY_INDEX_FILE}"
        COMMAND ${CMAKE_COMMAND} --version
        COMMAND ls -la "${_DL_SCRIPT}" || true
        COMMAND head -5 "${_DL_SCRIPT}" || true
        COMMAND ${CMAKE_COMMAND}
            "-DURL=${ASTROMETRY_INDEX_URL}"
            "-DDEST=${ASTROMETRY_INDEX_FILE}"
            -DOPTIONAL=1
            -P "${_DL_SCRIPT}"
        COMMAND ${CMAKE_COMMAND} -E touch "${ASTROMETRY_INDEX_STAMP}"
        COMMENT "Downloading astrometry.net index file index-5206-30.fits"
        VERBATIM
    )

    add_custom_target(download_astrometry_index
        DEPENDS "${ASTROMETRY_INDEX_STAMP}"
        COMMENT "Astrometry.net index file ready"
    )

    add_definitions(-DASTROMETRY_NET_INDEX_FILE="${ASTROMETRY_INDEX_FILE}")
    message("${Green}-- Astrometry index will be downloaded to:${ColourReset} ${ASTROMETRY_INDEX_FILE}")
else()
    # Leave ASTROMETRY_NET_DEFAULT_PATH unset so #cmakedefine comments it out
    message("${Yellow}-- astrometry.net solve-field not found; ASTROMETRY_NET_DEFAULT_PATH will be undefined${ColourReset}")
endif()

# scamp install path (matches the DESTINATION in the install() command below)
get_filename_component(SCAMP_DEFAULT_PATH "${CMAKE_INSTALL_PREFIX}/exe" ABSOLUTE)
add_definitions(-DSCAMP_TEST_DATA_DIR="${CMAKE_SOURCE_DIR}/scamp/tests")

configure_file(
    ${CMAKE_SOURCE_DIR}/modules/DSTastrometry_config.h.in
    ${CMAKE_CURRENT_BINARY_DIR}/include/DSTools/DSTastrometry_config.h
)

#-------------------- scamp (Autotools, submodule)
include(ExternalProject)

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
    NO_DEFAULT_PATH
)
find_library(OPENBLAS_LIBRARY NAMES openblas
    HINTS $ENV{CONDA_PREFIX}/lib
        /opt/homebrew/opt/openblas/lib
          /usr/local/opt/openblas/lib
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
    )
    find_path(_LAPACKE_INCLUDE_DIR lapacke.h
      HINTS $ENV{CONDA_PREFIX}/include
          $ENV{CONDA_PREFIX}/include/openblas
    )
    if(_LAPACKE_LIB)
        set(_OB_WRAPPER_DIR "${CMAKE_BINARY_DIR}/_openblas_wrapper")

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
            # Symlink cblas.h from OpenBLAS include dir
            execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
                "${OPENBLAS_INCLUDE_DIR}/cblas.h"
                "${_OB_WRAPPER_DIR}/include/cblas.h")
            # Symlink lapacke.h from lapacke-dev include dir
            execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
                "${_LAPACKE_INCLUDE_DIR}/lapacke.h"
                "${_OB_WRAPPER_DIR}/include/lapacke.h")
            # Also symlink lapacke_mangling.h if present (pulled in by lapacke.h)
            if(EXISTS "${_LAPACKE_INCLUDE_DIR}/lapacke_mangling.h")
                execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
                    "${_LAPACKE_INCLUDE_DIR}/lapacke_mangling.h"
                    "${_OB_WRAPPER_DIR}/include/lapacke_mangling.h")
            endif()
            # Also symlink lapacke_config.h if present
            if(EXISTS "${_LAPACKE_INCLUDE_DIR}/lapacke_config.h")
                execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
                    "${_LAPACKE_INCLUDE_DIR}/lapacke_config.h"
                    "${_OB_WRAPPER_DIR}/include/lapacke_config.h")
            endif()
            set(OPENBLAS_INCLUDE_DIR "${_OB_WRAPPER_DIR}/include")
            message("${Green}-- Created OpenBLAS+LAPACKE include wrapper: ${_OB_WRAPPER_DIR}/include/${ColourReset}")
        endif()
    endif()
endif()

# Step 1: run autogen.sh to generate the configure script from configure.ac.
# Prepend the automake bin dir to PATH so autogen.sh can locate automake/aclocal
# even when the package is keg-only and not linked into /opt/homebrew/bin.
add_custom_command(
    OUTPUT  ${CMAKE_SOURCE_DIR}/scamp/configure
    COMMAND ${CMAKE_COMMAND} -E env PATH=${AUTOMAKE_BIN_DIR}:$ENV{PATH} bash -c "cd ${CMAKE_SOURCE_DIR}/scamp && ./autogen.sh"
    DEPENDS ${CMAKE_SOURCE_DIR}/scamp/configure.ac
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/scamp
    COMMENT "Running scamp/autogen.sh (generating configure script)"
    VERBATIM
)
add_custom_target(autogen_scamp DEPENDS ${CMAKE_SOURCE_DIR}/scamp/configure)

# Step 2: configure, build and install scamp.
# Use bash -c to ensure env vars and all flags are forwarded correctly
# (same pattern as autogen_scamp above).
ExternalProject_Add(
    scamp_external
    SOURCE_DIR      ${CMAKE_SOURCE_DIR}/scamp
    BINARY_DIR      ${CMAKE_SOURCE_DIR}/build/scamp
    CONFIGURE_COMMAND ${CMAKE_COMMAND} -E env
        CC=${CMAKE_C_COMPILER}
        CFLAGS=-fPIC
        CPPFLAGS=-I${OPENBLAS_INCLUDE_DIR}
        ${CMAKE_SOURCE_DIR}/scamp/configure
        --prefix=${CMAKE_SOURCE_DIR}/build
        --with-fftw-incdir=${FFTW_INCLUDE_DIR}
        --with-fftw-libdir=${FFTW_LIB_DIR}
        --enable-openblas
        --with-openblas-incdir=${OPENBLAS_INCLUDE_DIR}
        --with-openblas-libdir=${OPENBLAS_LIB_DIR}
    BUILD_COMMAND     ${CMAKE_COMMAND} -E env MAKEFLAGS= ${CMAKE_MAKE_PROGRAM} -j ${CMAKE_BUILD_PARALLEL_LEVEL} -C src
    INSTALL_COMMAND   ${CMAKE_COMMAND} -E env MAKEFLAGS= ${CMAKE_MAKE_PROGRAM} -C src install
    BUILD_BYPRODUCTS  ${CMAKE_SOURCE_DIR}/build/bin/scamp
)

add_dependencies(scamp_external autogen_scamp)
message("${Green}-- scamp will be built from submodule${ColourReset}")

# Helper target: run 'scamp --help' once to observe its options and output format.
# Usage: make scamp_probe
add_custom_target(scamp_probe
    COMMAND ${CMAKE_SOURCE_DIR}/build/bin/scamp --help || true
    DEPENDS scamp_external
    COMMENT "Running scamp --help so you can observe its output format"
    VERBATIM
)


GET_DIRECTORY_PROPERTY(include_dir_list INCLUDE_DIRECTORIES)
MESSAGE("${Green}INCLUDE DIRS${ColourReset} ${include_dir_list}")

#-------------------- LIST SOURCE FILE

file(GLOB CORE_SRC "src/*.cxx")
file(GLOB GEN_SRC "${CMAKE_CURRENT_BINARY_DIR}/src/*.cxx")
file(GLOB CORE_HEADER "include/DSTools/*.h")
file(GLOB GEN_HEADER "${CMAKE_CURRENT_BINARY_DIR}/test/*.h")
file(GLOB TEST_SRC "test/*.cxx")

set(HEADER_FILES ${CORE_HEADER})


#-------------------- GPU backends (optional)

# Backend.cpp is always compiled (contains the Backend::create() factory)
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include/DSTools/gpu")
set(GPU_SRC "${CMAKE_CURRENT_SOURCE_DIR}/src/gpu/Backend.cpp")
set(GPU_LIBS "")

# ---- Vulkan Compute backend ----
option(USE_VULKAN "Enable Vulkan Compute GPU backend" OFF)
if(USE_VULKAN)
    find_package(Vulkan REQUIRED)
    message("${Green}-- Vulkan FOUND:${ColourReset} ${Vulkan_LIBRARIES}")

    # Compile GLSL compute shader to SPIR-V at build time
    set(SHADER_SRC "${CMAKE_CURRENT_SOURCE_DIR}/src/gpu/shaders/morpho.comp")
    set(SHADER_SPV "${CMAKE_CURRENT_BINARY_DIR}/morpho.spv")

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

    list(APPEND GPU_SRC "${CMAKE_CURRENT_SOURCE_DIR}/src/gpu/VulkanBackend.cpp")
    include_directories(${Vulkan_INCLUDE_DIRS})
    list(APPEND GPU_LIBS Vulkan::Vulkan)
    add_definitions(-DUSE_VULKAN)
    add_definitions(-DDST_VULKAN_SPIRV_PATH="${SHADER_SPV}")
    message("${Green}-- USE_VULKAN enabled.  SPIR-V -> ${SHADER_SPV}${ColourReset}")
endif()

# ---- MPI backend ----
option(USE_MPI "Enable MPI domain-decomposition backend" OFF)
if(USE_MPI)
    find_package(MPI REQUIRED)
    message("${Green}-- MPI FOUND:${ColourReset} ${MPI_CXX_LIBRARIES}")
    list(APPEND GPU_SRC "${CMAKE_CURRENT_SOURCE_DIR}/src/gpu/MPIBackend.cpp")
    include_directories(${MPI_CXX_INCLUDE_DIRS})
    list(APPEND GPU_LIBS MPI::MPI_CXX)
    add_definitions(-DUSE_MPI)
    message("${Green}-- USE_MPI enabled${ColourReset}")
endif()

#-------------------- Configuration de la bibliothèque

add_library(
    DSTools
    SHARED
    ${CORE_SRC}
    ${GPU_SRC}
    ${CORE_HEADER}
    ${version_file}
)

if(USE_VULKAN)
    add_dependencies(DSTools morpho_spv)
endif()

if(COMPILER_SUPPORTS_CXX23)
    target_link_libraries( DSTools DSTfits DSTCore DSTdisplay DSTmath Minuit2 ${QT_LIBS} ${GPU_LIBS})
else(COMPILER_SUPPORTS_CXX23)
    target_link_libraries( DSTools DSTfits DSTCore DSTdisplay DSTmath Minuit2 ${QT_LIBS} ${Boost_LIBRARIES} ${GPU_LIBS})
endif(COMPILER_SUPPORTS_CXX23)

#-------------------- Astrometry test reference FITS (M11, DSS2 Red, 512x512, 0.5° FOV)
# Downloaded from the CDS HiPS2FITS service — the REST API that ESA Sky uses internally.
# The file is fetched once at configure time; re-run cmake to retry a failed download.
# Use the 'download_astrometry_ref' target to force a fresh download at build time.
set(ASTROMETRY_REF_FITS "${CMAKE_SOURCE_DIR}/test/data/astrometry_reference.fits")
set(_HIPS2FITS_URL
    "https://alasky.cds.unistra.fr/hips-image-services/hips2fits?hips=CDS%2FP%2FDSS2%2Fred&width=512&height=512&fov=0.5&projection=TAN&coordsys=icrs&ra=282.7633&dec=-6.2717&format=fits"
)

file(MAKE_DIRECTORY "${CMAKE_SOURCE_DIR}/test/data")

if(NOT EXISTS "${ASTROMETRY_REF_FITS}")
    message("${Yellow}-- Downloading M11 reference FITS from HiPS2FITS (ESA Sky backend)...${ColourReset}")
    file(DOWNLOAD
        "${_HIPS2FITS_URL}"
        "${ASTROMETRY_REF_FITS}"
        STATUS _fits_download_status
        SHOW_PROGRESS
        TLS_VERIFY ON
    )
    list(GET _fits_download_status 0 _fits_status_code)
    if(NOT _fits_status_code EQUAL 0)
        list(GET _fits_download_status 1 _fits_error_msg)
        message(WARNING
            "Could not download M11 reference FITS: ${_fits_error_msg}\n"
            "   Astrometry WCS accuracy tests (Suite 5) will be skipped.\n"
            "   Run 'make download_astrometry_ref' to retry after network is available.")
        file(REMOVE "${ASTROMETRY_REF_FITS}")
    else()
        message("${Green}-- M11 reference FITS saved: ${ASTROMETRY_REF_FITS}${ColourReset}")
    endif()
else()
    message("${Green}-- M11 reference FITS already present: ${ASTROMETRY_REF_FITS}${ColourReset}")
endif()

# Build-time guard using a binary-tree stamp file.
# The OUTPUT is inside CMAKE_BINARY_DIR so CMake can reliably track it.
# SKIP_IF_EXISTS makes the script a no-op when the configure-time download
# already saved the file — no redundant network request on every `make`.
set(_ASTROMETRY_REF_STAMP "${CMAKE_BINARY_DIR}/astrometry_ref_download.stamp")
set(_DL_SCRIPT_REF "${CMAKE_SOURCE_DIR}/modules/DownloadFile.cmake")
add_custom_command(
    OUTPUT  "${_ASTROMETRY_REF_STAMP}"
    COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_SOURCE_DIR}/test/data"
    COMMAND ${CMAKE_COMMAND} -E echo "[download_astrometry_ref] cmake = ${CMAKE_COMMAND}"
    COMMAND ${CMAKE_COMMAND} -E echo "[download_astrometry_ref] script = ${_DL_SCRIPT_REF}"
    COMMAND ${CMAKE_COMMAND} -E echo "[download_astrometry_ref] URL = ${_HIPS2FITS_URL}"
    COMMAND ${CMAKE_COMMAND} -E echo "[download_astrometry_ref] DEST = ${ASTROMETRY_REF_FITS}"
    COMMAND ${CMAKE_COMMAND} --version
    COMMAND ls -la "${_DL_SCRIPT_REF}" || true
    COMMAND head -5 "${_DL_SCRIPT_REF}" || true
    COMMAND ${CMAKE_COMMAND}
        "-DURL=${_HIPS2FITS_URL}"
        "-DDEST=${ASTROMETRY_REF_FITS}"
        -DSKIP_IF_EXISTS=1
        -P "${_DL_SCRIPT_REF}"
    COMMAND ${CMAKE_COMMAND} -E touch "${_ASTROMETRY_REF_STAMP}"
    COMMENT "Ensuring M11 reference FITS is present"
    VERBATIM
)
add_custom_target(download_astrometry_ref
    DEPENDS "${_ASTROMETRY_REF_STAMP}"
)

#-------------------- Configuration des test
add_subdirectory(googletest   EXCLUDE_FROM_ALL)
# Make googletest headers available and enable GoogleTest helper which
# provides the gtest_discover_tests() function.
include_directories(${CMAKE_SOURCE_DIR}/googletest/googletest/include)
include(GoogleTest)
enable_testing()

foreach( test_source_file ${TEST_SRC} )
	get_filename_component( tSRC ${test_source_file} NAME_WE)
	add_executable(
		${tSRC}
		${test_source_file}
	)

    if(COMPILER_SUPPORTS_CXX23)
        target_link_libraries( ${tSRC} DSTools DSTfits DSTCore DSTdisplay DSTmath Minuit2 ${QT_LIBS} gtest gtest_main)
    else(COMPILER_SUPPORTS_CXX23)
        target_link_libraries( ${tSRC} DSTools DSTfits DSTCore DSTdisplay DSTmath Minuit2 ${QT_LIBS} ${Boost_LIBRARIES} gtest gtest_main)
    endif(COMPILER_SUPPORTS_CXX23)

    set_target_properties(${tSRC} PROPERTIES
        BUILD_RPATH   "${CMAKE_CURRENT_BINARY_DIR}/build/lib;${DeepSkyROOT}/lib"
        INSTALL_RPATH "${DeepSkyROOT}/lib"
    )

    add_dependencies(${tSRC} download_astrometry_ref)

    if(ASTROMETRY_NET_SOLVE_FIELD)
        add_dependencies(${tSRC} download_astrometry_index)
    endif()

endforeach( test_source_file ${TEST_SRC} )

foreach( test_source_file ${TEST_SRC} )
    get_filename_component( tSRC ${test_source_file} NAME_WE)
    gtest_discover_tests(${tSRC})
endforeach( test_source_file ${TEST_SRC} )

#-------------------- CONFIGURATION OF DOXYGEN DOC GENERATION

find_package(Doxygen)
if(DOXYGEN_FOUND)
    configure_file(${CMAKE_CURRENT_SOURCE_DIR}/DST.Doxyfile ${CMAKE_CURRENT_BINARY_DIR}/DST.Doxyfile @ONLY)
    add_custom_target(doc
    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/DST.Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM
)
endif(DOXYGEN_FOUND)

#######################
add_executable(
    FitsViewer
    utilities/FitsViewer/DSTfitsViewer.cxx
)

target_link_libraries(FitsViewer DSTfits DSTCore DSTdisplay DSTools DSTmath Minuit2 ${QT_LIBS})

# Embed the install RPATH at build time so that install_name_tool does not
# attempt to rewrite it on every subsequent "make install", which causes
# duplicate-path errors on macOS when the binary is already up to date.
set_target_properties(FitsViewer PROPERTIES
    BUILD_RPATH   "${CMAKE_CURRENT_BINARY_DIR}/build/lib;${DeepSkyROOT}/lib"
    INSTALL_RPATH "${CMAKE_INSTALL_RPATH}"
)

#-------------------- ADD DISTCLEAN CUSTOM TARGET

add_custom_target(distclean)
set(DISTCLEANED
cmake.depends
cmake.check_depends
CMakeCache.txt
cmake.check_cache
CMakeTmp
${version_file}
Makefile
CMakeFiles
*.cmake
*_autogen
./build/lib
./build/exe
./build/test
./doc
./build/scamp
./scamp/configure
./scamp/Makefile.in
./scamp/aclocal.m4
./googletest/CMakeFiles
./googletest/Makefile
./googletest/*.cmake
./lib
./bin
)

add_custom_target(
    distclean_run
    COMMENT "distribution clean"
    COMMAND ${CMAKE_COMMAND} -E echo "distribution clean"
    COMMAND rm -Rf ${DISTCLEANED}
    COMMAND /usr/bin/env bash -c "test -e install_manifest.txt && xargs rm -f < install_manifest.txt || true"
    COMMAND /usr/bin/env bash -c "test -e install_manifest.txt && rm -f install_manifest.txt || true"
    VERBATIM
)
add_dependencies(distclean distclean_run)


#-------------------- Installation proceedure
if(CMAKE_INSTALL_PREFIX)
    INSTALL(FILES ${CORE_HEADER} DESTINATION include/DSTools)
    INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/include/DSTools/DSTastrometry_config.h
            DESTINATION include/DSTools)
    INSTALL(FILES
        include/DSTools/gpu/Backend.h
        include/DSTools/gpu/CPUBackend.h
        include/DSTools/gpu/GPUContext.h
        include/DSTools/gpu/VulkanBackend.h
        include/DSTools/gpu/CUDABackend.h
        include/DSTools/gpu/HIPBackend.h
        include/DSTools/gpu/MPIBackend.h
        DESTINATION include/DSTools/gpu
    )
    # nlohmann/json headers (header-only, installed for downstream consumers)
    INSTALL(DIRECTORY ${CMAKE_SOURCE_DIR}/json/include/nlohmann
        DESTINATION include
        FILES_MATCHING PATTERN "*.hpp"
    )
    # scamp executable (built by ExternalProject, not a CMake target)
    install(PROGRAMS ${CMAKE_SOURCE_DIR}/build/bin/scamp DESTINATION exe)
    INSTALL(TARGETS DSTools FitsViewer
		RUNTIME DESTINATION exe
		LIBRARY DESTINATION lib
		ARCHIVE DESTINATION lib
)
endif(CMAKE_INSTALL_PREFIX)
