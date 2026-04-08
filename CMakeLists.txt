#-------------------- Configuration du projet
#-------------------- CMAKE Configuration
CMAKE_MINIMUM_REQUIRED(VERSION 3.5)

project(DeepSkyProject CXX)
message("${Green}-- Project source directory :" ${CMAKE_SOURCE_DIR} "${ColourReset}")

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
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/test ./include ./include/DSTools)

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
    target_link_libraries( DSTools DSTfits DSTCore DSTdisplay DSTmath cfitsio Minuit2 ${QT_LIBS} ${GPU_LIBS})
else(COMPILER_SUPPORTS_CXX23)
    target_link_libraries( DSTools DSTfits DSTCore DSTdisplay DSTmath cfitsio Minuit2 ${QT_LIBS} ${Boost_LIBRARIES} ${GPU_LIBS})
endif(COMPILER_SUPPORTS_CXX23)

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
        target_link_libraries( ${tSRC} DSTools DSTfits DSTCore DSTdisplay DSTmath cfitsio Minuit2 ${QT_LIBS} gtest gtest_main)
    else(COMPILER_SUPPORTS_CXX23)
        target_link_libraries( ${tSRC} DSTools DSTfits DSTCore DSTdisplay DSTmath cfitsio Minuit2 ${QT_LIBS} ${Boost_LIBRARIES} gtest gtest_main)
    endif(COMPILER_SUPPORTS_CXX23)

    set_target_properties(${tSRC} PROPERTIES
        BUILD_RPATH   "${CMAKE_CURRENT_BINARY_DIR}/build/lib;${DeepSkyROOT}/lib"
        INSTALL_RPATH "${DeepSkyROOT}/lib"
    )

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

target_link_libraries(FitsViewer DSTfits DSTCore DSTdisplay DSTools DSTmath cfitsio Minuit2 ${QT_LIBS})

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
    INSTALL(TARGETS DSTools FitsViewer
		RUNTIME DESTINATION exe
		LIBRARY DESTINATION lib
		ARCHIVE DESTINATION lib
)
endif(CMAKE_INSTALL_PREFIX)
