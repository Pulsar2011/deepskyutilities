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
	# Linux specific code
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pipe -fexceptions -fPIC  -m64 -Wno-unknown-pragmas")

	set(OperatingSystem "Linux")
	ADD_DEFINITIONS(-DLinuxx86_64)
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Linux")

#-------------------- DARWIN specifics
IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
	# Mac OS X specific code
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -m64 -fPIC")

	set(OperatingSystem "Mac OS X")
	ADD_DEFINITIONS(-DDarwinx86_64)
ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")

option(WITH_DEBUG "Activate debug mode" OFF)
if(WITH_DEBUG)
	message("${Blue}-- DEBUG MODE ACTIVATED${ColourReset}")
	add_definitions(-D__DEBUG__)
	set(CMAKE_BUILD_TYPE DEBUG)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0")
else(WITH_DEBUG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
endif(WITH_DEBUG)

message("${Green}-- CXX COMPILER:${ColourReset} ${CMAKE_CXX_COMPILER}")
message("${Green}-- DEFAULT CXX FLAGS:${ColourReset} ${CMAKE_CXX_FLAGS}")

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



#-------------------- LIST SOURCE FILE

file(GLOB CORE_SRC "src/*.cxx")
file(GLOB GEN_SRC "${CMAKE_CURRENT_BINARY_DIR}/src/*.cxx")
file(GLOB CORE_HEADER "include/DSTools/*.h")
file(GLOB GEN_HEADER "${CMAKE_CURRENT_BINARY_DIR}/test/*.h")
file(GLOB TEST_SRC "test/*.cxx")

set(HEADER_FILES ${CORE_HEADER})

#-------------------- Configuration de la bibliothèque

add_library(
    DSTools
    SHARED
    ${CORE_SRC}
    ${CORE_HEADER}
    ${version_file}
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

    if(COMPILER_SUPPORTS_CXX17)
        target_link_libraries( ${tSRC} DSTools DSTFits DSTCore DSTdisplay DSTools DSTmath cfitsio Minuit2 ${QT_LIBS} gtest gtest_main)
    else(COMPILER_SUPPORTS_CXX17)
        target_link_libraries( ${tSRC} DSTools DSTFits DSTCore DSTdisplay DSTools DSTmath cfitsio Minuit2 ${QT_LIBS} ${Boost_LIBRARIES} gtest gtest_main)
    endif(COMPILER_SUPPORTS_CXX17)


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
    utilities/FitsViewer/DSTFitsViewer.cxx
)

target_link_libraries(FitsViewer DSTFits DSTCore DSTdisplay DSTools DSTmath cfitsio Minuit2 ${QT_LIBS})

# Ensure RPATH is embedded into the executable
set_target_properties(FitsViewer PROPERTIES
    BUILD_RPATH "${CMAKE_BUILD_RPATH}"
    INSTALL_RPATH "${CMAKE_INSTALL_RPATH}"
)

#-------------------- ADD DISTCLEAN CUSTOM TARGET

ADD_CUSTOM_TARGET (distclean @echo cleaning for source distribution)
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

ADD_CUSTOM_COMMAND(
DEPENDS clean
COMMENT "distribution clean"
COMMAND rm ARGS    -Rf ${DISTCLEANED}
COMMAND test -e install_manifest.txt || xargs ARGS rm -f < install_manifest.txt
COMMAND test -e install_manifest.txt || rm ARGS -f install_manifest.txt
TARGET  distclean
)


#-------------------- Installation proceedure
if(CMAKE_INSTALL_PREFIX)
    INSTALL(FILES ${CORE_HEADER} DESTINATION include/DSTools)
    INSTALL(TARGETS DSTools FitsViewer
		RUNTIME DESTINATION exe
		LIBRARY DESTINATION lib
		ARCHIVE DESTINATION lib
)
endif(CMAKE_INSTALL_PREFIX)
