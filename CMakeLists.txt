cmake_minimum_required (VERSION 3.22)

project(expos
  LANGUAGES CXX
  VERSION 0.0.0
)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

option(ENABLE_DEBUG_FLAGS "Enable extra debug compile options" OFF)

option(MAKE_MAIN "Compile main expos binary" ON)
option(MAKE_TEST "Compile tests binary" OFF)
option(MAKE_DAUGHTER "Compile daughter tool binaries" OFF)

# user paths to htslib if needed
set(HTSLIB_INCLUDE_DIR "" CACHE PATH "Path to directory that contains htslib/hts.h")
set(HTSLIB_LIBRARY "" CACHE FILEPATH "Absolute path to libhts.{so,dylib,a}")

# //--- find deps ---//
find_package(PkgConfig QUIET)

# //--- cxxopts ---//
include(FetchContent)

FetchContent_Declare(
  cxxopts
  GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git
  GIT_TAG        v3.3.1  # what I developed against, but any 3.x.x should be fine
)
FetchContent_MakeAvailable(cxxopts)

# ------- htslib -------
set(HTSLIB_TARGET "")
set(_htslib_inc "")
set(_htslib_lib "")

if(HTSLIB_INCLUDE_DIR AND HTSLIB_LIBRARY)
  message(STATUS "Using user-specified htslib install. Must be >=1.14, but this is not enforced.")
  set(_htslib_inc "${HTSLIB_INCLUDE_DIR}")
  set(_htslib_lib "${HTSLIB_LIBRARY}")
endif()

if(_htslib_inc AND _htslib_lib)
  add_library(htslib UNKNOWN IMPORTED GLOBAL)
  add_library(htslib::htslib ALIAS htslib)
  set_target_properties(htslib PROPERTIES
    IMPORTED_LOCATION "${_htslib_lib}"
    INTERFACE_INCLUDE_DIRECTORIES "${_htslib_inc}"
  )
  set(HTSLIB_TARGET htslib::htslib)
endif()

if(NOT HTSLIB_TARGET AND PkgConfig_FOUND)
  # assuming the pkg-config module is called "htslib"
  pkg_check_modules(HTSLIB QUIET IMPORTED_TARGET htslib>=1.14)
  if(HTSLIB_FOUND)
    message(STATUS "Located htslib>=1.14 install via pkg-config.")
    set(HTSLIB_TARGET PkgConfig::HTSLIB)
  endif()
endif()

if(NOT HTSLIB_TARGET)
  message(FATAL_ERROR
    "Could not find htslib directly or via pkg-config. "
    "Set both HTSLIB_INCLUDE_DIR and HTSLIB_LIBRARY to use user install."
  )
endif()
# //--- end deps ---//

set(COMMON_INCLUDES
  ${CMAKE_CURRENT_SOURCE_DIR}/src
)

set(COMMON_WARNINGS
  -Wall -Wextra -Wpedantic -Wshadow
  -Wconversion -Wsign-conversion
  -Wfloat-equal -Wnon-virtual-dtor
  -Woverloaded-virtual -Wnull-dereference
  -Wdouble-promotion -Wformat=2
)

# //--- internal libs ---//
add_library(lib_stats INTERFACE)

target_include_directories(lib_stats INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# //--- testing ---//
if(MAKE_TEST)
  message(STATUS "Finding catch2 for test binary")
  Include(FetchContent)

  FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.8.1
  )

  FetchContent_MakeAvailable(Catch2)

  add_executable(test-expos tests/test.cpp)

  target_compile_options(test-expos PRIVATE
    ${COMMON_WARNINGS}
  )
  target_include_directories(test-expos PRIVATE
    ${COMMON_INCLUDES}
  )
  target_link_libraries(test-expos PRIVATE
    Catch2::Catch2WithMain
  )
endif(MAKE_TEST)
# //--- end testing ---//


if(MAKE_MAIN)
  add_executable(expos
    src/expos/main.cpp
    src/expos/pileup.cpp
  )
  target_compile_options(expos PRIVATE
    ${COMMON_WARNINGS}
  )
  target_include_directories(expos PRIVATE
    ${COMMON_INCLUDES}
  )
  target_link_libraries(expos
    lib_stats
    cxxopts::cxxopts
    ${HTSLIB_TARGET}
  )
endif(MAKE_MAIN)


if(MAKE_DAUGHTER)
  add_executable(estimate-entropy
    src/estimate-entropy/main.cpp
  )
  target_compile_options(estimate-entropy PRIVATE
    ${COMMON_WARNINGS}
  )
  target_include_directories(estimate-entropy PRIVATE
    ${COMMON_INCLUDES}
  )
  target_link_libraries(estimate-entropy
    lib_stats
  )
endif(MAKE_DAUGHTER)
