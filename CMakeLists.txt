cmake_minimum_required (VERSION 3.22)

project(expos
  LANGUAGES CXX
  VERSION 0.0.0
)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE)
  message(STATUS "No CMAKE_BUILD_TYPE specified; defaulting to RelWithDebInfo")
  set(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING "Build type" FORCE)
endif()

option(ENABLE_DEBUG_FLAGS "Enable sanitizers and full debug info" OFF)
option(ENABLE_VEC_REPORT "Report loop vectorisation decisions at compile time" OFF)

option(MAKE_MAIN "Compile main expos binary" ON)
option(MAKE_TEST "Compile tests binary" OFF)
option(MAKE_DAUGHTER "Compile daughter tool binaries" OFF)

# user paths to htslib if needed
set(HTSLIB_INCLUDE_DIR "" CACHE PATH "Path to directory that contains htslib/hts.h")
set(HTSLIB_LIBRARY "" CACHE FILEPATH "Absolute path to libhts.{so,dylib,a}")

# //--- find deps ---//
find_package(PkgConfig REQUIRED)

add_library(argparse INTERFACE)
target_include_directories(argparse SYSTEM INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}/vendor/argparse/include
)

add_library(fmt INTERFACE)
target_include_directories(fmt SYSTEM INTERFACE
  ${CMAKE_CURRENT_SOURCE_DIR}/vendor/fmt/include
)
target_compile_definitions(fmt INTERFACE FMT_HEADER_ONLY)
add_library(fmt::fmt-header-only ALIAS fmt)

set(HTSLIB_TARGET "")
if(HTSLIB_INCLUDE_DIR AND HTSLIB_LIBRARY)
  message(STATUS "Using user-specified htslib install. Must be >=1.14, but this is not enforced by the installer.")
  add_library(htslib UNKNOWN IMPORTED GLOBAL)
  add_library(htslib::htslib ALIAS htslib)
  set_target_properties(htslib PROPERTIES
    IMPORTED_LOCATION "${HTSLIB_LIBRARY}"
    INTERFACE_INCLUDE_DIRECTORIES "${HTSLIB_INCLUDE_DIR}"
  )
  set(HTSLIB_TARGET htslib::htslib)
endif()

if(NOT HTSLIB_TARGET)
  message(STATUS "Looking for htslib via PkgConfig")
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
  -Wimplicit-fallthrough
)

# Honour `#pragma omp simd`
set(COMMON_SIMD -fopenmp-simd)

set(COMMON_VEC_REPORT "")
if(ENABLE_VEC_REPORT)
  if(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
    set(COMMON_VEC_REPORT
      -Rpass=loop-vectorize
      -Rpass-missed=loop-vectorize
      -Rpass-analysis=loop-vectorize
    )
  elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(COMMON_VEC_REPORT
      -fopt-info-vec
      -fopt-info-vec-missed
    )
  else()
    message(WARNING "ENABLE_VEC_REPORT: no known flags for ${CMAKE_CXX_COMPILER_ID}")
  endif()
  message(STATUS "ENABLE_VEC_REPORT: reporting vectorisation decisions")
endif()

set(COMMON_DEBUG_FLAGS "")
if(ENABLE_DEBUG_FLAGS)
  set(COMMON_DEBUG_FLAGS
    -g3 -fno-omit-frame-pointer
    -fsanitize=address,undefined
    -fno-sanitize-recover=all
    -UNDEBUG  # optimised build types define NDEBUG, stripping assert()
  )
  message(STATUS "ENABLE_DEBUG_FLAGS: building with ASan + UBSan")
endif()

# //--- testing ---//
if(MAKE_TEST)
  message(STATUS "Fetching Catch2 for test binary")
  include(FetchContent)
  set(FETCHCONTENT_QUIET OFF)

  FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.8.1
    GIT_SHALLOW    TRUE
    SYSTEM
  )
  FetchContent_MakeAvailable(Catch2)

  enable_testing()

  add_executable(test-expos
    tests/test_stats.cpp
    tests/test_variant_stats.cpp
    tests/test_guards.cpp
    tests/test_skip.cpp
    tests/test_rng.cpp
    src/shared/stats.cpp
    src/expos/compute_info_field.cpp
    src/expos/encode_info_field.cpp
    src/expos/skip.cpp
    src/expos/guards.cpp
  )
  # Keep assert() preconditions live regardless of build type -- every
  # optimised config defines NDEBUG, which would otherwise strip them.
  target_compile_options(test-expos PRIVATE
    ${COMMON_WARNINGS}
    ${COMMON_SIMD}
    ${COMMON_VEC_REPORT}
    ${COMMON_DEBUG_FLAGS}
    -UNDEBUG
  )
  target_link_options(test-expos PRIVATE
    ${COMMON_DEBUG_FLAGS}
  )
  target_include_directories(test-expos PRIVATE
    ${COMMON_INCLUDES}
    ${CMAKE_CURRENT_SOURCE_DIR}/tests
  )
  target_link_libraries(test-expos PRIVATE
    Catch2::Catch2WithMain
    ${HTSLIB_TARGET}
  )

  add_test(NAME test-expos COMMAND test-expos)
endif(MAKE_TEST)
# //--- end testing ---//


if(MAKE_MAIN)
  add_executable(expos
    src/expos/main.cpp
    src/expos/annotate.cpp
    src/expos/guards.cpp
    src/expos/extract_pileup.cpp
    src/expos/pileup_features.cpp
    src/expos/compute_info_field.cpp
    src/expos/encode_info_field.cpp
    src/expos/skip.cpp
    src/expos/vcf_record.cpp
    src/hts/hts_types.cpp
    src/shared/stats.cpp
  )
  target_compile_options(expos PRIVATE
    ${COMMON_WARNINGS}
    ${COMMON_SIMD}
    ${COMMON_VEC_REPORT}
    ${COMMON_DEBUG_FLAGS}
  )
  target_link_options(expos PRIVATE
    ${COMMON_DEBUG_FLAGS}
  )
  target_compile_definitions(expos PRIVATE
    EXPOS_VERSION="${PROJECT_VERSION}"
  )
  target_include_directories(expos PRIVATE
    ${COMMON_INCLUDES}
  )
  target_link_libraries(expos PRIVATE
    argparse
    fmt::fmt-header-only
    ${HTSLIB_TARGET}
  )
endif(MAKE_MAIN)


if(MAKE_DAUGHTER)
  add_executable(estimate-entropy
    src/estimate-entropy/main.cpp
    src/shared/stats.cpp
  )
  target_compile_options(estimate-entropy PRIVATE
    ${COMMON_WARNINGS}
    ${COMMON_SIMD}
    ${COMMON_VEC_REPORT}
    ${COMMON_DEBUG_FLAGS}
  )
  target_link_options(estimate-entropy PRIVATE
    ${COMMON_DEBUG_FLAGS}
  )
  target_include_directories(estimate-entropy PRIVATE
    ${COMMON_INCLUDES}
  )
  target_link_libraries(estimate-entropy PRIVATE
    argparse
  )
endif(MAKE_DAUGHTER)
