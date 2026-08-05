cmake_minimum_required (VERSION 3.22)

project(expos
  LANGUAGES CXX
  VERSION 0.0.0
)

set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_POLICY_VERSION_MINIMUM 3.5)  # for argparse

option(ENABLE_DEBUG_FLAGS "Enable extra debug compile options" OFF)

option(MAKE_MAIN "Compile main expos binary" ON)
option(MAKE_TEST "Compile tests binary" OFF)
option(MAKE_DAUGHTER "Compile daughter tool binaries" OFF)

# user paths to htslib if needed
set(HTSLIB_INCLUDE_DIR "" CACHE PATH "Path to directory that contains htslib/hts.h")
set(HTSLIB_LIBRARY "" CACHE FILEPATH "Absolute path to libhts.{so,dylib,a}")

# //--- find deps ---//
find_package(PkgConfig REQUIRED)

include(FetchContent)
set(FETCHCONTENT_QUIET OFF)

FetchContent_Declare(
    plog
    GIT_REPOSITORY https://github.com/SergiusTheBest/plog
    GIT_TAG        1.1.10
    GIT_SHALLOW    true
    SYSTEM
)
FetchContent_MakeAvailable(plog)
FetchContent_Declare(
    argparse
    GIT_REPOSITORY https://github.com/p-ranav/argparse.git
    GIT_TAG v3.2
    GIT_SHALLOW TRUE
    SYSTEM
)
FetchContent_MakeAvailable(argparse)
FetchContent_Declare(
    fmt
    GIT_REPOSITORY https://github.com/fmtlib/fmt.git
    GIT_TAG        12.2.0
    GIT_SHALLOW    TRUE
    SYSTEM
)
FetchContent_MakeAvailable(fmt)

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
)

# //--- internal libs ---//
# add_library(lib_stats INTERFACE)

# target_include_directories(lib_stats INTERFACE
#   ${CMAKE_CURRENT_SOURCE_DIR}/src
# )

# if(MAKE_MAIN OR MAKE_TEST)
#   add_library(lib_expos STATIC
#     src/expos/pileup.cpp
#     src/expos/support.cpp
#     src/expos/variant_record.cpp
#     src/expos/cli.cpp
#   )
#   target_compile_options(lib_expos PRIVATE
#     ${COMMON_WARNINGS}
#   )
#   target_include_directories(lib_expos PUBLIC
#     ${CMAKE_CURRENT_SOURCE_DIR}/src
#     ${CMAKE_CURRENT_SOURCE_DIR}/src/expos
#   )
#   target_link_libraries(lib_expos PUBLIC
#     plog::plog
#     cxxopts::cxxopts
#     ${HTSLIB_TARGET}
#   )
# endif()

# //--- testing ---//
if(MAKE_TEST)
  message(STATUS "Fetching Catch2 for test binary")
  include(FetchContent)

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
    tests/test_background_guard.cpp
    src/shared/stats.cpp
    src/expos/variant_stats.cpp
    src/expos/background_guard.cpp
  )
  target_compile_options(test-expos PRIVATE
    ${COMMON_WARNINGS}
  )
  target_include_directories(test-expos PRIVATE
    ${COMMON_INCLUDES}
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
    src/expos/background_guard.cpp
    src/expos/extract_pileup.cpp
    src/expos/pileup_features.cpp
    src/expos/variant_stats.cpp
    src/expos/vcf_record.cpp
    src/hts/hts_types.cpp
    src/shared/stats.cpp
  )
  target_compile_options(expos PRIVATE
    ${COMMON_WARNINGS}
  )
  target_compile_definitions(expos PRIVATE
    EXPOS_VERSION="${PROJECT_VERSION}"
  )
  target_include_directories(expos PRIVATE
    ${COMMON_INCLUDES}
  )
  target_link_libraries(expos PRIVATE
    plog::plog
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
  )
  target_include_directories(estimate-entropy PRIVATE
    ${COMMON_INCLUDES}
  )
  target_link_libraries(estimate-entropy PRIVATE
    argparse
  )
endif(MAKE_DAUGHTER)
