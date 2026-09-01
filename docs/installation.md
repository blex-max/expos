# Installation

## **Requirements**

- C++23-capable compiler (GCC 13+ or Clang 16+)
- CMake 3.22+
- [htslib](https://github.com/samtools/htslib) ≥ 1.14

[argparse](https://github.com/p-ranav/argparse) and [fmt](https://github.com/fmtlib/fmt) are vendored directly under `vendor/`. Tests additionally fetch [Catch2](https://github.com/catchorg/Catch2) via CMake FetchContent when `MAKE_TEST=ON` - making the test binary therefore requires network access for at least the configure step.

## **Build from Source**

```bash
# from the cloned repo
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build .
./expos --help
```

If htslib is not discoverable via `pkg-config`, point CMake to it directly:

```bash
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DHTSLIB_INCLUDE_DIR=/path/to/htslib/include \
  -DHTSLIB_LIBRARY=/path/to/libhts.so
```

## **Containers**

A Dockerfile is provided with the repository. No user is set; the image produced is therefore compatible with
`singularity pull`.  

Prebuilt docker images for linux/arm64 and linux/amd64 are also provided via the GitHub container registry.
They can be found in the `packages` section of the GitHub page for the repo.

## **Optional: Build Daughter Tool**

To also compile `estimate-entropy` (see [Extras](extras.md)):

```bash
cmake .. -DMAKE_DAUGHTER=ON -DCMAKE_BUILD_TYPE=Release
cmake --build .
```
