# Installation

## **Requirements**

- C++23-capable compiler (GCC 13+ or Clang 16+)
- CMake 3.22+
- [htslib](https://github.com/samtools/htslib) ≥ 1.14

[cxxopts](https://github.com/jarro2783/cxxopts) is fetched automatically by CMake — no manual install needed.

## **Build from Source**

```bash
# from the cloned repo
mkdir build
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build .
./expos --help
```

## **Non-Standard htslib Location**

If htslib is not discoverable via `pkg-config`, point CMake to it directly:

```bash
cmake .. \
  -DCMAKE_BUILD_TYPE=Release \
  -DHTSLIB_INCLUDE_DIR=/path/to/htslib/include \
  -DHTSLIB_LIBRARY=/path/to/libhts.so
```

## **Optional: Build Daughter Tool**

To also compile `estimate-entropy` (see [Extras](extras.md)):

```bash
cmake .. -DMAKE_DAUGHTER=ON -DCMAKE_BUILD_TYPE=Release
cmake --build .
```
