# Build

## Requirements

* [CMake](https://cmake.org/)
* C++11 compiler toolsuite
* [CUDA](https://developer.nvidia.com/cuda-downloads) is required for the
`TGMM` executable, but not `ProcessStack_woGPU`.
* (Optionally) [Git](https://git-scm.com/) to download the software.

## Quick Start

Download the source code:

```sh
git clone https://nclack@bitbucket.org/fernandoamat/tgmm-paper.git
cd tgmm-paper
git submodule update --init --recursive
```

Make a build folder and configure the build.  Here, we're telling to put the
final build products into an "install" folder located next to the build 
folder.

```sh
mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=`pwd`/../install ..
```

And build!

```sh
cmake --build . --config Release --target INSTALL
```