# TGMM 2.0 Fork

This fork is confirmed to build and run properly with Ubuntu 20.04 LTS and 
NVIDIA CUDA toolkit 11.1.  There are improvements that allow user-set
temporal window for logical rules, and bug fixes to the 3D Haar GPU division
classifier.


## Build Requirements

* [CMake](https://cmake.org/) builds with v3.15.3
* [CUDA](https://developer.nvidia.com/cuda-downloads) install CUDA toolkit 11.1
* [Git](https://git-scm.com/) to download the software


## Quick Start

Download the source code:

```sh
git clone https://github.com/mhdominguez/tgmm-fork.git
cd tgmm-fork
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

Build...

```sh
cmake --build . --config Release --target install

```

And install...

```sh
sudo mkdir /opt/tgmm
sudo cp -rf ../install/* /opt/tgmm

```







## Resources

[Original TGMM Repository](https://bitbucket.org/fernandoamat/tgmm-paper)

## References 

McDole K, Guignard L, Amat F, Berger A, Malandain G, Royer LA, Turaga SC, Branson K, Keller PJ
Cell. 2018 Oct 10;175(3):859-876. doi: [10.1016/j.cell.2018.09.031](http://doi.org/10.1016/j.cell.2018.09.031)

Amat F, Höckendorf B, Wan Y, Lemon WC, McDole K, Keller PJ
Nature Protocols. 2015 Oct 2;10(11):1679-96. doi: [10.1038/nprot.2015.111](http://doi.org/10.1038/nprot.2015.111)

Amat F, Lemon W, Mossing DP, McDole K, Wan Y, Branson K, Myers EW, Keller PJ.
Nature Methods. 2014 Jul 20;11(9):951-8. doi: [10.1038/nmeth.3036](http://doi.org/10.1038/nmeth.3036)
