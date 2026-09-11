```
███    ██ ███████ ██   ██ ██████  ███████
████   ██ ██      ██  ██  ██   ██ ██     
██ ██  ██ █████   █████   ██████  ███████
██  ██ ██ ██      ██  ██  ██   ██      ██
██   ████ ███████ ██   ██ ██   ██ ███████ 
(c) 2019-2025 UCHICAGO ARGONNE, LLC
```

[![Build Status](https://travis-ci.com/Nek5000/nekRS.svg?branch=master)](https://travis-ci.com/Nek5000/nekRS)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-orange.svg)](https://opensource.org/licenses/BSD-3-Clause)

**nekRS** is a fast and scalable computational fluid dynamics (CFD) solver targeting HPC applications. The code started as an early fork of [libParanumal](https://github.com/paranumal/libparanumal) in 2019, with the intention of supplanting its precursor, [nek5000](https://github.com/Nek5000/Nek5000).

Capabilities:

* Incompressible and low Mach-number Navier-Stokes + scalar transport 
* High-order curvilinear conformal Hex spectral elements in space 
* Variable time step 2nd/3rd order semi-implicit time integration
* MPI + [OCCA](https://github.com/libocca/occa) supporting CUDA, HIP, DPC++, SERIAL (C++)
* LES and RANS turbulence models
* Arbitrary-Lagrangian-Eulerian moving mesh
* Lagrangian phase model
* Overlapping overset grids
* Conjugate fluid-solid heat transfer
* Various boundary conditions
* VisIt & Paraview for data analysis and visualization including in-situ support through Ascent
* Legacy interface

## Build Instructions

Requirements:
* Linux, Mac OS X (Microsoft WSL and Windows is not supported) 
* GNU/oneAPI/NVHPC/ROCm compilers (C++17/C99 compatible)
* MPI-3.1 or later
* CMake version 3.21 or later 

Download the latest [release](https://github.com/Nek5000/nekRS/releases) available under

```sh
https://github.com/Nek5000/nekRS/archive/refs/tags/v26.0.tar.gz 
```

or clone our GitHub repository:

```sh
https://github.com/Nek5000/nekRS.git
```
The [master](https://github.com/Nek5000/nekRS) branch always points to the latest stable release while [next](https://github.com/Nek5000/nekRS/tree/next) 
provides an early preview of the next upcoming release (do not use in a production environment).

#
If you're on an HPC system, ensure you log in to a compute node. You can find installation instructions and job submission scripts for common HPC systems [here](https://github.com/Nek5000/nekRS_HPCsupport).
Now, just run:

```sh
CC=mpicc CXX=mpic++ FC=mpif77 ./build.sh [-DCMAKE_INSTALL_PREFIX=$HOME/.local/nekrs] [<options>]
```
Adjust the compilers as necessary. Make sure to remove the previous build and installation directory if updating.

### Apple Silicon (MacBook M3)

AppleClang is not currently supported on arm64. The following CPU-only setup uses
Homebrew GCC 15 and builds MPICH with the same GNU toolchain. Keeping MPI and
nekRS on one compiler toolchain avoids mixing AppleClang and GNU objects.

Install the prerequisites (Homebrew uses `/opt/homebrew` on Apple Silicon):

```sh
xcode-select --install
brew install cmake gcc wget
```

Download MPICH 4.3.2, then configure and install it in a user-owned prefix:

```sh
mkdir -p "$HOME/src"
cd "$HOME/src"
wget https://www.mpich.org/static/downloads/4.3.2/mpich-4.3.2.tar.gz
tar -xzf mpich-4.3.2.tar.gz

export MPICH_SRC=$HOME/src/mpich-4.3.2
export MPICH_HOME=$HOME/.local/mpich-gcc15

cd "$MPICH_SRC"
CC=/opt/homebrew/bin/gcc-15 \
CXX=/opt/homebrew/bin/g++-15 \
FC=/opt/homebrew/bin/gfortran-15 \
F77=/opt/homebrew/bin/gfortran-15 \
./configure --prefix="$MPICH_HOME"
make -j$(sysctl -n hw.logicalcpu)
make install
export PATH="$MPICH_HOME/bin:$PATH"
```

From the nekRS source directory, build the SERIAL C++ backend in a separate
build directory:

```sh
export NEKRS_HOME=$HOME/.local/nekrs

cmake -S . -B build-m3 \
  -DCMAKE_BUILD_TYPE=RelWithDebInfo \
  -DCMAKE_INSTALL_PREFIX="$NEKRS_HOME" \
  -DCMAKE_C_COMPILER="$MPICH_HOME/bin/mpicc" \
  -DCMAKE_CXX_COMPILER="$MPICH_HOME/bin/mpicxx" \
  -DCMAKE_Fortran_COMPILER="$MPICH_HOME/bin/mpifort" \
  -DOCCA_CXX=/opt/homebrew/bin/g++-15 \
  -DOCCA_ENABLE_CUDA=OFF \
  -DOCCA_ENABLE_HIP=OFF \
  -DOCCA_ENABLE_DPCPP=OFF \
  -DOCCA_ENABLE_OPENCL=OFF \
  -DOCCA_ENABLE_METAL=OFF \
  -DENABLE_HYPRE_GPU=OFF \
  -DENABLE_ADIOS=OFF
cmake --build build-m3 --target install -j$(sysctl -n hw.logicalcpu)
```

Add nekRS and the matching MPICH installation to the shell environment:

```sh
export MPICH_HOME=$HOME/.local/mpich-gcc15
export NEKRS_HOME=$HOME/.local/nekrs
export PATH="$NEKRS_HOME/bin:$MPICH_HOME/bin:$PATH"
```

Set `CPUONLY=1` when running a case so OCCA selects the SERIAL backend:

```sh
cp -a "$NEKRS_HOME/examples/ethier" "$HOME/ethier-test"
cd "$HOME/ethier-test"
CPUONLY=1 mpirun -np 1 nekrs --setup ethier.par
```

## Setting the Environment

Assuming you run `bash` and your install directory is $HOME/.local/nekrs, 
add the following line to your $HOME/.bash_profile:

```sh
export NEKRS_HOME=$HOME/.local/nekrs
export PATH=$NEKRS_HOME/bin:$PATH
```
then type `source $HOME/.bash_profile` in the current terminal window. 

## Run the Code

We try hard not to break userland but the code is evolving quickly so things might change from one version to another without being backward compatible. Please consult `RELEASE.md` *before* using the code.  

```sh
cd <directory outside of installation/source folder>
cp -a $NEKRS_HOME/examples/turbPipePeriodic .
mpirun -np 2 nekrs --setup turbPipe.par
```
For convenience we provide various launch scripts in the `bin` directory.

## Documentation 
For documentation, see [readthedocs page](https://nekrs.readthedocs.io/en/latest/). 
The manual pages for the `par` file and environment variables `env` can be accessed through `nrsman`

## Discussion Group
Please visit [GitHub Discussions](https://github.com/Nek5000/nekRS/discussions). Here we help, find solutions, share ideas, and follow discussions.

## Contributing
Our project is hosted on [GitHub](https://github.com/Nek5000/nekRS). To learn how to contribute, see `CONTRIBUTING.md`.

## Reporting Bugs
All bugs are reported and tracked through [Issues](https://github.com/Nek5000/nekRS/issues). If you are having trouble installing the code or getting your case to run properly, you should first visit our discussion group.

## License
nekRS is released under the BSD 3-clause license (see `LICENSE` file). 
All new contributions must be made under the BSD 3-clause license.

## Citing
If you find our project useful, please cite [NekRS, a GPU-Accelerated Spectral Element Navier-Stokes Solver](https://www.sciencedirect.com/science/article/abs/pii/S0167819122000710) 

## Acknowledgment
This research was supported by the Exascale Computing Project (17-SC-20-SC), 
a joint project of the U.S. Department of Energy's Office of Science and National Nuclear Security 
Administration, responsible for delivering a capable exascale ecosystem, including software, 
applications, and hardware technology, to support the nation's exascale computing imperative.
