#!/bin/bash
set -e

# Default installation directory
INSTALL_DIR="$HOME/.local/nekrs"

echo "=========================================================="
echo " Starting Automated nekRS Linux Installation (with CUDA)  "
echo "=========================================================="

# Ensure MPI compilers are set
export CC=mpicc
export CXX=mpic++
export FC=mpif77

# OpenMPI's mpif77 wrapper looks for a plain `gfortran` on PATH.
# Only versioned gfortran-9 is installed here, so tell the wrapper which to use.
export OMPI_FC=gfortran-9
export OMPI_F77=gfortran-9

echo "[1/3] Cleaning up previous builds..."
rm -rf build
# Optionally clean previous install. Uncomment if you want to wipe it cleanly every time.
# rm -rf "$INSTALL_DIR"

echo "[2/3] Configuring CMake..."
cmake -S . -B build \
      -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR" \
      -DOCCA_ENABLE_CUDA=ON

echo "[3/3] Building and Installing (using all available CPU cores)..."
# Use all available cores for a faster build (fallback to 4 if cannot detect)
CORES=$(nproc 2>/dev/null || sysctl -n hw.logicalcpu 2>/dev/null || echo 4)
cmake --build ./build --target install -j${CORES}

echo "=========================================================="
echo "✅ Installation Complete!"
echo " nekRS is installed at: $INSTALL_DIR"
echo ""
echo " To use nekRS, please run or add the following to your ~/.bashrc:"
echo "   export NEKRS_HOME=$INSTALL_DIR"
echo "   export PATH=\$NEKRS_HOME/bin:\$PATH"
echo "=========================================================="
