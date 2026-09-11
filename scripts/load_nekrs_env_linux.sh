export NEKRS_HOME=$HOME/.local/nekrs
export PATH=$NEKRS_HOME/bin:$PATH

# nekRS JIT-builds Fortran at runtime. OpenMPI's mpif77/mpif90 wrappers
# look for a plain `gfortran` on PATH; only versioned gfortran-9 is
# installed on this box, so point the wrappers at it.
export OMPI_FC=gfortran-9
export OMPI_F77=gfortran-9

# Force OpenMPI's one-sided communication (OSC) backend to pt2pt.
# nekRS calls MPI_Win_create in the Nek5000 mesh reader; OpenMPI's default
# OSC component (rdma/ucx) fails with `MPI_ERR_WIN: invalid window` on this
# workstation (no proper RDMA/UCX fabric). The pt2pt backend implements
# one-sided ops over regular point-to-point messages and works everywhere.
export OMPI_MCA_osc=pt2pt
