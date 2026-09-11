# For Apple Silicon (MacBook M3) CPU-only setup
export MPICH_HOME=$HOME/.local/mpich-gcc15
export NEKRS_HOME=$HOME/.local/nekrs
export PATH="$NEKRS_HOME/bin:$MPICH_HOME/bin:$PATH"
