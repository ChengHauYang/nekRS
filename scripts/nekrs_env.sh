#!/usr/bin/env bash

# Source this file to configure the GNU/MPICH CPU-only nekRS environment.
# Environment variables set before sourcing this file take precedence.

if [[ -n "${ZSH_VERSION:-}" ]]; then
  _nekrs_env_script="${(%):-%x}"
  case "${ZSH_EVAL_CONTEXT}" in
    *:file) ;;
    *)
      printf '%s\n' "This script must be sourced: source ${_nekrs_env_script}" >&2
      exit 1
      ;;
  esac
elif [[ -n "${BASH_VERSION:-}" ]]; then
  _nekrs_env_script="${BASH_SOURCE[0]}"
  if [[ "${_nekrs_env_script}" == "$0" ]]; then
    printf '%s\n' "This script must be sourced: source ${_nekrs_env_script}" >&2
    exit 1
  fi
else
  printf '%s\n' "nekrs_env.sh requires bash or zsh." >&2
  return 1 2>/dev/null || exit 1
fi

_nekrs_env_script_dir="$(cd -- "$(dirname -- "${_nekrs_env_script}")" && pwd)"
_nekrs_env_temp_root="/var/folders/9l/dhxn3q9n12z63pgxgf5nttg40000gq/T/opencode"

export NEKRS_SOURCE_ROOT="${NEKRS_SOURCE_ROOT:-$(cd -- "${_nekrs_env_script_dir}/.." && pwd)}"
export NEKRS_BUILD_DIR="${NEKRS_BUILD_DIR:-${NEKRS_SOURCE_ROOT}/build-m3}"

if [[ -z "${MPICH_HOME:-}" ]]; then
  if [[ -x "${HOME}/.local/mpich-gcc15/bin/mpicc" ]]; then
    export MPICH_HOME="${HOME}/.local/mpich-gcc15"
  else
    export MPICH_HOME="${_nekrs_env_temp_root}/mpich-gcc15"
    printf 'nekRS environment: using temporary MPICH installation: %s\n' "${MPICH_HOME}" >&2
  fi
fi

export NEKRS_HOME="${NEKRS_HOME:-${HOME}/.local/nekrs}"

export CC="${CC:-${MPICH_HOME}/bin/mpicc}"
export CXX="${CXX:-${MPICH_HOME}/bin/mpicxx}"
export FC="${FC:-${MPICH_HOME}/bin/mpifort}"
export OCCA_CXX="${OCCA_CXX:-/opt/homebrew/bin/g++-15}"

# Force OCCA's SERIAL backend on Apple Silicon.
export CPUONLY="${CPUONLY:-1}"

case ":${PATH}:" in
  *":${NEKRS_HOME}/bin:"*) ;;
  *) export PATH="${NEKRS_HOME}/bin:${PATH}" ;;
esac

case ":${PATH}:" in
  *":${MPICH_HOME}/bin:"*) ;;
  *) export PATH="${MPICH_HOME}/bin:${PATH}" ;;
esac

_nekrs_env_missing=0
for _nekrs_env_program in "${CC}" "${CXX}" "${FC}" "${OCCA_CXX}"; do
  if [[ ! -x "${_nekrs_env_program}" ]]; then
    printf 'nekRS environment: required compiler not found: %s\n' \
      "${_nekrs_env_program}" >&2
    _nekrs_env_missing=1
  fi
done

if [[ -x "${NEKRS_HOME}/bin/nekrs" ]]; then
  printf 'nekRS environment ready: %s\n' "${NEKRS_HOME}/bin/nekrs"
elif [[ "${_nekrs_env_missing}" -eq 0 ]]; then
  printf 'nekRS environment ready for compilation. nekRS is not installed at %s yet.\n' \
    "${NEKRS_HOME}" >&2
  printf 'Build with: cmake -S "%s" -B "%s" [options], then cmake --build "%s" --target install\n' \
    "${NEKRS_SOURCE_ROOT}" "${NEKRS_BUILD_DIR}" "${NEKRS_BUILD_DIR}" >&2
fi

unset _nekrs_env_program _nekrs_env_missing _nekrs_env_script_dir
unset _nekrs_env_script _nekrs_env_temp_root
