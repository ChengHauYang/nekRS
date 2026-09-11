#!/usr/bin/env python3
"""Convergence study driver for the staticIBM example.

Automates multi-parameter convergence sweeps for the nekRS immersed boundary
method implementation.  Each sweep varies a single parameter while holding all
others at their default values from staticIBM.par.

Supported sweeps
----------------
- polynomialOrder : N = 3, 5, 7, 9
- dt              : 4e-3, 2e-3, 1e-3, 5e-4
- markerSpacing   : 0.80, 0.40, 0.20, 0.10
- gaussianWidth   : 0.50, 0.25, 0.125, 0.0625
- supportRadius   : 1.0, 0.50, 0.25, 0.125
- shellThickness  : 0.20, 0.10, 0.05, 0.025

Usage
-----
    python convergence_study.py --nekrs-home /path/to/nekRS --np 4
    python convergence_study.py --parameter markerSpacing --dry-run
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
from dataclasses import dataclass, fields
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

CASE_NAME = "staticIBM"
CASE_DIR = Path(__file__).resolve().parent

# Default parameter values matching staticIBM.par
DEFAULTS: Dict[str, float] = {
    "polynomialOrder": 5,
    "dt": 2e-3,
    "markerSpacing": 0.40,
    "gaussianWidth": 0.25,
    "supportRadius": 0.50,
    "shellThickness": 0.05,
}

# Sweep values for each parameter
SWEEPS: Dict[str, List[float]] = {
    "polynomialOrder": [3, 5, 7, 9],
    "dt": [4e-3, 2e-3, 1e-3, 5e-4],
    "markerSpacing": [0.80, 0.40, 0.20, 0.10],
    "gaussianWidth": [0.50, 0.25, 0.125, 0.0625],
    "supportRadius": [1.0, 0.50, 0.25, 0.125],
    "shellThickness": [0.20, 0.10, 0.05, 0.025],
}

# Mapping from Python parameter names to .par file keys/sections
PAR_KEY_MAP: Dict[str, Tuple[str, str]] = {
    "polynomialOrder": ("GENERAL", "polynomialOrder"),
    "dt": ("GENERAL", "dt"),
    "markerSpacing": ("CASEDATA", "marker_spacing"),
    "gaussianWidth": ("CASEDATA", "gaussian_width"),
    "supportRadius": ("CASEDATA", "support_radius"),
    "shellThickness": ("CASEDATA", "shell_thickness"),
}


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------

@dataclass
class RunResult:
    """Metrics extracted from a single nekRS run."""

    parameter_value: float = 0.0
    marker_count: int = 0
    max_slip: float = float("nan")
    rms_slip: float = float("nan")
    force_x: float = float("nan")
    force_y: float = float("nan")
    force_z: float = float("nan")


# ---------------------------------------------------------------------------
# .par file generation
# ---------------------------------------------------------------------------

def generate_par_content(overrides: Dict[str, float]) -> str:
    """Return the full text of a staticIBM.par with *overrides* applied.

    Parameters
    ----------
    overrides : dict
        Maps Python parameter names (e.g. ``"markerSpacing"``) to their
        desired numeric values.  Only the keys present in *overrides* will
        differ from the defaults.
    """
    poly_order = int(overrides.get("polynomialOrder", DEFAULTS["polynomialOrder"]))
    dt = overrides.get("dt", DEFAULTS["dt"])
    marker_spacing = overrides.get("markerSpacing", DEFAULTS["markerSpacing"])
    gaussian_width = overrides.get("gaussianWidth", DEFAULTS["gaussianWidth"])
    support_radius = overrides.get("supportRadius", DEFAULTS["supportRadius"])
    shell_thickness = overrides.get("shellThickness", DEFAULTS["shellThickness"])

    return f"""\
userSections = CASEDATA

[GENERAL]
polynomialOrder = {poly_order}
stopAt = numSteps
numSteps = 1
dt = {dt:.6g}
timeStepper = tombo3
checkpointInterval = 0

[PROBLEMTYPE]
equation = stokes

[FLUID PRESSURE]
residualTol = 1e-08
solver = cg+flexible
smootherType = Jac+Cheby

[FLUID VELOCITY]
boundaryTypeMap = udfDirichlet
residualTol = 1e-10
rho = 1.0
viscosity = 0.02

[CASEDATA]
stl_file = geometry/reference.stl
marker_spacing = {marker_spacing:.8g}
shell_thickness = {shell_thickness:.8g}
support_radius = {support_radius:.8g}
gaussian_width = {gaussian_width:.8g}
max_subdivision_depth = 30
periodic_dimensions = none
"""


# ---------------------------------------------------------------------------
# Output parsing
# ---------------------------------------------------------------------------

# Pattern for the IBM geometry summary printed during setup:
#   IBM geometry: triangles=... markers=... ...
_RE_GEOMETRY = re.compile(
    r"IBM geometry:.*?markers=(\d+)"
)

# Pattern for the diagnostics line printed by UDF_ExecuteStep.
_RE_DIAGNOSTICS = re.compile(
    r"slipMax=([\d.eE+\-]+)\s+slipRms=([\d.eE+\-]+)\s+"
    r"force=\(([\d.eE+\-]+),\s*([\d.eE+\-]+),\s*([\d.eE+\-]+)\)"
)

# Divergence metric (if printed)
_RE_DIVERGENCE = re.compile(r"div[Ee]rr\s*[:=]\s*([\d.eE+\-]+)")


def parse_stdout(stdout: str) -> RunResult:
    """Extract metrics from nekRS stdout."""
    result = RunResult()

    # Marker count
    match = _RE_GEOMETRY.search(stdout)
    if match:
        result.marker_count = int(match.group(1))

    # Diagnostics (last occurrence wins – captures the final timestep)
    for match in _RE_DIAGNOSTICS.finditer(stdout):
        result.max_slip = float(match.group(1))
        result.rms_slip = float(match.group(2))
        result.force_x = float(match.group(3))
        result.force_y = float(match.group(4))
        result.force_z = float(match.group(5))

    return result


# ---------------------------------------------------------------------------
# Runner
# ---------------------------------------------------------------------------

def run_case(
    nekrs_home: Path,
    np: int,
    overrides: Dict[str, float],
    dry_run: bool = False,
) -> RunResult:
    """Set up a temporary directory, run nekRS, and return parsed results.

    Parameters
    ----------
    nekrs_home : Path
        Root of the nekRS installation (contains ``bin/nrsmpi``).
    np : int
        Number of MPI ranks.
    overrides : dict
        Parameter overrides for the .par file.
    dry_run : bool
        If *True*, print the command instead of executing it.
    """
    nrsmpi = nekrs_home / "bin" / "nrsmpi"
    if not dry_run and not nrsmpi.exists():
        raise FileNotFoundError(f"nrsmpi not found at {nrsmpi}")

    work_dir = Path(tempfile.mkdtemp(prefix=f"{CASE_NAME}_conv_"))
    par_path = work_dir / f"{CASE_NAME}.par"
    par_path.write_text(generate_par_content(overrides))

    # Symlink required case files into the working directory. IBM sources are
    # provided by the nekRS installation include and kernel paths.
    for name in [f"{CASE_NAME}.re2", f"{CASE_NAME}.udf", "geometry"]:
        src = CASE_DIR / name
        dst = work_dir / name
        if not src.exists():
            raise FileNotFoundError(f"required case input not found: {src}")
        dst.symlink_to(src)

    cmd = [str(nrsmpi), str(par_path.stem), str(np)]
    print(f"  cmd: {' '.join(cmd)}  (cwd={work_dir})")

    if dry_run:
        shutil.rmtree(work_dir, ignore_errors=True)
        return RunResult()

    try:
        proc = subprocess.run(
            cmd,
            cwd=str(work_dir),
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            timeout=3600,
        )
        result = parse_stdout(proc.stdout)

        # Dump raw log before validating the run so failures remain inspectable.
        log_path = work_dir / "stdout.log"
        log_path.write_text(proc.stdout)
        print(f"  log: {log_path}")

        if proc.returncode != 0:
            raise RuntimeError(
                f"nekRS exited with code {proc.returncode}; see {log_path}"
            )
        metrics = (result.max_slip, result.rms_slip,
                   result.force_x, result.force_y, result.force_z)
        if result.marker_count <= 0 or not all(math.isfinite(value) for value in metrics):
            raise RuntimeError(f"nekRS output is missing IBM diagnostics; see {log_path}")
    except subprocess.TimeoutExpired as error:
        raise RuntimeError("nekRS timed out after 3600 s") from error
    finally:
        # Keep logs but remove heavy checkpoint/field data
        for pattern in ("*.fld", "*.f0*"):
            for fld in work_dir.glob(pattern):
                fld.unlink(missing_ok=True)

    return result


# ---------------------------------------------------------------------------
# Sweep driver
# ---------------------------------------------------------------------------

def run_sweep(
    parameter: str,
    nekrs_home: Path,
    np: int,
    output_dir: Path,
    dry_run: bool = False,
) -> List[RunResult]:
    """Run a convergence sweep over *parameter*.

    Returns a list of :class:`RunResult` objects, one per sweep value.
    """
    values = SWEEPS[parameter]
    results: List[RunResult] = []

    print(f"\n{'=' * 60}")
    print(f"Sweep: {parameter}  values={values}")
    print(f"{'=' * 60}")

    for value in values:
        overrides: Dict[str, float] = {parameter: value}
        print(f"\n--- {parameter} = {value} ---")
        result = run_case(nekrs_home, np, overrides, dry_run=dry_run)
        result.parameter_value = value
        results.append(result)

    # Save CSV
    output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = output_dir / f"sweep_{parameter}.csv"
    _write_csv(csv_path, results)
    print(f"\nResults written to {csv_path}")

    # Print summary table
    _print_table(parameter, results)

    return results


def _write_csv(path: Path, results: Sequence[RunResult]) -> None:
    """Write *results* to a CSV file at *path*."""
    fieldnames = [f.name for f in fields(RunResult)]
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        for r in results:
            row = {f.name: getattr(r, f.name) for f in fields(RunResult)}
            writer.writerow(row)


def _print_table(parameter: str, results: Sequence[RunResult]) -> None:
    """Print a human-readable summary table to stdout."""
    header = (
        f"{'value':>12s} {'markers':>8s} {'max_slip':>12s} "
        f"{'rms_slip':>12s} {'force_x':>12s} {'force_y':>12s} {'force_z':>12s}"
    )
    print(f"\nSummary: {parameter}")
    print(header)
    print("-" * len(header))
    for r in results:
        print(
            f"{r.parameter_value:12.6g} {r.marker_count:8d} "
            f"{r.max_slip:12.4e} {r.rms_slip:12.4e} "
            f"{r.force_x:12.4e} {r.force_y:12.4e} {r.force_z:12.4e}"
        )


# ---------------------------------------------------------------------------
# Plotting (optional)
# ---------------------------------------------------------------------------

def generate_plots(
    parameter: str,
    results: Sequence[RunResult],
    output_dir: Path,
) -> None:
    """Generate log-log convergence plots if matplotlib is available."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  matplotlib not available – skipping plots")
        return

    values = [r.parameter_value for r in results]
    metrics = {
        "max_slip": [r.max_slip for r in results],
        "rms_slip": [r.rms_slip for r in results],
        "force_x": [abs(r.force_x) for r in results],
        "force_y": [abs(r.force_y) for r in results],
        "force_z": [abs(r.force_z) for r in results],
    }

    fig, axes = plt.subplots(1, len(metrics), figsize=(4 * len(metrics), 4),
                             squeeze=False)
    axes = axes.ravel()

    for ax, (metric_name, metric_values) in zip(axes, metrics.items()):
        # Filter out NaN / zero values for log-log
        pairs = [
            (v, m) for v, m in zip(values, metric_values)
            if m == m and m > 0  # NaN check: NaN != NaN
        ]
        if not pairs:
            ax.set_title(metric_name)
            continue
        xp, yp = zip(*pairs)
        ax.loglog(xp, yp, "o-", linewidth=1.5, markersize=6)
        ax.set_xlabel(parameter)
        ax.set_ylabel(metric_name)
        ax.set_title(metric_name)
        ax.grid(True, which="both", linestyle="--", alpha=0.5)

    fig.suptitle(f"Convergence: {parameter}", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.95])

    plot_path = output_dir / f"convergence_{parameter}.png"
    fig.savefig(str(plot_path), dpi=150)
    plt.close(fig)
    print(f"  Plot saved to {plot_path}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Run IBM convergence sweeps for staticIBM.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    parser.add_argument(
        "--nekrs-home",
        type=Path,
        default=Path(os.environ.get("NEKRS_HOME", "")),
        help="Path to nekRS installation (default: $NEKRS_HOME)",
    )
    parser.add_argument(
        "--np",
        type=int,
        default=1,
        help="Number of MPI ranks (default: 1)",
    )
    parser.add_argument(
        "--parameter",
        choices=list(SWEEPS.keys()) + ["all"],
        default="all",
        help="Which parameter to sweep (default: all)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("convergence_results"),
        help="Directory for output files (default: convergence_results)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without executing",
    )
    return parser.parse_args(argv)


def main(argv: Optional[Sequence[str]] = None) -> None:
    """Entry point."""
    args = parse_args(argv)

    if not args.dry_run and not args.nekrs_home:
        print("ERROR: --nekrs-home is required (or set $NEKRS_HOME)",
              file=sys.stderr)
        sys.exit(1)

    parameters = list(SWEEPS.keys()) if args.parameter == "all" else [args.parameter]

    all_results: Dict[str, List[RunResult]] = {}
    for param in parameters:
        results = run_sweep(
            parameter=param,
            nekrs_home=args.nekrs_home,
            np=args.np,
            output_dir=args.output_dir,
            dry_run=args.dry_run,
        )
        all_results[param] = results

        if not args.dry_run:
            generate_plots(param, results, args.output_dir)

    print(f"\n{'=' * 60}")
    print("All sweeps complete.")
    print(f"CSV files and plots saved in: {args.output_dir.resolve()}")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
