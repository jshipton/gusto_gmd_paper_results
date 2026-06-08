"""
Computes error data for Williamson 6 temporal convergence test
"""

import os
from firedrake import Function, errornorm, CheckpointFile
from netCDF4 import Dataset

# -----------------------------------------
# Configuration
# -----------------------------------------

raw_data_dir = '/data/users/alex.brown/data/isambard/w6_convergence/'
base_results_dir = f'{os.path.abspath(os.path.dirname(__file__))}/../results/'
base_results_name = 'w6'

NCELLS = 128
element = "RTCF1"

# Implicit schemes and their dt values
IMPLICIT_SCHEMES = ["backward_euler", "siqn", "imex_ark5", "imex_ssp3"]
IMPLICIT_DTS     = [240.0, 120.0, 60.0, 30.0]
IMPLICIT_DT_STRS = ["240", "120", "60", "30"]

# Explicit schemes and their dt values
EXPLICIT_SCHEMES = ["heun", "ssprk3"]
EXPLICIT_DTS     = [60.0, 30.0, 15.0, 7.5, 3.75, 1.875]
EXPLICIT_DT_STRS = ["60", "30", "15", "7p5", "3p75", "1p875"]


# Colours and markers (one per scheme)
ALL_SCHEMES = IMPLICIT_SCHEMES + EXPLICIT_SCHEMES

# -----------------------------------------
# Path helpers
# -----------------------------------------

def chkpt_path(base, scheme, dt_str):
    dirname = f"williamson6_dt{dt_str}_C{NCELLS}_{scheme}"
    return os.path.join(base, scheme, dirname, "chkpt.h5")

def true_chkpt_path(base):
    return os.path.join(base, "true", "chkpt.h5")

# -----------------------------------------
# Load checkpoint
# -----------------------------------------

def load_chkpt(path):
    with CheckpointFile(path, 'r') as f:
        mesh = f.load_mesh("firedrake_default")
        D    = f.load_function(mesh, "D")
        u    = f.load_function(mesh, "u")
        t    = f.get_attr("/", "time")
        step = f.get_attr("/", "step")
    return D, u, mesh, t, step

# -----------------------------------------
# Compute errors for one scheme
# -----------------------------------------

def compute_errors(base, scheme, dts, dt_strs, D_true, t_true):
    for dt, dt_str in zip(dts, dt_strs):
        path = chkpt_path(base, scheme, dt_str)
        if not os.path.exists(path):
            print(f"  WARNING: not found, skipping: {path}")
            continue
        try:
            # Make directory if it doesn't exist
            results_dir = f'{base_results_dir}/{base_results_name}_{scheme}_{element}_n{NCELLS}_dt{dt_str}'
            os.makedirs(results_dir, exist_ok=True)

            D, _, mesh, t, step = load_chkpt(path)
            print(f"  [{scheme}] dt={dt}: t={t}, step={step}")
            if abs(t - t_true) > 1e-6 * t_true:
                print(f"  WARNING: t={t} does not match true t={t_true}, skipping.")
                continue

            # Calculate error
            D_ref = Function(D.function_space()).interpolate(D_true)
            l2_error = errornorm(D_ref, D, mesh=mesh)

            # Write error to new netcdf file
            error_file = os.path.join(results_dir, 'error_diagnostic.nc')
            with Dataset(error_file, 'w') as ncfile:
                ncfile.createDimension('time', None)
                error_var = ncfile.createVariable('D_error', 'f8', ('time',))
                error_var.units = 'm'
                error_var.long_name = 'L2 error in D against reference solution'
                error_var[0] = l2_error

        except Exception as e:
            print(f"  WARNING: failed to process {path}: {e}")

# -----------------------------------------
# Main
# -----------------------------------------

def main():
    # Load true solution
    true_path = true_chkpt_path(raw_data_dir)
    print(f"Loading true solution: {true_path}")
    D_true, _, _, t_true, step_true = load_chkpt(true_path)
    print(f"  True solution: t={t_true}, step={step_true}")

    # Compute errors for all schemes
    scheme_results = {}

    for scheme in IMPLICIT_SCHEMES:
        print(f"\n=== {scheme} (implicit) ===")
        scheme_results[scheme] = compute_errors(
            raw_data_dir, scheme, IMPLICIT_DTS, IMPLICIT_DT_STRS,
            D_true, t_true
        )

    for scheme in EXPLICIT_SCHEMES:
        print(f"\n=== {scheme} (explicit) ===")
        scheme_results[scheme] = compute_errors(
            raw_data_dir, scheme, EXPLICIT_DTS, EXPLICIT_DT_STRS,
            D_true, t_true
        )

if __name__ == "__main__":
    main()