"""
A linearised form of Test Case 2 (solid-body rotation) of Williamson et al 1992:
``A standard test set for numerical approximations to the shallow water
equations in spherical geometry'', JCP, with a perturbation to trigger a series
of gravity waves.
"""

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from firedrake import (
  Function, SpatialCoordinate, as_vector, pi, exp, errornorm, CheckpointFile
)
from firedrake.fml import keep
from gusto import (
    Domain, IO, OutputParameters, SemiImplicitQuasiNewton, DefaultTransport,
    ForwardEuler, ZonalComponent, MeridionalComponent, ShallowWaterParameters,
    LinearShallowWaterEquations, GeneralIcosahedralSphereMesh, Difference,
    GeneralCubedSphereMesh, SSPRK3, RK4, BDF2, Leapfrog, IMEX_Euler,
    Timestepper, TrapeziumRule, great_arc_angle, lonlatr_from_xyz
)
from gusto.core.labels import (
    explicit, implicit, pressure_gradient, transport, coriolis
)
from netCDF4 import Dataset
from os import path

linear_sw_gw_defaults = {
    'ncells_per_edge': 16,        # number of cells per icosahedron edge
    'dt': 900.0,                  # 15 minutes
    'tmax': 3.*60.*60.,           # 3 hours
    'element_family': 'BDM',      # family for velocity space
    'element_order': 1,           # order of the elements
    'dumpfreq': 12,               # once with default options
    'stepper': 'trapezium_rule',  # default is TrapeziumRule
    'reference_error': False,     # whether to measure against reference solution
    'reference_dirname': None,    # name of directory with reference solution
    'dirname': 'linear_sw_gw'
}


def linear_sw_gw(
        ncells_per_edge=linear_sw_gw_defaults['ncells_per_edge'],
        dt=linear_sw_gw_defaults['dt'],
        tmax=linear_sw_gw_defaults['tmax'],
        element_family=linear_sw_gw_defaults['element_family'],
        element_order=linear_sw_gw_defaults['element_order'],
        dumpfreq=linear_sw_gw_defaults['dumpfreq'],
        stepper=linear_sw_gw_defaults['stepper'],
        reference_error=linear_sw_gw_defaults['reference_error'],
        reference_dirname=linear_sw_gw_defaults['reference_dirname'],
        dirname=linear_sw_gw_defaults['dirname']
):

    # ------------------------------------------------------------------------ #
    # Parameters for test case
    # ------------------------------------------------------------------------ #

    R0 = pi/9.                         # radius of perturbation (rad)
    lamda_c = -pi/2.                   # longitudinal centre of perturbation (rad)
    phi_c = pi/3.                      # latitudinal centre of perturbation (rad)
    Dpert_mag = 100.0                  # magnitude of perturbation (m)
    radius = 6371220.                  # planetary radius (m)
    mean_depth = 2000.                 # reference depth (m)
    u_max = 2*pi*radius/(12*24*60*60)  # Max amplitude of the zonal wind (m/s)

    # ------------------------------------------------------------------------ #
    # Set up model objects
    # ------------------------------------------------------------------------ #

    # Domain
    if element_family in ['RT', 'BDM']:
        mesh = GeneralIcosahedralSphereMesh(radius, ncells_per_edge, degree=4)
    elif element_family in ['RTCF', 'BDMCF']:
        mesh = GeneralCubedSphereMesh(radius, ncells_per_edge, degree=4)
    else:
        raise ValueError(f"Unknown element family: {element_family}")

    x, y, z = SpatialCoordinate(mesh)
    domain = Domain(mesh, dt, element_family, element_order)

    # Equation
    parameters = ShallowWaterParameters(mesh, H=mean_depth)
    eqns = LinearShallowWaterEquations(domain, parameters)

    # I/O
    # Adjust default directory name
    dirname = f'{dirname}_{stepper}_{element_family}{element_order}_n{ncells_per_edge}_dt{int(dt)}'

    chkptfreq = dumpfreq if not reference_error else None
    to_checkpoint = not reference_error
    output = OutputParameters(
        dirname=dirname, dumpfreq=dumpfreq, dump_nc=True, dump_vtus=False,
        checkpoint=to_checkpoint, chkptfreq=chkptfreq, dumplist=['D']
    )
    diagnostic_fields = [
        ZonalComponent('u'), MeridionalComponent('u'),
        Difference('D', 'Dbalanced')
    ]
    io = IO(domain, output, diagnostic_fields=diagnostic_fields)

    if stepper == 'siqn':
        # Transport schemes
        transport_schemes = [ForwardEuler(domain, "D")]
        transport_methods = [DefaultTransport(eqns, "D")]

        # Time stepper
        timestepper = SemiImplicitQuasiNewton(
            eqns, io, transport_schemes, transport_methods
        )
    elif stepper == 'trapezium_rule':
        timestepper = Timestepper(eqns, TrapeziumRule(domain), io)
    elif stepper == 'ssprk3':
        timestepper = Timestepper(eqns, SSPRK3(domain), io)
    elif stepper == 'rk4':
        timestepper = Timestepper(eqns, RK4(domain), io)
    elif stepper == 'bdf2':
        timestepper = Timestepper(eqns, BDF2(domain), io)
    elif stepper == 'leapfrog':
        timestepper = Timestepper(eqns, Leapfrog(domain), io)
    elif stepper == 'imex_euler':
        eqns.residual = eqns.residual.label_map(
            lambda t: t.has_label(pressure_gradient, coriolis),
            map_if_true=lambda t: implicit(t),
            map_if_false=keep
        )
        eqns.residual = eqns.residual.label_map(
            lambda t: t.has_label(transport),
            map_if_true=lambda t: explicit(t),
            map_if_false=keep
        )
        timestepper = Timestepper(eqns, IMEX_Euler(domain), io)
    else:
        raise ValueError(f"Unknown stepper: {stepper}")

    # ------------------------------------------------------------------------ #
    # Initial conditions
    # ------------------------------------------------------------------------ #

    g = parameters.g
    Omega = parameters.Omega

    u0 = timestepper.fields("u")
    D0 = timestepper.fields("D")

    uexpr = as_vector([-u_max*y/radius, u_max*x/radius, 0.0])
    Dexpr = - ((radius*Omega*u_max) * (z/radius)**2) / g

    Dbalanced = timestepper.fields("Dbalanced", space=D0.function_space())
    Dbalanced.interpolate(Dexpr)

    # Add perturbation to depth to trigger gravity wave
    lamda, phi, _ = lonlatr_from_xyz(x, y, z)
    great_arc = great_arc_angle(lamda, phi, lamda_c, phi_c)
    Dpert = Dpert_mag*exp(-(great_arc/R0)**2)
    Dexpr += Dpert

    u0.project(uexpr)
    D0.interpolate(Dexpr)

    Dbar = Function(D0.function_space()).assign(mean_depth)
    timestepper.set_reference_profiles([('D', Dbar)])

    # ------------------------------------------------------------------------ #
    # Run
    # ------------------------------------------------------------------------ #

    timestepper.run(t=0, tmax=tmax)

    # ------------------------------------------------------------------------ #
    # Compare against reference solution, if requested
    # ------------------------------------------------------------------------ #

    if reference_error:
        # Open reference solution
        ref_file = path.join(reference_dirname, 'chkpt.h5')
        with CheckpointFile(ref_file, 'r', domain.mesh.comm) as ref_chk:
            high_res_mesh = ref_chk.load_mesh('firedrake_default')
            D_high_res = ref_chk.load_function(high_res_mesh, 'D')
            # Interpolate final state onto high res mesh
            D_low_res = Function(D_high_res.function_space()).interpolate(timestepper.fields("D"))

            # Compute error against reference solution
            l2_error = errornorm(D_high_res, D_low_res)

            # Write error to new netcdf file
            error_file = path.join('results', dirname, 'error_diagnostic.nc')
            with Dataset(error_file, 'w') as ncfile:
                ncfile.createDimension('time', None)
                error_var = ncfile.createVariable('D_error', 'f8', ('time',))
                error_var.units = 'm'
                error_var.long_name = 'L2 error in D against reference solution'
                error_var[0] = l2_error

# ---------------------------------------------------------------------------- #
# MAIN
# ---------------------------------------------------------------------------- #


if __name__ == "__main__":

    parser = ArgumentParser(
        description=__doc__,
        formatter_class=ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        '--ncells_per_edge',
        help="The number of cells per edge of icosahedron",
        type=int,
        default=linear_sw_gw_defaults['ncells_per_edge']
    )
    parser.add_argument(
        '--dt',
        help="The time step in seconds.",
        type=float,
        default=linear_sw_gw_defaults['dt']
    )
    parser.add_argument(
        "--tmax",
        help="The end time for the simulation in seconds.",
        type=float,
        default=linear_sw_gw_defaults['tmax']
    )
    parser.add_argument(
        '--element_family',
        help="The element family to use for the velocity space.",
        type=str,
        default=linear_sw_gw_defaults['element_family']
    )
    parser.add_argument(
        '--element_order',
        help="The order of the elements to use.",
        type=int,
        default=linear_sw_gw_defaults['element_order']
    )
    parser.add_argument(
        '--dumpfreq',
        help="The frequency at which to dump field output.",
        type=int,
        default=linear_sw_gw_defaults['dumpfreq']
    )
    parser.add_argument(
        '--reference_error',
        help="Whether to test against the reference solution.",
        action='store_true',
        default=linear_sw_gw_defaults['reference_error']
    )
    parser.add_argument(
        '--reference_dirname',
        help="Name of directory with reference solution.",
        type=str,
        default=linear_sw_gw_defaults['reference_dirname']
    )
    parser.add_argument(
        '--dirname',
        help="The name of the directory to write to.",
        type=str,
        default=linear_sw_gw_defaults['dirname']
    )
    parser.add_argument(
        '--stepper',
        help=(
            "Which timestepper to use."
        ),
        type=str,
        default=linear_sw_gw_defaults['stepper']
    )
    args, unknown = parser.parse_known_args()

    linear_sw_gw(**vars(args))
