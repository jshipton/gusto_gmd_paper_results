"""
A linearised form of Test Case 2 (solid-body rotation) of Williamson et al 1992:
``A standard test set for numerical approximations to the shallow water
equations in spherical geometry'', JCP.
"""

from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from firedrake import (
  Function, SpatialCoordinate, pi, errornorm, sin, grad
)
from firedrake.fml import keep
from gusto import (
    Domain, IO, OutputParameters, SemiImplicitQuasiNewton, DefaultTransport,
    ForwardEuler, ZonalComponent, MeridionalComponent, ShallowWaterParameters,
    LinearShallowWaterEquations, GeneralIcosahedralSphereMesh, Difference,
    GeneralCubedSphereMesh, SSPRK3, RK4, BDF2, Leapfrog, IMEX_Euler,
    Timestepper, TrapeziumRule, lonlatr_from_xyz
)
from gusto.core.labels import (
    explicit, implicit, pressure_gradient, transport, coriolis
)
from netCDF4 import Dataset
from os import path

linear_sw_w2_defaults = {
    'ncells_per_edge': 16,        # number of cells per icosahedron edge
    'dt': 900.0,                  # 15 minutes
    'tmax': 3.*60.*60.,           # 3 hours
    'element_family': 'BDM',      # family for velocity space
    'element_order': 1,           # order of the elements
    'dumpfreq': 12,               # once with default options
    'stepper': 'trapezium_rule',  # default is TrapeziumRule
    'dirname': 'linear_sw_w2'
}


def linear_sw_w2(
        ncells_per_edge=linear_sw_w2_defaults['ncells_per_edge'],
        dt=linear_sw_w2_defaults['dt'],
        tmax=linear_sw_w2_defaults['tmax'],
        element_family=linear_sw_w2_defaults['element_family'],
        element_order=linear_sw_w2_defaults['element_order'],
        dumpfreq=linear_sw_w2_defaults['dumpfreq'],
        stepper=linear_sw_w2_defaults['stepper'],
        dirname=linear_sw_w2_defaults['dirname']
):

    # ------------------------------------------------------------------------ #
    # Parameters for test case
    # ------------------------------------------------------------------------ #

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

    chkptfreq = None
    to_checkpoint = False
    output = OutputParameters(
        dirname=dirname, dumpfreq=dumpfreq, dump_nc=True, dump_vtus=True,
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

    _, lat, _ = lonlatr_from_xyz(x, y, z)
    g = parameters.g
    Omega = parameters.Omega
    H1 = domain.spaces('H1')

    u0 = timestepper.fields("u")
    D0 = timestepper.fields("D")
    psi = Function(H1)

    psi.interpolate(-radius*u_max*sin(lat))
    uexpr = domain.perp(grad(psi))
    Dexpr = - ((radius*Omega*u_max) * (z/radius)**2) / g

    Dbalanced = timestepper.fields("Dbalanced", space=D0.function_space())
    Dbalanced.project(Dexpr)

    u0.project(uexpr)
    D0.project(Dexpr)

    Dbar = Function(D0.function_space()).assign(mean_depth)
    timestepper.set_reference_profiles([('D', Dbar)])

    # ------------------------------------------------------------------------ #
    # Run
    # ------------------------------------------------------------------------ #

    timestepper.run(t=0, tmax=tmax)

    # ------------------------------------------------------------------------ #
    # Compare against reference solution, if requested
    # ------------------------------------------------------------------------ #

    # Compute error against reference solution
    l2_error = errornorm(Dbalanced, D0)

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
        default=linear_sw_w2_defaults['ncells_per_edge']
    )
    parser.add_argument(
        '--dt',
        help="The time step in seconds.",
        type=float,
        default=linear_sw_w2_defaults['dt']
    )
    parser.add_argument(
        "--tmax",
        help="The end time for the simulation in seconds.",
        type=float,
        default=linear_sw_w2_defaults['tmax']
    )
    parser.add_argument(
        '--element_family',
        help="The element family to use for the velocity space.",
        type=str,
        default=linear_sw_w2_defaults['element_family']
    )
    parser.add_argument(
        '--element_order',
        help="The order of the elements to use.",
        type=int,
        default=linear_sw_w2_defaults['element_order']
    )
    parser.add_argument(
        '--dumpfreq',
        help="The frequency at which to dump field output.",
        type=int,
        default=linear_sw_w2_defaults['dumpfreq']
    )
    parser.add_argument(
        '--dirname',
        help="The name of the directory to write to.",
        type=str,
        default=linear_sw_w2_defaults['dirname']
    )
    parser.add_argument(
        '--stepper',
        help=(
            "Which timestepper to use."
        ),
        type=str,
        default=linear_sw_w2_defaults['stepper']
    )
    args, unknown = parser.parse_known_args()

    linear_sw_w2(**vars(args))
