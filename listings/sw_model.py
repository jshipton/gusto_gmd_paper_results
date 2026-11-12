"""
Test Case 2 (solid-body rotation with geostrophically-balanced flow) of
Williamson et al, 1992:
``A standard test set for numerical approximations to the shallow water
equations in spherical geometry'', JCP.

The example here uses the icosahedral sphere mesh and degree 1 spaces.
"""

from firedrake import *
from gusto import *

# Parameters for this case
radius = 6371220                   # planetary radius (m)
mean_depth = 5960                  # reference depth (m)
day = 24*60*60                     # length of day (s)
u_max = 2*pi*radius/(12*day)       # Max amplitude of the zonal wind (m/s)
dt = 900.                          # timestep (s)

# Set up model domain
# Use an icosahedral mesh with 16 cells along each panel edge and quadratic (degree 2) coordinates
mesh = GeneralIcosahedralSphereMesh(radius, num_cells_per_edge_of_panel=16, degree=2)

# Set up I/O to write out every 96 timesteps (1 day) to results/sw_geostrophic_balance directory
output = OutputParameters(dirname='sw_geostrophic_balance', dumpfreq=96)

# Set up parameters and equation
parameters = ShallowWaterParameters(mesh, H=mean_depth)
eqns = ShallowWaterEquations

# Set up model
model = SIQNModel(mesh, dt, parameters, eqns, family='BDM')
model.setup(output)

# Velocity and depth fields for initialisation
u0 = model.stepper.fields('u')
D0 = model.stepper.fields('D')

# Set parameters and coordinates
g = parameters.g                  # use default value of gravity for Earth
Omega = parameters.Omega          # use default value of planetary rotation rate for Earth
xyz = SpatialCoordinate(mesh)     # Cartesian coordinates of mesh
_, lat, _ = lonlatr_from_xyz(xyz[0], xyz[1], xyz[2])   # spherical coordinates, only need latitude
e_lon = xyz_vector_from_lonlatr(1.0, 0.0, 0.0, xyz)    # unit longitudinal basis vector

# Initialise zonal flow velocity (equations 90-91 of Williamson et al, (1992))
u0.project(u_max*cos(lat)*e_lon)
# Initialise balanced fluid depth (equation 95 of Williamson et al, (1992))
D0.interpolate(mean_depth - (radius * Omega * u_max + 0.5*u_max**2)*(sin(lat))**2/g)

# Set reference profiles
Dbar = Function(D0.function_space()).assign(D0)
model.stepper.set_reference_profiles([('D', Dbar)])
    
# Run for 5 days
model.stepper.run(t=0, tmax=5*day)
