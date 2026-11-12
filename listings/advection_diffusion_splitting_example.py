from firedrake import (PeriodicIntervalMesh, SpatialCoordinate, Constant,
                       FunctionSpace, Function, TestFunction, TrialFunction,
                       LinearVariationalProblem, LinearVariationalSolver,
                       exp, dx)
from firedrake.fml import all_terms, drop, keep, replace_subject, subject
from firedrake.pyplot import plot
from gusto.core import time_derivative, transport, diffusion
import matplotlib.pyplot as plt

# define mesh
L = 1
nx = 100
mesh = PeriodicIntervalMesh(nx, L)
x = SpatialCoordinate(mesh)[0]

# define function space, initial condition and parameters
V = FunctionSpace(mesh, "Lagrange", 1)
d = 0.05
qn = Function(V).interpolate(exp(-((x-0.5*L)/d)**2))
c = Constant(0.5)
kappa = Constant(0.01)

# define equation residual
q = Function(V)
v = TestFunction(V)
residual = (
    time_derivative(q * v * dx) +
    transport(c * q.dx(0) * v * dx) +
    diffusion(kappa * q.dx(0) * v.dx(0) * dx)
    )
residual = subject(residual, q)    # label q as the subject

# set up solvers
q_trial = TrialFunction(V)
qstar = Function(V)
qnp1 = Function(V)
dt = 0.01

# extract the time derivative, or mass, term from the residual
mass = residual.label_map(lambda t: t.has_label(time_derivative),
                          map_if_true=keep,
                          map_if_false=drop)
# the lhs of the explicit problem is the mass term with unknown q_trial
lhs_explicit = mass.label_map(all_terms,
                              map_if_true=replace_subject(q_trial))
# extract the transport term from the residual and multiply by -dt
transport_rhs = residual.label_map(lambda t: t.has_label(transport),
                                   map_if_true=lambda t: -dt*t,
                                   map_if_false=drop)
# the rhs of the explicit problem is mass + transport_rhs evaluated at qn
rhs_explicit = mass + transport_rhs
rhs_explicit = rhs_explicit.label_map(all_terms,
                                      map_if_true=replace_subject(qn))
# setup the problem and solver, storing the result in qstar
prob_explicit = LinearVariationalProblem(lhs_explicit.form,
                                         rhs_explicit.form,
                                         qstar)
solver_explicit = LinearVariationalSolver(prob_explicit)

# extract the diffusion term from the residual and multiply by dt
diffusion_lhs = residual.label_map(lambda t: t.has_label(diffusion),
                                   map_if_true=lambda t: dt*t,
                                   map_if_false=drop)
# the lhs of the implicit problem is mass + diffusion_lhs with unknown q_trial
lhs_implicit = mass + diffusion_lhs
lhs_implicit = lhs_implicit.label_map(all_terms,
                                      map_if_true=replace_subject(q_trial))
# the rhs of the implicit problem is the mass term evaluated at qstar
rhs_implicit = mass.label_map(all_terms,
                              map_if_true=replace_subject(qstar))
# setup the problem and solver, storing the result in qnp1
prob_implicit = LinearVariationalProblem(lhs_implicit.form,
                                         rhs_implicit.form,
                                         qnp1)
solver_implicit = LinearVariationalSolver(prob_implicit)

# Timestepping loop from t=0 to t=0.5
tmax = 0.5
t = 0
while t < tmax:
    t += dt
    solver_explicit.solve()
    solver_implicit.solve()
    qn.assign(qnp1)
    plot(qn)
    plt.show()
