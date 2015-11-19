"""Incompressible flow past a periodic array of cylinders. (42 hours)


See Ellero and Adams, International Journal for Numerical Methods in
Engineering, 2011, vol 86, pp 1027-1040 for the detailed parameters for this
problem and also  Adami, Hu and Adams, JCP, 2013, vol 241, pp 292-307.

In particular, we note that we set c0 from Ellero and Adams as using the
value from Adami et al. will cause the solution to blow up.

If one sets c0=10*Umax and sets pb=300*p0, that will cause particles to void
at the rear of the cylinder.

"""

# PySPH imports
from pysph.base.nnps import DomainManager
from pysph.base.utils import get_particle_array
from pysph.base.kernels import QuinticSpline
from pysph.solver.solver import Solver
from pysph.solver.application import Application
from pysph.sph.integrator import PECIntegrator
from pysph.sph.integrator_step import TransportVelocityStep

# the eqations
from pysph.sph.equation import Group
from pysph.sph.wc.transport_velocity import (SummationDensity,
    SetWallVelocity, StateEquation,
    MomentumEquationPressureGradient, MomentumEquationViscosity,
    MomentumEquationArtificialStress,
    SolidWallPressureBC, SolidWallNoSlipBC, VolumeSummation)

# numpy
import os
import numpy as np

# domain and reference values
L = 0.12; Umax = 1.2e-4
a = 0.02; H = 4*a
fx = 2.5e-4

# c0 is set from Ellero and Adams.
# Note that setting this to 0.1*np.sqrt(a*fx) as per Adami Hu and Adams is
# incorrect and will actually cause a blow up of the solution.
c0 = 0.02
rho0 = 1000.0
p0 = c0*c0*rho0
pb = p0

# Reynolds number and kinematic viscosity
nu = 0.1/rho0; Re = a*Umax/nu

# Numerical setup
nx = 144; dx = L/nx
ghost_extent = 5 * 1.5 * dx
hdx = 1.2

# adaptive time steps
h0 = hdx * dx
dt_cfl = 0.25 * h0/( c0 + Umax )
dt_viscous = 0.125 * h0**2/nu
dt_force = 0.25 * np.sqrt(h0/abs(fx))

T = a/Umax

tf = 2.5*T
dt = min(dt_cfl, dt_viscous, dt_force)


class PeriodicCylinders(Application):

    def create_domain(self):
        # domain for periodicity
        domain = DomainManager(xmin=0, xmax=L, periodic_in_x=True)
        return domain

    def create_particles(self):
        # create all the particles
        _x = np.arange( dx/2, L, dx )
        _y = np.arange( -ghost_extent, H+ghost_extent, dx )
        x, y = np.meshgrid(_x, _y); x = x.ravel(); y = y.ravel()

        # sort out the fluid and the solid
        indices = []
        cx = 0.5 * L; cy = 0.5 * H
        for i in range(x.size):
            xi = x[i]; yi = y[i]
            if ( np.sqrt( (xi-cx)**2 + (yi-cy)**2 ) > a ):
                if ( (yi > 0) and (yi < H) ):
                    indices.append(i)

        # create the arrays
        solid = get_particle_array(name='solid', x=x, y=y)

        # remove the fluid particles from the solid
        fluid = solid.extract_particles(indices); fluid.set_name('fluid')
        solid.remove_particles(indices)

        print("Periodic cylinders :: Re = %g, nfluid = %d, nsolid=%d, dt = %g"%(
            Re, fluid.get_number_of_particles(),
            solid.get_number_of_particles(), dt))
        print("tf = %f"%tf)

        # add requisite properties to the arrays:

        # volume from number density for fluid and solid
        fluid.add_property('V')
        solid.add_property('V' )

        # extrapolated velocities for the solid
        for name in ['uf', 'vf', 'wf']:
            solid.add_property(name)

        # advection velocities and accelerations for the fluid
        for name in ('uhat', 'vhat', 'what', 'auhat', 'avhat', 'awhat'):
            fluid.add_property(name)

        # kernel summation correction for the solid
        solid.add_property('wij')

        # dummy velocities for the solid walls
        # required for the no-slip BC
        for name in ['ug', 'vg', 'wg']:
            solid.add_property(name)

        # imposed accelerations on the solid
        solid.add_property('ax')
        solid.add_property('ay')
        solid.add_property('az')
        solid.add_output_arrays(['p'])
        # magnitude of velocity
        fluid.add_property('vmag2')
        fluid.add_output_arrays(['vmag2', 'p'])

        # setup the particle properties
        volume = dx * dx

        # mass is set to get the reference density of rho0
        fluid.m[:] = volume * rho0
        solid.m[:] = volume * rho0

        # initial particle density
        fluid.rho[:] = rho0
        solid.rho[:] = rho0

        # volume is set as dx^2. V is the number density form of the
        # particle volume and will be computed in the equations for the
        # fluid phase. The initial values are used for the solid phase
        fluid.V[:] = 1./volume
        solid.V[:] = 1./volume

        # particle smoothing lengths
        fluid.h[:] = hdx * dx
        solid.h[:] = hdx * dx

        # return the particle list
        return [fluid, solid]

    def create_solver(self):
        kernel = QuinticSpline(dim=2)

        # The predictor corrector integrator for the TV formulation. As per
        # the paper, the integrator is defined for PEC mode with only one
        # function evaluation per time step.
        integrator = PECIntegrator(fluid=TransportVelocityStep())

        # Create a solver. Damping is performed for 100 iterations.
        solver = Solver(
            kernel=kernel, dim=2, integrator=integrator,
            adaptive_timestep=False, tf=tf, dt=dt, n_damp=100, pfreq=500)
        return solver

    def create_equations(self):
        equations = [

            # Summation density along with volume summation for the fluid
            # phase. This is done for all local and remote particles. At the
            # end of this group, the fluid phase has the correct density
            # taking into consideration the fluid and solid
            # particles.
            Group(
                equations=[
                    VolumeSummation(
                        dest='solid', sources=['fluid', 'solid']
                    ),
                    SummationDensity(dest='fluid', sources=['fluid','solid']),
                    ], real=False),

            # Once the fluid density is computed, we can use the EOS to set
            # the fluid pressure. Additionally, the dummy velocity for the
            # channel is set, which is later used in the no-slip wall BC.
            Group(
                equations=[
                    StateEquation(dest='fluid', sources=None, p0=p0, rho0=rho0, b=1.0),
                    SetWallVelocity(dest='solid', sources=['fluid']),
                    ], real=False),

            # Once the pressure for the fluid phase has been updated, we can
            # extrapolate the pressure to the ghost particles. After this
            # group, the fluid density, pressure and the boundary pressure has
            # been updated and can be used in the integration equations.
            Group(
                equations=[
                    SolidWallPressureBC(dest='solid', sources=['fluid'],
                                        gx=fx, b=1.0, rho0=rho0, p0=p0),
                    ], real=False),

            # The main accelerations block. The acceleration arrays for the
            # fluid phase are updated in this stage for all local particles.
            Group(
                equations=[
                    # Pressure gradient terms
                    MomentumEquationPressureGradient(
                        dest='fluid', sources=['fluid', 'solid'], gx=fx, pb=pb),

                    # fluid viscosity
                    MomentumEquationViscosity(
                        dest='fluid', sources=['fluid'], nu=nu),

                    # No-slip boundary condition. This is effectively a
                    # viscous interaction of the fluid with the ghost
                    # particles.
                    SolidWallNoSlipBC(
                        dest='fluid', sources=['solid'], nu=nu),

                    # Artificial stress for the fluid phase
                    MomentumEquationArtificialStress(dest='fluid', sources=['fluid']),

                    ], real=True),
        ]
        return equations

    def post_process(self, info_fname):
        info = self.read_info(info_fname)
        if len(self.output_files) == 0:
            return

        t, cd = self._plot_cd_vs_t()
        res = os.path.join(self.output_dir, 'results.npz')
        np.savez(res, t=t, cd=cd)

    def _plot_cd_vs_t(self):
        from pysph.solver.utils import iter_output, load
        from pysph.tools.sph_evaluator import SPHEvaluator
        from pysph.sph.rigid_body import (NumberDensity, ViscosityRigidBody,
                                          PressureRigidBody)
        from pysph.sph.equation import Group

        data = load(self.output_files[0])
        solid = data['arrays']['solid']
        fluid = data['arrays']['fluid']
        x, y = solid.x.copy(), solid.y.copy()
        cx = 0.5 * L; cy = 0.5 * H
        inside = np.sqrt((x-cx)**2 + (y-cy)**2) <= a
        dest = solid.extract_particles(inside.nonzero()[0])
        # We use the same equations for this as the simulation, except that we
        # do not include the acceleration terms as these are externally
        # imposed.  The goal of these is to find the force of the fluid on the
        # cylinder, thus, gx=0.0 is used in the following.
        equations = [
            Group(
                equations=[
                    VolumeSummation(
                        dest='fluid', sources=['fluid', 'solid']
                    ),
                    VolumeSummation(
                        dest='solid', sources=['fluid', 'solid']
                    ),
                    ], real=False),

            Group(
                equations=[
                    SetWallVelocity(dest='solid', sources=['fluid']),
                    ], real=False),

            Group(
                equations=[
                    SolidWallPressureBC(dest='solid', sources=['fluid'],
                                        gx=0.0, b=1.0, rho0=rho0, p0=p0),
                    ], real=False),

            Group(
                equations=[
                    # Pressure gradient terms
                    MomentumEquationPressureGradient(
                        dest='fluid', sources=['solid'], gx=0.0, pb=pb),
                    SolidWallNoSlipBC(
                        dest='fluid', sources=['solid'], nu=nu),
                    ], real=True),
        ]

        sph_eval = SPHEvaluator(
            arrays=[dest, fluid], equations=equations, dim=2,
            kernel=QuinticSpline(dim=2)
        )

        t, cd = [], []
        for sd, fluid in iter_output(self.output_files, 'fluid'):
            t.append(sd['t'])
            sph_eval.update_particle_arrays([dest, fluid])
            sph_eval.evaluate()
            Fx = np.sum(-fluid.au*fluid.m)
            cd.append(Fx/(nu*rho0*Umax))

        t, cd = list(map(np.asarray, (t, cd)))

        # Now plot the results.
        from matplotlib import pyplot as plt
        f = plt.figure()
        plt.plot(t, cd)
        plt.xlabel('$t$'); plt.ylabel(r'$C_D$')
        fig = os.path.join(self.output_dir, "cd_vs_t.png")
        plt.savefig(fig, dpi=300)
        plt.close()

        return t, cd

if __name__ == '__main__':
    app = PeriodicCylinders()
    app.run()
    app.post_process(app.info_filename)
