"""Rayleigh-Taylor instability problem. (16 hours)
"""

# numpy
import numpy as np

# PySPH imports
from pysph.base.utils import get_particle_array
from pysph.solver.solver import Solver
from pysph.solver.application import Application

# Diffusion equation(s) I am working on
from pysph.sph.misc.diffusion import DiffuseMaterial
# Material class I am working on
from pysph.sph.misc.mixture_particle import MixtureParticle
# Modified TaitEOSHGCorrection to account for different rho0 values of Mixture Particles
from pysph.sph.misc.mixture_particle_equations import TaitEOSHGCorrectionVariableRho, StateEquationVariableRho0, TVFSchemeMixtureParticle
from pysph.sph.wc.transport_velocity import StateEquation


# domain and reference values
gy = -1.0
Lx = 1.0; Ly = 2.0
Re = 420; Vmax = np.sqrt(0.5*Ly*abs(gy))
nu = Vmax*Ly/Re

# density for the two phases
rho1 = 1.8; rho2 = 1.0

# speed of sound and reference pressure
Fr = 0.01
c0 = Vmax/Fr
p1 = c0**2 * rho1
p2 = c0**2 * rho2

# Numerical setup
nx = 50; dx = Lx/nx
ghost_extent = 5 * dx
hdx = 1.2

# adaptive time steps
h0 = hdx * dx
dt_cfl = 0.25 * h0/( c0 + Vmax )
dt_viscous = 0.125 * h0**2/nu
dt_force = 0.25 * np.sqrt(h0/abs(gy))

tf = 25
dt = 0.5 * min(dt_cfl, dt_viscous, dt_force)

class RayleighTaylorTemperature(Application):
    def create_particles(self):
        # create all the particles
        _x = np.arange( -ghost_extent - dx/2, Lx + ghost_extent + dx/2, dx )
        _y = np.arange( -ghost_extent - dx/2, Ly + ghost_extent + dx/2, dx )
        x, y = np.meshgrid(_x, _y); x = x.ravel(); y = y.ravel()

        # sort out the fluid and the solid
        indices = []
        for i in range(x.size):
            if ( (x[i] > 0.0) and (x[i] < Lx) ):
                if ( (y[i] > 0.0)  and (y[i] < Ly) ):
                    indices.append(i)

        # create the arrays
        solid = get_particle_array(name='solid', x=x, y=y)

        # remove the fluid particles from the solid
        fluid = solid.extract_particles(indices); fluid.set_name('fluid')
        solid.remove_particles(indices)

        # sort out the two fluid phases
        indices = []
        for i in range(fluid.get_number_of_particles()):
            if fluid.y[i] > 1 - 0.15*np.sin(2*np.pi*fluid.x[i]):
                indices.append(i)

        fluid1 = fluid.extract_particles(indices); fluid1.set_name('fluid1')
        fluid2 = fluid
        fluid2.set_name('fluid2')
        fluid2.remove_particles(indices)

        materials = MixtureParticle(baseProperties={'rho0': rho1, 'p0': p1, 't0': 270, 't': 270})

        print('props', materials.generateFullParticleProperties())

        # add requisite properties to the arrays:
        self.scheme.setup_properties([fluid1, fluid2, solid])

        for propertyName in materials.generateFullParticleProperties():
            fluid1.add_property(propertyName, type='double', default=materials.getPropertyValue(propertyName))
            fluid2.add_property(propertyName, type='double', default=materials.getPropertyValue(propertyName))

        fluid1.rho[:] = rho1
        fluid1.rho0[:] = rho1
        fluid1.p[:] = p1
        fluid1.p0[:] = p1
        fluid2.rho[:] = rho2
        fluid2.rho0[:] = rho2
        fluid2.p[:] = p2
        fluid2.p0[:] = p2

        print("Rayleigh Taylor Instability problem :: Re = %d, nfluid = %d, nsolid=%d, dt = %g"%(
            Re, fluid1.get_number_of_particles() + fluid2.get_number_of_particles(),
            solid.get_number_of_particles(), dt))

        

        # setup the particle properties
        volume = dx * dx

        # mass is set to get the reference density of each phase
        fluid1.m[:] = volume * rho1
        fluid2.m[:] = volume * rho2

        # volume is set as dx^2
        fluid1.V[:] = 1./volume
        fluid2.V[:] = 1./volume
        solid.V[:] = 1./volume

        # smoothing lengths
        fluid1.h[:] = hdx * dx
        fluid2.h[:] = hdx * dx
        solid.h[:] = hdx * dx

        # return the arrays
        return [fluid1, fluid2, solid]

    def create_scheme(self):
        s = TVFSchemeMixtureParticle(
            ['fluid1', 'fluid2'], ['solid'], dim=2, rho0=rho1, c0=c0, nu=nu,
            p0=p1, pb=p1, h0=dx*hdx, gy=gy
        )
        s.configure_solver(tf=tf, dt=dt, pfreq=500)
        return s

    def create_equations(self):
        equations = super(RayleighTaylorTemperature, self).create_equations()
        from pysph.sph.equation import Group

        return equations

if __name__ == '__main__':
    app = RayleighTaylorTemperature()
    app.run()
