"""Hydrostatic tank example with materials

"""

import os.path

import numpy as np 

# PyZoltan imports
from pyzoltan.core.carray import LongArray

# PySPH imports
from pysph.base.utils import get_particle_array_wcsph as gpa 
from pysph.base.kernels import Gaussian, WendlandQuintic, CubicSpline, QuinticSpline
from pysph.solver.solver import Solver
from pysph.solver.application import Application
from pysph.sph.integrator import PECIntegrator
from pysph.sph.integrator_step import WCSPHStep

# the equations
from pysph.sph.equation import Group

from pysph.sph.basic_equations import SummationDensity
# Diffusion equation(s) I am working on
from pysph.sph.misc.diffusion import DiffuseMaterial

# Equations for REF1
from pysph.sph.wc.transport_velocity import VolumeFromMassDensity, \
    ContinuityEquation, \
    MomentumEquationPressureGradient, \
    MomentumEquationArtificialViscosity, \
    SolidWallPressureBC

# Monaghan type repulsive boundary forces used in REF(2)
from pysph.sph.boundary_equations import MonaghanBoundaryForce, \
    MonaghanKajtarBoundaryForce

# Equations for the standard WCSPH formulation and dynamic boundary
# conditions defined in REF3
from pysph.sph.wc.basic import TaitEOS, TaitEOSHGCorrection, MomentumEquation
from pysph.sph.basic_equations import XSPHCorrection, \
    MonaghanArtificialViscosity

# domain and reference values
length_x = 2.0; length_y = 1.0; liquid_height = 0.5; air_height = 0.5
gravity_y = -1.0
Vmax = np.sqrt(abs(gravity_y) * liquid_height)
c0 = 10 * Vmax; 
rho0 = 1000.0
p0 = c0*c0*rho0
gamma = 1.0

# Reynolds number and kinematic viscosity
ReynoldsNumber = 100; nu = Vmax * length_y/ReynoldsNumber

# Numerical setup
number_x = 100; dx = length_x/number_x
ghost_extent = 5.5 * dx
hdx = 1.2

# adaptive time steps
h0 = hdx * dx
dt_cfl = 0.25 * h0/( c0 + Vmax )
dt_viscous = 0.125 * h0**2/nu
dt_force = 0.25 * np.sqrt(h0/abs(gravity_y))

tdamp = 1.0
tf = 4.0
dt = 0.75 * min(dt_cfl, dt_viscous, dt_force)
output_at_times = np.arange(0.25, 2.1, 0.25)

def damping_factor(t, tdamp):
    if t < tdamp:
        return 0.5 * ( np.sin((-0.5 + t/damp)*np.pi) + 1.0 )
    else:
        return 1.0

class HydrostaticTankMaterials(Application):
    def add_user_options(self, group):
        group.add_argument(
            '--bc-type', action='store', type=int,
            dest='bc_type', default=1,
            help="Specify the implementation type one of (1, 2, 3)"
        )
    
    def create_particles(self):
        # create all the particles
        _x = np.arange( -ghost_extent, length_x + ghost_extent, dx )
        _y = np.arange( -ghost_extent, length_y, dx )
        x, y = np.meshgrid(_x, _y); x = x.ravel(); y = y.ravel()

        # sort out the fluid and the solid
        fluid_indices = []
        air_indices = []
        for i in range(x.size):
            if ( (x[i] > 0.0 ) and (x[i] < length_x) ):
                if (y[i] > 0.0):
                    if (y[i] < liquid_height):
                        fluid_indices.append(i)
                    elif (y[i] < liquid_height + air_height):
                        air_indices.append(i)
        
        # create the arrays
        solid = gpa(name='solid', x=x, y=y)

        

        # remove the fluid particles from the solid
        fluid = solid.extract_particles(fluid_indices); fluid.set_name('fluid')
        indices_to_remove = fluid_indices + air_indices
        solid.remove_particles(indices_to_remove)

        # remove the lid to generate an open tank
        
            # actually, i am skipping this
        
        ''' print("Hydrostatic tank :: nfluid = %d, nsolid = %d, dt = %g"%(
            fluid.get_number_of_particles(),
            solid.get_number_of_particles(), dt)) '''
        
        ###### ADD PARTICLE PROPS FOR MULTI-PHASE SPH ######

        # particle volume
        fluid.add_property('V')
        solid.add_property('V')

        # material stuff I am working on
        fluid.add_property('material_amount')
        fluid.add_property('material_velocity')
        solid.add_property('material_amount')
        solid.add_property('material_velocity')

        # kernel sum term for boundary particles
        solid.add_property('wij')

        # advection velocities and accelerations
        for name in ('auhat', 'avhat', 'awhat'):
            fluid.add_property(name)
        
        ##### INITIALIZE PARTICLE PROPS #####
        fluid.rho[:] = rho0
        solid.rho[:] = rho0

        fluid.rho0[:] = rho0
        solid.rho0[:] = rho0

        # mass is set to get the reference density of rho0
        volume = dx * dx

        # volume is set as dx^2
        fluid.V[:] = 1./volume
        solid.V[:] = 1./volume

        fluid.m[:] = volume * rho0
        solid.m[:] = volume * rho0

        # smoothing lengths
        fluid.h[:] = hdx * dx
        solid.h[:] = hdx * dx

        # Giving fluid particles properties 'material_amount' and 'material_velocity'
        fluid.material_amount[:] = 5.0
        fluid.material_amount[50] = 50
        fluid.material_velocity[:] = 0.0
        solid.material_amount[:] = 50.0
        solid.material_velocity[:] = 0.0
        
        properties_to_save = ['pid', 'tag', 'material_velocity', 'gid', 'material_amount', 'rho', 'V', 'h', 'm', 'p', 'u', 'w', 'v', 'y', 'x', 'div', 'z']

        solid.set_output_arrays(properties_to_save)
        fluid.set_output_arrays(properties_to_save)
        

        #print('fluid material amount: ', fluid.material_amount, 'fluid props', fluid.properties.keys())

        # return the particle list
        return [fluid, solid]

    def create_solver(self):
        # Create the kernel
        #kernel = Gaussian(dim=2)
        kernel = QuinticSpline(dim=2)

        integrator = PECIntegrator(fluid=WCSPHStep())

        # Create a solver
        solver = Solver(kernel = kernel, dim=2, integrator=integrator,
                        tf=tf, dt=dt, output_at_times=output_at_times)
        return solver

    def create_equations(self):
        # Formulation for REF1
        # (using only first set of equations for simplicity)
        equations = [

            Group(equations=[
                DiffuseMaterial(dest='fluid', sources=['fluid'], diffusion_speed=0.1)
                ], ),

            # For the multi-phase formulation, we require an estimate of the
            # particle volume. This can be either defined from the particle
            # number density or simply as the ratio of mass to density
            Group(equations=[
                VolumeFromMassDensity(dest='fluid', sources=None)
                ], ),
            
            # Equation of state is typically the Tait EOS with a suitable
            # exponent gamma
            Group(equations=[
                TaitEOS(dest='fluid', sources=None, rho0=rho0, c0=c0, gamma=gamma),
                ], ),
            
            # The boundary conditions are imposed by extrapolating the fluid
            # pressure, taking into consideration the boundary acceleration
            Group(equations=[
                SolidWallPressureBC(dest='solid', sources=['fluid'], b=1.0, gy=gravity_y,
                                    rho0=rho0, p0=p0)
            ], ),

            # Main acceleration block
            Group(equations=[

                # Continuity equation
                ContinuityEquation(dest='fluid', sources=['fluid', 'solid']),

                # Pressure gradient with acceleration damping
                MomentumEquationPressureGradient(
                    dest='fluid', sources=['fluid', 'solid'], pb=0.0, gy=gravity_y,
                    tdamp=tdamp),

                # artificial viscosity for stability
                MomentumEquationArtificialViscosity(
                    dest='fluid', sources=['fluid', 'solid'], alpha=0.24, c0=c0),

                # Position step with XSPH
                XSPHCorrection(dest='fluid', sources=['fluid'], eps=0.0)
                ]),
        ]

        return equations

if __name__ == '__main__':
    app = HydrostaticTankMaterials()
    app.run()