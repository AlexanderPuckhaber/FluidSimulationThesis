"""Hydrostatic tank example with materials that mix together,
   but they are not balanced as one is denser than the other

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
# Material class I am working on
from pysph.sph.misc.mixture_particle import MixtureParticle

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
liquid1_length = length_x / 2.0; liquid2_length = length_x / 2.0

fluid_name = 'water'
solid_name = ['solid']

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

        # create mixture particle
        materials = MixtureParticle()
        materials.addMaterial('water', {'_frac':0.8, 'm':0.10})
        materials.addMaterial('oil', {'_frac':0.2, 'm':0.09})

        print(materials.getMaterialList())

        print('FULL PARTICLE PROPERTIES:', materials.generateFullParticleProperties())

        materials.addMaterialProperties('oil', {'rho':.008})
        materials.addMaterialDriftVelocities(dim=2)

        print(materials.getMaterialList())

        print(materials.generateFullParticleProperties())

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
        print('len', fluid.get_number_of_particles())

        # make half of the fluid particles water and the other half oil
        water_indicies = []
        oil_indicies = []
        for i in range(fluid.get_number_of_particles()):
            if (i < fluid.get_number_of_particles() / 2.0):
                oil_indicies.append(i)
            else:
                water_indicies.append(i)

        water = fluid.extract_particles(water_indicies); water.set_name('fluid')
        oil = fluid.extract_particles(oil_indicies); oil.set_name('fluid')

        fluid.remove_particles(fluid_indices)

        indices_to_remove = fluid_indices + air_indices
        solid.remove_particles(indices_to_remove)

        # remove the lid to generate an open tank
        
            # actually, i am skipping this
        
        ''' print("Hydrostatic tank :: nfluid = %d, nsolid = %d, dt = %g"%(
            fluid.get_number_of_particles(),
            solid.get_number_of_particles(), dt)) '''
        
        ###### ADD PARTICLE PROPS FOR MULTI-PHASE SPH ######

        # particle volume
        water.add_property('V')
        oil.add_property('V')
        solid.add_property('V')

        # kernel sum term for boundary particles
        solid.add_property('wij')

        # advection velocities and accelerations
        for name in ('auhat', 'avhat', 'awhat'):
            water.add_property(name)
            oil.add_property(name)

        for propertyName in materials.generateFullParticleProperties():
            water.add_property(propertyName, type='double', default=materials.getPropertyValue(propertyName))
            oil.add_property(propertyName, type='double', default=materials.getPropertyValue(propertyName))

        print('water default', water.default_values)
        print('oil default', oil.default_values)

        # water.copy_over_properties(materials.generateFullParticleProperties())
        # oil.copy_over_properties(materials.generateFullParticleProperties())
        
        ##### INITIALIZE PARTICLE PROPS #####
        water.rho[:] = rho0 * 1.00
        water.rho0[:] = rho0 * 1.00
        waterRho = rho0 * 1.00
        oil.rho[:] = rho0
        oil.rho0[:] = rho0
        oilRho = rho0
        solid.rho[:] = rho0
        solid.rho0[:] = rho0
        solidRho = rho0


        # mass is set to get the reference density of rho0
        volume = dx * dx

        VTarget = 1./volume

        waterMass = volume * waterRho
        oilMass = volume * oilRho * 0.5
        solidMass = volume * solidRho

        water.m[:] = waterMass
        oil.m[:] = oilMass
        solid.m[:] = solidMass

        print('volume', volume, 'VTarget', VTarget, 'waterRho', waterRho, 'waterMass', waterMass, 'water.V', waterRho/waterMass)
        print('volume', volume, 'VTarget', VTarget, 'oilRho', oilRho, 'oilMass', oilMass, 'oil.V', oilRho/oilMass)

        # volume is set to density/mass
        water.V[:] = waterRho/waterMass
        oil.V[:] = oilRho/oilMass
        solid.V[:] = solidRho/solidMass

        # smoothing lengths
        water.h[:] = hdx * dx
        oil.h[:] = hdx * dx
        solid.h[:] = hdx * dx

        
        properties_to_save = ['pid', 'tag', 'gid', 'rho', 'V', 'h', 'm', 'p', 'u', 'w', 'v', 'y', 'x', 'div', 'z']

        for propertyName in materials.generateFullParticleProperties():
            properties_to_save.append(propertyName)

        print('properties to save', properties_to_save)

        water.align_particles()
        oil.align_particles()

        water.set_output_arrays(properties_to_save)
        oil.set_output_arrays(properties_to_save)

        # generate one particle array containing both water and oil
        water_and_oil = water
        water_and_oil.append_parray(oil)
        

        #print('fluid material amount: ', fluid.material_amount, 'fluid props', fluid.properties.keys())

        # return the particle list
        return [water_and_oil, solid]

    def create_solver(self):
        # Create the kernel
        # kernel = Gaussian(dim=2)
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


            # For the multi-phase formulation, we require an estimate of the
            # particle volume. This can be either defined from the particle
            # number density or simply as the ratio of mass to density
            Group(equations=[
                VolumeFromMassDensity(dest='fluid', sources=None)
                ], ),
            
            # Equation of state is typically the Tait EOS with a suitable
            # exponent gamma
            Group(equations=[
                TaitEOSHGCorrection(dest='fluid', sources=None, rho0=rho0, c0=c0, gamma=gamma),
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