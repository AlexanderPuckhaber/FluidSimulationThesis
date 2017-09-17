"""
Diffusion equations
###################
"""

from pysph.sph.equation import Equation
from textwrap import dedent

class DiffuseMaterial(Equation):
    r"""Diffusion of materials. Yeah
    """

    def __init__(self, dest, sources, diffusion_speed, dt):
        r"""
        Parameters
        ----------
        DiffusionSpeed: float. Need to find units/actual equations
        """

        self.diffusion_speed = diffusion_speed
        self.dt = dt

        super(DiffuseMaterial, self).__init__(dest, sources)

    def initialize(self, d_idx, d_material_velocity):
        d_material_velocity[d_idx] = 0.0
    
    def loop(self, d_idx, s_idx, d_material_velocity, s_material_velocity, d_material_amount, s_material_amount):
        material_gradient = (s_material_amount[s_idx] - d_material_amount[d_idx])
        d_material_velocity[d_idx] += material_gradient*self.diffusion_speed
        s_material_velocity[s_idx] -= material_gradient*self.diffusion_speed

    def post_loop(self, d_idx, d_material_amount, d_material_velocity):
        d_material_amount[d_idx] += d_material_velocity[d_idx] * self.dt

