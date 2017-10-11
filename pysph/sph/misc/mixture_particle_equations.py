# from pysph.sph.basic_equations import Equation
from pysph.sph.wc.basic import TaitEOSHGCorrection as TaitEOSHGCorrection
from pysph.sph.wc.transport_velocity import StateEquation
from pysph.sph.scheme import Scheme, TVFScheme
from pysph.sph.equation import Equation

class TaitEOSHGCorrectionVariableRho(TaitEOSHGCorrection):
    def __init__(self, dest, sources, c0, gamma):
        rho0 = 1000.0
        self.gamma = gamma
        self.rho0 = rho0
        self.rho01 = 1.0/rho0
        self.c0 = c0
        self.gamma = gamma
        self.gamma1 = 0.5*(gamma - 1.0)
        self.B = rho0*c0*c0/gamma
        super(TaitEOSHGCorrection, self).__init__(dest, sources)
        #TaitEOSHGCorrection.__init__(self, dest, sources, rho0, c0, gamma)
    
    def loop(self, d_idx, d_rho, d_rho00, d_p, d_cs):
        self.rho0 = d_rho00[d_idx]
        self.rho01 = 1.0/self.rho0
        #self.rho0 = 1000 + (d_rho0[d_idx] - 1000)*0.5
        #self.rho01 = 1.0/self.rho0
        self.B = self.rho0*self.c0*self.c0/self.gamma
        # Can't figure out how to use super, so copy-pasting for now...

        #super(self, d_idx, d_rho, d_rho0, d_p, d_cs)
        #super(TaitEOSHGCorrectionVariableRho, self).loop(self, d_idx, d_rho, d_p, d_cs)
        #self.loop(d_idx, d_rho, d_rho0, d_p, d_cs)
        #super(TaitEOSHGCorrection, self).loop()
        #super().loop(self, d_idx, d_rho, d_p, d_cs)

        if d_rho[d_idx] < self.rho0:
            d_rho[d_idx] = self.rho0

        ratio = d_rho[d_idx] * self.rho01
        tmp = pow(ratio, self.gamma)

        d_p[d_idx] = self.B * (tmp - 1.0)
        d_cs[d_idx] = self.c0 * pow( ratio, self.gamma1 )

class StateEquationVariableRho0(Equation): ## STATE EQUATION IS THE CULPRIT
    def __init__(self, dest, sources, b=1.0):
        self.b = b
        super(StateEquationVariableRho0, self).__init__(dest, sources)
    
    def loop(self, d_idx, d_p, d_p0, d_rho, d_rho0):
        d_p[d_idx] = d_p0[d_idx] * (d_rho[d_idx]/d_rho0[d_idx] - self.b)

class TVFSchemeMixtureParticle(TVFScheme):
    def __init__(self, fluids, solids, dim, rho0, c0, nu, p0, pb, h0, gx=0.0, gy=0.0, gz=0.0, alpha=0.0, tdamp=0.0):
        self.fluids = fluids
        self.solids = solids
        self.solver = None
        self.rho0 = rho0
        self.c0 = c0
        self.pb = pb
        self.p0 = p0
        self.nu = nu
        self.dim = dim
        self.h0 = h0
        self.gx = gx
        self.gy = gy
        self.gz = gz
        self.alpha = alpha
        self.tdamp = 0.0
    
    def get_equations(self): # fix copy pasta later
        from pysph.sph.equation import Group
        from pysph.sph.wc.transport_velocity import (
            SummationDensity, StateEquation, MomentumEquationPressureGradient,
            MomentumEquationArtificialViscosity,
            MomentumEquationViscosity, MomentumEquationArtificialStress,
            SolidWallPressureBC, SolidWallNoSlipBC, SetWallVelocity
        )
        equations = []
        all = self.fluids + self.solids
        g1 = []
        for fluid in self.fluids:
            g1.append(SummationDensity(dest=fluid, sources=all))

        equations.append(Group(equations=g1, real=False))

        g2 = []
        for fluid in self.fluids:
            g2.append(StateEquationVariableRho0(
                dest=fluid, sources=None, b=1.0
            ))
        for solid in self.solids:
            g2.append(SetWallVelocity(dest=solid, sources=self.fluids))

        equations.append(Group(equations=g2, real=False))

        g3 = []
        for solid in self.solids:
            g3.append(SolidWallPressureBC(
                dest=solid, sources=self.fluids, b=1.0, rho0=self.rho0,
                p0=self.p0, gx=self.gx, gy=self.gy, gz=self.gz
            ))

        equations.append(Group(equations=g3, real=False))

        g4 = []
        for fluid in self.fluids:
            g4.append(
                MomentumEquationPressureGradient(
                    dest=fluid, sources=all, pb=self.pb, gx=self.gx,
                    gy=self.gy, gz=self.gz, tdamp=self.tdamp
                )
            )
            if self.alpha > 0.0:
                g4.append(
                    MomentumEquationArtificialViscosity(
                        dest=fluid, sources=all, c0=self.c0,
                        alpha=self.alpha
                    )
                )
            if self.nu > 0.0:
                g4.append(
                    MomentumEquationViscosity(
                        dest=fluid, sources=self.fluids, nu=self.nu
                    )
                )
                if len(self.solids) > 0:
                    g4.append(
                        SolidWallNoSlipBC(
                            dest=fluid, sources=self.solids, nu=self.nu
                        )
                    )

            g4.append(
                MomentumEquationArtificialStress(
                    dest=fluid, sources=self.fluids)
            )

        equations.append(Group(equations=g4))
        return equations
