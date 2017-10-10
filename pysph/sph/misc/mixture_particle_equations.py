# from pysph.sph.basic_equations import Equation
from pysph.sph.wc.basic import TaitEOSHGCorrection as TaitEOSHGCorrection

class TaitEOSHGCorrectionVariableRho(TaitEOSHGCorrection):
    def __init__(self, dest, sources, c0, gamma):
        rho0 = 1000.0
        self.gamma = gamma
        TaitEOSHGCorrection.__init__(self, dest, sources, rho0, c0, gamma)
    
    def loop(self, d_idx, d_rho, d_rho0, d_p, d_cs):
        #self.rho0 = d_rho0[d_idx]
        #self.rho01 = 1.0/self.rho0
        #self.rho0 = 1000 + (d_rho0[d_idx] - 1000)*0.5
        #self.rho01 = 1.0/self.rho0
        #self.B = self.rho0*self.c0*self.c0/self.gamma
        # Can't figure out how to use super, so copy-pasting for now...

        #super(self, d_idx, d_rho, d_rho0, d_p, d_cs)
        #super(TaitEOSHGCorrectionVariableRho, self).loop(self, d_idx, d_rho, d_p, d_cs)
        #self.loop(d_idx, d_rho, d_rho0, d_p, d_cs)
        #super(TaitEOSHGCorrection, self).loop()
        #super().loop(self, d_idx, d_rho, d_p, d_cs)

        if d_rho[d_idx] < d_rho0[d_idx]:
            d_rho[d_idx] = d_rho[d_idx]

        ratio = d_rho[d_idx] * 1.0/1000.0
        tmp = pow(ratio, self.gamma)

        d_p[d_idx] = self.B * (tmp - 1.0)
        d_cs[d_idx] = self.c0 * pow( ratio, self.gamma1 )