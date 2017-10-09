from pysph.sph.basic_equations import Equation

class TaitEOSHGCorrectionVariableRho(object):
    def __init__(self, dest, sources, rho0, c0, gamma):
        super(equation, self).__init__()
    
    def loop(self, d_idx, d_rho, d_rho0, d_p, d_cs):
        self.rho0 = d_rho0[d_idx]
        self.rho1 = 1/self.rho0
        super(self, d_idx, d_rho, d_rho0, d_p, d_cs)