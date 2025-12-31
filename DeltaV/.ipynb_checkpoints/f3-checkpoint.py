import numpy as np

class F_DeltaV:
    def fitter(self, TmTc_MB, Mc, a0, b0, gamma):
        TmTc = TmTc_MB[0]
        MB = TmTc_MB[1]
        return (MB-Mc)**gamma*( a0*TmTc + b0*TmTc**2 )

    def d0_fitter(self, TmTc_MB, Mc, a0, b0, gamma):
        TmTc = TmTc_MB[0]
        MB = TmTc_MB[1]
        return (MB-Mc)**gamma*( a0 + 2.0*b0*TmTc )

    def d1_fitter(self, TmTc_MB, Mc, a0, b0, gamma):
        TmTc = TmTc_MB[0]
        MB = TmTc_MB[1]
        return gamma*(MB-Mc)**(gamma-1.0)*( a0*TmTc + b0*TmTc**2 )

    def __init__( self ):
        self.fp = np.loadtxt( "../DeltaV/coeffs_DeltaV_TmTc_MB.dat" )
        self.fcp = np.loadtxt( "../DeltaV/cov_DeltaV_TmTc_MB.dat" )
        self.Dfs = [self.Df0, self.Df1, self.Df2, self.Df3]

    def __call__( self, TmTc, MB ):
        return self.fitter( [TmTc, MB], self.fp[0], self.fp[1], self.fp[2], self.fp[3] )

    def Df0( self, TmTc, MB ):
        Mc = self.fp[0]
        a0 = self.fp[1]
        b0 = self.fp[2]
        gamma = self.fp[3]
        return -gamma*(MB-Mc)**(gamma-1.0)*( a0*TmTc + b0*TmTc**2 )
    
    def Df1( self, TmTc, MB ):
        Mc = self.fp[0]
        a0 = self.fp[1]
        b0 = self.fp[2]
        gamma = self.fp[3]
        return (MB-Mc)**gamma*TmTc

    def Df2( self, TmTc, MB ):
        Mc = self.fp[0]
        a0 = self.fp[1]
        b0 = self.fp[2]
        gamma = self.fp[3]
        return (MB-Mc)**gamma*TmTc**2
    
    def Df3( self, TmTc, MB ):
        Mc = self.fp[0]
        a0 = self.fp[1]
        b0 = self.fp[2]
        gamma = self.fp[3]
        return np.log(MB-Mc)*(MB-Mc)**gamma*( a0*TmTc + b0*TmTc**2 )

    def D( self, TmTc, MB ):
        df = np.array([ f(TmTc, MB) for f in self.Dfs ])
        return np.sqrt( df@self.fcp@df )