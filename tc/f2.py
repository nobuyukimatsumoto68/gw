import numpy as np

class F_Tc:
    def fitter( self, MB, b0, b1, b2 ):
        return b0 + b1*MB + b2*MB**2

    def dfitter( self, MB, b0, b1, b2 ):
        return b1 + 2.0*b2*MB

    def __init__( self, is_w_global, is_variate=False, eps=1.0 ):
        if is_w_global:
            self.fp = np.loadtxt( "../tc/coeffs_Tc_MB_w.dat" )
            self.fcp = np.loadtxt( "../tc/cov_Tc_MB_w.dat" )
        else:
            self.fp = np.loadtxt( "../tc/coeffs_Tc_MB_wo.dat" )
            self.fcp = np.loadtxt( "../tc/cov_Tc_MB_wo.dat" )

        if is_variate:
            self.fp = np.random.multivariate_normal(self.fp, eps**2 * self.fcp, size=None, check_valid='warn', tol=1e-8)
        self.Dfs = [self.Df0, self.Df1, self.Df2]

    def __call__( self, MB ):
        return self.fitter( MB, self.fp[0], self.fp[1], self.fp[2] )
    
    def Df0( self, MB ):
        return 1.0
    
    def Df1( self, MB ):
        return MB

    def Df2( self, MB ):
        return MB**2

    def d( self, MB ):
        return self.dfitter( MB, self.fp[0], self.fp[1], self.fp[2] )
    
    def D( self, MB ):
        df = np.array([ f(MB) for f in self.Dfs ])
        return np.sqrt( df@self.fcp@df )

# f_Tc = F_Tc()

# eps=1.0e-6

# i=1

# fpp = np.copy(f_Tc.fp)
# fpm = np.copy(f_Tc.fp)

# fpp[i]+=eps
# fpm[i]-=eps

# MB=2.4

# f_Tc( MB )

# fp = f_Tc.fitter( MB, fpp[0], fpp[1] )
# fm = f_Tc.fitter( MB, fpm[0], fpm[1] )
# (fp-fm)/(2.0*eps), f_Tc.Dfs[i](MB)