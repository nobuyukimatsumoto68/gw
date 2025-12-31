import numpy as np

# class F_S3hat:
#     def fitter(self, TmTc_MB, Mc, c1, c2, gam ):
#         TmTc = TmTc_MB[0]
#         MB = TmTc_MB[1]
#         return (MB-Mc)**gam  * ( c1/TmTc + c2/TmTc**2 )

#     def d0_fitter(self, TmTc_MB, Mc, c1, c2, gam ):
#         TmTc = TmTc_MB[0]
#         MB = TmTc_MB[1]
#         return (MB-Mc)**gam  * ( -1.0*c1/TmTc**2 -2.0 * c2/TmTc**3 )

#     def d1_fitter(self, TmTc_MB, Mc, c1, c2, gam ):
#         TmTc = TmTc_MB[0]
#         MB = TmTc_MB[1]
#         return gam*(MB-Mc)**(gam-1.0)  * ( c1/TmTc + c2/TmTc**2 )

#     def __init__( self ):
#         self.fp = np.loadtxt( "../S3/coeffs_S3hat_TmTc_MB.dat" )
#         self.fcp = np.loadtxt( "../S3/cov_S3hat_TmTc_MB.dat" )
#         self.Dfs = [self.Df0, self.Df1, self.Df2, self.Df3]

#     def __call__( self, TmTc, MB ):
#         return self.fitter( [TmTc, MB], self.fp[0], self.fp[1], self.fp[2], self.fp[3] )

#     def Df0( self, TmTc, MB ):
#         Mc = self.fp[0]
#         a0 = self.fp[1]
#         b0 = self.fp[2]
#         gamma = self.fp[3]
#         return -gamma*(MB-Mc)**(gamma-1.0)*( a0/TmTc + b0/TmTc**2 )
    
#     def Df1( self, TmTc, MB ):
#         Mc = self.fp[0]
#         a0 = self.fp[1]
#         b0 = self.fp[2]
#         gamma = self.fp[3]
#         return (MB-Mc)**gamma/TmTc

#     def Df2( self, TmTc, MB ):
#         Mc = self.fp[0]
#         a0 = self.fp[1]
#         b0 = self.fp[2]
#         gamma = self.fp[3]
#         return (MB-Mc)**gamma/TmTc**2
    
#     def Df3( self, TmTc, MB ):
#         Mc = self.fp[0]
#         a0 = self.fp[1]
#         b0 = self.fp[2]
#         gamma = self.fp[3]
#         return np.log(MB-Mc)*(MB-Mc)**gamma*( a0/TmTc + b0/TmTc**2 )

#     def D( self, TmTc, MB ):
#         df = np.array([ f(TmTc, MB) for f in self.Dfs ])
#         return np.sqrt( df@self.fcp@df )


# class F_S3hat:
#     def fitter(self, TmTc_MB, Mc, c1, c2, gam ):
#         TmTc = TmTc_MB[0]
#         MB = TmTc_MB[1]
#         return (MB-Mc)**gam  * ( c1/TmTc + c2/TmTc**2 )

#     def d0_fitter(self, TmTc_MB, Mc, c1, c2, gam ):
#         TmTc = TmTc_MB[0]
#         MB = TmTc_MB[1]
#         return (MB-Mc)**gam  * ( -1.0*c1/TmTc**2 -2.0 * c2/TmTc**3 )

#     def d1_fitter(self, TmTc_MB, Mc, c1, c2, gam ):
#         TmTc = TmTc_MB[0]
#         MB = TmTc_MB[1]
#         return gam*(MB-Mc)**(gam-1.0)  * ( c1/TmTc + c2/TmTc**2 )

#     def __init__( self ):
#         self.fp = np.loadtxt( "../S3/coeffs_S3hat_TmTc_MB.dat" )
#         self.fcp = np.loadtxt( "../S3/cov_S3hat_TmTc_MB.dat" )
#         self.Dfs = [self.Df0, self.Df1, self.Df2, self.Df3]

#     def __call__( self, TmTc, MB ):
#         return self.fitter( [TmTc, MB], self.fp[0], self.fp[1], self.fp[2], self.fp[3] )

#     def Df0( self, TmTc, MB ):
#         Mc = self.fp[0]
#         a0 = self.fp[1]
#         b0 = self.fp[2]
#         gamma = self.fp[3]
#         return -gamma*(MB-Mc)**(gamma-1.0)*( a0/TmTc + b0/TmTc**2 )
    
#     def Df1( self, TmTc, MB ):
#         Mc = self.fp[0]
#         a0 = self.fp[1]
#         b0 = self.fp[2]
#         gamma = self.fp[3]
#         return (MB-Mc)**gamma/TmTc

#     def Df2( self, TmTc, MB ):
#         Mc = self.fp[0]
#         a0 = self.fp[1]
#         b0 = self.fp[2]
#         gamma = self.fp[3]
#         return (MB-Mc)**gamma/TmTc**2
    
#     def Df3( self, TmTc, MB ):
#         Mc = self.fp[0]
#         a0 = self.fp[1]
#         b0 = self.fp[2]
#         gamma = self.fp[3]
#         return np.log(MB-Mc)*(MB-Mc)**gamma*( a0/TmTc + b0/TmTc**2 )

#     def D( self, TmTc, MB ):
#         df = np.array([ f(TmTc, MB) for f in self.Dfs ])
#         return np.sqrt( df@self.fcp@df )

class F_S3hat:
    def fitter(self, TmTc_MB, Mc, c0, c1, c2, gam ):
        TmTc = TmTc_MB[0]
        MB = TmTc_MB[1]
        return (MB-Mc)**gam  * ( c0 + c1/TmTc + c2/TmTc**2 )

    def d0_fitter(self, TmTc_MB, Mc, c0, c1, c2, gam ):
        TmTc = TmTc_MB[0]
        MB = TmTc_MB[1]
        return (MB-Mc)**gam  * ( -1.0*c1/TmTc**2 -2.0 * c2/TmTc**3 )

    def d1_fitter(self, TmTc_MB, Mc, c0, c1, c2, gam ):
        TmTc = TmTc_MB[0]
        MB = TmTc_MB[1]
        return gam*(MB-Mc)**(gam-1.0)  * (c0 + c1/TmTc + c2/TmTc**2 )

    def __init__( self, is_variate=False, eps=1.0 ):
        self.fp = np.loadtxt( "../S3/coeffs_S3hat_TmTc_MB.dat" )
        self.fcp = np.loadtxt( "../S3/cov_S3hat_TmTc_MB.dat" )
        if is_variate:
            self.fp = np.random.multivariate_normal(self.fp, eps**2 * self.fcp, size=None, check_valid='warn', tol=1e-8)
        self.Dfs = [self.Df0, self.Df1, self.Df2, self.Df3, self.Df4]

    def __call__( self, TmTc, MB ):
        return self.fitter( [TmTc, MB], self.fp[0], self.fp[1], self.fp[2], self.fp[3], self.fp[4] )

    def Df0( self, TmTc, MB ):
        Mc = self.fp[0]
        c0 = self.fp[1]
        c1 = self.fp[2]
        c2 = self.fp[3]
        gam = self.fp[4]
        return -gam*(MB-Mc)**(gam-1.0)*( c0+ c1/TmTc + c2/TmTc**2 )

    def Df1( self, TmTc, MB ):
        Mc = self.fp[0]
        c0 = self.fp[1]
        c1 = self.fp[2]
        c2 = self.fp[3]
        gam = self.fp[4]
        return (MB-Mc)**gam

    def Df2( self, TmTc, MB ):
        Mc = self.fp[0]
        c0 = self.fp[1]
        c1 = self.fp[2]
        c2 = self.fp[3]
        gam = self.fp[4]
        return (MB-Mc)**gam/TmTc

    def Df3( self, TmTc, MB ):
        Mc = self.fp[0]
        c0 = self.fp[1]
        c1 = self.fp[2]
        c2 = self.fp[3]
        gam = self.fp[4]
        return (MB-Mc)**gam/TmTc**2
    
    def Df4( self, TmTc, MB ):
        Mc = self.fp[0]
        c0 = self.fp[1]
        c1 = self.fp[2]
        c2 = self.fp[3]
        gam = self.fp[4]
        return np.log(MB-Mc)*(MB-Mc)**gam*( c0 + c1/TmTc + c2/TmTc**2 )

    def D( self, TmTc, MB ):
        df = np.array([ f(TmTc, MB) for f in self.Dfs ])
        return np.sqrt( df@self.fcp@df )

    def dT( self, TmTc, MB ):
        Mc = self.fp[0]
        c0 = self.fp[1]
        c1 = self.fp[2]
        c2 = self.fp[3]
        gam = self.fp[4]
        return (MB-Mc)**gam  * ( -1.0*c1/TmTc**2 -2.0 * c2/TmTc**3 )