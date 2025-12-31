# fp_ainv = np.loadtxt( "../global_fit/coeffs_ainv_beta_mqhat.dat" )
# fcp_ainv = np.loadtxt( "../global_fit/cov_ainv_beta_mqhat.dat" )
import numpy as np

class F_ainv:
    # def fitter( self, x, c0, c1_0, c1_1, c2_00, c2_01, c2_11 ):
    #     ainv = x[0]
    #     M = x[1]
    #     return c0 + c1_0*ainv + c1_1*M + c2_00*ainv*ainv + c2_01*ainv*M + c2_11*M*M
    def fitter( self, x, c0, c1_0, c1_1, c2_01, c2_11 ):
        ainv = x[0]
        M = x[1]
        return c0 + c1_0*ainv + c1_1*M + c2_01*ainv*M + c2_11*M*M

    # fp_ainv = np.loadtxt( "../global_fit/coeffs_ainv_beta_mqhat.dat" )
    # fcp_ainv = np.loadtxt( "../global_fit/cov_ainv_beta_mqhat.dat" )

    def __init__( self ):
        self.fp_ainv = np.loadtxt( "../global_fit/coeffs_ainv_beta_mqhat.dat" )
        self.fcp_ainv = np.loadtxt( "../global_fit/cov_ainv_beta_mqhat.dat" )
        self.Df_ainvs = [self.Df0_ainv, self.Df1_ainv, self.Df2_ainv, self.Df3_ainv, self.Df4_ainv]

    def __call__( self, beta_, mq_ ):
        return self.fitter( [beta_, mq_], self.fp_ainv[0], self.fp_ainv[1], self.fp_ainv[2], self.fp_ainv[3], self.fp_ainv[4] )

    # check
    # f_ainv( 10.95, 0.19 )
    
    def Df0_ainv( self, beta_, mq_ ):
        return 1.0
    
    def Df1_ainv( self, beta_, mq_ ):
        return beta_
    
    def Df2_ainv( self, beta_, mq_ ):
        return mq_
    
    def Df3_ainv( self, beta_, mq_ ):
        return beta_*mq_
    
    def Df4_ainv( self, beta_, mq_ ):
        return mq_**2

    def D( self, beta_, mq_ ):
        df = np.array([ f(beta_,mq_) for f in self.Df_ainvs ])
        return np.sqrt( df@self.fcp_ainv@df )

class F_TmTc:
    def fitter( self, x, c1_0, c2_01 ):
        dbeta = x[0]
        mq = x[1]
        return c1_0*dbeta + c2_01*dbeta*mq

    # fp_ainv = np.loadtxt( "../global_fit/coeffs_ainv_beta_mqhat.dat" )
    # fcp_ainv = np.loadtxt( "../global_fit/cov_ainv_beta_mqhat.dat" )

    def __init__( self, Nt ):
        self.fp_ainv = np.loadtxt( "../global_fit/coeffs_ainv_beta_mqhat.dat" )
        self.fcp_ainv = np.loadtxt( "../global_fit/cov_ainv_beta_mqhat.dat" )
        self.Df_ainvs = [self.Df1_ainv, self.Df3_ainv]
        self.Nt = Nt

    def __call__( self, beta_, mq_ ):
        return self.fitter( [beta_, mq_], self.fp_ainv[1], self.fp_ainv[3] )/self.Nt

    # check
    # f_ainv( 10.95, 0.19 )
    
    def Df1_ainv( self, dbeta, mq ):
        return dbeta
        
    def Df3_ainv( self, dbeta, mq ):
        return dbeta*mq
    

    def D( self, dbeta, mq ):
        df = np.array([ f(dbeta, mq) for f in self.Df_ainvs ])
        # print( df.shape )
        cov = self.fcp_ainv[ [1,3] ].T[ [1,3] ]
        # print( cov.shape )
        return np.sqrt( df@ cov @df )/self.Nt

# check
# ii=3

# eps=1.0e-7
# beta_=10.95
# mq_=0.19

# a = Df_ainvs[ii]( beta_, mq_ )
# fp_ainvP = np.copy(fp_ainv)
# fp_ainvM = np.copy(fp_ainv)
# fp_ainvP[ii] += eps
# fp_ainvM[ii] -= eps

# b = ( fitter( np.array([beta_, mq_]), fp_ainvP[0], fp_ainvP[1], fp_ainvP[2], fp_ainvP[3], fp_ainvP[4], fp_ainvP[5] )-fitter( np.array([beta_, mq_]), fp_ainvM[0], fp_ainvM[1], fp_ainvM[2], fp_ainvM[3], fp_ainvM[4], fp_ainvM[5] ))/(2.0*eps)
# a, b



#####################################################


class F_MB:
    def fitter( self, x, c0, c1_0, c1_1, c2_01, c2_11 ):
        ainv = x[0]
        M = x[1]
        return c0 + c1_0*ainv + c1_1*M + c2_01*ainv*M + c2_11*M*M
    # fp_ainv = np.loadtxt( "../global_fit/coeffs_ainv_beta_mqhat.dat" )
    # fcp_ainv = np.loadtxt( "../global_fit/cov_ainv_beta_mqhat.dat" )

    def __init__( self ):
        self.fp_MB = np.loadtxt( "../global_fit/coeffs_MB_beta_mqhat.dat" )
        self.fcp_MB = np.loadtxt( "../global_fit/cov_MB_beta_mqhat.dat" )
        self.Df_MBs = [self.Df0_MB, self.Df1_MB, self.Df2_MB, self.Df4_MB, self.Df5_MB]

    def __call__( self, beta_, mq_ ):
        return self.fitter( [beta_, mq_], self.fp_MB[0], self.fp_MB[1], self.fp_MB[2], self.fp_MB[3], self.fp_MB[4] )

    # check
    # f_MB( 10.95, 0.19 )
    
    def Df0_MB( self, beta_, mq_ ):
        return 1.0
    
    def Df1_MB( self, beta_, mq_ ):
        return beta_
    
    def Df2_MB( self, beta_, mq_ ):
        return mq_
    
    # def Df3_MB( self, beta_, mq_ ):
    #     return beta_**2
    
    def Df4_MB( self, beta_, mq_ ):
        return beta_*mq_
    
    def Df5_MB( self, beta_, mq_ ):
        return mq_**2

    def D( self, beta_, mq_ ):
        df = np.array([ f(beta_,mq_) for f in self.Df_MBs ])
        return np.sqrt( df@self.fcp_MB@df )



# # check
# ii=2

# eps=1.0e-7
# beta_=10.95
# mq_=0.19

# a = Df_MBs[ii]( beta_, mq_ )
# fp_MBP = np.copy(fp_MB)
# fp_MBM = np.copy(fp_MB)
# fp_MBP[ii] += eps
# fp_MBM[ii] -= eps

# b = ( fitter( np.array([beta_, mq_]), fp_MBP[0], fp_MBP[1], fp_MBP[2], fp_MBP[3], fp_MBP[4], fp_MBP[5] )-fitter( np.array([beta_, mq_]), fp_MBM[0], fp_MBM[1], fp_MBM[2], fp_MBM[3], fp_MBM[4], fp_MBM[5] ))/(2.0*eps)
# a, b


# def Df_MB( beta_, mq_ ):
#     df = np.array([ f(beta_,mq_) for f in Df_MBs ])
#     return np.sqrt( df@fcp_MB@df )