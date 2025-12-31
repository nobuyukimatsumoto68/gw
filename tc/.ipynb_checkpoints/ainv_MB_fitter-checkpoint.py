# fp_ainv = np.loadtxt( "../global_fit/coeffs_ainv_beta_mqhat.dat" )
# fcp_ainv = np.loadtxt( "../global_fit/cov_ainv_beta_mqhat.dat" )

class f_ainv:
    def fitter( x, c0, c1_0, c1_1, c2_00, c2_01, c2_11 ):
        ainv = x[0]
        M = x[1]
        return c0 + c1_0*ainv + c1_1*M + c2_00*ainv*ainv + c2_01*ainv*M + c2_11*M*M
    # fp_ainv = np.loadtxt( "../global_fit/coeffs_ainv_beta_mqhat.dat" )
    # fcp_ainv = np.loadtxt( "../global_fit/cov_ainv_beta_mqhat.dat" )

    def __init__( self ):
        self.fp_ainv = np.loadtxt( "../global_fit/coeffs_ainv_beta_mqhat.dat" )
        self.fcp_ainv = np.loadtxt( "../global_fit/cov_ainv_beta_mqhat.dat" )
        self.Df_ainvs = [self.Df0_ainv, self.Df1_ainv, self.Df2_ainv, self.Df3_ainv, self.Df4_ainv, self.Df5_ainv]

    def __call__( beta_, mq_ ):
        return self.fitter( [beta_, mq_], self.fp_ainv[0], self.fp_ainv[1], self.fp_ainv[2], self.fp_ainv[3], self.fp_ainv[4], self.fp_ainv[5] )

    # check
    # f_ainv( 10.95, 0.19 )
    
    def Df0_ainv( beta_, mq_ ):
        return 1.0
    
    def Df1_ainv( beta_, mq_ ):
        return beta_
    
    def Df2_ainv( beta_, mq_ ):
        return mq_
    
    def Df3_ainv( beta_, mq_ ):
        return beta_**2
    
    def Df4_ainv( beta_, mq_ ):
        return beta_*mq_
    
    def Df5_ainv( beta_, mq_ ):
        return mq_**2

    def D( beta_, mq_ ):
        df = np.array([ f(beta_,mq_) for f in self.Df_ainvs ])
        return np.sqrt( df@self.fcp_ainv@df )

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


class f_MB:
    def fitter( x, c0, c1_0, c1_1, c2_00, c2_01, c2_11 ):
        ainv = x[0]
        M = x[1]
        return c0 + c1_0*ainv + c1_1*M + c2_00*ainv*ainv + c2_01*ainv*M + c2_11*M*M
    # fp_ainv = np.loadtxt( "../global_fit/coeffs_ainv_beta_mqhat.dat" )
    # fcp_ainv = np.loadtxt( "../global_fit/cov_ainv_beta_mqhat.dat" )

    def __init__( self ):
        self.fp_MB = np.loadtxt( "../global_fit/coeffs_MB_beta_mqhat.dat" )
        self.fcp_MB = np.loadtxt( "../global_fit/cov_MB_beta_mqhat.dat" )
        self.Df_MBs = [self.Df0_MB, self.Df1_MB, self.Df2_MB, self.Df3_ainv, self.Df4_ainv, self.Df5_ainv]

    def __call__( beta_, mq_ ):
        return self.fitter( [beta_, mq_], self.fp_MB[0], self.fp_MB[1], self.fp_MB[2], self.fp_MB[3], self.fp_MB[4], self.fp_MB[5] )

    # check
    # f_MB( 10.95, 0.19 )
    
    def Df0_MB( beta_, mq_ ):
        return 1.0
    
    def Df1_MB( beta_, mq_ ):
        return beta_
    
    def Df2_MB( beta_, mq_ ):
        return mq_
    
    def Df3_MB( beta_, mq_ ):
        return beta_**2
    
    def Df4_MB( beta_, mq_ ):
        return beta_*mq_
    
    def Df5_MB( beta_, mq_ ):
        return mq_**2

    def D( beta_, mq_ ):
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