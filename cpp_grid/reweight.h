#pragma once

#include <Grid/Grid.h>
#include <algorithm>
#include <random>
#include <filesystem>
#include <limits>


using namespace std;
using namespace Grid;


std::string obsinfo="polyakov";

auto get_prefix = [](const std::string& beta, const std::string& mass, const int type=1){
                    if(type==1) return "ckpoint_lat";
                    else assert(false);
                  };

// type=2 : sungwoo, 24c, cold
auto get_configname = [](const std::string& beta, const std::string& mass, const int type=1, const std::string& lat="248"){
                        if(type==1){
                          char id[200];
                          std::sprintf(id,
                                       "beta%sm%s",
                                       beta.data(), mass.data());
                          std::string f_id = id;
                          return f_id;
                        }
                        else if(type==2) {
                          // return "conf_nc4nf1_"+lat+"_b"+beta+"c0_m"+mass;
                          return "b"+beta+"c0_m"+mass;
                        }
                        else if(type==3) {
                          // return "conf_nc4nf1_"+lat+"_b"+beta+"h0_m"+mass;
                          return "b"+beta+"h0_m"+mass;
                        }
                        else if(type==4) {
                          // return "conf_nc4nf1_"+lat+"_b"+beta+"c0_m"+mass;
                          return "b"+beta+"c_m"+mass;
                        }
                        else if(type==5) {
                          // return "conf_nc4nf1_"+lat+"_b"+beta+"h0_m"+mass;
                          return "b"+beta+"_m"+mass;
                        }
                        else assert(false);
                      };
// auto get_configname = [](const Real beta, const Real mass){
//                         if(type==1){
//                           char id[100];
//                           if(mass_digits==1){
//                             std::sprintf(id,
//                                          "beta%.2fm%.1f",
//                                          beta, mass);
//                           }
//                           else if(mass_digits==2){
//                             std::sprintf(id,
//                                          "beta%.2fm%.2f",
//                                          beta, mass);
//                           }
//                           std::string f_id = id;
//                           return f_id;
//                         }
//                         else assert(false);
//                       };






// std::string base_dir = "/p/lustre1/matsumoto5/16c/";
// std::string basedir="/usr/workspace/lsd/matsumoto5/su4_dane_analysis/reweight/data/16c/";
// std::string basedir2 = basedir;
// std::string base_dir = "/p/lustre1/matsumoto5/16c_h/";
// std::string basedir="/usr/workspace/lsd/matsumoto5/su4_dane_analysis/reweight/data/16c_h/";

// Real mass=0.4; // 0.4
// const int interval=20;
// std::vector<Real> betas({11.01, 11.02, 11.03, 11.04, 11.05, 11.06, 11.07, 11.08, 11.09});
// std::string obsinfo="polyakov";

// Real mass=0.05;
// const int interval=20;
// std::vector<Real> betas({10.88, 10.89, 10.9, 10.91, 10.92, 10.93, 10.94, 10.95, 10.96});
// std::string obsinfo="polyakov";
// int mass_digits=2;


// Real mass=0.15;
// const int interval=20;
// std::vector<Real> betas({10.94, 10.95, 10.96, 10.97, 10.98, 10.99, 11.00, 11.01, 11.02});
// std::string obsinfo="polyakov";
// int mass_digits=2;



// Real mass=0.3;
// const int interval=20;
// std::vector<Real> betas({11.00, 11.01, 11.02, 11.03, 11.04, 11.05, 11.06, 11.07, 11.08});
// std::string obsinfo="polyakov";

// Real mass=0.2;
// const int interval=20;
// std::vector<Real> betas({10.96, 10.97, 10.98, 10.99, 11.00, 11.01, 11.02, 11.03, 11.04});
// std::string obsinfo="polyakov";


// auto get_prefix = [](const Real beta, const Real mass){
//                     return "ckpoint_lat";
//                   };

// auto get_configname = [](const Real beta, const Real mass){
//                         char id[200];
//                         std::sprintf(id,
//                                      "beta%.2fm%.1f",
//                                      beta, mass);
//                         std::string f_id = id;
//                         return f_id;
//                       };
// auto get_configname = [](const Real beta, const Real mass){
//                         char id[100];
//                         if(mass_digits==1){
//                           std::sprintf(id,
//                                        "beta%.2fm%.1f",
//                                        beta, mass);
//                         }
//                         else if(mass_digits==2){
//                           std::sprintf(id,
//                                        "beta%.2fm%.2f",
//                                        beta, mass);
//                         }
//                         std::string f_id = id;
//                         return f_id;
//                       };

// --------------------------------


// std::string base_dir = "/p/lustre1/matsumoto5/16c/";
// std::string basedir="/usr/workspace/lsd/matsumoto5/su4_dane_analysis/reweight/data/16c/";
// // std::string base_dir = "/p/lustre1/matsumoto5/16c_h/";
// // std::string basedir="/usr/workspace/lsd/matsumoto5/su4_dane_analysis/reweight/data/16c_h/";

// // Real mass=0.4; // 0.4
// // const int interval=20;
// // std::vector<Real> betas({11.01, 11.02, 11.03, 11.04, 11.05, 11.06, 11.07, 11.08, 11.09});
// // std::string obsinfo="polyakov";

// Real mass=0.35;
// const int interval=20;
// std::vector<Real> betas({11.00, 11.01, 11.02, 11.03, 11.04, 11.05, 11.06, 11.07, 11.08});
// std::string obsinfo="polyakov";

// // Real mass=0.3;
// // const int interval=20;
// // std::vector<Real> betas({11.00, 11.01, 11.02, 11.03, 11.04, 11.05, 11.06, 11.07, 11.08});
// // std::string obsinfo="polyakov";

// // Real mass=0.2;
// // const int interval=20;
// // std::vector<Real> betas({10.96, 10.97, 10.98, 10.99, 11.00, 11.01, 11.02, 11.03, 11.04});
// // std::string obsinfo="polyakov";


// auto get_prefix = [](const Real beta, const Real mass){
//                     return "ckpoint_lat";
//                   };

// // auto get_configname = [](const Real beta, const Real mass){
// //                         char id[200];
// //                         std::sprintf(id,
// //                                      "beta%.2fm%.1f",
// //                                      beta, mass);
// //                         std::string f_id = id;
// //                         return f_id;
// //                       };
// auto get_configname = [](const Real beta, const Real mass){
//                         char id[100];
//                         std::sprintf(id,
//                                      "beta%.2fm%.2f",
//                                      beta, mass);
//                         std::string f_id = id;
//                         return f_id;
//                       };




// ------------------------------

// Real mass=0.4;
// const int interval=10;
// std::vector<Real> betas({
//                          // 11.000, 11.005,
//                          11.010, 11.015,
//                          11.020, 11.025,
//                          11.030, 11.035,
//                          11.040
//                          , 11.050
//                          // , 11.1, 11.2
//   });
// std::string obsinfo="pl";

////

// std::string base_dir = "/p/lustre2/matsumoto5/24c/";
// std::string basedir="/usr/workspace/lsd/matsumoto5/su4_dane_analysis/reweight/data/24c/";

// auto get_prefix = [](const Real beta, const Real mass){
//                     char id[200];
//                     std::sprintf(id,
//                                  "conf_nc4nf1_248_b%.3fc0_m%.4f_lat",
//                                  beta, mass);
//                     std::string f_id = id;
//                     std::replace( f_id.begin(), f_id.end(), '.', 'p');
//                     return f_id;
//                   };

// auto get_configname = [](const Real beta, const Real mass){
//                         char id[200];
//                         std::sprintf(id,
//                                      "conf_nc4nf1_248_b%.3fc0_m%.4f",
//                                      beta, mass);
//                         std::string f_id = id;
//                         std::replace( f_id.begin(), f_id.end(), '.', 'p');
//                         return f_id;
//                       };

////

// std::string base_dir = "/p/lustre2/matsumoto5/24h/";
// std::string basedir="/usr/workspace/lsd/matsumoto5/su4_dane_analysis/reweight/data/24c_h/";

// auto get_prefix = [](const Real beta, const Real mass){
//                     char id[200];
//                     std::sprintf(id,
//                                  "conf_nc4nf1_248_b%.3fh0_m%.4f_lat",
//                                  beta, mass);
//                     std::string f_id = id;
//                     std::replace( f_id.begin(), f_id.end(), '.', 'p');
//                     return f_id;
//                   };

// auto get_configname = [](const Real beta, const Real mass){
//                         char id[200];
//                         std::sprintf(id,
//                                      "conf_nc4nf1_248_b%.3fh0_m%.4f",
//                                      beta, mass);
//                         std::string f_id = id;
//                         std::replace( f_id.begin(), f_id.end(), '.', 'p');
//                         return f_id;
//                       };


// // ------------------------------

// Real mass=0.4;
// const int interval=10;
// std::vector<Real> betas({11.000,
//                          11.010, 11.015,
//                          11.020, 11.025,
//                          11.030, 11.035,
//                          11.040,
//                          11.050
//                          // ,
//                          // 11.1, 11.2
//   });
// std::string obsinfo="polyakov";

// ////

// std::string base_dir = "/p/lustre2/matsumoto5/32c/";
// std::string basedir="/usr/workspace/lsd/matsumoto5/su4_dane_analysis/reweight/data/32c/";

// auto get_prefix = [](const Real beta, const Real mass){
//                     char id[200];
//                     std::sprintf(id,
//                                  "conf_nc4nf1_328_b%.3fc0_m%.4f_lat",
//                                  beta, mass);
//                     std::string f_id = id;
//                     std::replace( f_id.begin(), f_id.end(), '.', 'p');
//                     return f_id;
//                   };

// auto get_configname = [](const Real beta, const Real mass){
//                         char id[200];
//                         std::sprintf(id,
//                                      "conf_nc4nf1_328_b%.3fc0_m%.4f",
//                                      beta, mass);
//                         std::string f_id = id;
//                         std::replace( f_id.begin(), f_id.end(), '.', 'p');
//                         return f_id;
//                       };

// ////

// // std::string base_dir = "/p/lustre2/matsumoto5/32h/";
// // std::string basedir="/usr/workspace/lsd/matsumoto5/su4_dane_analysis/reweight/data/32c_h/";

// // auto get_prefix = [](const Real beta, const Real mass){
// //                     char id[200];
// //                     std::sprintf(id,
// //                                  "conf_nc4nf1_328_b%.3fh0_m%.4f_lat",
// //                                  beta, mass);
// //                     std::string f_id = id;
// //                     std::replace( f_id.begin(), f_id.end(), '.', 'p');
// //                     return f_id;
// //                   };

// // auto get_configname = [](const Real beta, const Real mass){
// //                         char id[200];
// //                         std::sprintf(id,
// //                                      "conf_nc4nf1_328_b%.3fh0_m%.4f",
// //                                      beta, mass);
// //                         std::string f_id = id;
// //                         std::replace( f_id.begin(), f_id.end(), '.', 'p');
// //                         return f_id;
// //                       };



// // ------------------------------

// Real mass=0.2;
// const int interval=20;
// std::vector<Real> betas({// 10.90,
//                          // 10.95,
//                          10.960,
//                          10.965,
//                          10.970,
//                          10.975,
//                          10.980,
//                          10.985,
//                          10.990// ,
//                          // 10.99
//   });
// std::string obsinfo="polyakov";

// ////

// std::string base_dir = "/usr/WS2/lsd/sungwoo/SU4_sdm/run_getobs/";
// std::string basedir="/usr/workspace/lsd/matsumoto5/su4_dane_analysis/reweight/data/24c/";
// std::string basedir2 = base_dir;

// auto get_prefix = [](const Real beta, const Real mass){
//                     char id[200];
//                     std::sprintf(id,
//                                  "conf_nc4nf1_248_b%.3fc_m%.4f_lat",
//                                  beta, mass);
//                     std::string f_id = id;
//                     std::replace( f_id.begin(), f_id.end(), '.', 'p');
//                     return f_id;
//                   };

// auto get_configname = [](const Real beta, const Real mass){
//                         char id[200];
//                         std::sprintf(id,
//                                      "conf_nc4nf1_248_b%.3fc_m%.4f",
//                                      beta, mass);
//                         std::string f_id = id;
//                         std::replace( f_id.begin(), f_id.end(), '.', 'p');
//                         return f_id;
//                       };

// ////



// std::string base_dir = "/p/lustre2/matsumoto5/32h/";
// std::string basedir="/usr/workspace/lsd/matsumoto5/su4_dane_analysis/reweight/data/32c_h/";

// auto get_prefix = [](const Real beta, const Real mass){
//                     char id[200];
//                     std::sprintf(id,
//                                  "conf_nc4nf1_328_b%.3fh0_m%.4f_lat",
//                                  beta, mass);
//                     std::string f_id = id;
//                     std::replace( f_id.begin(), f_id.end(), '.', 'p');
//                     return f_id;
//                   };

// auto get_configname = [](const Real beta, const Real mass){
//                         char id[200];
//                         std::sprintf(id,
//                                      "conf_nc4nf1_328_b%.3fh0_m%.4f",
//                                      beta, mass);
//                         std::string f_id = id;
//                         std::replace( f_id.begin(), f_id.end(), '.', 'p');
//                         return f_id;
//                       };


// ######################################################################


// std::string base_dir = "/p/lustre1/matsumoto5/16c/";
// std::string basedir="/usr/workspace/lsd/matsumoto5/su4_dane_analysis/reweight/data/16c/";
// // std::string base_dir = "/p/lustre1/matsumoto5/16c_h/";
// // std::string basedir="/usr/workspace/lsd/matsumoto5/su4_dane_analysis/reweight/data/16c_h/";

// Real mass=0.35;
// const int interval=20;
// std::vector<Real> betas({11.00, 11.01, 11.02, 11.03, 11.04, 11.05, 11.06, 11.07, 11.08});
// std::string obsinfo="polyakov";

// auto get_prefix = [](const Real beta, const Real mass){
//                     return "ckpoint_lat";
//                   };

// auto get_configname = [](const Real beta, const Real mass){
//                         char id[200];
//                         std::sprintf(id,
//                                      "beta%.2fm%.2f",
//                                      beta, mass);
//                         std::string f_id = id;
//                         return f_id;
//                       };





class Obs : public Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(
                                  Obs,
                                  std::string, beta,
                                  //
                                  int, conf,
                                  //
                                  std::string, base_dir,
                                  std::string, id,
                                  //
                                  RealD, energy,
                                  RealD, O
                                  );
};




class EnsembleOfObs : public Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(
                                  EnsembleOfObs,
                                  std::string, beta,
                                  //
                                  // int, conf_min,
                                  int, nconf,
                                  // int, interval,
                                  int, binsize,
                                  int, nbin,
                                  //
                                  std::string, base_dir,
                                  std::string, id,
                                  //
                                  std::vector<RealD>, energies,
                                  std::vector<RealD>, Os,
                                  std::vector<RealD>, Osqs
                                  );

  EnsembleOfObs& operator=(const EnsembleOfObs& ens){
    beta = ens.beta;
    // conf_min = ens.conf_min;
    nconf = ens.nconf;
    // interval = ens.interval;
    binsize = ens.binsize;
    nbin = ens.nbin;
    base_dir = ens.base_dir;
    id = ens.id;
    energies = ens.energies;
    Os = ens.Os;
    Osqs = ens.Osqs;

    return *this;
  }

  int size() const {
    // return (conf_max-conf_min)/interval;
    assert( energies.size()==Os.size() );
    return energies.size();
  }

  void jackknife_drop(const int nbin, const int ibin) {
    assert( nbin * (int)(size()/nbin) - size() == 0 );
    const int start = size()/nbin * ibin;
    const int end = size()/nbin * (ibin+1);

    energies.erase( energies.begin()+start, energies.begin()+end );
    Os.erase( Os.begin()+start, Os.begin()+end );
    Osqs.erase( Osqs.begin()+start, Osqs.begin()+end );
  }
};



RealD regularized_delta( const Real O, const Real x, const Real delta ){
  const Real expn = 0.5/delta/delta * (O-x)*(O-x);
  return std::exp( -expn ) / delta / std::sqrt(2.0 * M_PI);
}


RealD mean( const std::vector<RealD>& Os ){
  RealD res=0.0;
  for(RealD elem : Os) res += elem;
  res /= Os.size();
  return res;
}





// Geoge's email on 2/18/2025
class FreeEnergies : public Serializable {
public:
  GRID_SERIALIZABLE_CLASS_MEMBERS(
                                  FreeEnergies,
                                  std::vector<RealD>, f
                                  );
  std::vector<RealD> fnew;
  std::vector<RealD> fold;

  int nparam;

  FreeEnergies(const int nparam)
    : nparam(nparam)
    , f(nparam, 0.0)
    , fnew(nparam, 0.0)
    , fold(nparam, 0.0)
  {}

  FreeEnergies& operator=(const FreeEnergies& fs){
    f = fs.f;
    fnew = fs.fnew;
    fold = fs.fold;
    nparam = fs.nparam;

    return *this;
  }


  // @@@ TODO: fwrite, fread, forecasting

  void recenter_fnew(){
    const RealD fmin = *std::min_element(fnew.begin(), fnew.end());
    const RealD fmax = *std::max_element(fnew.begin(), fnew.end());
    for(RealD& elem : fnew) elem = elem - 0.5 * (fmin + fmax) ;
  }

  void update(){
    for (int j=0; j<nparam; j++) {
      f[j]    = fnew[j] ;
      fold[j] = fnew[j] ;
    }
  }

  RealD operator[](const int i) const { return f[i]; }
  RealD o(const int i) const { return fold[i]; }
  RealD n(const int i) const { return fnew[i]; }

  RealD& operator[](const int i) { return f[i]; }
  RealD& o(const int i) { return fold[i]; }
  RealD& n(const int i) { return fnew[i]; }

  RealD diff_norm_sq() const {
    double res = 0.0;
    for(int i=0; i<nparam; i++) res += (fold[i]-fnew[i])*(fold[i]-fnew[i]);
    return res;
  }

  RealD norm_sq() const {
    double res = 0.0;
    for(int i=0; i<nparam; i++) res += fold[i]*fold[i];
    return res;
  }
};



// lse_max from George
RealD lse_max(const std::vector<RealD>& a) {
  const RealD max = *std::max_element(a.begin(), a.end());
  RealD sum = 0.0;
  for(const RealD& elem : a) {
    // if( elem!=-std::numeric_limits<double>::infinity() ) sum += std::exp( elem - max );
    sum += std::exp( elem - max );
  }
  return max + std::log(sum);
}


// // lse_max from George
// RealD lse_max(const std::vector<RealD>& a) {
//   const RealD max = *std::max_element(a.begin(), a.end());
//   RealD sum = 0.0;
//   for(const RealD& elem : a) {
//     if( elem!=std::numeric_limits<double>::infinity() ) sum += std::exp( elem - max );
//   }
//   return max + std::log(sum);
// }


// RealD lse_max(const std::vector<std::vector<RealD>>& a) {
//   std::vector<RealD> tmp;
//   for(auto& v : a) tmp.push_back( lse_max(v) );
//   return lse_max(tmp);
// }



class LogRs {
public:
  std::vector<RealD> logRs;

  // reweighting betas
  Real beta_i;
  Real beta_j;
  Real beta_k;

  LogRs(){}

  void set(const Real beta_i_,
           const Real beta_j_,
           const Real beta_k_,
           const std::vector<RealD>& energies_k){
    beta_i = beta_i_;
    beta_j = beta_j_;
    beta_k = beta_k_;

    for(const RealD& H : energies_k){
      logRs.push_back( -(beta_j-beta_i)*H );
    }
  }

  RealD operator[](const int i) const { return logRs[i]; }
  int size() const { return logRs.size(); } // Nk
};



class LogGs {
public:
  std::vector<RealD> logGs;
  std::vector<RealD> logGPs;

  Real beta_i;
  Real beta_k;
  int size; // sample size of ensemble k
  int nparam;

  LogGs(){};

  LogGs(const std::vector<LogRs>& vj_logRs, const FreeEnergies& fs, const std::vector<Real>& logNs){
    set(vj_logRs, fs, logNs);
  }

  void set( const std::vector<LogRs>& vj_logRs, const FreeEnergies& fs, const std::vector<Real>& logNs ){
    beta_i = vj_logRs[0].beta_i;
    beta_k = vj_logRs[0].beta_k;

    size = vj_logRs[0].size();
    nparam = vj_logRs.size();

    for(const LogRs& logRs : vj_logRs) {
      assert( std::abs(logRs.beta_i-beta_i)<1.0e-14 ); // same shift
      assert( std::abs(logRs.beta_k-beta_k)<1.0e-14 ); // same ensemble
      assert( size==logRs.size() ); // same size
    }
    assert(fs.nparam==nparam); // same nparam
    assert(logNs.size()==nparam); // same nparam

    for(int ik=0; ik<size; ik++) {
      std::vector<RealD> b(fs.nparam);
      for(int j=0; j<nparam; j++) {
        b[j] = logNs[j] + fs[j] + vj_logRs[j][ik];
      }
      logGs.push_back( -lse_max(b) );
    }
  }

  void add_logGPs( const std::vector<RealD>& Os ){
    assert( size>0 );
    logGPs.resize( size );

    for(int ik=0; ik<size; ik++) {
      // assert( Os[ik]>0.0);
      if( Os[ik]>0.0 ){
        logGPs[ik] = std::log( Os[ik] ) + logGs[ik];
      }
      else if( std::abs(Os[ik])<1.0e-14 ){
        // assert( std::numeric_limits<double>::is_iec559 );
        logGPs[ik] = -std::numeric_limits<double>::max();
      }
      // else assert(false);
    }
  }


  // void add_OGs( const std::vector<RealD>& Os ){
  //   assert( size>0 );
  //   OGs.resize( size );

  //   for(int ik=0; ik<size; ik++) {
  //     OGs[ik] = Os[ik] * std::exp(logGs[ik]);
  //     else assert(false);
  //   }
  // }

  RealD operator[](const int i) const { return logGs[i]; }

  RealD get_B() const {
    return lse_max( logGs );
  }
  RealD get_BP() const {
    return lse_max( logGPs );
  }
};







void run_multihistogram( const std::vector<Real>& betas,
                         const std::vector<std::vector<RealD>>& v_energies,
                         FreeEnergies& fs,
                         const std::string& filename){
  const int nparam = betas.size();
  // assert( fs.nparam==nparam );
  if( fs.nparam!=nparam ){
    fs = FreeEnergies(nparam);
  }

  std::vector<Real> logNs( nparam );
  for(int j=0; j<nparam; j++){
    logNs[j] = std::log(v_energies[j].size());
  }

  // -------------------------
  // pre calc
  // -------------------------
  std::vector<std::vector<std::vector<LogRs>>> logRs_ikj;
  // just resizing
  logRs_ikj.resize( nparam );
  for(auto& logRs_kj : logRs_ikj){
    logRs_kj.resize( nparam );
    for(auto& logRs_j : logRs_kj){
      logRs_j.resize( nparam );
    }
  }

#ifdef _OPENMP
#pragma omp parallel for collapse(3)
#endif
  for(int i=0; i<nparam; i++) {
    for(int k=0; k<nparam; k++){
      for(int j=0; j<nparam; j++){
        logRs_ikj[i][k][j].set( betas[i], betas[j], betas[k], v_energies[k] );
      }}}

  std::cout << GridLogMessage << "logR set" << std::endl;

  // -------------------------
  // main loop
  // -------------------------

  do{
    fs.update();

    std::vector<std::vector<RealD>> B_ik(nparam, std::vector<RealD>(nparam));
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
    for(int i=0; i<nparam; i++) {
      for(int k=0; k<nparam; k++) {
        const LogGs logGs( logRs_ikj[i][k], fs, logNs );
        B_ik[i][k] = logGs.get_B();
      }
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int i=0; i<nparam; i++) {
      fs.n(i) = - lse_max( B_ik[i] );
    }

    fs.recenter_fnew();
    // @@@ TODO:forecasting

    std::cout << GridLogMessage << "norm = " << fs.diff_norm_sq() << std::endl;
    BinaryWriter WR( filename );
    write(WR, "f", fs.f );

  } while( fs.diff_norm_sq() > 1.0e-17 * fs.norm_sq() );

}





class ReweightedObs : public Serializable {
  const Real beta_min;
  const Real beta_max;

public:
  const int nmeas;
  const int nparam;

  std::vector<std::vector<RealD>> v_energies;
  std::vector<std::vector<RealD>> v_Os;
  std::vector<std::vector<RealD>> v_Osqs;
  std::vector<Real> logNs;

  GRID_SERIALIZABLE_CLASS_MEMBERS(
                                  ReweightedObs,
                                  //
                                  std::vector<Real>, betas,
                                  std::vector<RealD>, fs,
                                  std::vector<RealD>, fPs,
                                  std::vector<RealD>, fPs_sq
                                  );

  ReweightedObs( const Real beta_min_,
                 const Real beta_max_,
                 const int nmeas_,
                 const std::vector<EnsembleOfObs> ensembles
                 )
    : beta_min(beta_min_)
    , beta_max(beta_max_)
    , nmeas(nmeas_)
    , nparam(ensembles.size())
    , v_energies(nparam)
    , v_Os(nparam)
    , v_Osqs(nparam)
    , logNs(nparam)
    , betas( nmeas+1 )
    , fs( nmeas+1 )
    , fPs( nmeas+1 )
    , fPs_sq( nmeas+1 )
  {
    for(int k=0; k<nparam; k++){
      v_energies[k] = ensembles[k].energies;
      v_Os[k] = ensembles[k].Os;
      v_Osqs[k] = ensembles[k].Osqs;
      logNs[k] = std::log(v_energies[k].size());
    }
    for(int i=0; i<=nmeas; i++) betas[i] = beta_min + i*(beta_max-beta_min)/nmeas;
  }

  RealD est( const int i ) const {
    assert( i<=nmeas );
    return std::exp( fs[i] - fPs[i] );
  }

  RealD est_sq( const int i ) const {
    assert( i<=nmeas );
    return std::exp( fs[i] - fPs_sq[i] );
  }
};





void run_reweighting( const std::vector<Real>& betas,
                      const std::vector<EnsembleOfObs>& ensembles,
                      const FreeEnergies& fs,
                      ReweightedObs& meas,
                      const bool calculate_square=true,
                      const int nparallel=1){
  const int nparam=meas.nparam;
  const int nmeas=meas.nmeas;

  std::vector<std::vector<std::vector<LogRs>>> logRs_ikj;
  logRs_ikj.resize( nmeas+1 );
  for(auto& logRs_kj : logRs_ikj){
    logRs_kj.resize( nparam );
    for(auto& logRs_j : logRs_kj){
      logRs_j.resize( nparam );
    }
  }

#ifdef _OPENMP
#pragma omp parallel for collapse(3) num_threads(nparallel)
#endif
  for(int i=0; i<=nmeas; i++) {
    for(int k=0; k<nparam; k++){
      for(int j=0; j<nparam; j++){
        logRs_ikj[i][k][j].set( meas.betas[i], betas[j], betas[k], meas.v_energies[k] );
      }}}

  std::cout << GridLogMessage << "logR set" << std::endl;


  // =================================


  std::vector<std::vector<LogGs>> logGs_ik( nmeas+1, std::vector<LogGs>(nparam) );
#ifdef _OPENMP
#pragma omp parallel for collapse(2) num_threads(nparallel)
#endif
  for(int i=0; i<=nmeas; i++) {
    for(int k=0; k<nparam; k++) {
      logGs_ik[i][k].set( logRs_ikj[i][k], fs, meas.logNs );
    }
  }

  std::cout << GridLogMessage << "logGs set" << std::endl;

  std::vector<std::vector<RealD>> Bs_ik(nmeas+1, std::vector<RealD>(nparam) );
#ifdef _OPENMP
#pragma omp parallel for collapse(2) num_threads(nparallel)
#endif
  for(int i=0; i<=nmeas; i++){
    for(int k=0; k<nparam; k++) Bs_ik[i][k] = logGs_ik[i][k].get_B();
  }

  std::cout << GridLogMessage << "Bs set" << std::endl;


#ifdef _OPENMP
#pragma omp parallel for num_threads(nparallel)
#endif
  for(int i=0; i<=nmeas; i++){
    meas.fs[i] = -lse_max( Bs_ik[i] );
  }


  std::vector<std::vector<RealD>> BPs_ik( nmeas+1, std::vector<RealD>(nparam) );
#ifdef _OPENMP
#pragma omp parallel for collapse(2) num_threads(nparallel)
#endif
  for(int i=0; i<=nmeas; i++) {
    for(int k=0; k<nparam; k++){
      logGs_ik[i][k].add_logGPs( meas.v_Os[k] );
      BPs_ik[i][k] = logGs_ik[i][k].get_BP();
    }
  }

#ifdef _OPENMP
#pragma omp parallel for num_threads(nparallel)
#endif
  for(int i=0; i<=nmeas; i++) {
    meas.fPs[i] = -lse_max( BPs_ik[i] );
  }

  if(calculate_square){
    std::vector<std::vector<RealD>> BPs_sq_ik( nmeas+1, std::vector<RealD>(nparam) );
#ifdef _OPENMP
#pragma omp parallel for collapse(2) num_threads(nparallel)
#endif
    for(int i=0; i<=nmeas; i++) {
      for(int k=0; k<nparam; k++){
        logGs_ik[i][k].add_logGPs( meas.v_Osqs[k] );
        BPs_sq_ik[i][k] = logGs_ik[i][k].get_BP();
      }
    }

#ifdef _OPENMP
#pragma omp parallel for num_threads(nparallel)
#endif
    for(int i=0; i<=nmeas; i++) {
      meas.fPs_sq[i] = -lse_max( BPs_sq_ik[i] );
    }
  }
}







