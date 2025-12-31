// #include "reweight.h"
#include "reweight2.h"

#include <regex>

// #define NOSIGN

int main(int argc, char **argv) {
  Grid_init(&argc, &argv);
  int threads = GridThread::GetThreads();
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

  // -------------------------------------
  // for reading data
  const int conf_min0=atoi(argv[1]);
  const std::string base_dir(argv[2]);
  const std::string basedir(argv[3]);
  const std::string basedir2c(argv[4]);
  const std::string basedir2h(argv[5]);
  const std::string mass(argv[6]);
  const int conf_max0=atoi(argv[7]);
  // const int binsize=atoi(argv[9]);
  const int nbin_global=atoi(argv[8]);
  int nbeta=atoi(argv[9]);
  int runtype=atoi(argv[10]);
  // obsinfo = argv[12]; // dummy

  std::cout << "nbeta = " << nbeta << std::endl;
  std::cout << "basedir = " << basedir << std::endl;

  // std::string obsinfoR="flow_polyakov_re";
  // std::string obsinfoI="flow_polyakov_im";

  int type=1;
  if(basedir2c==basedir2h) {
    if(mass=="0p1000" || mass=="0p4000") type=2;
    else if(mass=="0p0500" || mass=="0p0100" || mass=="0p2000" || mass=="0p3000") type=4;
    else if(mass=="0") type=6;
  }
  else if(mass=="0") type = 0;

  // -----------------------

  std::vector<std::string> betas;
  std::vector<double> betas_double;
  {
    for(int i=11; i<11+nbeta; i++) {
      std::string str(argv[i]);
      betas.push_back(str);

      std::regex to_replace("p");
      str = std::regex_replace(str, to_replace, ".");
      std::cout << "str replaced = " << str << std::endl;

      // std::replace( str, "p", "." );
      betas_double.push_back(std::stod(str));
    }
  }

  using Impl = PeriodicGimplR;

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(std::string beta : betas){
    std::cout << GridLogMessage << "beta = " << beta << std::endl;
    ComplexEnsembleOfObs ens;
    ens.beta = beta;
    ens.base_dir=base_dir;
    if(type!=1) ens.id=get_configname(beta, mass, type+1);
    else ens.id=get_configname(beta, mass, type);

    std::string path_ens = basedir+get_configname(beta, mass, type)+"_Pcplx.h5";
    if( std::filesystem::exists( path_ens ) && (runtype==0) ) {
      continue;
    }

    // ---------------------------
    // cold start
    // ---------------------------

    ens.id=get_configname(beta, mass, type);

    for(int conf=conf_min0; conf<conf_max0; conf++){
      std::string pathO;
      if(type>1) pathO = basedir2c+ens.id+std::to_string(conf)+".h5";
      else pathO = basedir2c+ens.id+std::to_string(conf)+".bin";
      if( !std::filesystem::exists( pathO ) ) continue; // assert(false);

#pragma omp critical
      {
        std::unique_ptr<Hdf5Reader> WR;
        WR = std::make_unique<Hdf5Reader>( pathO );
        Obs obs;
        read(*WR, "obs", obs );
        std::vector<ComplexD> all_polyakov;
        read(*WR, "polyakov", all_polyakov );
        ens.Os.push_back( all_polyakov.back() );
        ens.energies.push_back(obs.energy);
      }
    } // end for conf

    // ---------------------------
    // hot start
    // ---------------------------

    if(type!=1) {
      ens.id=get_configname(beta, mass, type+1);
    }

    for(int conf=conf_min0; conf<conf_max0; conf++){
      // std::string pathO = basedir2h+ens.id+std::to_string(conf)+".bin";
      std::string pathO;
      if(type>1) pathO = basedir2c+ens.id+std::to_string(conf)+".h5";
      else pathO = basedir2c+ens.id+std::to_string(conf)+".bin";
      if( !std::filesystem::exists( pathO ) ) continue; // assert(false);

#pragma omp critical
      {
        std::unique_ptr<Hdf5Reader> WR;
        WR = std::make_unique<Hdf5Reader>( pathO );
        Obs obs;
        read(*WR, "obs", obs );
        std::vector<ComplexD> all_polyakov;
        read(*WR, "polyakov", all_polyakov );
        ens.Os.push_back( all_polyakov.back() );
        ens.energies.push_back(obs.energy);
      }
    } // end for conf

    ens.nbin=nbin_global;
    ens.binsize = ens.size()/nbin_global;
    ens.nconf = ens.nbin*ens.binsize;

    ens.compute_sq();

    { // trimming
      const int start = ens.nconf;
      ens.energies.erase( ens.energies.begin()+start, ens.energies.end() );
      ens.Os.erase( ens.Os.begin()+start, ens.Os.end() );
      ens.Osqs.erase( ens.Osqs.begin()+start, ens.Osqs.end() );
    }

    std::cout << GridLogMessage << "beta = " << beta << " size = " << ens.size() << " nbin = " << ens.nbin << " nconf = " << ens.nconf << std::endl;

#pragma omp critical
    {
      if( std::filesystem::exists( path_ens ) ) std::filesystem::remove( path_ens );
      std::unique_ptr<Hdf5Writer> WR;
      WR = std::make_unique<Hdf5Writer>( path_ens );
      // std::cout << "debug. write" << std::endl;
      // std::cout << "debug. --e" << std::endl;
      // write(*WR, "energies", ens.energies );
      // std::cout << "debug. --O" << std::endl;
      // write(*WR, "Os", ens.Os );
      // std::cout << "debug. --Osq" << std::endl;
      // write(*WR, "Osqs", ens.Osqs );
      write(*WR, "ens", ens );
      // std::cout << "debug. write, done" << std::endl;
    }
  } // for end beta



  // ################################################################ //
  std::cout << GridLogMessage << "run3 part" << std::endl;

  {
    std::vector<ComplexEnsembleOfObs> ensembles0( nbeta );

// #ifdef _OPENMP
// #pragma omp parallel for
// #endif
    for(int j=0; j<nbeta; j++){
      const std::string beta = betas[j];
      std::string id_str = get_configname(beta, mass, type);
      // std::string path = basedir+id_str+obsinfo+".bin";
      std::string path_ens = basedir+id_str+"_Pcplx.h5";

#pragma omp critical
      {
        std::cout << GridLogMessage << "Pcplx" << std::endl;
        std::unique_ptr<Hdf5Reader> WR;
        WR = std::make_unique<Hdf5Reader>( path_ens );
        read(*WR, "ens", ensembles0[j] );
        // auto& ens = ensembles0[j];
        // XmlReader WR( path_ens );
        // read(*WR, "energies", ens.energies );
        // read(*WR, "Os", ens.Os );
        // read(*WR, "Osqs", ens.Osqs );
        // real part only
        // std::unique_ptr<Hdf5Writer> WR;
        // WR = std::make_unique<Hdf5Writer>( path_ens );
      }
    }

    char id[200];
    std::sprintf(id, "m%s", mass.data());

    // -------------------

    FreeEnergies fs(nbeta);

    std::string id_str = id;
    id_str = id_str+"_nojk";

    // @@@ -> xml?
    // BinaryReader WR(basedir+id_str+"fs.bin");
    // read(WR, "f", fs.f );
#pragma omp critical
    {
      std::cout << GridLogMessage << "f" << std::endl;
      std::unique_ptr<Hdf5Reader> WR;
      WR = std::make_unique<Hdf5Reader>( basedir+id_str+"fs.bin" );
      read(*WR, "f", fs.f );
      // read(WR, "ens", ensembles0[j] );
    }

    // =================================

    std::vector<EnsembleOfObs> ensembles_re( nbeta );
    std::vector<EnsembleOfObs> ensembles_im( nbeta );
    std::vector<EnsembleOfObs> ensembles_abs( nbeta );

    for(int ibeta=0; ibeta<nbeta; ibeta++){
      const auto& ens0 = ensembles0[ibeta];
      auto& ens_re = ensembles_re[ibeta];
      auto& ens_im = ensembles_im[ibeta];
      auto& ens_abs = ensembles_abs[ibeta];

      ens_re.energies = ens0.energies;
      ens_im.energies = ens0.energies;
      ens_abs.energies = ens0.energies;

      for(const ComplexD elem : ens0.Os){
        ens_re.Os.push_back( real(elem) );
        ens_im.Os.push_back( imag(elem) );
        ens_abs.Os.push_back( abs(elem) );
      }

      ens_re.compute_sq();
      ens_im.compute_sq();
      ens_abs.compute_sq();
    }

    ComplexReweightedObs meas( betas_double.front(), betas_double.back(), nbeta*100, ensembles0 );
    ReweightedObs meas_re( betas_double.front(), betas_double.back(), nbeta*100, ensembles_re );
    ReweightedObs meas_im( betas_double.front(), betas_double.back(), nbeta*100, ensembles_im );
    ReweightedObs meas_abs( betas_double.front(), betas_double.back(), nbeta*100, ensembles_abs );
    const int nmeas=meas.nmeas;

    run_reweighting( betas_double, fs, meas );
    run_reweighting( betas_double, fs, meas_re );
    run_reweighting( betas_double, fs, meas_im );
    run_reweighting( betas_double, fs, meas_abs );

//     std::cout << "fs = ";
//     for(auto elem : meas.fPs) std::cout << elem << " ";
//     std::cout << std::endl;

//     std::cout << "fPs = ";
//     for(auto elem : meas.fPs) std::cout << elem << " ";
//     std::cout << std::endl;

// // #ifndef NOSIGN
// //     std::cout << "sign_P = ";
// //     for(auto elem : meas.sign_Ps) std::cout << elem << " ";
// //     std::cout << std::endl;
// // #endif

//     std::cout << "fP_sq = ";
//     for(auto elem : meas.fPs_sq) std::cout << elem << " ";
//     std::cout << std::endl;

    std::string path_res = basedir+id_str+"meas_Pcplx.h5";
#pragma omp critical
    {
      if( std::filesystem::exists( path_res ) ) std::filesystem::remove( path_res );
      std::cout << "path_res = " << path_res << std::endl;
      std::unique_ptr<Hdf5Writer> WR;
      WR = std::make_unique<Hdf5Writer>( path_res );
      write(*WR, "beta", meas.betas );
      write(*WR, "f", meas.fs );
      write(*WR, "fP", meas.fPs );
      write(*WR, "fPs_sq", meas.fPs_sq );
      //
      write(*WR, "fP_re", meas_re.fPs );
      write(*WR, "fPs_re_sq", meas_re.fPs_sq );
      write(*WR, "fP_im", meas_im.fPs );
      write(*WR, "fPs_im_sq", meas_im.fPs_sq );
      write(*WR, "fP_abs", meas_abs.fPs );
      write(*WR, "fPs_abs_sq", meas_abs.fPs_sq );
    }
  }

  Grid_finalize();
}
