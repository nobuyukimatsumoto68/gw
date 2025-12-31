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
  const int nbin_global=atoi(argv[8]);
  int nbeta=atoi(argv[9]);
  int runtype=atoi(argv[10]);
  const double beta_meas_min=atof(argv[11]);
  const double beta_meas_max=atof(argv[12]);
  const int nbeta_meas=atoi(argv[13]);

  std::cout << "nbeta = " << nbeta << std::endl;
  std::cout << "basedir = " << basedir << std::endl;

  int type=1;
  if(basedir2c==basedir2h) {
    if(mass=="0p1000" || mass=="0p4000") type=2;
    else if(mass=="0p0500" || mass=="0p0100" || mass=="0p2000" || mass=="0p3000") type=4;
    else if(mass=="0") type=6;
  }
  else if(mass=="0") type = 0;

  if(runtype>=2) type+=1;
  runtype -= 2;

  // -----------------------

  std::vector<std::string> betas;
  std::vector<double> betas_double;
  {
    for(int i=14; i<14+nbeta; i++) {
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

    EnsembleOfObs ens;
    ens.beta = beta;
    ens.base_dir=base_dir;
    if(type!=1) ens.id=get_configname(beta, mass, type+1);
    else ens.id=get_configname(beta, mass, type);

    std::string path_ens = basedir+get_configname(beta, mass, type)+"_chcond.h5";
    if( std::filesystem::exists( path_ens ) && (runtype==0) ) continue;

    // ---------------------------
    // cold start
    // ---------------------------

    ens.id=get_configname(beta, mass, type);

    for(int conf=conf_min0; conf<conf_max0; conf++){
      std::string pathO;
      if(type>1) pathO = basedir2c+ens.id+std::to_string(conf)+".h5";
      else pathO = basedir2c+ens.id+std::to_string(conf)+".bin";
      std::string pathO3 = basedir2c+ens.id+"chcond"+std::to_string(conf)+".xml";
      if( !std::filesystem::exists( pathO3 ) ) continue; // assert(false);
      if( !std::filesystem::exists( pathO ) ) continue; // assert(false);
      // if( !std::filesystem::exists( pathO3 ) ) {
      //   std::cout << "pathO3 = " << pathO3 << std::endl;
      //   assert(false);
      // }


#pragma omp critical
      {
        std::unique_ptr<Hdf5Reader> WR;
        WR = std::make_unique<Hdf5Reader>( pathO );
        Obs obs;
        read(*WR, "obs", obs );
        ens.energies.push_back(obs.energy);
      }
      {
        XmlReader WR( pathO3 );

        if(type>1){
          Complex chcond;
          read(WR, "chcond", chcond );
          ens.Os.push_back( real(chcond) );
        }
        else{
          Real chcond;
          read(WR, "chcond", chcond );
          ens.Os.push_back( chcond );
        }
      }

    } // end for conf

    // ---------------------------
    // hot start
    // ---------------------------

    if(type!=1) {
      ens.id=get_configname(beta, mass, type+1);
    }

    for(int conf=conf_min0; conf<conf_max0; conf++){
      std::string pathO;
      if(type>1) pathO = basedir2h+ens.id+std::to_string(conf)+".h5";
      else pathO = basedir2c+ens.id+std::to_string(conf)+".bin";
      std::string pathO3 = basedir2h+ens.id+"chcond"+std::to_string(conf)+".xml";
      if( !std::filesystem::exists( pathO ) ) continue; // assert(false);
      if( !std::filesystem::exists( pathO3 ) ) continue; // assert(false);
      // if( !std::filesystem::exists( pathO3 ) ) {
      //   std::cout << "pathO3 = " << pathO3 << std::endl;
      //   assert(false);
      // }

#pragma omp critical
      {
        std::unique_ptr<Hdf5Reader> WR;
        WR = std::make_unique<Hdf5Reader>( pathO );
        Obs obs;
        read(*WR, "obs", obs );
        ens.energies.push_back(obs.energy);
      }
      {
        XmlReader WR( pathO3 );
        // Real chcond;
        // read(WR, "chcond", chcond );
        // ens.Os.push_back( chcond );
        if(type>1){
          Complex chcond;
          read(WR, "chcond", chcond );
          ens.Os.push_back( real(chcond) );
        }
        else{
          Real chcond;
          read(WR, "chcond", chcond );
          ens.Os.push_back( chcond );
        }

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
      write(*WR, "ens", ens );
    }
  } // for end beta



  // ################################################################ //
  std::cout << GridLogMessage << "run3 part" << std::endl;

  {
    std::vector<EnsembleOfObs> ensembles0( nbeta );

    for(int j=0; j<nbeta; j++){
      const std::string beta = betas[j];
      std::string id_str = get_configname(beta, mass, type);
      std::string path_ens = basedir+id_str+"_chcond.h5";

#pragma omp critical
      {
        std::unique_ptr<Hdf5Reader> WR;
        WR = std::make_unique<Hdf5Reader>( path_ens );
        read(*WR, "ens", ensembles0[j] );
      }
    }

    char id[200];
    std::sprintf(id, "m%s", mass.data());

    // -------------------

    FreeEnergies fs(nbeta);

    std::string id_str = id;
    id_str = id_str+"_nojk";

#pragma omp critical
    {
      std::cout << GridLogMessage << "f" << std::endl;
      std::unique_ptr<Hdf5Reader> WR;
      WR = std::make_unique<Hdf5Reader>( basedir+id_str+"fs.bin" );
      read(*WR, "f", fs.f );
    }

    // =================================

    ReweightedObs meas( beta_meas_min, beta_meas_max, nbeta_meas, ensembles0 );
    const int nmeas=meas.nmeas;

    run_reweighting( betas_double, fs, meas );

    std::string path_res = basedir+id_str+"meas_chcond.h5";
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
    }
  }

  Grid_finalize();
}
