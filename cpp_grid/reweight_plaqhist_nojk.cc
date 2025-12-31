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
  // const std::string basedir2(argv[4]);
  const std::string mass(argv[6]);
  const int conf_max0=atoi(argv[7]);
  const int nbin_global=atoi(argv[8]);
  int nbeta=atoi(argv[9]);
  int runtype=atoi(argv[10]);
  const int nbins_histo_x=atoi(argv[11]);
  double min_x = atof(argv[12]);
  double max_x = atof(argv[13]);

  std::cout << "nbeta = " << nbeta << std::endl;
  std::cout << "basedir = " << basedir << std::endl;
  std::cout << "nbins_histo_x = " << nbins_histo_x << std::endl;
  std::cout << "min_x = " << min_x << std::endl;
  std::cout << "max_x = " << max_x << std::endl;

  // int type=1;
  // if(basedir2c==basedir2h) {
  //   if(mass=="0p1000") type=2;
  //   else if(mass=="0p0500" || mass=="0p0100" || mass=="0p2000") type=4;
  // }

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

      betas_double.push_back(std::stod(str));
    }
  }

  using Impl = PeriodicGimplR;


  // ################################################################ //
  std::cout << GridLogMessage << "run histogram part" << std::endl;

  char id[200];
  std::sprintf(id, "m%s", mass.data());

  // no openmp loop so far
  {
    // read fs
    FreeEnergies fs(nbeta);
    {
      std::string id_str = id;
      id_str = id_str+"_nojk";

      std::unique_ptr<Hdf5Reader> WR;
      WR = std::make_unique<Hdf5Reader>( basedir+id_str+"fs.bin" );
      read(*WR, "f", fs.f );
    }

    std::vector<std::vector<EnsembleOfObs>> array_of_hist_ens( nbins_histo_x,
                                                               std::vector<EnsembleOfObs>(nbeta)
                                                               );

#pragma omp parallel for
    for(int j=0; j<nbeta; j++){
      // EnsembleOfObs ens = array_of_hist_ens;
      const std::string obs_id = get_configname(betas[j], mass);

      for(std::string basedir2 : std::vector<std::string>{basedir2c, basedir2h} ){
        for(int conf=conf_min0; conf<conf_max0; conf++){
          // std::cout << "debug. conf = " << conf << std::endl;

          std::string pathO = basedir2+obs_id+std::to_string(conf)+".bin";
          std::string pathO3 = basedir2+obs_id+"plaqhist"+std::to_string(conf)+".bin";
          if( !std::filesystem::exists( pathO3 ) || !std::filesystem::exists( pathO ) ) continue; // assert(false);

          VectorObs hist;
          Obs obs;
#pragma omp critical
          {
            std::unique_ptr<Hdf5Reader> WR;
            WR = std::make_unique<Hdf5Reader>( pathO3 );
            read(*WR, "obs", hist );
          }
#pragma omp critical
          {
            std::unique_ptr<Hdf5Reader> WR;
            WR = std::make_unique<Hdf5Reader>( pathO );
            read(*WR, "obs", obs );
          }

          assert( obs.beta == betas[j] );

          // push_back
          for(int ibx=0; ibx<nbins_histo_x; ibx++){ // histogram binning for
            array_of_hist_ens[ibx][j].Os.push_back( hist.vO[ibx] );
            array_of_hist_ens[ibx][j].energies.push_back( obs.energy );
          }
        } // conf
      } // cold, hot
      for(int ibx=0; ibx<nbins_histo_x; ibx++){ // histogram binning for
        assert( array_of_hist_ens[ibx][j].size()!=0 );
      }
      // ensembles.push_back(ens);
    } // for j
    // assert( ensembles.size() == nbeta );


    std::cout << "debug. hist" << std::endl;


#pragma omp parallel for num_threads(omp_get_max_threads()/4)
    for(int ibx=0; ibx<nbins_histo_x; ibx++){ // histogram binning for
        std::string id_str = id;
        id_str = id_str+"plaqhist_ibx"+std::to_string(ibx)+"_nojk";
        std::string path_res = basedir+id_str+"meas.bin";
        if( std::filesystem::exists( path_res ) && (runtype==0) ) continue; // assert(false);

        ReweightedObs meas( betas_double.front(), betas_double.back(), nbeta*100, array_of_hist_ens[ibx] );
        const int nmeas=meas.nmeas;
        run_reweighting( betas_double, fs, meas, false, 4 );

#pragma omp critical
        {
          std::unique_ptr<Hdf5Writer> WR;
          WR = std::make_unique<Hdf5Writer>( path_res );
          write(*WR, "beta", meas.betas );
          write(*WR, "f", meas.fs );
          write(*WR, "fP", meas.fPs );
        }

    } // end for ib
  }


  Grid_finalize();
}
