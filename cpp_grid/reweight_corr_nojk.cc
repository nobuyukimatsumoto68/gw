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

  std::cout << "nbeta = " << nbeta << std::endl;
  std::cout << "basedir = " << basedir << std::endl;
  // std::cout << "nbins_histo_x = " << nbins_histo_x << std::endl;
  // std::cout << "min_x = " << min_x << std::endl;
  // std::cout << "max_x = " << max_x << std::endl;
  // std::cout << "nbins_histo_y = " << nbins_histo_y << std::endl;
  // std::cout << "min_y = " << min_y << std::endl;
  // std::cout << "max_y = " << max_y << std::endl;

  int type=1;
  if(basedir2c==basedir2h) {
    if(mass=="0p1000" || mass=="0p4000") type=2;
    else if(mass=="0p0500" || mass=="0p0100" || mass=="0p2000" || mass=="0p3000") type=4;
    else if(mass=="0") type=6;
  }
  else if(mass=="0") type = 0;

  if(runtype>=2) {
    type+=1;
    runtype -= 2;
  }


  // int type=1;
  // if(basedir2c==basedir2h) {
  //   if(mass=="0p1000") type=2;
  //   else if(mass=="0p0500" || mass=="0p0100" || mass=="0p2000") type=4;
  // }

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

      betas_double.push_back(std::stod(str));
    }
  }

  using Impl = PeriodicGimplR;

  GridCartesian         *UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
                                                                  GridDefaultSimd(Nd,vComplex::Nsimd()),
                                                                  GridDefaultMpi());



  // ################################################################ //
  std::cout << GridLogMessage << "run histogram part" << std::endl;

  char id[200];
  std::sprintf(id, "m%s", mass.data());

  std::cout << "debug. f" << std::endl;

  // no openmp loop so far
  {
    // read fs
    FreeEnergies fs(nbeta);
    {
      std::string id_str = id;
      id_str = id_str+"_nojk";

      std::unique_ptr<Hdf5Reader> WR;
      std::cout << "debug. fpath = " << basedir+id_str+"fs.bin" << std::endl;
      WR = std::make_unique<Hdf5Reader>( basedir+id_str+"fs.bin" );
      read(*WR, "f", fs.f );
    }

    int X = UGrid->GlobalDimensions()[0];
    std::vector<std::vector<EnsembleOfObs>> array_of_corr_ens( X,
                                                               std::vector<EnsembleOfObs>(nbeta)
                                                               );

    std::cout << "debug. read" << std::endl;

#pragma omp parallel for
    for(int j=0; j<nbeta; j++){
      // EnsembleOfObs ens = array_of_hist_ens;
      const std::string obs_id = get_configname(betas[j], mass, type);

      for(std::string basedir2 : std::vector<std::string>{basedir2c, basedir2h} ){
        for(int conf=conf_min0; conf<conf_max0; conf++){
          // std::cout << "debug. conf = " << conf << std::endl;

          std::string pathO; // = basedir2+obs_id+std::to_string(conf)+".bin";
          if(type>1) pathO = basedir2+get_configname(betas[j], mass, type)+std::to_string(conf)+".h5";
          else pathO = basedir2+get_configname(betas[j], mass, type)+std::to_string(conf)+".bin";
          std::string pathO3 = basedir2+obs_id+"P_corr"+std::to_string(conf)+".h5";
          if( !std::filesystem::exists( pathO3 ) || !std::filesystem::exists( pathO ) ) continue; // assert(false);

          // HistObs hist;
          std::vector<RealD> xcorr(X);
          Obs obs;
#pragma omp critical
          {
            std::unique_ptr<Hdf5Reader> WR;
            WR = std::make_unique<Hdf5Reader>( pathO3 );
            read(*WR, "xcorr", xcorr );
          }
#pragma omp critical
          {
            std::unique_ptr<Hdf5Reader> WR;
            WR = std::make_unique<Hdf5Reader>( pathO );
            read(*WR, "obs", obs );
          }

          // assert( obs.beta == betas[j] );

          // push_back
          for(int x=0; x<X; x++) {
            array_of_corr_ens[x][j].Os.push_back( xcorr[x] );
            array_of_corr_ens[x][j].energies.push_back( obs.energy );
          }
        } // conf
      } // cold, hot
    } // for j

    // std::cout << "debug. hist" << std::endl;

#pragma omp parallel for num_threads(omp_get_max_threads()/4)
    for(int x=0; x<X; x++){ // histogram binning for
      std::string id_str = id;
      id_str = id_str+"corr_x"+std::to_string(x)+"_nojk";
      std::string path_res = basedir+id_str+"meas.bin";
      if( std::filesystem::exists( path_res ) && (runtype==0) ) continue; // assert(false);

      ReweightedObs meas( betas_double.front(), betas_double.back(), nbeta*100, array_of_corr_ens[x] );
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
    } // end for x
  }


  Grid_finalize();
}
