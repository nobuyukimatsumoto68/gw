// #define TEST
// #undef _OPENMP

#include "reweight.h"

int main(int argc, char **argv) {
  Grid_init(&argc, &argv);


  int threads = GridThread::GetThreads();
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

  // -------------------------------------

  const std::string base_dir(argv[1]);
  const std::string basedir(argv[2]);
  const std::string basedir2(argv[3]);
  const std::string mass(argv[4]);
  std::cout << "mass = " << mass << std::endl;
  const int nbeta=atoi(argv[5]);
  std::cout << "nbeta = " << nbeta << std::endl;

  std::vector<std::string> betas;
  std::vector<double> betas_double;
  {
    for(int i=6; i<6+nbeta; i++) {
      std::string str(argv[i]);
      betas.push_back(str);
      betas_double.push_back(std::stod(str));
    }
  }


  std::cout << "betas = " << std::endl;
  for(int i=0; i<nbeta; i++){
    std::cout << betas[i] << " ";
  }
  std::cout << std::endl;


  std::vector<EnsembleOfObs> ensembles0( nbeta );

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int j=0; j<nbeta; j++){
    const std::string beta = betas[j];
    std::string id_str = get_configname(beta, mass);
    std::string path = basedir+id_str+obsinfo+".bin";
    std::cout << "path = " << path << std::endl;

    BinaryReader WR( path );
    read(WR, "ens", ensembles0[j] );
  }

  char id[200];
  std::sprintf(id, "m%s", mass.data());

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(int jdrop=0; jdrop<nbeta; jdrop++){
    const int nbin=ensembles0[jdrop].nbin;
    for(int ibin=0; ibin<nbin; ibin++){
      std::vector<EnsembleOfObs> ensembles(nbeta);
      for(int i=0; i<nbeta; i++) ensembles[i] = ensembles0[i];
      ensembles[jdrop].jackknife_drop( nbin, ibin );

      // -------------------

      FreeEnergies fs(nbeta);

      std::string id_str = id;
      id_str = id_str+"_jk_"+std::to_string(jdrop)+"_"+std::to_string(nbin)+"_"+std::to_string(ibin);

      BinaryReader WR(basedir+id_str+"fs.bin");
      read(WR, "f", fs.f );

      // =================================

      ReweightedObs meas( betas_double.front(), betas_double.back(), nbeta*100, ensembles );
      const int nmeas=meas.nmeas;
      run_reweighting( betas_double, ensembles, fs, meas );

      {
        BinaryWriter WR(basedir+id_str+"meas_beta.bin");
        write(WR, "beta", meas.betas );
      }
      {
        BinaryWriter WR(basedir+id_str+"meas_f.bin");
        write(WR, "f", meas.fs );
      }
      {
        BinaryWriter WR(basedir+id_str+"meas_fP_"+obsinfo+".bin");
        write(WR, "fP", meas.fPs );
      }
      {
        BinaryWriter WR(basedir+id_str+"meas_fP_sq_"+obsinfo+".bin");
        write(WR, "fP_sq", meas.fPs_sq );
      }
    }
  }

  Grid_finalize();
}
