#include "reweight.h"


// #define TEST
// #undef _OPENMP


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

  std::cout << GridLogMessage << "Ensembles Ready?" << std::endl;

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

      std::vector<std::vector<RealD>> v_energies(nbeta);
      for(int k=0; k<nbeta; k++) v_energies[k] = ensembles[k].energies;
      run_multihistogram( betas_double, v_energies, fs );

      char id[200];
      std::sprintf(id, "m%s", mass.data());
      std::string id_str = id;
      id_str = id_str+"_jk_"+std::to_string(jdrop)+"_"+std::to_string(nbin)+"_"+std::to_string(ibin);

      BinaryWriter WR(basedir+id_str+"fs.bin");
      write(WR, "f", fs.f );
    }
  }



  Grid_finalize();
}
