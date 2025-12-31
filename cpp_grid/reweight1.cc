#include "reweight.h"


int main(int argc, char **argv) {
  Grid_init(&argc, &argv);
  int threads = GridThread::GetThreads();
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

  GridCartesian         *UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()), GridDefaultMpi());

  // -------------------------------------
  // for reading data
  const int conf_min=atoi(argv[1]);
  // const int conf_max=atoi(argv[2]);
  const std::string base_dir(argv[2]);
  const std::string basedir(argv[3]);
  const std::string basedir2(argv[4]);
  const std::string mass(argv[5]);
  const int interval=atoi(argv[6]);
  const int binsize=atoi(argv[7]);
  const int nbeta=atoi(argv[8]);

  // -----------------------

  std::cout << "nbeta = " << nbeta << std::endl;

  std::vector<std::string> betas;
  {
    for(int i=9; i<9+nbeta; i++) {
      std::string str(argv[i]);
      betas.push_back(str);
    }
  }

  for(auto elem : betas) std::cout << elem << std::endl;

  using Impl = PeriodicGimplR;

  for(std::string beta : betas){
    std::cout << "beta = " << beta << std::endl;
    WilsonGaugeActionR Waction(std::stod(beta));

    int conf_max=conf_min;
    for(int conf=conf_min; conf<conf_min+10000000*interval; conf+=interval){
      char f[200];
      std::string prefix = get_prefix(beta, mass);
      std::sprintf(f,
                   "%s.%d",
                   prefix.data(), conf);

      const std::string f_str = f;
      const std::string path = base_dir+get_configname(beta, mass)+"/"+f_str;
      std::cout << "path = " << path << std::endl;

      conf_max = conf-interval;
      if( !std::filesystem::exists( path ) ) {
        break;
      }
    }

    const int nbin = (conf_max-conf_min)/(interval*binsize);
    conf_max = conf_min + nbin*interval*binsize;

    EnsembleOfObs ens;
    ens.beta = beta;
    ens.conf_min = conf_min;
    ens.conf_max = conf_max;
    ens.interval = interval;
    ens.binsize = binsize;
    ens.nbin = nbin;

    ens.base_dir=base_dir;
    ens.id=get_configname(beta, mass);

    for(int conf=conf_min; conf<conf_max; conf+=interval){
      std::string pathO = basedir2+ens.id+obsinfo+std::to_string(conf)+".bin";
      std::cout << "pathO = " << pathO << std::endl;
      if( !std::filesystem::exists( pathO ) ) {
        assert(false);
      }

      Obs obs;
      BinaryReader WR( pathO );
      read(WR, "obs", obs );

      ens.Os.push_back(obs.O);
      ens.Osqs.push_back(obs.O*obs.O);
      ens.energies.push_back(obs.energy);
    } // end for conf

    std::string path = basedir+ens.id+obsinfo+".bin";
    std::cout << "path = " << path << std::endl;
    BinaryWriter WR( path );
    write(WR, "ens", ens );
  }

  Grid_finalize();
}
