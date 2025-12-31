#include "reweight.h"


int main(int argc, char **argv) {
  Grid_init(&argc, &argv);
  int threads = GridThread::GetThreads();
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

  GridCartesian         *UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(), GridDefaultSimd(Nd,vComplex::Nsimd()), GridDefaultMpi());

  // -------------------------------------
  // for reading data
  // const Real mass=0.4;
  // std::vector<Real> betas({11.01, 11.02, 11.03, 11.04, 11.05, 11.06, 11.07, 11.08, 11.09});
  // std::vector<Real> betas({10.98, 10.99, 11.00});
  // @@@@ should be arrays ---------------
  // const int conf_min=1000;
  // const int conf_max=3000;
  const int conf_min=atoi(argv[1]);
  const int conf_max=atoi(argv[2]);
  // -----------------------

  using Impl = PeriodicGimplR;

  for(Real beta : betas){
    std::cout << "beta = " << beta << std::endl;
    WilsonGaugeActionR Waction(beta);

    EnsembleOfObs ens;
    ens.beta = beta;
    ens.conf_min = conf_min;
    ens.conf_max = conf_max;
    ens.interval = interval;

    ens.base_dir=base_dir;
    ens.id=get_configname(beta, mass);

    for(int conf=conf_min; conf<conf_max; conf+=interval){

      char f[200];
      std::string prefix = get_prefix(beta, mass);
      std::sprintf(f,
                   "%s.%d",
                   prefix.data(), conf);

      const std::string f_str = f;
      const std::string path = ens.base_dir+ens.id+"/"+f_str;
      std::cout << "path = " << path << std::endl;
      FieldMetaData header;
      LatticeGaugeField U(UGrid);
      NerscIO::readConfiguration(U, header, path);
      std::cout << "read. " << std::endl;

      const ComplexD p    = WilsonLoops<Impl>::avgPolyakovLoop(U);
      // ens.Os.push_back(abs(p));
      ens.Os.push_back(abs(p));
      ens.Osqs.push_back(abs(p)*abs(p));
      const ComplexD S    = Waction.S(U);
      ens.energies.push_back(real(S)/beta);
    }

    std::string path = basedir+ens.id+obsinfo+".bin";
    std::cout << "path = " << path << std::endl;
    BinaryWriter WR( path );
    write(WR, "ens", ens );
  }

  Grid_finalize();
}
