#include "reweight.h"
#include <sstream>


int main(int argc, char **argv) {
  Grid_init(&argc, &argv);
  int threads = GridThread::GetThreads();
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

  GridCartesian         *UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
                                                                  GridDefaultSimd(Nd,vComplex::Nsimd()),
                                                                  GridDefaultMpi());

  // -------------------------------------
  // for reading data
  const int conf_min=atoi(argv[1]);

  const std::string base_dir(argv[2]); // directory of lattice config
  // -> string path

  const std::string basedir(argv[3]); // not used

  const std::string basedir2(argv[4]); // directory for output .bin
  // -> string pathO

  const std::string mass(argv[5]);
  const int interval=atoi(argv[6]);
  // -> 10, 20

  const int nbeta=atoi(argv[7]);
  // 

  // -----------------------

  std::vector<std::string> betas;
  {
    for(int i=8; i<8+nbeta; i++) {
      std::string str(argv[i]);
      betas.push_back(str);
    }
    // std::cout << "str = " << str << std::endl;
    // std::istringstream s( str );
    // std::string d;
    // while (s >> d) betas.push_back(d);
  }

  for(auto elem : betas) std::cout << elem << std::endl;

  // const int binsize = atoi(argv[8]);
  // const int nbin = atoi(argv[8]);
  // const int ibin = atoi(argv[9]);

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
      if( !std::filesystem::exists( path ) ) break;
      // conf_max = conf_min-interval;
    }

    std::vector<std::string> obsinfolist={"polyakov_abs","polyakov_re","polyakov_im","plaquette"};

    for(int conf=conf_min; conf<conf_max; conf+=interval){
      std::cout << "obs" << std::endl;
      Obs obs;
      obs.beta = beta;
      obs.base_dir=base_dir;
      obs.id=get_configname(beta, mass);

      // std::cout << "pathO" << std::endl;
      // std::cout << "pathO = " << pathO << std::endl;
      bool exists = true;
      for(std::string oi : obsinfolist){
        std::string pathO = basedir2+obs.id+oi+std::to_string(conf)+".bin";
        if( !std::filesystem::exists( pathO ) ) exists = false;
      }
      if(exists) continue;

      std::cout << "prefix" << std::endl;
      char f[200];
      std::string prefix = get_prefix(beta, mass);
      std::sprintf(f,
                   "%s.%d",
                   prefix.data(), conf);

      std::cout << "path" << std::endl;
      const std::string f_str = f;
      const std::string path = obs.base_dir+obs.id+"/"+f_str;
      std::cout << "path2 = " << path << std::endl;
      FieldMetaData header;
      LatticeGaugeField U(UGrid);
      NerscIO::readConfiguration(U, header, path);
      std::cout << "read2. " << std::endl;

      {
        obsinfo = "polyakov_abs";
        std::string pathO = basedir2+obs.id+obsinfo+std::to_string(conf)+".bin";
        const ComplexD p    = WilsonLoops<Impl>::avgPolyakovLoop(U);
        obs.O = abs(p);
        const ComplexD S    = Waction.S(U);
        obs.energy = real(S)/std::stod(beta);

        BinaryWriter WR( pathO );
        write(WR, "obs", obs );
      }

      {
        obsinfo = "polyakov_re";
        std::string pathO = basedir2+obs.id+obsinfo+std::to_string(conf)+".bin";
        const ComplexD p    = WilsonLoops<Impl>::avgPolyakovLoop(U);
        obs.O = real(p);
        const ComplexD S    = Waction.S(U);
        obs.energy = real(S)/std::stod(beta);

        BinaryWriter WR( pathO );
        write(WR, "obs", obs );
      }

      {
        obsinfo = "polyakov_im";
        std::string pathO = basedir2+obs.id+obsinfo+std::to_string(conf)+".bin";
        const ComplexD p    = WilsonLoops<Impl>::avgPolyakovLoop(U);
        obs.O = real(p);
        const ComplexD S    = Waction.S(U);
        obs.energy = real(S)/std::stod(beta);

        BinaryWriter WR( pathO );
        write(WR, "obs", obs );
      }

      {
        obsinfo = "plaquette";
        std::string pathO = basedir2+obs.id+obsinfo+std::to_string(conf)+".bin";
        const ComplexD S    = Waction.S(U);
        obs.energy = real(S)/std::stod(beta);
        double volume = GridDefaultLatt()[0]*GridDefaultLatt()[1]*GridDefaultLatt()[2]*GridDefaultLatt()[3];
        obs.O = real(S)/std::stod(beta)/volume/6;

        BinaryWriter WR( pathO );
        write(WR, "obs", obs );
      }
    }
  }

  Grid_finalize();
}
