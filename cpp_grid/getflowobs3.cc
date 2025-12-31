#include "reweight2.h"
#include <sstream>


template <class Impl>
class FlowObs {
  int nsteps;
  double step_size;
  int meas_interval;
  double t0;

public:
  // here forces the Impl to be of gauge fields
  // if not the compiler will complain
  INHERIT_GIMPL_TYPES(Impl);

  // necessary for HmcObservable compatibility
  typedef typename Impl::Field Field;

  FlowObs(const int nsteps,
          const double step_size,
          const int meas_interval )
    : nsteps(nsteps)
    , step_size(step_size)
    , meas_interval(meas_interval)
    , t0(step_size*nsteps)
  {
  }

  ComplexD operator()( std::vector<RealD>& all_polyakov_re,
                       std::vector<RealD>& all_polyakov_im,
                       std::vector<RealD>& all_T0,
                       std::vector<RealD>& all_Q,
                       Field& Usmear, const Field &U ){
    int def_prec = std::cout.precision();

    all_polyakov_re.clear();
    all_polyakov_im.clear();
    all_T0.clear();
    all_Q.clear();

    // Real q    = WilsonLoops<Impl>::TopologicalCharge(U);
    // std::cout << GridLogMessage
    //           << std::setprecision(std::numeric_limits<Real>::digits10 + 1)
    //           << "Topological Charge (before flow): "<< q << std::endl;
    ComplexD p    = WilsonLoops<Impl>::avgPolyakovLoop(U);
    std::cout << GridLogMessage
              << std::setprecision(std::numeric_limits<Real>::digits10 + 1)
              << "Polyakov Loop (before flow): "<< p << std::endl;

    WilsonFlow<PeriodicGimplR> WF(step_size, nsteps, meas_interval);

    WF.addMeasurement(meas_interval,
                      [&](int step, RealD t, const typename Impl::GaugeField &U){
                        const ComplexD p1 = WilsonLoops<Impl>::avgPolyakovLoop(U);
                        all_polyakov_re.push_back( real(p1) );
                        all_polyakov_im.push_back( imag(p1) );
                        std::cout << GridLogMessage
                                  << std::setprecision(std::numeric_limits<Real>::digits10 + 1)
                                  << "[WilsonFlow] Polyakov loop           : "  << step << "  "
                                  << p1
                                  << std::endl;
                      });
    WF.addMeasurement(meas_interval,
                      [&](int step, RealD t, const typename Impl::GaugeField &U){
                        const RealD p1 = WF.energyDensityPlaquette(t, U);
                        all_T0.push_back( p1 );
                        std::cout << GridLogMessage
                                  << std::setprecision(std::numeric_limits<Real>::digits10 + 1)
                                  << "[WilsonFlow] T0           : "  << step << "  "
                                  << p1
                                  << std::endl;
                      });
    WF.addMeasurement(meas_interval,
                      [&](int step, RealD t, const typename Impl::GaugeField &U){
                        const RealD p1 = WilsonLoops<Impl>::TopologicalCharge(U);
                        all_Q.push_back( p1 );
                        std::cout << GridLogMessage
                                  << std::setprecision(std::numeric_limits<Real>::digits10 + 1)
                                  << "[WilsonFlow] topo charge        : "  << step << "  "
                                  << p1
                                  << std::endl;
                      });

    WF.smear(Usmear, U);
    Real T0   = WF.energyDensityPlaquette(t0, Usmear);
    std::cout << GridLogMessage << std::setprecision(std::numeric_limits<Real>::digits10 + 1)
              << "T0                : " << T0 << std::endl;

    // q    = WilsonLoops<Impl>::TopologicalCharge(Usmear);
    // std::cout << GridLogMessage
    //           << std::setprecision(std::numeric_limits<Real>::digits10 + 1)
    //           << "Topological Charge: "<< q << std::endl;
    p    = WilsonLoops<Impl>::avgPolyakovLoop(Usmear);
    std::cout << GridLogMessage
              << std::setprecision(std::numeric_limits<Real>::digits10 + 1)
              << "Polyakov Loop: "<< p << std::endl;

    std::cout.precision(def_prec);
    return p;
  }

};




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
  const int Nt=8;
  const int cinv=2;//atoi(argv[8]);


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

  const double step_size = 0.01;
  const int nsteps = 300; // (int)(Nt*Nt/cinv/cinv/8/step_size);
  const int meas_interval = 1;

  FlowObs<PeriodicGimplD> flowobs(nsteps, step_size, meas_interval);


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

    std::vector<std::string> obsinfolist={"polyakov_abs","polyakov_re","polyakov_im","plaquette","flow_polyakov_re","flow_polyakov_im", "flow_T0", "flow_Q"};

    for(int conf=conf_min; conf<conf_max; conf+=interval){
      // std::cout << "obs" << std::endl;
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

      // std::cout << "prefix" << std::endl;
      char f[200];
      std::string prefix = get_prefix(beta, mass);
      std::sprintf(f,
                   "%s.%d",
                   prefix.data(), conf);

      // std::cout << "path" << std::endl;
      const std::string f_str = f;
      const std::string path = obs.base_dir+obs.id+"/"+f_str;
      std::cout << "path2 = " << path << std::endl;
      FieldMetaData header;
      LatticeGaugeField U(UGrid);
      NerscIO::readConfiguration(U, header, path);
      // std::cout << "read2. " << std::endl;

      const ComplexD S    = Waction.S(U);

      {
        obsinfo = "polyakov_abs";
        std::string pathO = basedir2+obs.id+obsinfo+std::to_string(conf)+".bin";
        const ComplexD p    = WilsonLoops<Impl>::avgPolyakovLoop(U);
        obs.O = abs(p);
        obs.energy = real(S)/std::stod(beta);

        BinaryWriter WR( pathO );
        write(WR, "obs", obs );
      }

      {
        obsinfo = "polyakov_re";
        std::string pathO = basedir2+obs.id+obsinfo+std::to_string(conf)+".bin";
        const ComplexD p    = WilsonLoops<Impl>::avgPolyakovLoop(U);
        obs.O = real(p);
        obs.energy = real(S)/std::stod(beta);

        BinaryWriter WR( pathO );
        write(WR, "obs", obs );
      }

      {
        obsinfo = "polyakov_im";
        std::string pathO = basedir2+obs.id+obsinfo+std::to_string(conf)+".bin";
        const ComplexD p    = WilsonLoops<Impl>::avgPolyakovLoop(U);
        obs.O = real(p);
        obs.energy = real(S)/std::stod(beta);

        BinaryWriter WR( pathO );
        write(WR, "obs", obs );
      }

      {
        obsinfo = "plaquette";
        std::string pathO = basedir2+obs.id+obsinfo+std::to_string(conf)+".bin";
        obs.energy = real(S)/std::stod(beta);
        double volume = GridDefaultLatt()[0]*GridDefaultLatt()[1]*GridDefaultLatt()[2]*GridDefaultLatt()[3];
        obs.O = real(S)/std::stod(beta)/volume/6;

        BinaryWriter WR( pathO );
        write(WR, "obs", obs );
      }

      Impl::Field Usmear = U;
      std::vector<RealD> all_polyakov_re;
      std::vector<RealD> all_polyakov_im;
      std::vector<RealD> all_T0;
      std::vector<RealD> all_Q;
      const ComplexD p = flowobs( all_polyakov_re, all_polyakov_im, all_T0, all_Q,
                                  Usmear, U );

      VectorObs vobs;
      vobs.beta = beta;
      vobs.base_dir=base_dir;
      vobs.id=get_configname(beta, mass);

      {
        obsinfo = "flow_polyakov_re";
        std::string pathO = basedir2+vobs.id+obsinfo+std::to_string(conf)+".bin";
        vobs.vO = all_polyakov_re;
        vobs.energy = real(S)/std::stod(beta);

        BinaryWriter WR( pathO );
        write(WR, "obs", vobs );
      }
      {
        obsinfo = "flow_polyakov_im";
        std::string pathO = basedir2+vobs.id+obsinfo+std::to_string(conf)+".bin";
        vobs.vO = all_polyakov_im;
        vobs.energy = real(S)/std::stod(beta);

        BinaryWriter WR( pathO );
        write(WR, "obs", vobs );
      }
      {
        obsinfo = "flow_T0";
        std::string pathO = basedir2+vobs.id+obsinfo+std::to_string(conf)+".bin";
        vobs.vO = all_T0;
        vobs.energy = real(S)/std::stod(beta);

        BinaryWriter WR( pathO );
        write(WR, "obs", vobs );
      }
      {
        obsinfo = "flow_Q";
        std::string pathO = basedir2+vobs.id+obsinfo+std::to_string(conf)+".bin";
        vobs.vO = all_Q;
        vobs.energy = real(S)/std::stod(beta);

        BinaryWriter WR( pathO );
        write(WR, "obs", vobs );
      }
    }
  }

  Grid_finalize();
}
