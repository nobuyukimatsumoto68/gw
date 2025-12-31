#include <Grid/Grid.h>
// #include "/mnt/hdd_barracuda/hmc/header.hpp"

using namespace std;
using namespace Grid;


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

  void operator()( const Field &U ){
    Field Usmear = U;
    int def_prec = std::cout.precision();

    Real q    = WilsonLoops<Impl>::TopologicalCharge(U);
    std::cout << GridLogMessage
              << std::setprecision(std::numeric_limits<Real>::digits10 + 1)
              << "Topological Charge (before flow): "<< q << std::endl;
    ComplexD p    = WilsonLoops<Impl>::avgPolyakovLoop(U);
    std::cout << GridLogMessage
              << std::setprecision(std::numeric_limits<Real>::digits10 + 1)
              << "Polyakov Loop (before flow): "<< p << std::endl;

    WilsonFlow<PeriodicGimplR> WF(step_size, nsteps, meas_interval);

    WF.addMeasurement(meas_interval,
                      [](int step, RealD t, const typename Impl::GaugeField &U){
                        std::cout << GridLogMessage
                                  << std::setprecision(std::numeric_limits<Real>::digits10 + 1)
                                  << "[WilsonFlow] Polyakov loop           : "  << step << "  "
                                  << WilsonLoops<Impl>::avgPolyakovLoop(U)
                                  << std::endl;
                      });

    WF.smear(Usmear, U);
    Real T0   = WF.energyDensityPlaquette(t0, Usmear);
    std::cout << GridLogMessage << std::setprecision(std::numeric_limits<Real>::digits10 + 1)
              << "T0                : " << T0 << std::endl;

    q    = WilsonLoops<Impl>::TopologicalCharge(Usmear);
    std::cout << GridLogMessage
              << std::setprecision(std::numeric_limits<Real>::digits10 + 1)
              << "Topological Charge: "<< q << std::endl;
    p    = WilsonLoops<Impl>::avgPolyakovLoop(Usmear);
    std::cout << GridLogMessage
              << std::setprecision(std::numeric_limits<Real>::digits10 + 1)
              << "Polyakov Loop: "<< p << std::endl;

    std::cout.precision(def_prec);
  }

};


int main (int argc, char ** argv)
{
  Grid_init(&argc,&argv);

  GridCartesian* UGrid = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
                                                        GridDefaultSimd(Nd,vComplex::Nsimd()),
                                                        GridDefaultMpi());

  LatticeGaugeField U(UGrid);

  // std::string path_base = "/mnt/hdd_barracuda/hmc/configs_grid";
  // int traj = 0;
  // std::string config;
  if( argc > 1 && argv[1][0] != '-' ){
    // std::cout<<GridLogMessage <<"Loading configuration from "<<argv[1]<<std::endl;
    // FieldMetaData header;
    // traj = std::atoi(argv[1]);
    // std::string path = path_base+"/ckpoint_"+std::to_string(traj)+".lat";
    // NerscIO::readConfiguration(U, header, path);
    std::cout << GridLogMessage << "Loading configuration from " << argv[1] << std::endl;
    FieldMetaData header;
    NerscIO::readConfiguration(U, header, argv[1]);
    // config=argv[1];
  }
  else {
    SU<Nc>::ColdConfiguration(U);
    // config="cold";
  }

  // RealD plaq = WilsonLoops<Impl>::avgPlaquette(U);

  const int interval = 1;
  const int nsteps = 200;
  const double step_size = 0.01;
  const int meas_interval = 10;

  FlowObs<PeriodicGimplD> obs(nsteps, step_size, meas_interval);
  obs( U );

  Grid_finalize();
}


