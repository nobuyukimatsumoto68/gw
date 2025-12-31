// #include "reweight.h"
#include "reweight2.h"

#include <regex>

// #define NOSIGN



// template <class Gimpl, class GaugeField>
// void PolyakovField( LatticeComplexD& out, const GaugeField& Umu) {  //assume Nd=4
//   typedef typename Gimpl::GaugeLinkField GaugeMat;
//   GaugeMat Ut(Umu.Grid()), P(Umu.Grid());
//   int T = Umu.Grid()->GlobalDimensions()[3];
//   int X = Umu.Grid()->GlobalDimensions()[0];
//   int Y = Umu.Grid()->GlobalDimensions()[1];
//   int Z = Umu.Grid()->GlobalDimensions()[2];

//   Ut = peekLorentz(Umu,3); //Select temporal direction
//   P = Ut;
//   for (int t=1;t<T;t++){
//     P = Gimpl::CovShiftForward(Ut,3,P);
//   }
//   RealD norm = 1.0/(Nc*X*Y*Z*T);
//   out = trace(P)*norm;
// }





int main(int argc, char **argv) {
  Grid_init(&argc, &argv);
  int threads = GridThread::GetThreads();
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

  // -------------------------------------
  // for reading data
  const int conf_min0=atoi(argv[1]);
  const std::string base_dir(argv[2]);
  const std::string basedir(argv[3]);
  const std::string basedir2(argv[4]);
  // const std::string basedir2h(argv[5]);
  const std::string mass(argv[5]);
  const int interval=atoi(argv[6]); // !!!!!!!!! @@@@@@@@@@@
  const int conf_max0=atoi(argv[7]);
  // const int binsize=atoi(argv[9]);
  const int nbin_global=atoi(argv[8]); // dummy
  int nbeta=atoi(argv[9]);
  int runtype=atoi(argv[10]);
  // obsinfo = argv[12]; // dummy
  const int nbins_histo_x=atoi(argv[11]);
  double min_x = atof(argv[12]);
  double max_x = atof(argv[13]);

  std::cout << "nbeta = " << nbeta << std::endl;
  std::cout << "basedir = " << basedir << std::endl;
  std::cout << "nbins_histo_x = " << nbins_histo_x << std::endl;
  std::cout << "min_x = " << min_x << std::endl;
  std::cout << "max_x = " << max_x << std::endl;

  // int type=1;
  // if(basedir2c==basedir2h) { !!!!!!! @@@@@@@
  //   if(mass=="0p1000") type=2;
  //   else if(mass=="0p0500" || mass=="0p0100" || mass=="0p2000") type=4;
  // }

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


  // -----------------------




  std::vector<double> xs(nbins_histo_x);
  const Real delta_x = (max_x-min_x)/nbins_histo_x;
  for(int i=0; i<nbins_histo_x; i++) xs[i] = min_x + delta_x*(0.5+i);
  std::cout << "xs = ";
  for(auto elem : xs) std::cout << elem << " ";
  std::cout << std::endl;


  GridCartesian         *UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
                                                                  GridDefaultSimd(Nd,vComplex::Nsimd()),
                                                                  GridDefaultMpi());

  for(std::string beta : betas){
    std::cout << "beta = " << beta << std::endl;

    int conf_max=conf_min0;
    for(int conf=conf_min0; conf<conf_max0; conf+=interval){
      char f[200];
      std::string prefix = get_prefix(beta, mass);
      std::sprintf(f,
                   "%s.%d",
                   prefix.data(), conf);

      std::string pathO2 = basedir2+get_configname(beta, mass)+"chcond"+std::to_string(conf)+".bin";
      if( !std::filesystem::exists( pathO2 ) ) {
        std::cout << "missing pathO2 = " << pathO2 << std::endl;
        break;
      }
      conf_max = conf;
    }

    std::cout << "conf_max = " << conf_max << std::endl;

    for(int conf=conf_min0; conf<conf_max; conf+=interval){
      std::string obs_id=get_configname(beta, mass);

      std::string pathO3 = basedir2+obs_id+"chcond_hist"+std::to_string(conf)+".bin";
      if( std::filesystem::exists( pathO3 ) && (runtype==0) ) {
        // std::cout << GridLogMessage << "skip." << std::endl;
        continue;
      }

      Complex chcond;
      std::string pathO2 = basedir2+obs_id+"chcond"+std::to_string(conf)+".bin";
      {
        XmlReader WR( pathO2 );
        read(WR, "chcond", chcond );
      }
      // std::cout << "debug. chcond = " << chcond << std::endl;

      std::vector<double> hist( nbins_histo_x, 0.0 );
      for(int ibx=0; ibx<nbins_histo_x; ibx++){ // histogram binning for
        hist[ibx] = regularized_delta( real(chcond), xs[ibx], delta_x );
      } // end for ibx

      {
        XmlWriter WR( pathO3 );
        write(WR, "hist", hist );
        write(WR, "xs", xs );
      }
    } // end conf
  } // end beta



  // -----------------------





    Grid_finalize();
}
