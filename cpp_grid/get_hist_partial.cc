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
  const int nbins_histo_y=atoi(argv[14]);
  double min_y = atof(argv[15]);
  double max_y = atof(argv[16]);

  const int blocksize=8;

  std::cout << "nbeta = " << nbeta << std::endl;
  std::cout << "basedir = " << basedir << std::endl;
  std::cout << "nbins_histo_x = " << nbins_histo_x << std::endl;
  std::cout << "min_x = " << min_x << std::endl;
  std::cout << "max_x = " << max_x << std::endl;
  std::cout << "nbins_histo_y = " << nbins_histo_y << std::endl;
  std::cout << "min_y = " << min_y << std::endl;
  std::cout << "max_y = " << max_y << std::endl;

  // std::string obsinfoR="flow_polyakov_re";
  // std::string obsinfoI="flow_polyakov_im";

  int type=1;
  if(basedir2==basedir) {
    if(mass=="0p1000" || mass=="0p4000") type=2;
    else if(mass=="0p0500" || mass=="0p0100" || mass=="0p2000" || mass=="0p3000") type=4;
    else if(mass=="0") type=6;
  }
  else if(mass=="0") type = 0;

  if(runtype>=2) {
    type+=1;
    runtype -= 2;
  }


  std::vector<std::string> betas;
  std::vector<double> betas_double;
  {
    for(int i=17; i<17+nbeta; i++) {
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

  std::vector<double> ys(nbins_histo_y);
  const Real delta_y = (max_y-min_y)/nbins_histo_y;
  for(int i=0; i<nbins_histo_y; i++) ys[i] = min_y + delta_y*(0.5+i);
  std::cout << "ys = ";
  for(auto elem : ys) std::cout << elem << " ";
  std::cout << std::endl;



  std::cout << "conf_max0 = " << conf_max0 << std::endl;

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

      std::string pathO;
      if(type>1) pathO = basedir2+get_configname(beta, mass, type)+std::to_string(conf)+".h5";
      else pathO = basedir2+get_configname(beta, mass, type)+std::to_string(conf)+".bin";
      if( !std::filesystem::exists( pathO ) ) {
        std::cout << "missing pathO = " << pathO << std::endl;
        break;
      }
      conf_max = conf;
    }

    std::cout << "conf_max = " << conf_max << std::endl;

    for(int conf=conf_min0; conf<conf_max; conf+=interval){
      std::string obs_id=get_configname(beta, mass, type);

      LatticeComplexD pfield( UGrid );
      std::string pathO2 = basedir2+obs_id+"P"+std::to_string(conf)+".bin";
      {
        emptyUserRecord record;
        Grid::ScidacReader wr; // (UGrid->IsBoss());
        wr.open( pathO2 );
        wr.readScidacFieldRecord( pfield, record );
        wr.close();
      }

      const int T = UGrid->GlobalDimensions()[3];
      const int X = UGrid->GlobalDimensions()[0];
      const int Y = UGrid->GlobalDimensions()[1];
      const int Z = UGrid->GlobalDimensions()[2];

      std::vector<std::vector<double>> hist( nbins_histo_x,
                                             std::vector<double>(nbins_histo_y, 0.0));

      const int nblocks = X/blocksize;
      for(int bx=0; bx<nblocks; bx++){
        for(int by=0; by<nblocks; by++){
          for(int bz=0; bz<nblocks; bz++){
            //
            ComplexD blockavg = 0.0;
            for(int ix=bx*blocksize; ix<(bx+1)*blocksize; ix++){
              for(int iy=by*blocksize; iy<(by+1)*blocksize; iy++){
                for(int iz=bz*blocksize; iz<(bz+1)*blocksize; iz++){
                  ComplexD val = 0.0;
                  const Coordinate p({ix,iy,iz,0});
                  peekSite(val, pfield, p );
                  blockavg += val;
                }}}
            blockavg /= blocksize*blocksize*blocksize;

            for(int ibx=0; ibx<nbins_histo_x; ibx++){ // histogram binning for
              for(int iby=0; iby<nbins_histo_y; iby++){ // histogram binning for
                hist[ibx][iby] += regularized_delta( real(blockavg), imag(blockavg),
                                                     xs[ibx], ys[iby],
                                                     delta_x, delta_y );
              } // end for iby
            } // end for ibx
          }}} // end block

      for(int ibx=0; ibx<nbins_histo_x; ibx++){ // histogram binning for
        for(int iby=0; iby<nbins_histo_y; iby++){ // histogram binning for
          hist[ibx][iby] /= nblocks*nblocks*nblocks;
        }}

      std::string pathO3 = basedir2+obs_id+"P_partial"+std::to_string(conf)+".xml";
      // std::cout << "pathO3 = " << pathO3 << std::endl;
      if( std::filesystem::exists( pathO3 ) && (runtype==0) ) {
        std::cout << "pathO3 = " << pathO3 << std::endl;
        std::cout << GridLogMessage << "skip." << std::endl;
        continue;
      }
      {
        XmlWriter WR( pathO3 );
        WR.scientificFormat( true ); // @@@ double check
        WR.setPrecision( std::numeric_limits<Real>::digits10 + 1 );
        write(WR, "hist", hist );
        write(WR, "xs", xs );
        write(WR, "ys", ys );
      }
    } // end conf
  } // end beta



  // -----------------------


    Grid_finalize();
}
