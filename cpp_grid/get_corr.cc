// #include "reweight.h"
#include "reweight2.h"

#include <regex>

#include "fftw3.h"
// #define NOSIGN


template <class Field>
void fft3d( Field& Fv, const Field& v, const int sign=FFT::forward ){
  std::vector<int> dim_mask({1,1,1,0}); // 3d FFT
  FFT theFFT((GridCartesian *) v.Grid());
  theFFT.FFT_dim_mask( Fv, v, dim_mask, sign );
}



// template <class Field>
void fft3d( const std::vector<std::vector<std::vector<std::complex<double>>>>& sin,
            std::vector<std::vector<std::vector<std::complex<double>>>>& sout,
            const int X, const int Y, const int Z,
            const int sign=FFTW_FORWARD,
            const double factor = 1.0 ){
  fftw_complex *in, *out;
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * X*Y*Z);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * X*Y*Z);

  {
    int counter=0;
    for (int x=0;x<X;x++){
      for (int y=0;y<Y;y++){
        for (int z=0;z<Z;z++){
          in[counter][0] = sin[x][y][z].real();
          in[counter][1] = sin[x][y][z].imag();
          // std::cout << in[counter][0] << " " << in[counter][1] << std::endl;
          counter++;
          // double phx = 2.0*M_PI * x / X;
          // double phy = 2.0*M_PI * y / Y;
          // double phz = 2.0*M_PI * z / Z;
          // in[counter][0] = std::cos( phx + phy + phz );
          // in[counter][1] = std::sin( phx + phy + phz );
        }}}
  }

  fftw_plan p;
  p = fftw_plan_dft_3d(X, Y, Z,
                       in, out,
                       sign, FFTW_ESTIMATE);

  fftw_execute(p); /* repeat as needed */
  {
    int counter=0;
    for (int x=0;x<X;x++){
      for (int y=0;y<Y;y++){
        for (int z=0;z<Z;z++){
          sout[x][y][z] = std::complex<double>(out[counter][0], out[counter][1]) * factor;
          // std::cout << out[counter][0] << " " << out[counter][1] << std::endl;
          counter++;
        }}}
  }

  fftw_destroy_plan(p);
  fftw_free(in); fftw_free(out);
}




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

  std::cout << "nbeta = " << nbeta << std::endl;
  std::cout << "basedir = " << basedir << std::endl;

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


  // int type=1;
  // if(basedir2c==basedir2h) { !!!!!!! @@@@@@@
  //   if(mass=="0p1000") type=2;
  //   else if(mass=="0p0500" || mass=="0p0100" || mass=="0p2000") type=4;
  // }

  std::vector<std::string> betas;
  std::vector<double> betas_double;
  {
    for(int i=11; i<11+nbeta; i++) {
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



  GridCartesian         *UGrid   = SpaceTimeGrid::makeFourDimGrid(GridDefaultLatt(),
                                                                  GridDefaultSimd(Nd,vComplex::Nsimd()),
                                                                  GridDefaultMpi());



  for(std::string beta : betas){
    std::cout << "beta = " << beta << std::endl;
    std::string obs_id=get_configname(beta, mass, type);

    int conf_max=conf_min0;
    for(int conf=conf_min0; conf<conf_max0; conf+=interval){
      char f[200];
      std::string prefix = get_prefix(beta, mass);
      std::sprintf(f,
                   "%s.%d",
                   prefix.data(), conf);


      std::string pathO2 = basedir2+get_configname(beta, mass, type)+"P"+std::to_string(conf)+".bin";
      if( !std::filesystem::exists( pathO2 ) ) break;
      conf_max = conf;
    }

    std::cout << "conf_max = " << conf_max << std::endl;

    for(int conf=conf_min0; conf<conf_max; conf+=interval){

      std::string pathO3 = basedir2+obs_id+"P_corr"+std::to_string(conf)+".h5";
      if( std::filesystem::exists( pathO3 ) && (runtype==0) ) {
        continue;
      }

      LatticeComplexD pfield( UGrid );
      std::string pathO2 = basedir2+obs_id+"P"+std::to_string(conf)+".bin";

      {
        emptyUserRecord record;
        Grid::ScidacReader wr; // (UGrid->IsBoss());
        wr.open( pathO2 );
        wr.readScidacFieldRecord( pfield, record );
        wr.close();
      }
      int T = UGrid->GlobalDimensions()[3];
      int X = UGrid->GlobalDimensions()[0];
      int Y = UGrid->GlobalDimensions()[1];
      int Z = UGrid->GlobalDimensions()[2];

      std::vector<std::vector<std::vector<std::complex<double>>>> p
        (X, std::vector<std::vector<std::complex<double>>>
         (Y, std::vector<std::complex<double>>(Z, 0.0)));
      std::vector<std::vector<std::vector<std::complex<double>>>> pf
        (X, std::vector<std::vector<std::complex<double>>>
         (Y, std::vector<std::complex<double>>(Z, 0.0)));
      std::vector<std::vector<std::vector<std::complex<double>>>> pf2
        (X, std::vector<std::vector<std::complex<double>>>
         (Y, std::vector<std::complex<double>>(Z, 0.0)));


#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (int x=0;x<X;x++){
        for (int y=0;y<Y;y++){
          for (int z=0;z<Z;z++){
            Coordinate r({x,y,z,0});
            peekLocalSite( p[x][y][z], pfield, r );
          }}}

      fft3d( p, pf, X, Y, Z, FFTW_FORWARD, 1.0/std::sqrt(X*Y*Z) );

#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (int x=0;x<X;x++){
        for (int y=0;y<Y;y++){
          for (int z=0;z<Z;z++){
            pf2[x][y][z] = pf[x][y][z] * std::conj(pf[x][y][z]);
          }}}

      fft3d( pf2, p, X, Y, Z, FFTW_BACKWARD, 1.0/std::sqrt(X*Y*Z) );

      // pf[0][0][0] = 0.0;
      // fft3d( pf, pf, X, Y, Z, FFTW_BACKWARD, 1.0/std::sqrt(X*Y*Z) );
      std::vector<RealD> xcorr(X, 0.0);
      for (int x=0;x<X;x++) xcorr[x] = real(p[x][0][0]);

      {
        if (UGrid->IsBoss()){
          std::unique_ptr<Hdf5Writer> WR;
          WR = std::make_unique<Hdf5Writer>( pathO3 );
          write(*WR, "xcorr", xcorr );
          // write(*WR, "p", p );
          // write(*WR, "pf", pf2 );
          // write(*WR, "pf0", pf );
        }
      }
    } // end conf
  } // for beta // end data prep


  // -----------------------





    Grid_finalize();
}
