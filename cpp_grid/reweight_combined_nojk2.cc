// #include "reweight.h"
#include "reweight2.h"

#include <regex>

// #define NOSIGN

int main(int argc, char **argv) {
  Grid_init(&argc, &argv);
  int threads = GridThread::GetThreads();
  std::cout << GridLogMessage << "Grid is setup to use " << threads << " threads" << std::endl;

  // -------------------------------------
  // for reading data
  const int conf_min0=atoi(argv[1]);
  const std::string base_dir(argv[2]);
  const std::string basedir(argv[3]);
  const std::string basedir2c(argv[4]);
  const std::string basedir2h(argv[5]);
  const std::string mass(argv[6]);
  const int conf_max0=atoi(argv[7]);
  // const int binsize=atoi(argv[9]);
  const int nbin_global=atoi(argv[8]);
  int nbeta=atoi(argv[9]);
  int runtype=atoi(argv[10]);
  // obsinfo = argv[12]; // dummy
  const int nbins_histo_x=atoi(argv[11]);
  double min_x = atof(argv[12]);
  double max_x = atof(argv[13]);
  const int nbins_histo_y=atoi(argv[14]);
  double min_y = atof(argv[15]);
  double max_y = atof(argv[16]);

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
  if(basedir2c==basedir2h) {
    if(mass=="0p1000") type=2;
    else if(mass=="0p0500" || mass=="0p0100" || mass=="0p2000") type=4;
  }

  // -----------------------

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

  using Impl = PeriodicGimplR;

#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(std::string beta : betas){
    std::cout << GridLogMessage << "beta = " << beta << std::endl;
    EnsembleOfObs ensR;
    ensR.beta = beta;
    ensR.base_dir=base_dir;
    if(type!=1) ensR.id=get_configname(beta, mass, type+1);
    else ensR.id=get_configname(beta, mass, type);

    EnsembleOfObs ensI;
    ensI.beta = beta;
    ensI.base_dir=base_dir;
    if(type!=1) ensI.id=get_configname(beta, mass, type+1);
    else ensI.id=get_configname(beta, mass, type);

    std::string path_ens = basedir+ensR.id+".bin";
    // std::cout << "path_ens = " << path_ens << std::endl;
    if( std::filesystem::exists( path_ens ) && (runtype==0) ) {
      // std::cout << GridLogMessage << "skip." << std::endl;
      continue;
    }

    // ---------------------------
    // cold start
    // ---------------------------

    ensR.id=get_configname(beta, mass, type);
    ensI.id=get_configname(beta, mass, type);

    for(int conf=conf_min0; conf<conf_max0; conf++){
      // std::string pathO = basedir2c+ensR.id+obsinfoR+std::to_string(conf)+".bin";
      std::string pathO = basedir2c+ensR.id+std::to_string(conf)+".bin";
      // std::cout << "pathO = " << pathO << std::endl;
      if( !std::filesystem::exists( pathO ) ) continue; // assert(false);
      // std::cout << "debug. pt1" << std::endl;

      Obs obsR, obsI;
      // if (UGrid->IsBoss()){
#pragma omp critical
      {
        std::unique_ptr<Hdf5Reader> WR;
        WR = std::make_unique<Hdf5Reader>( pathO );
        // Hdf5Reader WR( pathO );
        read(*WR, "obs", obsR );
        // std::cout << "debug. pt2" << std::endl;
        // read(WR, "obs", obsR );
        obsI = obsR;
        // std::cout << "debug. pt3" << std::endl;
        std::vector<ComplexD> all_polyakov;
        read(*WR, "polyakov", all_polyakov );
        // read(WR, "polyakov", all_polyakov );
        // std::cout << "debug. pt4" << std::endl;
        obsR.O = real( all_polyakov.back() );
        obsI.O = imag( all_polyakov.back() );
        // std::cout << "debug. pt5" << std::endl;
      }
      ensR.Os.push_back(obsR.O);
      ensR.Osqs.push_back(obsR.O*obsR.O);
      ensR.energies.push_back(obsR.energy);
      //
      ensI.Os.push_back(obsI.O);
      ensI.Osqs.push_back(obsI.O*obsI.O);
      ensI.energies.push_back(obsI.energy);
    } // end for conf

    // ---------------------------
    // hot start
    // ---------------------------

    if(type!=1) {
      ensR.id=get_configname(beta, mass, type+1);
      ensI.id=get_configname(beta, mass, type+1);
    }

    for(int conf=conf_min0; conf<conf_max0; conf++){
      // std::string pathO = basedir2h+ensR.id+obsinfoR+std::to_string(conf)+".bin";
      std::string pathO = basedir2h+ensR.id+std::to_string(conf)+".bin";
      if( !std::filesystem::exists( pathO ) ) continue; // assert(false);

      Obs obsR, obsI;
      // std::unique_ptr<Hdf5Reader> WR;
      // if (UGrid->IsBoss()){
#pragma omp critical
      {
        // WR = std::make_unique<Hdf5Reader>( pathO );
        // Hdf5Reader WR( pathO );
        std::unique_ptr<Hdf5Reader> WR;
        WR = std::make_unique<Hdf5Reader>( pathO );
        // read(*WR, "obs", obsR );
        read(*WR, "obs", obsR );
        obsI = obsR;
        std::vector<ComplexD> all_polyakov;
        // read(*WR, "polyakov", all_polyakov );
        read(*WR, "polyakov", all_polyakov );
        obsR.O = real( all_polyakov.back() );
        obsI.O = imag( all_polyakov.back() );
      }
      // BinaryReader WR( pathO );
      // read(WR, "obs", obs );
      ensR.Os.push_back(obsR.O);
      ensR.Osqs.push_back(obsR.O*obsR.O);
      ensR.energies.push_back(obsR.energy);
      //
      ensI.Os.push_back(obsI.O);
      ensI.Osqs.push_back(obsI.O*obsI.O);
      ensI.energies.push_back(obsI.energy);
    } // end for conf

    ensR.nbin=nbin_global;
    ensR.binsize = ensR.size()/nbin_global;
    ensR.nconf = ensR.nbin*ensR.binsize;
    ensI.nbin=nbin_global;
    ensI.binsize = ensI.size()/nbin_global;
    ensI.nconf = ensI.nbin*ensI.binsize;

    { // trimming
      const int start = ensR.nconf;
      // std::cout << GridLogMessage << "start = " << start << " size = " << ensR.size() << std::endl;
      ensR.energies.erase( ensR.energies.begin()+start, ensR.energies.end() );
      ensR.Os.erase( ensR.Os.begin()+start, ensR.Os.end() );
      ensR.Osqs.erase( ensR.Osqs.begin()+start, ensR.Osqs.end() );
    }
    { // trimming
      const int start = ensI.nconf;
      // std::cout << GridLogMessage << "start = " << start << " size = " << ensI.size() << std::endl;
      ensI.energies.erase( ensI.energies.begin()+start, ensI.energies.end() );
      ensI.Os.erase( ensI.Os.begin()+start, ensI.Os.end() );
      ensI.Osqs.erase( ensI.Osqs.begin()+start, ensI.Osqs.end() );
    }

    std::cout << GridLogMessage << "beta = " << beta << " size = " << ensR.size() << " nbin = " << ensR.nbin << " nconf = " << ensR.nconf << std::endl;

#pragma omp critical
    {
      {
        std::unique_ptr<Hdf5Writer> WR;
        WR = std::make_unique<Hdf5Writer>( path_ens );
        write(*WR, "ensR", ensR );
        write(*WR, "ensI", ensI );
      }
    }
  } // for end beta



  // ################################################################ //
  std::cout << GridLogMessage << "run2 part v.2" << std::endl;

  if(type!=1) type+=1;

  {
    std::vector<EnsembleOfObs> ensembles0( nbeta );

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int j=0; j<nbeta; j++){
      // std::cout << GridLogMessage << "j = " << j << std::endl;
      const std::string beta = betas[j];
      std::string id_str = get_configname(beta, mass, type);
      // std::string path = basedir+id_str+obsinfoR+".bin";
      std::string path = basedir+id_str+".bin";
      // std::cout << GridLogMessage << "path = " << path << std::endl;

      // std::unique_ptr<Hdf5Reader> WR;
#pragma omp critical
      {
        std::unique_ptr<Hdf5Reader> WR;
        WR = std::make_unique<Hdf5Reader>( path );
        read(*WR, "ensR", ensembles0[j] );
      }
      // BinaryReader WR( path );
      // read(WR, "ens", ensembles0[j] );
    }

    // std::cout << GridLogMessage << "debug 1" << std::endl;

    {
      char id[200];
      std::sprintf(id, "m%s", mass.data());
      std::string id_str = id;
      std::string path = basedir+id_str+"_binsize.txt";

      std::ofstream file(path, std::ios::trunc);
      for(int j=0; j<nbeta; j++){
        file << ensembles0[j].binsize << std::endl;
      }
    }

    // std::cout << GridLogMessage << "debug 2" << std::endl;


    char id[200];
    std::sprintf(id, "m%s", mass.data());
    std::string id_str = id;
    id_str = id_str+"_nojk";

    // -------------------
    FreeEnergies fs(nbeta);

    {
      std::string id_str2 = id;
      id_str2 = id_str2+"_nojk";
      std::string filename2=basedir+id_str2+"fs.bin";
      if( std::filesystem::exists( filename2 ) ){
#pragma omp critical
        {
          std::unique_ptr<Hdf5Reader> WR;
          WR = std::make_unique<Hdf5Reader>( filename2 );
          read(*WR, "f", fs.f );
        }
        fs.fnew = fs.f;
        fs.fold = fs.f;
      }
    }

    // std::cout << GridLogMessage << "debug 3" << std::endl;
    std::vector<std::vector<RealD>> v_energies(nbeta);
    for(int k=0; k<nbeta; k++) v_energies[k] = ensembles0[k].energies;
    // std::cout << GridLogMessage << "debug 4" << std::endl;
    run_multihistogram( betas_double, v_energies, fs, basedir+id_str+"fs.bin" );
    // std::cout << GridLogMessage << "debug 5" << std::endl;

    // std::cout << GridLogxMessage << "debug. write " << std::endl;
    std::cout << "fs = ";
    for(auto elem : fs.f) std::cout << elem << " ";
    std::cout << std::endl;
    std::cout << GridLogMessage << "pathf = " << basedir+id_str+"fs.bin" << std::endl;

    // BinaryWriter WR(basedir+id_str+"fs.bin");
    // write(WR, "f", fs.f );
    #pragma omp critical
    // if (UGrid->IsBoss()){
    {
      std::unique_ptr<Hdf5Writer> WR;
      WR = std::make_unique<Hdf5Writer>( basedir+id_str+"fs.bin" );
      write(*WR, "f", fs.f );
    }
  }

//   // ################################################################ //
//   std::cout << GridLogMessage << "run3 part" << std::endl;

//   {
//     std::vector<EnsembleOfObs> ensembles0( nbeta );

// // #ifdef _OPENMP
// // #pragma omp parallel for
// // #endif
//     for(int j=0; j<nbeta; j++){
//       const std::string beta = betas[j];
//       std::string id_str = get_configname(beta, mass, type);
//       std::string path = basedir+id_str+obsinfoR+".bin";

//       BinaryReader WR( path );
//       read(WR, "ens", ensembles0[j] );
//     }

//     char id[200];
//     std::sprintf(id, "m%s", mass.data());

//     // -------------------

//     FreeEnergies fs(nbeta);

//     std::string id_str = id;
//     id_str = id_str+"_nojk";

//     BinaryReader WR(basedir+id_str+"fs.bin");
//     read(WR, "f", fs.f );

//     // =================================

//     ReweightedObs meas( betas_double.front(), betas_double.back(), nbeta*100, ensembles0 );
//     const int nmeas=meas.nmeas;
//     run_reweighting( betas_double, ensembles0, fs, meas );

//     std::cout << "fs = ";
//     for(auto elem : meas.fPs) std::cout << elem << " ";
//     std::cout << std::endl;

//     std::cout << "fPs = ";
//     for(auto elem : meas.fPs) std::cout << elem << " ";
//     std::cout << std::endl;

// #ifndef NOSIGN
//     std::cout << "sign_P = ";
//     for(auto elem : meas.sign_Ps) std::cout << elem << " ";
//     std::cout << std::endl;
// #endif

//     std::cout << "fP_sq = ";
//     for(auto elem : meas.fPs_sq) std::cout << elem << " ";
//     std::cout << std::endl;

//     {
//       std::cout << GridLogMessage << "path_beta = " << basedir+id_str+"meas_beta.bin" << std::endl;
//       BinaryWriter WR(basedir+id_str+"meas_beta.bin");
//       write(WR, "beta", meas.betas );
//     }
//     {
//       BinaryWriter WR(basedir+id_str+"meas_f.bin");
//       write(WR, "f", meas.fs );
//     }
//     {
//       BinaryWriter WR(basedir+id_str+"meas_fP_"+obsinfoR+".bin");
//       write(WR, "fP", meas.fPs );
//     }
// #ifndef NOSIGN
//     {
//       BinaryWriter WR(basedir+id_str+"meas_sign_P_"+obsinfoR+".bin");
//       write(WR, "sign_P", meas.sign_Ps );
//     }
// #endif
//     {
//       BinaryWriter WR(basedir+id_str+"meas_fP_sq_"+obsinfoR+".bin");
//       write(WR, "fP_sq", meas.fPs_sq );
//     }
//   }



  // ################################################################ //
  std::cout << GridLogMessage << "run histogram part" << std::endl;

  {

    std::vector<std::vector<std::vector<EnsembleOfObs>>> array_of_hist_ens( nbins_histo_x,
                                                                            std::vector<std::vector<EnsembleOfObs>>(nbins_histo_y,
                                                                                                                    std::vector<EnsembleOfObs>(nbeta)
                                                                                                                    )
                                                                            );
    char id[200];



    { //make hist
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


      // #ifdef _OPENMP
      // #pragma omp parallel for
      // #endif
      std::vector<EnsembleOfObs> ensemblesR( nbeta );
      std::vector<EnsembleOfObs> ensemblesI( nbeta );

      // obsinfo="flow_polyakov_re";
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(int j=0; j<nbeta; j++){
        // std::cout << GridLogMessage << "j = " << j << std::endl;
        const std::string beta = betas[j];
        std::string id_str = get_configname(beta, mass, type);
        std::string path = basedir+id_str+".bin";
        // std::cout << GridLogMessage << "path = " << path << std::endl;

#pragma omp critical
        {
          std::unique_ptr<Hdf5Reader> WR;
          WR = std::make_unique<Hdf5Reader>( path );
          read(*WR, "ensR", ensemblesR[j] );
          read(*WR, "ensI", ensemblesI[j] );
        }
        // BinaryReader WR( path );
        // read(WR, "ens", ensemblesR[j] );
      }

//       // obsinfo="flow_polyakov_im";
// #ifdef _OPENMP
// #pragma omp parallel for
// #endif
//       for(int j=0; j<nbeta; j++){
//         std::cout << GridLogMessage << "j = " << j << std::endl;
//         const std::string beta = betas[j];
//         std::string id_str = get_configname(beta, mass, type);
//         std::string path = basedir+id_str+obsinfoI+".bin";
//         std::cout << GridLogMessage << "path = " << path << std::endl;

//         BinaryReader WR( path );
//         read(WR, "ens", ensemblesI[j] );
//       }

      std::sprintf(id, "m%s", mass.data());


      for(int j=0; j<nbeta; j++){
        for(int ibx=0; ibx<nbins_histo_x; ibx++){ // histogram binning for
#ifdef _OPENMP
#pragma omp parallel for
#endif
          for(int iby=0; iby<nbins_histo_y; iby++){ // histogram binning for
            EnsembleOfObs& ibth_bin = array_of_hist_ens[ibx][iby][j];
            ibth_bin = ensemblesR[j];

            ibth_bin.Os.clear();
            //for(const RealD O : ensembles0[j].Os){
            for(int i=0; i<ensemblesR[j].Os.size(); i++){
              const RealD Ox = ensemblesR[j].Os[i];
              const RealD Oy = ensemblesI[j].Os[i];

              const RealD v = regularized_delta( Ox, Oy, xs[ibx], ys[iby], delta_x, delta_y );
              ibth_bin.Os.push_back( v );
            } // end for i
          } // end for ibx
        } // end for iby

        // std::cout << "histogram = " << std::endl;
        // for(int ib=0; ib<nbins_histo; ib++) {
        //         for(int ibx=0; ibx<nbins_histo_x; ibx++){ // histogram binning for
// #ifdef _OPENMP
// #pragma omp parallel for
// #endif
//           for(int iby=0; iby<nbins_histo_y; iby++){ // histogram binning for
//             double sum = 0.0;
//             EnsembleOfObs& ibth_bin = array_of_hist_ens[ibx][iby][j];
//             for(auto elem : ibth_bin.Os) sum += elem;
//             // std::cout << ibx << " " << iby << " " << sum << std::endl;
//           }
//         }
        // std::cout << std::endl;

      } // end for j
    } // close make hist

    // -------------------
    FreeEnergies fs(nbeta);

    {
      std::string id_str = id;
      id_str = id_str+"_nojk";

#pragma omp critical
      {
        std::unique_ptr<Hdf5Reader> WR;
        WR = std::make_unique<Hdf5Reader>( basedir+id_str+"fs.bin" );
        read(*WR, "f", fs.f );
      }
      // @@ BinaryReader WR(basedir+id_str+"fs.bin");
      // read(WR, "f", fs.f );
    }

    // =================================

// #ifdef _OPENMP
// #pragma omp parallel for
// #endif
    // for(int ib=0; ib<nbins_histo; ib++){
    for(int ibx=0; ibx<nbins_histo_x; ibx++){ // histogram binning for
#pragma omp parallel for num_threads(omp_get_max_threads()/20)
      for(int iby=0; iby<nbins_histo_y; iby++){ // histogram binning for

        ReweightedObs meas( betas_double.front(), betas_double.back(), nbeta*100, array_of_hist_ens[ibx][iby] );
        const int nmeas=meas.nmeas;
        // run_reweighting( betas_double, ensemblesR, fs, meas, false, 20 );
        run_reweighting( betas_double, fs, meas, false, 20 );

        std::string id_str = id;
        id_str = id_str+"hist_ibx"+std::to_string(ibx)+"_iby"+std::to_string(iby)+"_nojk";

        // if (UGrid->IsBoss()){
#pragma omp critical
        {
          std::unique_ptr<Hdf5Writer> WR;
          WR = std::make_unique<Hdf5Writer>( basedir+id_str+"meas.bin" );
          // write(WR, "ensR", ensR );
          // write(WR, "ensI", ensI );
          write(*WR, "beta", meas.betas );
          write(*WR, "f", meas.fs );
          write(*WR, "fP", meas.fPs );
        }
        // {
        //   BinaryWriter WR(basedir+id_str+"meas_xs.bin");
        //   write(WR, "beta", meas.betas );
        // }
        // {
        //   BinaryWriter WR(basedir+id_str+"meas_f.bin");
        //   write(WR, "f", meas.fs );
        // }
        // {
        //   std::cout << GridLogMessage << "path_fp = " << basedir+id_str+"meas_fP_2d.bin" << std::endl;
        //   BinaryWriter WR(basedir+id_str+"meas_fP_2d.bin");
        //   write(WR, "fP", meas.fPs );
        // }
      } // end for ib
    }
  }


  Grid_finalize();
}
