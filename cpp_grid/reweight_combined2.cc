#include "reweight2.h"

#include <regex>


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

  std::string obsinfoR="flow_polyakov_re";
  std::string obsinfoI="flow_polyakov_im";

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

  // ################################################################ //
  std::cout << GridLogMessage << "run1 part" << std::endl;
  // std::cout << GridLogMessage << "interval = " << interval << std::endl;
  std::cout << GridLogMessage << "conf_min0 = " << conf_min0 << std::endl;
  std::cout << GridLogMessage << "conf_max0 = " << conf_max0 << std::endl;

// #ifdef _OPENMP
// #pragma omp parallel for
// #endif
  // for(std::string beta : betas){
  //   std::cout << GridLogMessage << "beta = " << beta << std::endl;
  //   EnsembleOfObs ens;
  //   ens.beta = beta;
  //   // ens.conf_min = conf_min;
  //   // ens.interval = interval;
  //   // ens.binsize = binsize;

  //   ens.base_dir=base_dir;
  //   ens.id=get_configname(beta, mass, type);

  //   int conf_min=conf_min0, conf_max; // , conf_tot=0; // , nbin, nbin_tot=0;
  //   // int conf_max; // , conf_tot=0, nbin, nbin_tot=0;

  //   // ---------------------------
  //   // cold start
  //   // ---------------------------
  // for(std::string beta : betas){
  //   std::cout << GridLogMessage << "beta = " << beta << std::endl;
  //   EnsembleOfObs ens;
  //   ens.beta = beta;
  //   ens.base_dir=base_dir;
  //   ens.id=get_configname(beta, mass, type);

  //   // ---------------------------
  //   // cold start
  //   // ---------------------------

  //   for(int conf=conf_min0; conf<conf_max0; conf++){
  //     std::string pathO = basedir2c+ens.id+obsinfo+std::to_string(conf)+".bin";
  //     if( !std::filesystem::exists( pathO ) ) continue; // assert(false);

  //     Obs obs;
  //     BinaryReader WR( pathO );
  //     read(WR, "obs", obs );

  //     ens.Os.push_back(obs.O);
  //     ens.Osqs.push_back(obs.O*obs.O);
  //     ens.energies.push_back(obs.energy);
  //   } // end for conf

  //   // ---------------------------
  //   // hot start
  //   // ---------------------------

  //   if(type!=1) ens.id=get_configname(beta, mass, type+1);

  //   for(int conf=conf_min0; conf<conf_max0; conf++){
  //     std::string pathO = basedir2h+ens.id+obsinfo+std::to_string(conf)+".bin";
  //     if( !std::filesystem::exists( pathO ) ) continue; // assert(false);

  //     Obs obs;
  //     BinaryReader WR( pathO );
  //     read(WR, "obs", obs );

  //     ens.Os.push_back(obs.O);
  //     ens.Osqs.push_back(obs.O*obs.O);
  //     ens.energies.push_back(obs.energy);
  //   } // end for conf

  //   ens.nbin=nbin_global;
  //   ens.binsize = ens.size()/nbin_global;
  //   ens.nconf = ens.nbin*ens.binsize;

  //   { // trimming
  //     const int start = ens.nconf;
  //     ens.energies.erase( ens.energies.begin()+start, ens.energies.end() );
  //     ens.Os.erase( ens.Os.begin()+start, ens.Os.end() );
  //     ens.Osqs.erase( ens.Osqs.begin()+start, ens.Osqs.end() );
  //   }

  //   // std::cout << GridLogMessage << "beta = " << beta << " conf_tot = " << conf_tot << std::endl;

  //   std::string path_ens = basedir+ens.id+obsinfo+".bin";
  //   // std::cout << GridLogMessage << "path_ens = " << path_ens << std::endl;
  //   BinaryWriter WR( path_ens );
  //   write(WR, "ens", ens );
  // } // for end beta

  // ################################################################ //
  std::cout << GridLogMessage << "run2 part" << std::endl;

  if(type!=1) type+=1;

  {
    std::vector<EnsembleOfObs> ensembles0( nbeta );

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int j=0; j<nbeta; j++){
      std::cout << GridLogMessage << "j = " << j << std::endl;
      const std::string beta = betas[j];
      std::string id_str = get_configname(beta, mass, type);
      std::string path = basedir+id_str+obsinfoR+".bin";
      std::cout << GridLogMessage << "path = " << path << std::endl;

      BinaryReader WR( path );
      read(WR, "ens", ensembles0[j] );
    }

    // std::cout << GridLogMessage << "run2 part" << std::endl;



    // std::cout << GridLogMessage << "debug. pt1" << std::endl;

// #ifdef _OPENMP
// #pragma omp parallel for
// #endif
    for(int jdrop=0; jdrop<nbeta; jdrop++){
      // std::cout << GridLogMessage << "debug. jdrop = " << jdrop << std::endl;
      const int nbin=ensembles0[jdrop].nbin;
      // std::cout << GridLogMessage << "debug. nbin = " << nbin << std::endl;
      for(int ibin=0; ibin<nbin; ibin++){
        // std::cout << GridLogMessage << "debug. ibin = " << ibin << std::endl;

        std::vector<EnsembleOfObs> ensembles(nbeta);
        for(int i=0; i<nbeta; i++) ensembles[i] = ensembles0[i];

        // std::cout << GridLogMessage << "debug. pt1" << std::endl;
        ensembles[jdrop].jackknife_drop( nbin, ibin );
        // std::cout << GridLogMessage << "debug. pt2" << std::endl;

        char id[200];
        std::sprintf(id, "m%s", mass.data());
        std::string id_str = id;
        id_str = id_str+"_jk_"+std::to_string(jdrop)+"_"+std::to_string(nbin)+"_"+std::to_string(ibin);

        // -------------------
        FreeEnergies fs(nbeta);
        // std::cout << GridLogMessage << "debug. pt3" << std::endl;

        {
          std::string id_str2 = id;
          id_str2 = id_str2+"_nojk";
          std::string filename2=basedir+id_str2+"fs.bin";
          std::string filename0=basedir+id_str+"fs.bin";
          if(runtype==0){ // rerunning after termination
            if( std::filesystem::exists( filename0 ) ){
              BinaryReader WR(filename0);
              read(WR, "f", fs.f );
              fs.fnew = fs.f;
              fs.fold = fs.f;
            }
            else if( std::filesystem::exists( filename2 ) ){
              BinaryReader WR(filename2);
              read(WR, "f", fs.f );
              fs.fnew = fs.f;
              fs.fold = fs.f;
            }
          }
          else if(runtype==1){ // recalculating from nojk
            if( std::filesystem::exists( filename2 ) ){
              BinaryReader WR(filename2);
              read(WR, "f", fs.f );
              fs.fnew = fs.f;
              fs.fold = fs.f;
            }
            // else if( std::filesystem::exists( filename0 ) ){
            //   BinaryReader WR(filename0);
            //   read(WR, "f", fs.f );
            //   fs.fnew = fs.f;
            //   fs.fold = fs.f;
            // }
          }
        }

        // std::cout << GridLogMessage << "debug. calcf" << std::endl;

        std::vector<std::vector<RealD>> v_energies(nbeta);
        for(int k=0; k<nbeta; k++) v_energies[k] = ensembles[k].energies;
        run_multihistogram( betas_double, v_energies, fs, basedir+id_str+"fs.bin" );

        // std::cout << GridLogMessage << "debug. write " << std::endl;
        // std::cout << GridLogMessage << "pathf = " << basedir+id_str+"fs.bin" << std::endl;

        BinaryWriter WR(basedir+id_str+"fs.bin");
        write(WR, "f", fs.f );
      }
    }
  }

//   // ################################################################ //
//   std::cout << GridLogMessage << "run3 part" << std::endl;

//   {
//     std::vector<EnsembleOfObs> ensembles0( nbeta );

// #ifdef _OPENMP
// #pragma omp parallel for
// #endif
//     for(int j=0; j<nbeta; j++){
//       const std::string beta = betas[j];
//       std::string id_str = get_configname(beta, mass, type);
//       std::string path = basedir+id_str+obsinfo+".bin";

//       BinaryReader WR( path );
//       read(WR, "ens", ensembles0[j] );
//     }

//     char id[200];
//     std::sprintf(id, "m%s", mass.data());

// // #ifdef _OPENMP
// // #pragma omp parallel for
// // #endif
//     for(int jdrop=0; jdrop<nbeta; jdrop++){
//       const int nbin=ensembles0[jdrop].nbin;
//       for(int ibin=0; ibin<nbin; ibin++){
//         std::string id_str = id;
//         id_str = id_str+"_jk_"+std::to_string(jdrop)+"_"+std::to_string(nbin)+"_"+std::to_string(ibin);

//         if(runtype==0 && std::filesystem::exists( basedir+id_str+"meas_fP_"+obsinfo+".bin" ) ) continue;

//         std::vector<EnsembleOfObs> ensembles(nbeta);
//         for(int i=0; i<nbeta; i++) ensembles[i] = ensembles0[i];
//         ensembles[jdrop].jackknife_drop( nbin, ibin );

//         // -------------------

//         FreeEnergies fs(nbeta);

//         BinaryReader WR(basedir+id_str+"fs.bin");
//         read(WR, "f", fs.f );

//         // =================================

//         ReweightedObs meas( betas_double.front(), betas_double.back(), nbeta*100, ensembles );

//         // {
//         //   std::string filename2=basedir+id_str+"meas_f.bin";
//         //   if( std::filesystem::exists( filename2 ) ){
//         //     BinaryReader WR(filename);
//         //     read(WR, "f", meas.fs );
//         //   }
//         // }
//         // {
//         //   std::string filename2=basedir+id_str+"meas_fP_"+obsinfo+".bin";
//         //   if( std::filesystem::exists( filename2 ) ){
//         //     BinaryReader WR(filename);
//         //     read(WR, "fP", meas.fPs );
//         //   }
//         // }
//         // {
//         //   std::string filename2=basedir+id_str+"meas_fP_sq_"+obsinfo+".bin";
//         //   if( std::filesystem::exists( filename2 ) ){
//         //     BinaryReader WR(filename);
//         //     read(WR, "fP_sq", meas.fPs_sq );
//         //   }
//         // }

//         const int nmeas=meas.nmeas;
//         run_reweighting( betas_double, ensembles, fs, meas );

//         {
//           // std::cout << GridLogMessage << "path_beta = " << basedir+id_str+"meas_beta.bin" << std::endl;
//           BinaryWriter WR(basedir+id_str+"meas_beta.bin");
//           write(WR, "beta", meas.betas );
//         }
//         {
//           BinaryWriter WR(basedir+id_str+"meas_f.bin");
//           write(WR, "f", meas.fs );
//         }
//         {
//           BinaryWriter WR(basedir+id_str+"meas_fP_"+obsinfo+".bin");
//           write(WR, "fP", meas.fPs );
//         }
//         {
//           BinaryWriter WR(basedir+id_str+"meas_sign_P_"+obsinfo+".bin");
//           write(WR, "sign_P", meas.sign_Ps );
//         }
//         {
//           BinaryWriter WR(basedir+id_str+"meas_fP_sq_"+obsinfo+".bin");
//           write(WR, "fP_sq", meas.fPs_sq );
//         }
//       }
//     }
//   }

  // ################################################################ //
  std::cout << GridLogMessage << "run histogram part" << std::endl;

  {
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


    std::vector<EnsembleOfObs> ensemblesR( nbeta );

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int j=0; j<nbeta; j++){
      const std::string beta = betas[j];
      std::string id_str = get_configname(beta, mass, type);
      std::string path = basedir+id_str+obsinfoR+".bin";
      std::cout << GridLogMessage << "path = " << path << std::endl;

      BinaryReader WR( path );
      read(WR, "ens", ensemblesR[j] );
    }

    std::vector<EnsembleOfObs> ensemblesI( nbeta );


#ifdef _OPENMP
#pragma omp parallel for
#endif
    for(int j=0; j<nbeta; j++){
      std::cout << GridLogMessage << "j = " << j << std::endl;
      const std::string beta = betas[j];
      std::string id_str = get_configname(beta, mass, type);
      std::string path = basedir+id_str+obsinfoI+".bin";
      std::cout << GridLogMessage << "path = " << path << std::endl;

      BinaryReader WR( path );
      read(WR, "ens", ensemblesI[j] );
    }


    char id[200];
    std::sprintf(id, "m%s", mass.data());




    for(int jdrop=0; jdrop<nbeta; jdrop++){
      const int nbin=ensemblesR[jdrop].nbin;
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for(int ibin=0; ibin<nbin; ibin++){ // jackknife for
        std::vector<std::vector<std::vector<EnsembleOfObs>>> array_of_hist_ens( nbins_histo_x,
                                                                                std::vector<std::vector<EnsembleOfObs>>(nbins_histo_y,
                                                                                                                        std::vector<EnsembleOfObs>(nbeta)
                                                                                                                        )
                                                                                );

        std::vector<EnsembleOfObs> ensemblesR1(nbeta);
        std::vector<EnsembleOfObs> ensemblesI1(nbeta);
        for(int i=0; i<nbeta; i++) {
          ensemblesR1[i] = ensemblesR[i];
          ensemblesI1[i] = ensemblesI[i];
        }
        ensemblesR1[jdrop].jackknife_drop( nbin, ibin );
        ensemblesI1[jdrop].jackknife_drop( nbin, ibin );

        for(int j=0; j<nbeta; j++){
          // for(int ib=0; ib<nbins_histo; ib++){ // histogram binning for
          for(int ibx=0; ibx<nbins_histo_x; ibx++){ // histogram binning for
            for(int iby=0; iby<nbins_histo_y; iby++){ // histogram binning for

              EnsembleOfObs& ibth_bin = array_of_hist_ens[ibx][iby][j];
              ibth_bin = ensemblesR1[j];

              ibth_bin.Os.clear();
              for(int i=0; i<ensemblesR1[j].Os.size(); i++){
                const RealD Ox = ensemblesR1[j].Os[i];
                const RealD Oy = ensemblesI1[j].Os[i];

                const RealD v = regularized_delta( Ox, Oy, xs[ibx], ys[iby], delta_x, delta_y );
                ibth_bin.Os.push_back( v );
              } // end for i
            }} // end for ib

          std::cout << "histogram = " << std::endl;
          // for(int ib=0; ib<nbins_histo; ib++) {
// #pragma omp parallel for collapse(2) num_threads(omp_get_max_threads()/10)
          for(int ibx=0; ibx<nbins_histo_x; ibx++){ // histogram binning for
            for(int iby=0; iby<nbins_histo_y; iby++){ // histogram binning for
              double sum = 0.0;
              EnsembleOfObs& ibth_bin = array_of_hist_ens[ibx][iby][j];
              for(auto elem : ibth_bin.Os) sum += elem;
              std::cout << ibx << " " << iby << " " << sum << std::endl;
            }
          }
          std::cout << std::endl;

        } // end for j

        // -------------------
        FreeEnergies fs(nbeta);

        {
          std::string id_str = id;
          id_str = id_str+"_jk_"+std::to_string(jdrop)+"_"+std::to_string(nbin)+"_"+std::to_string(ibin);
          BinaryReader WR(basedir+id_str+"fs.bin");
          read(WR, "f", fs.f );
        }

        // =================================

// #ifdef _OPENMP
// #pragma omp parallel for num_threads(omp_get_max_threads()/10)
// #endif
        // for(int ib=0; ib<nbins_histo; ib++){
        for(int ibx=0; ibx<nbins_histo_x; ibx++){ // histogram binning for
          for(int iby=0; iby<nbins_histo_y; iby++){ // histogram binning for

            std::string id_str = id;
            id_str = id_str+"hist_ibx"+std::to_string(ibx)+"_iby"+std::to_string(iby)+"_jk_"+std::to_string(jdrop)+"_"+std::to_string(nbin)+"_"+std::to_string(ibin);
            if(runtype==0 && std::filesystem::exists( basedir+id_str+"meas_fP_"+obsinfo+".bin" ) ) continue;


            // ReweightedObs meas( betas_double.front(), betas_double.back(), nbeta*100, array_of_hist_ens[ib] );
            ReweightedObs meas( betas_double.front(), betas_double.back(), nbeta*100, array_of_hist_ens[ibx][iby] );
            const int nmeas=meas.nmeas;
            // run_reweighting( betas_double, ensemblesR, fs, meas, false );
            run_reweighting( betas_double, fs, meas, false );


            {
              BinaryWriter WR(basedir+id_str+"meas_xs.bin");
              write(WR, "beta", meas.betas );
            }
            {
              BinaryWriter WR(basedir+id_str+"meas_f.bin");
              write(WR, "f", meas.fs );
            }
            {
              // std::cout << GridLogMessage << "path_fp = " << basedir+id_str+"meas_fP_"+obsinfo+".bin" << std::endl;
              BinaryWriter WR(basedir+id_str+"meas_fP_"+obsinfo+".bin");
              write(WR, "fP", meas.fPs );
            }
          }} // end for ib
      }} // end for jackknife
  }


  Grid_finalize();
}
