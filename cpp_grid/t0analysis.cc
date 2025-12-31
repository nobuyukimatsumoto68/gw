#include "reweight2.h"
#include "jk.h"

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
  // const int nbins_histo_x=atoi(argv[11]);
  // double min_x = atof(argv[12]);
  // double max_x = atof(argv[13]);
  // const int nbins_histo_y=atoi(argv[14]);
  // double min_y = atof(argv[15]);
  // double max_y = atof(argv[16]);

  std::cout << "nbeta = " << nbeta << std::endl;
  std::cout << "basedir = " << basedir << std::endl;
  // std::cout << "nbins_histo_x = " << nbins_histo_x << std::endl;
  // std::cout << "min_x = " << min_x << std::endl;
  // std::cout << "max_x = " << max_x << std::endl;
  // std::cout << "nbins_histo_y = " << nbins_histo_y << std::endl;
  // std::cout << "min_y = " << min_y << std::endl;
  // std::cout << "max_y = " << max_y << std::endl;

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

  using Impl = PeriodicGimplR;

  // ################################################################ //
  std::cout << GridLogMessage << "run1 part" << std::endl;
  // std::cout << GridLogMessage << "interval = " << interval << std::endl;
  std::cout << GridLogMessage << "conf_min0 = " << conf_min0 << std::endl;
  std::cout << GridLogMessage << "conf_max0 = " << conf_max0 << std::endl;


#ifdef _OPENMP
#pragma omp parallel for
#endif
  for(std::string beta : betas){
    std::string obsinfo = "flow_T0";

    std::string id_str = get_configname(beta, mass, type);
    std::cout << GridLogMessage << "beta = " << beta << std::endl;

    std::string path1 = basedir+id_str+"poly_Q.txt";
    std::ofstream file1(path1, std::ios::trunc);
    file1 << "# poly Q " << std::endl;
    // for(int j=0; j<nbeta; j++){

    Corr sum;
    int counter=0;

    //////
    // TO DO: Skipping

    // ---------------------------
    // cold start
    // ---------------------------

    // predef, reused
    VectorObs vobs;
    vobs.beta = beta;
    vobs.base_dir=base_dir;
    vobs.id=get_configname(beta, mass, type);

    for(int conf=conf_min0; conf<conf_max0; conf++){
      {
        obsinfo = "flow_T0";
        std::string pathO = basedir2c+vobs.id+obsinfo+std::to_string(conf)+".bin";
        if( !std::filesystem::exists( pathO ) ) continue; // assert(false);
        BinaryReader WR( pathO );
        read(WR, "vobs", vobs );
      }
      // std::cout << "read: pathO = " << pathO << std::endl;

      sum += Corr(vobs.vO);
      counter++;

      std::vector<std::string> obsinfolist={"flow_polyakov_re", "flow_polyakov_im", "flow_Q"};
      for(std::string elem : obsinfolist){
        std::string pathO = basedir2c+vobs.id+elem+std::to_string(conf)+".bin";
        if( !std::filesystem::exists( pathO ) ) continue;
      }

      {
        obsinfo = "flow_polyakov_re";
        std::string pathO = basedir2c+vobs.id+obsinfo+std::to_string(conf)+".bin";
        if( !std::filesystem::exists( pathO ) ) assert(false);
        BinaryReader WR( pathO );
        read(WR, "vobs", vobs );

        file1 << vobs.vO[200] << " ";
      }
      {
        obsinfo = "flow_polyakov_im";
        std::string pathO = basedir2c+vobs.id+obsinfo+std::to_string(conf)+".bin";
        if( !std::filesystem::exists( pathO ) ) assert(false);
        BinaryReader WR( pathO );
        read(WR, "vobs", vobs );

        file1 << vobs.vO[200] << " ";
      }
      {
        obsinfo = "flow_Q";
        std::string pathO = basedir2c+vobs.id+obsinfo+std::to_string(conf)+".bin";
        if( !std::filesystem::exists( pathO ) ) assert(false);
        BinaryReader WR( pathO );
        read(WR, "vobs", vobs );

        file1 << vobs.vO[200] << std::endl;
      }

    } // end for conf


    sum /= counter;


    // char id[200];
    // std::sprintf(id, "m%s", mass.data());
    // std::string id_str = id;
    // std::string id_str = get_configname(beta, mass, type);
    std::string path = basedir+id_str+obsinfo+".txt";
    std::ofstream file(path, std::ios::trunc);
    file << "# t T0 avg. for beta = " << beta << std::endl;
    // for(int j=0; j<nbeta; j++){
    counter = 0;
    double dt = 0.01;
    for(double elem : sum.data) {
      file << dt*counter << " " << elem << std::endl;
      counter++;
    }
    // }


  } // for end beta




  Grid_finalize();
}
