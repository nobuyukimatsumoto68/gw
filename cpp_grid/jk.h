struct Corr {
  std::vector<double> data;

  Corr(){}

  Corr( const std::vector<double>& data2 ){
    data.clear();
    data = data2;
  }

  Corr( const std::string& file ){
    data.clear();
    readfile( file );
  }

  std::size_t size() const { return data.size(); }

  void clear() {
    data.clear();
  }

  void readfile( const std::string& file ){
    std::string str;
    std::ifstream ifs(file);
    if(!ifs) assert(false);

    while (std::getline(ifs, str)){
      std::istringstream iss(str);
      double v1;
      iss >> v1;
      data.push_back( v1 );
    }
  }

  double operator[](const int i) const { return data[i]; }
  double& operator[](const int i) { return data[i]; }

  Corr& operator+=(const Corr& rhs){
    if( data.size()==0 ) {
      data = rhs.data;
    }
    else{
      // std::cout << "# debug. " << data.size() << " " << rhs.data.size() << std::endl;
      assert( data.size()==rhs.data.size() );
// #ifdef _OPENMP
// #pragma omp parallel for
// #endif
      for(int i=0; i<data.size(); i++) data[i] += rhs.data[i];
    }
    return *this;
  }

  Corr& operator-=(const Corr& rhs){
    if( data.size()==0 ) {
      data = rhs.data;
    }
    else{
      assert( data.size()==rhs.data.size() );
// #ifdef _OPENMP
// #pragma omp parallel for
// #endif
      for(int i=0; i<data.size(); i++) data[i] -= rhs.data[i];
    }
    return *this;
  }

  friend Corr operator-(Corr v, const Corr& w) { v -= w; return v; }

  Corr& operator*=(const Corr& rhs){
    assert( data.size()==rhs.data.size() );
// #ifdef _OPENMP
// #pragma omp parallel for
// #endif
    for(int i=0; i<data.size(); i++) data[i] *= rhs[i];
    return *this;
  }

  friend Corr operator*(Corr v, const Corr& w) { v *= w; return v; }

  Corr& operator*(const double& rhs){
// #ifdef _OPENMP
// #pragma omp parallel for
// #endif
    for(int i=0; i<data.size(); i++) data[i] *= rhs;
    return *this;
  }

  Corr& operator*=(const double& rhs){
// #ifdef _OPENMP
// #pragma omp parallel for
// #endif
    for(int i=0; i<data.size(); i++) data[i] *= rhs;
    return *this;
  }

  Corr& operator/=(const double& rhs){
// #ifdef _OPENMP
// #pragma omp parallel for
// #endif
    for(int i=0; i<data.size(); i++) data[i] /= rhs;
    return *this;
  }
};



template<typename T>
struct Jackknife {
  std::vector<T> binavgs;
  std::vector<T> jackavgs;
  T est;
  T var;

  int nbins;
  const int binsize;

  Jackknife( const int binsize )
    : binsize(binsize)
  {}

  void do_it(){
    nbins = binavgs.size();

    jackavgs.clear();
    jackavgs.resize( nbins );
    for(int i=0; i<nbins; i++){
      for(int j=0; j<nbins; j++){
        if(i==j) continue;
        jackavgs[i] += binavgs[j];
      }
      jackavgs[i] /= nbins-1.0;
    }

    est.clear();
    for(int i=0; i<nbins; i++){
      est += jackavgs[i];
    }
    est /= nbins;

    var.clear();
    for(int i=0; i<nbins; i++){
      T diff = jackavgs[i] - est;
      var += diff * diff;
    }
    var /= nbins;
    var *= nbins-1.0;
  }

};
