#include <iostream>

#include <Rcpp.h>
using namespace Rcpp;



//[[Rcpp::plugins(cpp11)]]

// Code shared with IsoSpec
//
//

#include <cmath>
#include <vector>
#include <cstdlib>
#include <limits>
#include <cstring>
#include <unordered_set>
#include <algorithm>

typedef int* Conf;


double* getLFactorials(int howmany)
{
  howmany++;
  double* ret = new double[howmany];
  for(int ii=0; ii<howmany; ii++)
    ret[ii] = lgamma(ii+1);
  return ret;
}


class Marginal
{
private:
  bool disowned;
  void setupInitialConf(const double* probs);
protected:
  const unsigned int isotopeNo;
  const unsigned int atomCnt;
  const double* atom_lProbs;
  const Conf mode_conf;
  double mode_lprob;
  const double* lfactorials;
  const double nom_lfact;
  
public:
  Marginal(
    const double* _probs,
    int _isotopeNo,                  // No of isotope configurations.
    int _atomCnt
  );
  Marginal(Marginal& other) = delete;
  Marginal& operator= (const Marginal& other) = delete;
  Marginal(Marginal&& other);
  virtual ~Marginal();
  
  inline int get_isotopeNo() const { return isotopeNo; };
  inline double getModeLProb() const { return mode_lprob; };
  inline double logProb(const int* conf)
  {
    double  res = 0.0;
    
    for(unsigned int i=0; i < isotopeNo; i++)
    {
      res -= lfactorials[conf[i]];
      res += conf[i] * atom_lProbs[i];
    }
    return res + nom_lfact;
  }
  
};



void Marginal::setupInitialConf(const double* probs)
{
  Conf res = mode_conf;
  
  for(unsigned int i = 0; i < isotopeNo; ++i )
  {
    res[i] = int( atomCnt * probs[i] ) + 1;
  }
  
  int s = 0;
  
  for(unsigned int i = 0; i < isotopeNo; ++i) s += res[i];
  
  int diff = atomCnt - s;
  
  // Too little: enlarging fist index.
  if( diff > 0 ){
    res[0] += diff;
  }
  // Too much: trying to redistribute the difference: hopefully the first element is the largest.
  if( diff < 0 ){
    diff = abs(diff);
    int i = 0, coordDiff = 0;
    
    while( diff > 0){
      coordDiff = res[i] - diff;
      
      if( coordDiff >= 0 ){
        res[i] -= diff;
        diff = 0;
      } else {
        res[i] = 0;
        i++;
        diff = abs(coordDiff);
      }
    }
  }
  
  // What we computed so far will be very close to the mode: hillclimb the rest of the way
  
  bool modified = true;
  double LP = logProb(res);
  double NLP;
  
  while(modified)
  {
    modified = false;
    for(unsigned int ii = 0; ii<isotopeNo; ii++)
      for(unsigned int jj = 0; jj<isotopeNo; jj++)
        if(ii != jj and res[ii] > 0)
        {
          res[ii]--;
          res[jj]++;
          NLP = logProb(res);
          if(NLP>LP or (NLP==LP and ii>jj))
          {
            modified = true;
            LP = NLP;
          }
          else
          {
            res[ii]++;
            res[jj]--;
          }
        }
        
        
  }
}



template <typename T> inline static T* array_copy(const T* A, int size)
{
  T* ret = new T[size];
  memcpy(ret, A, size*sizeof(T));
  return ret;
}

double* getMLogProbs(const double* probs, int isotopeNo)
{
  double* ret = new double[isotopeNo];
  for(int ii=0; ii<isotopeNo; ii++)
    ret[ii] = log(probs[ii]);
  return ret;
}

template <typename T> inline void copyConf(
    const T* source, T* destination,
    int dim
){
  memcpy(destination, source, dim*sizeof(T));
}

template <typename T> class Allocator{
private:
  T*      currentTab;
  int currentId;
  const int       dim, tabSize;
  std::vector<T*>  prevTabs;
public:
  Allocator(const int dim, const int tabSize = 10000);
  ~Allocator();
  
  void shiftTables();
  
  inline T* newConf()
  {
    currentId++;
    
    if (currentId >= tabSize)
      shiftTables();
    
    return &(currentTab[ currentId * dim ]);
  }
  
  inline T* makeCopy(const T* conf)
  {
    T* currentPlace = newConf();
    copyConf<T>( conf, currentPlace, dim );
    
    return currentPlace;
  }
  
  inline T* makeExternalCopy(const T* conf)
  {
    T* res = new T[dim];
    copyConf( conf, res, dim );
    
    return res;
  }
};

template <typename T>
Allocator<T>::Allocator(const int dim, const int tabSize): currentId(-1), dim(dim), tabSize(tabSize)
{
  currentTab = new T[dim * tabSize];
}

template <typename T>
Allocator<T>::~Allocator()
{
  for(unsigned int i = 0; i < prevTabs.size(); ++i)
  {
    delete [] prevTabs[i];
  }
  
  delete [] currentTab;
}

template <typename T>
void Allocator<T>::shiftTables()
{
  prevTabs.push_back(currentTab);
  currentTab      = new T[dim * tabSize];
  currentId       = 0;
}
template class Allocator<int>;


class KeyHasher
{
private:
  int dim;
public:
  KeyHasher(int dim);
  
  inline std::size_t operator()(const int* conf) const
  {
    // Following Boost...
    std::size_t seed = 0;
    for(int i = 0; i < dim; ++i )
      seed ^= conf[i] + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    return seed;
  };
};


class ConfEqual
{
private:
  int size;
public:
  ConfEqual(int dim);
  
  inline bool operator()(const int* conf1, const int* conf2) const
  {
    return not memcmp(conf1, conf2, size);
  }
};

class ConfOrderMarginalDescending
{
  //configurations comparator
  Marginal* marginal;
public:
  ConfOrderMarginalDescending(Marginal* marginal);
  
  inline bool operator()(const Conf conf1, const Conf conf2)
  {// Return true if conf1 is less probable than conf2.
    return marginal->logProb(conf1) > marginal->logProb(conf2);
  };
};


KeyHasher::KeyHasher(int dim)
  : dim(dim)
{}

ConfEqual::ConfEqual(int dim)
  : size( dim*sizeof(int) )
{}

ConfOrderMarginalDescending::ConfOrderMarginalDescending(Marginal* _marginal)
  : marginal(_marginal)
{}



Marginal::Marginal(
  const double* _probs,
  int _isotopeNo,
  int _atomCnt
) :
  disowned(false),
  isotopeNo(_isotopeNo),
  atomCnt(_atomCnt),
  atom_lProbs(getMLogProbs(_probs, isotopeNo)),
  mode_conf(new int[_isotopeNo]),
  lfactorials(getLFactorials(_atomCnt)),
  nom_lfact(lfactorials[_atomCnt])
{
  setupInitialConf(_probs);
  mode_lprob = logProb(mode_conf);
}

Marginal::Marginal(Marginal&& other) :
  disowned(other.disowned),
  isotopeNo(other.isotopeNo),
  atomCnt(other.atomCnt),
  atom_lProbs(other.atom_lProbs),
  mode_conf(other.mode_conf),
  mode_lprob(other.mode_lprob),
  lfactorials(other.lfactorials),
  nom_lfact(other.nom_lfact)
{
  other.disowned = true;
}

Marginal::~Marginal()
{
  if(not disowned)
  {
    delete[] atom_lProbs;
    delete[] mode_conf;
    delete[] lfactorials;
  }
}



class LayeredMarginal : public Marginal
{
  double current_threshold, new_threshold;
  std::vector<Conf> configurations;
  std::vector<Conf> fringe;
  Allocator<int> allocator;
  const ConfEqual equalizer;
  const KeyHasher keyHasher;
  const ConfOrderMarginalDescending orderMarginal;
  const int hashSize;
  std::unordered_set<Conf,KeyHasher,ConfEqual> visited;
  double opc;
  Conf currentConf;
  std::vector<Conf> new_fringe;
  
  
  
public:
  LayeredMarginal(Marginal&& m, int tabSize = 1000, int hashSize = 1000);
  bool extend(double new_threshold);
  bool next();
  inline double get_lProb() const { return opc; };
  inline const Conf& get_conf() const { return currentConf; };
  bool empty() { return fringe.empty() and new_fringe.empty(); };
  
};





LayeredMarginal::LayeredMarginal(Marginal&& m, int tabSize, int _hashSize)
  : Marginal(std::move(m)), new_threshold(1.0), allocator(isotopeNo, tabSize),
    equalizer(isotopeNo), keyHasher(isotopeNo), orderMarginal(this), hashSize(_hashSize),
    visited(hashSize,keyHasher,equalizer)
{
  new_fringe.push_back(mode_conf);
}

bool LayeredMarginal::extend(double _new_threshold)
{
  if(new_fringe.empty())
    return false;
  
  current_threshold = new_threshold;
  new_threshold = _new_threshold;
  
  visited.clear();
  
  for(unsigned int ii = 0; ii<new_fringe.size(); ii++)
    visited.insert(new_fringe[ii]);
  
  fringe.swap(new_fringe);
  
  new_fringe.clear();
  
  return true;
}

bool LayeredMarginal::next()
{
  double lpc;
  
  while(true)
  {
    if(fringe.empty())
      return false;
    
    currentConf = fringe.back();
    fringe.pop_back();
    
    opc = logProb(currentConf);
    if(opc < new_threshold)
      new_fringe.push_back(currentConf);
    else
      break;
  }
  
  for(unsigned int ii = 0; ii < isotopeNo; ii++ )
    for(unsigned int jj = 0; jj < isotopeNo; jj++ )
      if( ii != jj and currentConf[jj] > 0 )
      {
        currentConf[ii]++;
        currentConf[jj]--;
        
        lpc = logProb(currentConf);
        
        if (visited.count(currentConf) == 0 and lpc < current_threshold and
              (opc > lpc or (opc == lpc and ii > jj)))
        {
          Conf nc = allocator.makeCopy(currentConf);
          visited.insert(nc);
          if(lpc >= new_threshold)
            fringe.push_back(nc);
          else
            new_fringe.push_back(nc);
        }
        
        currentConf[ii]--;
        currentConf[jj]++;
        
      }
      
      return true;
}




class Summator{
  // Kahan algorithm
  double sum;
  double c;
  
public:
  inline Summator()
  { sum = 0.0; c = 0.0;}
  
  inline void add(double what)
  {
    double y = what - c;
    double t = sum + y;
    c = (t - sum) - y;
    sum = t;
  }
  
  inline double get()
  {
    return sum;
  }
};




//[[Rcpp::export]]
Rcpp::List jaccard_mca_rcpp(const double px,const double py, const int m, const double T_observed, const double accuracy = 1e-05 )
{
    double probs[] = {px*py, px*(1-py), py*(1-px),(1-px)*(1-py)}; // (1,1), (1,0), (0,1), (0,0)
    LayeredMarginal LM(Marginal(probs, 4, m));


    double current_threshold = -0.8 ;
    Summator s;
    Summator out;
    Conf current_conf;
    double current_statistic;
    double current_px,current_py,current_pxpy;
    while(s.get() < 1-accuracy and not LM.empty())
    {
        LM.extend(current_threshold);
        
        while(s.get() < 1-accuracy and LM.next())
        {
            current_conf = LM.get_conf();
            if(current_conf[3]==m ){
              current_statistic =0;
            } else {
              current_px = (double)(current_conf[0]+current_conf[1])/(double)m;
              current_py = (double)(current_conf[0]+current_conf[2])/(double)m;
              current_pxpy=current_px*current_py;
              current_statistic=(double)(current_conf[0])/(double)(m-current_conf[3])-current_pxpy/(current_px+current_py-current_pxpy);
            }
            //std::cout << exp(LM.get_lProb()) << "\t" << fabs(current_statistic) << std::endl;
            //std::cout << exp(LM.get_lProb()) << "\t" <<current_conf[0] << "\t" <<current_conf[1] << "\t" << current_conf[2] << "\t"<<current_conf[3] << "\t"<< fabs(current_statistic) << "\t"<< current_py << "\t"<< current_px << std::endl;
            if(fabs(current_statistic)>=fabs(T_observed)){
              
              out.add(exp(LM.get_lProb()));
              
            }
            s.add(exp(LM.get_lProb()));
        }
        current_threshold -= 0.8;
    }
    return Rcpp::List::create(Rcpp::Named("pvalue") = out.get(),
                              Rcpp::Named("accuracy") = s.get()
                                ); 
}


//[[Rcpp::export]]
Rcpp::List jaccard_mca_rcpp_known_p(const double px,const double py, const int m, const double T_observed, const double accuracy = 1e-05 )
{
  double probs[] = {px*py, px*(1-py), py*(1-px),(1-px)*(1-py)}; // (1,1), (1,0), (0,1), (0,0)
  LayeredMarginal LM(Marginal(probs, 4, m));
  double expectation =px*py/(px+py-px*py);
  
  double current_threshold = -0.8 ;
  Summator s;
  Summator out;
  Conf current_conf;
  double current_statistic;
  
  while(s.get() < 1-accuracy and not LM.empty())
  {
    LM.extend(current_threshold);
    
    while(s.get() < 1-accuracy and LM.next())
    {
      current_conf = LM.get_conf();
      if(current_conf[3]==m ){
        current_statistic =0;
      } else {
          current_statistic=(double)(current_conf[0])/(double)(m-current_conf[3])-expectation;
      }
      //std::cout << exp(LM.get_lProb()) << "\t" << fabs(current_statistic) << std::endl;
      //std::cout << exp(LM.get_lProb()) << "\t" <<current_conf[0] << "\t" <<current_conf[1] << "\t" << current_conf[2] << "\t"<<current_conf[3] << "\t"<< fabs(current_statistic) << "\t"<< current_py << "\t"<< current_px << std::endl;
      if(fabs(current_statistic)>=fabs(T_observed)){
        
        out.add(exp(LM.get_lProb()));
        
      }
      s.add(exp(LM.get_lProb()));
    }
    current_threshold -= 0.8;
  }
  return Rcpp::List::create(Rcpp::Named("pvalue") = out.get(),
                            Rcpp::Named("accuracy") = s.get()
  );  
}



