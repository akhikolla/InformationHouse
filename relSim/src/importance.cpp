#include <Rcpp.h>
#include <math.h>
#include <map>
#include <string>
#include <sstream>
#include <vector>
using namespace Rcpp;
using namespace std;

long factorial(long n){
  switch(n){
  case 0:
  case 1:
    return 1;
  case 2:
    return 2;
  case 3:
    return 6;
  case 4:
    return 24;
  case 5:
    return 120;
  case 6:
    return 720;
  case 7:
    return 5040;
  case 8:
    return 40320;
  case 9:
    return 362880;
  case 10:
    return 3628800;
  case 11:
    return 39916800;
  case 12:
    return 479001600;
  default:
    return n * factorial(n - 1);
  }
}

class ProfileGenerator;
class FreqInfo;

class Profile {
public:
  class Locus {
  protected:
    int nAlleles;
    int numTotalAlleles;
    long numTotalAllelesFactorial;
    vector<int> epg;
    double gx, hx, wx;
    
    vector<int> alleles;
    vector<double> freqs;
    map<int, int> mapCounts;
    

  public:
    Locus(const FreqInfo& fi, int numAlleles, int nTotalAlleles);
    
    Locus(const Locus& loc){
      nAlleles = loc.nAlleles;
      numTotalAlleles = loc.numTotalAlleles;
      numTotalAllelesFactorial = loc.numTotalAllelesFactorial;
      
      epg = loc.epg;
      gx = loc.gx;
      hx = loc.hx;
      wx = loc.wx;
      
      alleles = loc.alleles;
      mapCounts = loc.mapCounts;
      freqs = loc.freqs;
    }
    
    Locus& operator=(const Locus& loc){
      nAlleles = loc.nAlleles;
      numTotalAlleles = loc.numTotalAlleles;
      numTotalAllelesFactorial = loc.numTotalAllelesFactorial;
      
      epg = loc.epg;
      gx = loc.gx;
      hx = loc.hx;
      wx = loc.wx;
      
      alleles = loc.alleles;
      mapCounts = loc.mapCounts;
      freqs = loc.freqs;
      
      return *this;
    }

    int numAlleles(){
      return mapCounts.size();
    }

    // get
    const int operator[](int a) const{
      return mapCounts.at(a);
    }

    // set
    int& operator[](int a){
      return mapCounts[a];
    }

    NumericVector asNumericVector(){
      NumericVector result(nAlleles);
      map<int, int>::iterator i = mapCounts.begin();

      while(i != mapCounts.end()){
        result[i->first - 1] = i->second;
        i++;
      }

      return result;
    }

    void calcProb(const vector<double>& probs, bool bLog = true){
      map<int, int>::iterator i = mapCounts.begin();
      double dSum = 0;

      if(mapCounts.size() == 1)
        gx = (i->second) * std::log(probs[i->first]);

      double dFact = std::log(numTotalAllelesFactorial);

      while(i != mapCounts.end()){
        dSum += (i->second) * std::log(probs[i->first]);
        dFact -= std::log(factorial(i->second));
        i++;
      }

      gx = dSum + dFact;
    }
    
    double ISprob(const vector<NumericMatrix>& perms, bool bTail = false, bool bLog = true){
      double result = 0;
      
      NumericMatrix m;
      double correct = 0;
      
      if(bTail){
        m = as<NumericMatrix>(perms[nAlleles - 1]);
        correct = 0;//std::log(perms.size());
      }else{
        m = as<NumericMatrix>(perms[0]);
      }
      int numPerms = m.nrow();
      // Rprintf("%d, %d\n", nAlleles, numPerms);
      
      for(int i = 0; i < numPerms; i++){
        double p , s;
        p = s = freqs[m(i, 0) - 1];
        
        for(int j = 1; j < nAlleles; j++){
          double pj = freqs[m(i, j) - 1];
          p *=  pj / (1 - s);
          s += pj;
        }
        // Rprintf("p = %.7f\n", p);
        result += p;
      }
      
      // Rprintf("%.7f\n", result);
      // for(vector<double>::iterator i = freqs.begin(); i != freqs.end(); i++){
      //    Rprintf("f = %.7f\n", *i);
      // }
      double normConst = std::accumulate(freqs.begin(), freqs.end(), 0.0);
      // Rprintf("sum f = %.7f\n", normConst);
      
      int sumCounts = 0;
      double p2 = 0;
      
      for(int j = 0; j < nAlleles; j++){
        int a = alleles[j] - 1;
        int c = mapCounts[a];
        if(c > 1){
          p2 += (c - 1) * std::log(freqs[j] / normConst) - std::log(factorial(c - 1));
          sumCounts += (c - 1);
        }
        // Rprintf("%d %d %.7f %.7f\n", a, c, freqs[j], p2);
      }
      
      p2 += std::log(factorial(sumCounts));
      // Rprintf("%.7f %.7f %.7f\n", log(result), p2, correct);
      result = std::log(result) + p2 + correct;
      hx = result; // store in case you don't want to recalc
      // Rprintf("%.7f\n", result);
      return result;
    }
    
    double prob(void){
    		return gx;
    }

    void print(void){
      vector<int>::iterator i = epg.begin();
      while(i != epg.end()){
        Rprintf("%d ", *i);
        i++;
      }
      Rprintf("\n");
    }
    
    List toList(){
      List result;
      
      result["n"] = nAlleles;
      result["epg"] = Rcpp::wrap(epg);
      result["a"] = Rcpp::wrap(alleles);
      result["f"] = freqs;
      
      result["gx"] = gx;
      
      vector<int>::iterator i = alleles.begin();
      IntegerVector counts;
      
      while(i != alleles.end()){
        counts.push_back(mapCounts[*i -1]);
        i++;
      }
      result["c"] = counts;
      
      return result;
    }
    
  };

  vector<Locus> profile;
  
  Profile(const ProfileGenerator& pg, int numLoci, int numContributors, 
          const IntegerVector::iterator& numAllelesShowing);
  
  Profile(const ProfileGenerator& pg, int numLoci, int numContributors, int numAllelesShowing);
  
  Profile(const Profile& prof){
    profile = prof.profile;
  }
  
  Profile& operator=(const Profile& prof){
    profile = prof.profile;
    
    return *this;
  }
  
  NumericVector ISprob(const vector<NumericMatrix>& perms, bool bTail = false, bool bLog = true){
    
    int numLoci = profile.size();
    NumericVector locusResults(numLoci);
    
    for(int loc = 0; loc < numLoci; loc++){
      locusResults[loc] = profile[loc].ISprob(perms, bTail);
    }
    
    return locusResults;
  }
  
  NumericVector prob(){
	  int numLoci = profile.size();
	  NumericVector locusResults(numLoci);

	  for(int loc = 0; loc < numLoci; loc++){
	  	locusResults[loc] = profile[loc].prob();
	  }

	  return locusResults;
	}

  void print(void){
    vector<Locus>::iterator loc = profile.begin();
    while(loc != profile.end()){
      loc->print();
      loc++;
    }
  }
  
  List toList(void){
    List l;
    int numLoci = profile.size();
  
    for(int loc = 0; loc < numLoci; loc++){
      l.push_back(profile[loc].toList());
    }
    
    return l;
  }
  
  
};

class FreqInfo {
public:
  vector<double> probs;
  vector<double> cumProbs;
  int numAlleles;
  
public:
  FreqInfo(const vector<double>& lprobs){
    numAlleles = lprobs.size();
    
    vector<double>::const_iterator i = lprobs.begin();
    int a = 0;
    
    while(i != lprobs.end()){
      double f = *i;
      probs.push_back(f);
      
      if(a == 0){
        cumProbs.push_back(f);
      }else{
        cumProbs.push_back(cumProbs[a - 1] + f);
      }
      a++; 
      i++;
    }
  }
  
  FreqInfo& operator=(const vector<double>& lprobs){
    numAlleles = lprobs.size();
    
    vector<double>::const_iterator i = lprobs.begin();
    int a = 0;
    
    while(i != lprobs.end()){
      double f = *i;
      probs.push_back(f);
      
      if(a == 0){
        cumProbs.push_back(f);
      }else{
        cumProbs.push_back(cumProbs[a - 1] + f);
      }
      a++; 
      i++;
    }
    
    return *this;
  }
  
  void print(){
    for(int i =  0; i < numAlleles; i++){
      Rprintf("%d: %.7f %.7f\n", i + 1, probs[i], cumProbs[i]);
    }
  }
};
  

class ProfileGenerator {
public:
  vector<FreqInfo> freqs;
  int numLoci;

  ProfileGenerator(const List& freqs){
    numLoci = (int)freqs.size();
    List::const_iterator i = freqs.begin();
    
    while(i != freqs.end()){
      vector<double> f = as< vector<double> >((NumericVector)*i);
      this->freqs.push_back(f);
      i++;
    }
  };
  
  void print(void){
    vector<FreqInfo>::iterator i = freqs.begin();
    int loc = 1;
    
    while(i != freqs.end()){
      Rprintf("Locus %2d\n", loc);
      Rprintf("-----------\n");
      i->print();
      Rprintf("-----------\n\n");
      i++;
      loc++;
    }
  }
  
  // [[Rcpp::plugins("cpp11")]]
  
  Profile randProf(int numContributors, const IntegerVector::iterator& numAllelesShowing){
    Profile prof(*this, numLoci, numContributors, numAllelesShowing);
    return prof;
  }
  
  Profile randProf(int numContributors, int numAllelesShowing){
    Profile prof(*this, numLoci, numContributors, numAllelesShowing);
    return prof;
  }
};

Profile::Profile(const ProfileGenerator& pg, int numLoci, int numContributors, 
        const IntegerVector::iterator& numAllelesShowing){
  for(int loc = 0; loc < numLoci; loc++){
    Locus l(pg.freqs[loc], *(numAllelesShowing + loc), 2 * numContributors);
    profile.push_back(l);
  }
}

Profile::Profile(const ProfileGenerator& pg, int numLoci, int numContributors, int numAllelesShowing){
  for(int loc = 0; loc < numLoci; loc++){
    Locus l(pg.freqs[loc], numAllelesShowing, 2 * numContributors);
    profile.push_back(l);
  }
}

Profile::Locus::Locus(const FreqInfo& fi, int numAlleles, int nTotalAlleles){
  nAlleles = numAlleles;
  numTotalAlleles = nTotalAlleles;
  numTotalAllelesFactorial = factorial(numTotalAlleles);
  
  // choose numAllelesShowing alleles without replacement
  IntegerVector alleles = seq_len(fi.numAlleles);
 
  if(numAlleles != fi.numAlleles){
    alleles = sample(alleles, nAlleles, false, Rcpp::wrap(fi.probs));
    sort(alleles.begin(), alleles.end()); // this makes sure that the stored probabilities are in the right order
  }
  
  this->alleles = as< vector<int> >(alleles);

  double sum = 0;
  NumericVector newProbs;
  
  epg.resize(fi.numAlleles);

  for(IntegerVector::iterator a = alleles.begin(); a != alleles.end(); a++){
    int allele = *a - 1;
    epg[allele] = 1;
    double f = fi.probs[allele];
    sum += f;
    freqs.push_back(f); // NOTE: This should be in order since we ordered the alleles;
    mapCounts[allele] = 1;
 //   Rprintf("%d\n",*a);
  }

  // normalize
  newProbs = Rcpp::wrap(freqs);
  newProbs = newProbs / sum;
  alleles = sample(alleles, numTotalAlleles - numAlleles, true, newProbs);

  for(IntegerVector::iterator a = alleles.begin(); a != alleles.end(); a++){
    epg[*a - 1] += 1;
    mapCounts[*a - 1] += 1;
  }

  // calculate probability under target density
  calcProb(fi.probs);
}

NumericVector p1(List freqs, int numContributors){
  int numLoci = freqs.size();
  NumericVector probs(numLoci);
  
  for(int loc = 0; loc < numLoci; loc++){
    NumericVector locus = as<NumericVector>(freqs[loc]);
    probs[loc] = accumulate(locus.begin(), locus.end(), 0.0, 
                            [numContributors](double dSum, double const& pi){
                              return dSum + pow(pi, 2 * numContributors);
                            });
  }
  
  return probs;
}
        
NumericVector p2(List freqs, int numContributors){
  int numLoci = freqs.size();
  NumericVector probs(numLoci);
  
  IntegerVector m(2 * numContributors - 1);
  m[0] = 2 * numContributors;
  for(int x = 1; x < 2 * numContributors - 1; x++){
    m[x] = m[x - 1] * (2 * numContributors - x) / (x + 1);
  }
  
  for(int loc = 0; loc < numLoci; loc++){
    NumericVector px = as<NumericVector>(freqs[loc]);
    int numAlleles = px.size();
    
    for(int x = 1; x < 2 * numContributors; x++){
      for(int i = 0; i < numAlleles - 1; i++){
        for(int j = i + 1; j < numAlleles; j++){
          probs[loc] += m[x - 1] * pow(px[i], x) * pow(px[j], 2 * numContributors - x);
         // Rprintf("%d %.7f %.7f\n", m[x - 1], pow(px[i], x), pow(px[j], 2 * numContributors - x));
        }
      }
    }
  }
    
  return probs;
}

// [[Rcpp::export(".IS")]]
List IS(List freqs,int N, int numContributors, int maxAllelesShowing, List Perms, bool bTail = false){
  
  if(maxAllelesShowing == 1){
    List result;
    result["est"] = p1(freqs, numContributors);
    return result;
  }else if(maxAllelesShowing == 2 && !bTail){
    List result;
    result["est"] = p2(freqs, numContributors);
    return result;
  }
  //else
  // copy the elements of Perms into a vector for quick random access.
  vector<NumericMatrix> perms;
  List::iterator p = Perms.begin();
  while(p != Perms.end()){
    NumericMatrix m = as<NumericMatrix>(*p);
    perms.push_back(m);
    p++;
  }
  
  
  ProfileGenerator g(freqs);
  int numLoci = g.numLoci;

  List results;
  NumericMatrix denoms(N, numLoci);
  NumericMatrix numers(N, numLoci);
    
  
  if(bTail){
    if(maxAllelesShowing == 1){
      List result;
      result["est"] = p1(freqs, numContributors) + p2(freqs, numContributors);
      return result;
    }else if (maxAllelesShowing == 2){
      List result;
      result["est"] = p1(freqs, numContributors) + p2(freqs, numContributors);
      return result;
    }
    
    IntegerVector numAllelesShowing = sample(maxAllelesShowing - 2, N * numLoci, true) + 2;
 
     //vector<Profile> profiles;
    IntegerVector::iterator nA = numAllelesShowing.begin();

    for(int i = 0; i < N; i++){
      Profile p = g.randProf(numContributors, nA + i * numLoci);
      //p.print();
      numers(i,_) = p.prob();
      denoms(i,_) = p.ISprob(perms, true);
    }
    
    // store the number of peaks -- for debugging mostly.
    results["numPeaks"] = numAllelesShowing;
    
    results["p12"] = p1(freqs, numContributors) + p2(freqs, numContributors);
    
  }else{
    int numAllelesShowing = maxAllelesShowing;

    for(int i = 0; i < N; i++){
      Profile p = g.randProf(numContributors, numAllelesShowing);
      numers(i,_) = p.prob();
      denoms(i,_) = p.ISprob(perms, false);
    }
  }
  
  results["numerator"] = numers;
  results["denominator"] = denoms;
  return results;

  
  // NumericMatrix Alleles(N, freqs.size());
  // NumericVector probs(N);
  // IntegerVector numAllelesShowing(N);
  // 
  // numAllelesShowing = sample(maxAllelesShowing, N, TRUE);
  // 
  // for(int i = 0; i < N; i++){
  //   profileGenerator::Profile p = g.randProf(numContributors, numAllelesShowing[i]);
  //   Alleles(i, _) = p.asNumericVector();
  //   probs[i] = p.prob(freqs);
  // }
  // 
  // List results;
  // results["Alleles"] = Alleles;
  // results["probs"] = probs;
  // 
  // return results;
}

// // [[Rcpp::export]]
// NumericVector ISprob(const List& listCombs, const List& Perms){
//   int numCombs = listCombs.size();
//   NumericVector results(numCombs);
//   
//   for(int i = 0; i < numCombs; i++){
//     List lComb = as<List>(listCombs[i]);
//     NumericVector freqs = as<NumericVector>(lComb["f"]);
//     IntegerVector alleles = as<IntegerVector>(lComb["a"]);
//     IntegerVector counts = as<IntegerVector>(lComb["c"]);
//     int numAlleles = as<int>(lComb["n"]);
//     
//     NumericMatrix perms = as<NumericMatrix>(Perms[numAlleles - 1]);
//     int numPerms = perms.nrow();
//    
//     for(int j = 0; j < numPerms; j++){
//       double p , s;
//       p = s = freqs[perms(j, 0) - 1];
//       
//       for(int k = 1; k < numAlleles; k++){
//         double pk = freqs[perms(j, k) - 1];
//         p *=  pk / (1 - s);
//         s += pk;
//       }
//      // Rprintf("p = %.7f\n", p);
//       results[i] += p;
//     }
//     
//    // Rprintf("%.7f\n", results[i]);
//     freqs = freqs / sum(freqs);
// 
//     int sumCounts = 0;
//     double p2 = 0;
// 
//     for(int j = 0; j < numAlleles; j++){
//       // Rprintf("%d %d %.7f %.7f\n", alleles[j], counts[j], freqs[j], p2);
//       p2 += (counts[j] - 1) * std::log(freqs[j]) - std::log(factorial(counts[j] - 1));
//       sumCounts += counts[j] - 1;
//     }
// 
//     p2 += std::log(factorial(sumCounts));
//     // Rprintf("%.7f\n", p2);
//     results[i] = std::log(results[i]) + p2;
//   }
//   
//   return results;
// }
