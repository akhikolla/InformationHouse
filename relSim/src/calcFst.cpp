#include <vector>
#include <Rcpp.h>


using namespace Rcpp;
using namespace std;

void calculateAlleleFrequencies(const IntegerVector& Pop, const IntegerVector &SubpopIdx,
                                const int N, const int ns, const int nLoci,
                                const IntegerVector& NumLocusAlleles,
                                vector< vector < vector<double> > >& AlleleFreqs,
                                vector< vector < vector<double> > >& Hom,
                                vector<int>& SubPopSize){
  // Note - we can make this way more efficient
  // and if it's an issue we'll do it later.
  
  int r;
  int nLoc;
  
  
  for(r = 0; r < ns; r++){
    SubPopSize[r] = 0;
    
    // Clear the frequency array
    for(nLoc = 0; nLoc < nLoci; nLoc++){
      int nAlleles = NumLocusAlleles[nLoc];
      int nA;
      for(nA = 0; nA < nAlleles; nA++){
        AlleleFreqs[r][nLoc][nA] = 0;
        Hom[r][nLoc][nA] = 0;
      }
    }
  }
  
  // iterate through the profiles subpop by subpop
  
  int i;
  
  for(i = 0; i < N; i++){
    r = SubpopIdx[i] - 1;
    SubPopSize[r]++;
    
    // tally up the allele counts
    IntegerVector::const_iterator iProfile = Pop.begin() + i * 2 * nLoci;
    int nA1, nA2;
    int i1, i2;
    
    for(nLoc = 0; nLoc < nLoci; nLoc++){
      i1 = 2 * nLoc;
      i2 = 2 * nLoc + 1;
      
      nA1 = iProfile[i1] - 1; // Alleles assumed to be 1..nA
      nA2 = iProfile[i2] - 1;
      
      
       AlleleFreqs[r][nLoc][nA1] += 1.0;
       AlleleFreqs[r][nLoc][nA2] += 1.0;
     
       if(nA1 == nA2)
         Hom[r][nLoc][nA1] += 1.0;
    

    } 
  }
  for(r = 0; r < ns; r++){
    // divide each allele count by 2*nSubPopSize[r]
    for(nLoc = 0; nLoc < nLoci; nLoc++){
      int nA;
      for(nA = 0; nA < NumLocusAlleles[nLoc]; nA++){
        AlleleFreqs[r][nLoc][nA] /= 2.0 * SubPopSize[r];
        Hom[r][nLoc][nA] /= (double)SubPopSize[r];
        
        //calculate and store average value for each allele count
        //in AlleleFreqs[ns][nLoc][nA]
        //the last member of the AlleleFreqs 3-d vector
        
        double pi = (double)SubPopSize[r] / (double)N;
        if(r==0)
        AlleleFreqs[ns][nLoc][nA] = AlleleFreqs[r][nLoc][nA] * pi;
        else
        AlleleFreqs[ns][nLoc][nA] += AlleleFreqs[r][nLoc][nA] * pi;
      }
    }
  }
}


NumericVector  calcTheta(int nLoci, int nSubPop, const IntegerVector& NumLocusAlleles,
                         const vector<int>& SubPopSize,
                         const vector< vector < vector<double> > >& AlleleFreqs, 
                         const vector< vector < vector<double> > >& Hom){
  
  int i;
  double dSum_ni = 0;
  double dSum_niSq = 0;
  
  for(i = 0; i < nSubPop; i++){
    dSum_ni += SubPopSize[i];
    dSum_niSq += SubPopSize[i] * SubPopSize[i];
  }
  
  double nc = (dSum_ni - (dSum_niSq / (double)dSum_ni)) / (nSubPop - 1);
  double nbar = (double)dSum_ni / (double)nSubPop;
  double dNumerator = 0;
  double dDenominator = 0;
  
  int nLoc;
  int nA;
  int r = nSubPop;
  
  NumericVector result(nLoci + 1);
  
  for(nLoc = 0; nLoc < nLoci; nLoc++){
    double dLocusNumerator = 0;
    double dLocusDenominator = 0;
    
    int nAlleles = NumLocusAlleles[nLoc];
    
    for(nA = 0; nA < nAlleles; nA++){
      double dS1;
      double dS2,dS2a,dS2b,dS2c,dS2d,dS2e,dS2f,dS2g;
      double sASq=0;
      double pAbar = AlleleFreqs[r][nLoc][nA];//0.5*AlleleFreqs[r][nLoc][nA]/dSum_ni; // CHECK
      double HAbar = 0;
      double pAi;
      
      if(pAbar > 0){
        for(i = 0;i < r;i++){
          pAi = AlleleFreqs[i][nLoc][nA];
          sASq += SubPopSize[i] * (pAi - pAbar) * (pAi - pAbar);
          HAbar += 2 * (SubPopSize[i])*(pAi - Hom[i][nLoc][nA]);
        }
        
        sASq /= (r - 1) * nbar;
        HAbar /= dSum_ni;
        
        dS1 = sASq - ((pAbar * (1 - pAbar) - sASq * (r - 1) / r -
                       0.25 * HAbar) / (nbar - 1));
        dS2a = pAbar * (1 - pAbar);
        dS2b = nbar / (r * (nbar - 1));
        dS2c = r * (nbar - nc) / nbar;
        dS2d = nbar - 1;
        dS2e = (r - 1) * (nbar - nc);
        dS2f = sASq / nbar;
        dS2g = (nbar - nc) * HAbar / (4 * nc * nc);
        dS2 =  dS2a - dS2b * (dS2c * dS2a - dS2f * (dS2d + dS2e) - dS2g);
        
        dLocusNumerator += dS1;
        dLocusDenominator += dS2;
        result[nLoc] = dLocusNumerator /  dLocusDenominator;
        
        dNumerator += dS1;
        dDenominator += dS2;
        
      }
    }
  }
  
  result[nLoci] = dNumerator/dDenominator;
  return(result);
}

NumericVector  calcTheta2(int nLoci, int nSubPop, const IntegerVector& NumLocusAlleles,
                         const vector<int>& SubPopSize,
                         const vector< vector < vector<double> > >& AlleleFreqs, 
                         const vector< vector < vector<double> > >& Hom){
  
  int i;
  double dSum_ni = 0;
  double dSum_niSq = 0;
  
  for(i = 0; i < nSubPop; i++){
    dSum_ni += SubPopSize[i];
    dSum_niSq += SubPopSize[i] * SubPopSize[i];
  }
  
  double nc = (dSum_ni - (dSum_niSq / (double)dSum_ni)) / (nSubPop - 1);
  double nbar = (double)dSum_ni / (double)nSubPop;
  double dNumerator = 0;
  double dDenominator = 0;
  
  int nLoc;
  int nA;
  int r = nSubPop;
  
  NumericVector result(nLoci + 1);
  
  for(nLoc = 0; nLoc < nLoci; nLoc++){
    double dLocusNumerator = 0;
    double dLocusDenominator = 0;
    
    int nAlleles = NumLocusAlleles[nLoc];
    
    for(nA = 0; nA < nAlleles; nA++){
      double dS1;
      double dS2,dS2a,dS2b,dS2c,dS2d,dS2e,dS2f,dS2g;
      double dS3;
      double sASq=0;
      double pAbar = AlleleFreqs[r][nLoc][nA];//0.5*AlleleFreqs[r][nLoc][nA]/dSum_ni; // CHECK
      double HAbar = 0;
      double pAi;
      
      if(pAbar > 0){
        for(i = 0;i < r;i++){
          pAi = AlleleFreqs[i][nLoc][nA];
          sASq += SubPopSize[i] * (pAi - pAbar) * (pAi - pAbar);
          HAbar += 2 * (SubPopSize[i])*(pAi - Hom[i][nLoc][nA]);
        }
        
        sASq /= (r - 1) * nbar;
        HAbar /= dSum_ni;
        
        dS1 = sASq - ((pAbar * (1 - pAbar) - sASq * (r - 1) / r -
          0.25 * HAbar) / (nbar - 1));
        dS2a = pAbar * (1 - pAbar);
        dS2b = nbar / (r * (nbar - 1));
        dS2c = r * (nbar - nc) / nbar;
        dS2d = nbar - 1;
        dS2e = (r - 1) * (nbar - nc);
        dS2f = sASq / nbar;
        dS2g = (nbar - nc) * HAbar / (4 * nc * nc);
        dS2 =  dS2a - dS2b * (dS2c * dS2a - dS2f * (dS2d + dS2e) - dS2g);
        
        dS3 = 0.5 * nc * HAbar / nbar;
        
        Rcout << 1 - (dS3 / dS2) << endl;
        
        dLocusNumerator += dS1;
        dLocusDenominator += dS2;
        
        dNumerator += dS1;
        dDenominator += dS2;
        
      }
    }
    result[nLoc] = dLocusNumerator /  dLocusDenominator;
  }
  
  result[nLoci] = dNumerator/dDenominator;
  return(result);
}


// [[Rcpp::export(".calcFst")]]
NumericVector calcFst(const IntegerVector& Pop, IntegerVector SubPopIdx, int N, int ns,
               int nLoci,
               IntegerVector NumLocusAlleles){

    // Allocate memory for the allele frequencies
    // The dimension of this array is going to be (ns+1)*(nLoci)*(nAlleleMax)
    // way to create a 3-d array
    
    vector< vector < vector <double> > > AlleleFreqs(ns + 1, vector < vector <double> >(nLoci));
    vector< vector < vector <double> > > Hom(ns, vector < vector <double> >(nLoci));
    vector<int> SubPopSize(ns);

    int r, nLoc;

    for(r = 0; r <= ns; r++){
      for(nLoc = 0; nLoc < nLoci; nLoc++){
	      int nAlleles = NumLocusAlleles[nLoc];
	      AlleleFreqs[r][nLoc].resize(nAlleles);

      	if(r < ns)
	        Hom[r][nLoc].resize(nAlleles);
      }
    }
    
    // calculate the allele frequencies

    calculateAlleleFrequencies(Pop, SubPopIdx,
                               N, ns, nLoci,
                               NumLocusAlleles,
               		             AlleleFreqs,
                               Hom,
                               SubPopSize); 

    // calcFst
    return calcTheta2(nLoci, ns, NumLocusAlleles, SubPopSize, AlleleFreqs, Hom);
    
    return 0;

}

List  calcFStats(int nLoci, int nSubPop, const IntegerVector& NumLocusAlleles,
                          const vector<int>& SubPopSize,
                          const vector< vector < vector<double> > >& AlleleFreqs, 
                          const vector< vector < vector<double> > >& Hom){
  
  int i;
  double dSum_ni = 0;
  double dSum_niSq = 0;
  
  for(i = 0; i < nSubPop; i++){
    dSum_ni += SubPopSize[i];
    dSum_niSq += SubPopSize[i] * SubPopSize[i];
  }
  
  double nc = (dSum_ni - (dSum_niSq / (double)dSum_ni)) / (nSubPop - 1);
  double nbar = (double)dSum_ni / (double)nSubPop;
  double dThetaNumerator = 0;
  double dThetaDenominator = 0;
  double dFNumerator = 0;
  double dFDenominator = 0;

  int nLoc;
  int nA;
  int r = nSubPop;
  
  List result(nLoci + 1);
  
  for(nLoc = 0; nLoc < nLoci; nLoc++){
    double dThetaLocusNumerator = 0;
    double dThetaLocusDenominator = 0;
    
    double dFLocusNumerator = 0;
    double dFLocusDenominator = 0;
    
    int nAlleles = NumLocusAlleles[nLoc];
    
    for(nA = 0; nA < nAlleles; nA++){
      double dS1;
      double dS2,dS2a,dS2b,dS2c,dS2d,dS2e,dS2f,dS2g;
      double dS3;
      double sASq=0;
      double pAbar = AlleleFreqs[r][nLoc][nA];//0.5*AlleleFreqs[r][nLoc][nA]/dSum_ni; // CHECK
      double HAbar = 0;
      double pAi;
      
      if(pAbar > 0){
        for(i = 0;i < r;i++){
          pAi = AlleleFreqs[i][nLoc][nA];
          sASq += SubPopSize[i] * (pAi - pAbar) * (pAi - pAbar);
          HAbar += 2 * (SubPopSize[i])*(pAi - Hom[i][nLoc][nA]);
        }
        
        sASq /= (r - 1) * nbar;
        HAbar /= dSum_ni;
        
        dS1 = sASq - ((pAbar * (1 - pAbar) - sASq * (r - 1) / r -
          0.25 * HAbar) / (nbar - 1));
        dS2a = pAbar * (1 - pAbar);
        dS2b = nbar / (r * (nbar - 1));
        dS2c = r * (nbar - nc) / nbar;
        dS2d = nbar - 1;
        dS2e = (r - 1) * (nbar - nc);
        dS2f = sASq / nbar;
        dS2g = (nbar - nc) * HAbar / (4 * nc * nc);
        dS2 =  dS2a - dS2b * (dS2c * dS2a - dS2f * (dS2d + dS2e) - dS2g);
        
        dS3 = 0.5 * nc * HAbar / nbar;
        
        dThetaLocusNumerator += dS1;
        dThetaLocusDenominator += dS2;
        
        dFLocusNumerator += dS3;
        dFLocusDenominator += dS2;
        
        dThetaNumerator += dS1;
        dThetaDenominator += dS2;
        
        dFNumerator += dS3;
        dFDenominator += dS2;
        
      }
    }
    double dThetaLocus = dThetaLocusNumerator /  dThetaLocusDenominator;
    double dFLocus = 1 - (dFLocusNumerator /  dFLocusDenominator);
    double dfLocus = (dFLocus - dThetaLocus) / (1 - dThetaLocus);
    
    NumericVector locusResults {dThetaLocus, dFLocus, dfLocus};
    locusResults.names() = CharacterVector {"Theta (F_ST)", "F (F_IT)", "f (F_IS)"};
    result[nLoc] = locusResults;
  }
  
  double dTheta = dThetaNumerator / dThetaDenominator;
  double dF = 1 - (dFNumerator / dFDenominator);
  double df = (dF - dTheta) / (1 - dTheta);  
  
  NumericVector overallResults {dTheta, dF, df};
  overallResults.names() = CharacterVector {"Theta (F_ST)", "F (F_IT)", "f (F_IS)"};
  result[nLoci] = overallResults;
  
  return(result);
}

// [[Rcpp::export(".calcFStatistics")]]
List calcFStatistics(const IntegerVector& Pop, IntegerVector SubPopIdx, int N, int ns,
                     int nLoci,
                     IntegerVector NumLocusAlleles){
  
  // Allocate memory for the allele frequencies
  // The dimension of this array is going to be (ns+1)*(nLoci)*(nAlleleMax)
  // way to create a 3-d array
  
  vector< vector < vector <double> > > AlleleFreqs(ns + 1, vector < vector <double> >(nLoci));
  vector< vector < vector <double> > > Hom(ns, vector < vector <double> >(nLoci));
  vector<int> SubPopSize(ns);
  
  int r, nLoc;
  
  for(r = 0; r <= ns; r++){
    for(nLoc = 0; nLoc < nLoci; nLoc++){
      int nAlleles = NumLocusAlleles[nLoc];
      AlleleFreqs[r][nLoc].resize(nAlleles);
      
      if(r < ns)
        Hom[r][nLoc].resize(nAlleles);
    }
  }
  
  // calculate the allele frequencies
  
  calculateAlleleFrequencies(Pop, SubPopIdx,
                             N, ns, nLoci,
                             NumLocusAlleles,
                             AlleleFreqs,
                             Hom,
                             SubPopSize); 
  
  // calcFst
  return calcFStats(nLoci, ns, NumLocusAlleles, SubPopSize, AlleleFreqs, Hom);
  
  return 0;
  
}
