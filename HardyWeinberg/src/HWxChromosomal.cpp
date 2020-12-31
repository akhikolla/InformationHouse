#include <Rcpp.h>
#include "HWxChrom.h"
#include <math.h>
#include <time.h>
using namespace Rcpp;

#define lgammafn lgamma

/* #else
#include <R.h>
#include <Rmath.h>
#endif
*/

// Globals

COUNTTYPE * Rarray;
unsigned Rbytes;
int nAlleles;
unsigned ntotal;
time_t start;
int timeLimit;
int histobins, HN;
int statID;
int * mi;
double *lnFact;
double tableCount;
double pLLR, pU, pPr, pX2;  // P values
double maxLLR, maxlPr, minmaxU, minX2; // cutoff values
double statSpan, constProbTerm, constLLRterm, probSum, leftStat;
double *hProb;
double umean, uvariance;

// ---------------------------- NEW --------------------------
COUNTTYPE * alleleVect;  // vector of alleles + 1 (for males recursion)
int male;
int female;


static void homozygoteX(int r, double probl, double constMales, COUNTTYPE * R)
{
  // If the process takes longer than `timeLimit` seconds, set
  // `tableCount` negative to signify that the job is aborted
  //if (tableCount < 0) return;
  if (time(NULL) - start >= timeLimit) tableCount = -tableCount;
  
  COUNTTYPE * res, *resn;
  int lower, upper;//, exindix;
  int i, arr;
  double arrln2;
  COUNTTYPE * Rnew = R + nAlleles;
  memcpy(Rnew, R, Rbytes);
  //Find upper and lower limits for arr.
  res = R - 1;  // So res is a 1-based version of R
  resn = Rnew - 1; // resn is 1 based for Rnew
  lower = res[r];
  for (i = 1; i <= r - 1; i++) lower -= res[i];
  lower = lower < 2 ? 0 : lower / 2;
  upper = res[r] / 2;
  //For each possible value of arr, examine the heterozygote at r, r-1
  for (arr = lower; arr <= upper; arr++) {
    resn[r] = res[r] - 2 * arr;    // ------------------------------------ subtracting the homozygous a_rr
    arrln2 = arr * M_LN2;

    heterozygoteX(r,
                  r - 1,
                  probl + lnFact[arr] + arrln2,
                  constMales,
                  Rnew);
  }
}

static void heterozygoteX(int r, int c, double probl, double constMales, COUNTTYPE * R)
{
  
  COUNTTYPE *res, *resn;
  int lower, upper; 
  int i, arc, ar1, ar2, a31, a32, a11, a21, a22;
  int res1, res2, resTemp, dT;

  double problT, probl3, prob;
  COUNTTYPE * Rnew = R + nAlleles;
  
  res = R - 1; // to make res a 1-based version of R
  resn = Rnew - 1; // so resn is 1-based for Rnew
  lower = res[r];
  for (i = 1; i < c; i++) lower -= res[i];
  lower = fmax(0, lower);
  upper = fmin(res[r], res[c]);
  if (c > 2) for (arc = lower; arc <= upper; arc++) {
    memcpy(Rnew, R, Rbytes); // Put a fresh set of residuals from R into Rnew
    
    // decrement residuals for the current value of arc.
    resn[r] -= arc;
    resn[c] -= arc;

    heterozygoteX(r, c - 1,
                  probl + lnFact[arc],
                                constMales,
                                Rnew);
  } // for arc
  if (c == 2) {
    if (r > 3) for (ar2 = lower; ar2 <= upper; ar2++) {
      memcpy(Rnew, R, Rbytes); // Put a fresh set of residuals from R into Rnew
      // decrement residuals for the current value of arc.
      resn[r] -= ar2;
      resn[c] -= ar2;
      // The value of ar1 is now fixed, so no need for any more calls to heterozygote in this row
      ar1 = fmin(resn[r], resn[1]);
      resn[1] -= ar1;
      resn[r] -= ar1;

      homozygoteX(r - 1,
                  probl + lnFact[ar2] + lnFact[ar1], constMales,
                  Rnew);
    } // if r > 3
    if (r == 3) // and c = 2, then we can handle a series of two-allele cases with no deeper recursion
    {
      
      for (a32 = lower; a32 <= upper; a32++) {
        a31 = fmin(res[1], res[3] - a32); //Value of a31 is now fixed for each a32
        probl3 = probl + lnFact[a32] + lnFact[a31];

        // get residual allele counts for two-allele case
        res1 = res[1] - a31;
        res2 = res[2] - a32;

        if (res1 > res2) {            // make sure res1 <= res2. If they need swapping, then swap lookups too
          resTemp = res2;
          res2 = res1;
          res1 = resTemp;

        }
        
        // Now process two-allele case with allele counts res1 and res2
        tableCount += res1 / 2 + 1;
        for (a11 = 0; a11 <= res1 / 2; a11++) {
          a21 = res1 - a11 * 2; // integer arithmetic rounds down
          a22 = (res2 - a21) / 2;
          problT = probl3 + lnFact[a11] + lnFact[a21] + lnFact[a22];

          dT = a11 + a22;
          
          // Here come the actual probability and LLR and X2 and U values
          problT = constProbTerm - problT - dT * M_LN2 - constMales;   // ----------------- NEW -----------------
          prob = exp(problT);

          //Now process the new values of prob and stat
          probSum += prob;
          
          //if (problT <= maxlPr) pPr += prob;   // ----------------------------------------- Where we add the Probability
          if (nearlyEqual(problT,maxlPr)) pPr += prob; 
          else if(problT < maxlPr) pPr += prob;
          
        } // for a11
      } // for a32
    } // if r == 3
  } // if c == 2
}

static void twoAlleleSpecialCaseX() {
  int a11, a21, a22, res1, res2, resTemp;
  double problT, prob, constMales = 0.0; 
          int dT;
          int  nGenes = 0, nMales = 0, n;
          
          res1 = mi[2] - alleleVect[1]; // because they come ordered largest to smallest
          res2 = mi[1] - alleleVect[0];
          tableCount = res1 / 2 + 1;
          
          
          
          
          for (int i = 0; i < nAlleles; i++) {
            constMales += lgammafn(alleleVect[i] + 1);
            nMales += alleleVect[i];
            nGenes += Rarray[i];
          }
          n = (nGenes - nMales) / 2;
          constMales -= log(2.0)*n;    // add 2^n in the constant males
          
          if(res1 > res2) {            // make sure res1 <= res2. If they need swapping, then swap lookups too
            resTemp = res2;
            res2 = res1;
            res1 = resTemp;
            
          }
          
          for (a11 = 0; a11 <= res1 / 2; a11++) {
            a21 = res1 - a11 * 2; // integer arithmetic rounds down
            a22 = (res2 - a21) / 2;
            // -------------- OLD ------------------
            //problT = lgammafn(a11 +1) + lgammafn(a21+1) + lgammafn(a22+1);
            problT = lnFact[a11] + lnFact[a21] + lnFact[a22];

            // ---------- OLD ----------------
            dT = a11 + a22;   // 2^h   (h=n-d).. here only females
            //dT = a21;   // heterozygous
            // Here come the actual probability and LLR values
            problT = constProbTerm - problT - dT * M_LN2 - constMales;   //------------ formula (3) paper
            prob = exp(problT);
            
            //Now process the new values of prob and stat
            probSum += prob;
            
           // if (problT <= maxlPr) pPr += prob;
           if (nearlyEqual(problT,maxlPr)) pPr += prob; 
           else if(problT < maxlPr) pPr += prob;
            
          } // for a11
}
/* -------------------------- Auxiliaries functions -------------------*/

// starting with allele vector + 1 recursively add value to reach the sum required
static void recursiveEnumeration(int index) {
  //Rprintf("\n--------recursive enum --- %d --------\n", index);
  int value = alleleVect[index - 1];
  int nGenes;
  double constMales;
  while (alleleVect[index - 1] > 0) {
    constMales = 0.0;
    alleleVect[index - 1]--;
    
    if (alleleVect[index - 1] < 0)
      return;
    else {
      if (index == nAlleles) {  // last one
        if (sum() != male) {
          continue;
        }
        else {
          if (enoughFemale()) {
            
            if (nAlleles == 2) {
              twoAlleleSpecialCaseX();
            }
            else {
              nGenes = 0;
              for (int i = 0; i < nAlleles; i++) {
                constMales += lgammafn(alleleVect[i] + 1);
                Rarray[i] -= alleleVect[i];
                nGenes += Rarray[i];
              }
              constMales -= log(2.0)*(nGenes / 2);
              homozygoteX(nAlleles, 0, constMales, Rarray);
              // Reset Rarray
              for (int i = 0; i < nAlleles; i++) {
                Rarray[i] += alleleVect[i];
              }
            }
          }
          else {
            continue;
          }
        }
      }
      else {
        recursiveEnumeration(index + 1);
      }
      
    }
  }
  
  if (index != 1) {
    alleleVect[index - 1] = value;
  }
  
  return;
  
}


static int sum() {
  int s = 0, l = nAlleles;
  for (int i = 0; i < l; i++) {
    s += alleleVect[i];
  }
  return s;
}

static int enoughFemale() {
  int s = 0, l = nAlleles;
  for (int i = 0; i < l; i++) {
    s += Rarray[i] - alleleVect[i];
  }
  if (s / 2 == female)
    return 1;
  else
    return 0;
}

static bool nearlyEqual(double a, double b, double epsilon){
  double absA, absB, diff;
  bool y;
  absA = abs(a);
  absB = abs(b);
  diff = abs(a - b);
  if (a == b) {
    y = true;
  } else if (a == 0 || b == 0 || (absA + absB < std::numeric_limits<double>::epsilon())) {
    y = diff < (epsilon*std::numeric_limits<double>::epsilon());
  } else {
    y = diff/std::min((absA + absB), std::numeric_limits<double>::epsilon()) < epsilon;
  }
  return y;
}


void xChrom (int *rm,
             int *mf,
             int  *rk,
             double *robservedVals, // observed stats: LLR, Prob, U, X2
             double *rPvals, // computed P values: LLR, Prob, U, X2
             int  *rstatID, // which statistic to use for histogram (1-4)
             int  *rhistobins, // number of bins for histogram. (no histogram if 0)
             double *rhistobounds, // Two values indicating the range for histogram
             double  *rhistoData, // histogram data. length = histobounds.
             int  *rsafeSecs, // abort calculation after this many seconds
             double  *tables // the number of tables examined
)
{
  
  // Set up global variables used during recursion
  nAlleles = *rk;
  male = mf[0];     // ------------------- NEW-------------------------------------
  female = mf[1];

  pPr = 0;
  
  hProb = rhistoData;
  Rbytes = *rk * sizeof(COUNTTYPE);
  statID = *rstatID;
  timeLimit = *rsafeSecs;
  HN = *rhistobins;
  start = time(NULL);
  Rarray = R_Calloc(*rk * *rk * (*rk-1)/2, COUNTTYPE);
  alleleVect = R_Calloc(*rk * *rk * (*rk-1)/2, COUNTTYPE);
  for (int i = 0; i < nAlleles; i++) {
    Rarray[i] = rm[i];
    alleleVect[i] = rm[i]+1;
  }
  mi = rm-1; // 1-based list of allele counts------------------------------
  tableCount = 0;

  // Make lookup tables
  lnFact = R_Calloc(rm[0] + 1, double);
  
  
  lnFact[0] = 0;
  double lni=0;
  for (int i = 1; i <= rm[0]; i++) {
    lni = log((double) i);
    lnFact[i] = lnFact[i-1] + lni;
  }

  int  nGenes = 0;
  for(int i = 0; i < nAlleles; i++) nGenes += rm[i];
  ntotal = nGenes/2;

  // Get constant terms for LLR and Prob
  constProbTerm =  0; // constLLRterm = 0;
  for (int i = 0; i < nAlleles; i++) {
    //---------------OLD--------------------
    constProbTerm +=  lgammafn(rm[i] + 1); //lnFact[rm[i]];

  }
  //--------------------------- OLD----------------
  //constProbTerm += log(2)*ntotal + lgammafn(ntotal+1) - lgammafn(nGenes +1);
  int nt = male + 2*female;
  constProbTerm +=  lgammafn(male+1) + lgammafn(female+1) - lgammafn(nt+1);
  
  // Get cutoffs for the four test statistics
  double oneMinus = 0.9999999; // Guards against floating-point-equality-test errors
  if(robservedVals[0] > 0.000000000001) robservedVals[0] = 0; // positive values are rounding errors
  maxLLR = robservedVals[0] * oneMinus;
  maxlPr = log(robservedVals[1]) * oneMinus;
  minmaxU = robservedVals[2] * oneMinus;
  minX2 = robservedVals[3] * oneMinus;
  
  
  start = time(NULL);

  recursiveEnumeration(1);
  *tables = tableCount;
  rPvals[0] = pLLR;
  rPvals[1] = pPr;
  rPvals[2] = pU;
  rPvals[3] = pX2;

  
  //Free(xlnx);
  R_Free(lnFact);
  R_Free(Rarray);
  R_Free(alleleVect);
  //Free(exa); Free(uTerm1); Free(uTerm2);
  //Free(x211); Free(x221); Free(x222);
  
  
}

// [[Rcpp::export]]
double xChromosomal(IntegerVector rmV,
                       IntegerVector mfV,
                       int  rk,
                       NumericVector robservedValsV, // observed stats: LLR, Prob, U, X2
                       NumericVector rPvalsV, // computed P values: LLR, Prob, U, X2
                       int  rstatID, // which statistic to use for histogram (1-4)
                       int  rhistobins, // number of bins for histogram. (no histogram if 0)
                       NumericVector rhistoboundsV, // Two values indicating the range for histogram
                       double  rhistoData, // histogram data. length = histobounds.
                       int  rsafeSecs, // abort calculation after this many seconds
                       double  tables // the number of tables examined
)
{
  
  // from NumericVector to double *
  // .... call xChrom
  int i;
  int * rm;
  int * mf;
  double * robservedVals; // observed stats: LLR, Prob, U, X2
  double * rPvals; // computed P values: LLR, Prob, U, X2
  double * rhistobounds; // Two values indicating the range for histogram
  double result;
  
  rm = (int *)calloc(rmV.length(), sizeof(int));
  mf = (int *)calloc(mfV.length(), sizeof(int));
  robservedVals = (double *)calloc(robservedValsV.length(), sizeof(double));
  rPvals = (double *)calloc(rPvalsV.length(), sizeof(double));
  rhistobounds = (double *)calloc(rhistoboundsV.length(), sizeof(double));
  
  for(i=0; i<rmV.length();i++){
    rm[i] = rmV[i];
  }
  for(i=0; i<mfV.length();i++){
    mf[i] = mfV[i];
  }
  for(i=0; i<robservedValsV.length();i++){
    robservedVals[i] = robservedValsV[i];
  }
  for(i=0; i<rPvalsV.length();i++){
    rPvals[i] = rPvalsV[i];
  }
  for(i=0; i<rhistoboundsV.length();i++){
    rhistobounds[i] = rhistoboundsV[i];
  }
  
  xChrom (rm, mf, &rk, robservedVals, rPvals, &rstatID, &rhistobins, rhistobounds, &rhistoData, &rsafeSecs, &tables);

  for(i=0; i<rPvalsV.length();i++){
    rPvalsV[i] = rPvals[i];
  }

  free(rm);
  free(mf);
  free(robservedVals);
  free(rhistobounds);

  result = rPvals[1];
  free(rPvals);
  
  //return DataFrame::create(Named("rPvals")=rPvalsV);
  return result;
}

