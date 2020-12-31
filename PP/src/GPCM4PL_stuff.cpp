#include <Rcpp.h>
using namespace Rcpp;


// 4PL Stuff ***************************************************************************************


// P FUNCTION +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// [[Rcpp::export]]
NumericVector P_4pl(NumericVector delta, double alpha, double theta, double la, double ua) {
// resp kann man eigenltich weglassen
// hier wird nur die Wahrscheinlichkeit berechnet zu loesen

//int nthres = delta.size();
//double nenner = 0;
//double zae = 0;
NumericVector PP1I(3);

double beta = delta(1); // weil der erste muss 0 sein, weil ja auch in der thres matrix GPCM items drinstehen koennen

// probability
PP1I(0) = la + (ua - la) * exp(alpha*(theta - beta))/(1+exp(alpha*(theta - beta)));

// first derivate
PP1I(1) = alpha * (ua - PP1I(0)) * (PP1I(0) - la) / (ua - la);

// information
PP1I(2) = (PP1I(1)*PP1I(1))/(PP1I(0)*(1-PP1I(0)));

return PP1I;
}



// P FUNCTION +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// [[Rcpp::export]]
NumericVector P_4pl4wle(NumericVector delta, double alpha, double theta, double la, double ua) {
// resp kann man eigenltich weglassen
// hier wird nur die Wahrscheinlichkeit berechnet zu loesen

//int nthres = delta.size();
//double nenner = 0;
//double zae = 0;
NumericVector PP1I(6);

double beta = delta(1); // weil der erste muss 0 sein, weil ja auch in der thres matrix GPCM items drinstehen koennen

// probability
PP1I(0) = la + (ua - la) * exp(alpha*(theta - beta))/(1+exp(alpha*(theta - beta)));

// first derivate
PP1I(1) = alpha * (ua - PP1I(0)) * (PP1I(0) - la) / (ua - la);

// information
PP1I(2) = (PP1I(1)*PP1I(1))/(PP1I(0)*(1-PP1I(0)));

// second deriv
PP1I(3) = alpha / (ua-la) * (ua*PP1I(1) - 2*PP1I(0)*PP1I(1)  + la*PP1I(1));

// third deriv
PP1I(4) = alpha / (ua-la) * (ua*PP1I(3) - 2*PP1I(1)*PP1I(1) - 2*PP1I(0)*PP1I(3) + la*PP1I(3));

// information first deriv

double empp1i = 1-PP1I(0);

double oben  = (2*PP1I(1)*PP1I(3)*PP1I(0)*(1-PP1I(0)) - PP1I(1)*PP1I(1)*(PP1I(1)*(1-PP1I(0)) + PP1I(0)*(-PP1I(1))));
double unten = PP1I(0)*PP1I(0)*empp1i*empp1i;

double INF2 = oben/unten;

//PP1I(5) = (PP1I(2) * PP1I(1) * PP1I(3) * PP1I(1) - PP1I(0) * (PP1I(1) * (PP1I(2) * PP1I(4) - PP1I(3) * INF2) + PP1I(2) * PP1I(3) *PP1I(3)))/pow(PP1I(0),2);
double Qj = 1-PP1I(0);
double PP1IQj = PP1I(0) * Qj;
// J'
//PP1I(5) = ((PP1I(3)*PP1I(3) + PP1I(1)*PP1I(4))*PP1I(0)*Qj - PP1I(1)*PP1I(3)*(PP1I(1)*Qj - PP1I(0)*PP1I(1))) / pow(PP1I(0) * Qj,2);

// diese version stimmt mit dem alten PP package ueberein, hat aber mMn einen vorzeichenfehler

double J1 = ((PP1I(3)*PP1I(3) - PP1I(1)*PP1I(4))*PP1I(0)*Qj + PP1I(1)*PP1I(3)*(PP1I(1)*Qj - PP1I(0)*PP1I(1))) / (PP1IQj*PP1IQj);

double J = PP1I(1) * PP1I(3) / (PP1I(0)*Qj);

// this is the whole numerator of the first derivtae of the correction term 
PP1I(5) = J1 * PP1I(2) - J * INF2;


//PP1I(6) = INF2; // i hope we wont need this anymore

return PP1I;
}



// [[Rcpp::export]]
double r_huber_4pl(NumericVector delta, double alpha,
                          double theta, double la, double ua, double H) {

//double nenner = 0;
//double zae = 0;
double HU = 0;
double beta = delta(1); // weil der erste muss 0 sein, weil ja auch in der thres matrix GPCM items drinstehen koennen

// compute residuum
//double resid = log(la + (ua - la) * exp(alpha*(theta - beta)));
double resid = alpha*(theta - beta);
// compute Huber weight

if(resid < 0)
  {
  resid = resid * (-1);  
  }


if(resid <= H)
  {
  HU = 1;  
  } else 
    {
    HU = H/resid;
    }

return HU;
}


// LIKELIHOOD - L1, L2 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// [[Rcpp::export]]
NumericMatrix L4pl(IntegerMatrix awm, NumericMatrix DELTA, NumericVector ALPHA, 
                   NumericVector LOWA, NumericVector UPPA, NumericVector THETA, bool map, 
                   NumericVector mu, NumericVector sigma2) {

int npers = awm.nrow();
int nitem = awm.ncol();
int maxca = DELTA.nrow();
// ALPHA muss so lang sein wie nitem anzeigt
// THETA muss so lang sein wie npers anzeigt
// DELTA muss so viele Spalten haben wie nitem anzeigt
// maxca gibt an wieviele maximale kategorien


// 4 columns: 1st deriv of logL; 2nd deriv of logL; delta = 1st/2nd; theta - delta
NumericMatrix l1l2M(npers,4);

for(int it = 0; it < nitem; it++)
  {
    
  // Response vector of one item  
  IntegerVector respvec = awm(_,it);
  double alpha = ALPHA(it);
  NumericVector delta = DELTA(_,it);
  LogicalVector nas(maxca);
  double lowerA = LOWA(it);
  double upperA = UPPA(it);
  
  // find NA and kill them
  for (int fna = 0; fna < maxca; ++fna) {
    nas[fna] = NumericVector::is_na(delta[fna]);
  }
  // parameters without missing values. missing values should only be possible at the end of the matrix
  NumericVector delta1 = delta[!nas];
  
  //int kmax = delta1.size();
  
  for(int pe = 0; pe < npers; pe++)
    {
     
    int resp = respvec(pe);
    double theta = THETA(pe);
    NumericVector ergP(3);
    
    
    // NA handling 

    // if the i,j obs is NA, add nothing
    if(IntegerVector::is_na(resp))
    { // in case of missing value as response
    continue;
    } else 
      {
        
     ergP = P_4pl(delta1, alpha, theta, lowerA, upperA);  
      
      // l1
      double Qj = 1 - ergP(0);
      l1l2M(pe,0) += (resp - ergP(0))/(ergP(0)*Qj) * ergP(1);
      l1l2M(pe,1) +=  ergP(2);

      }

    }
    
  }

if(map)
  {
    NumericVector corrterm1(npers);
    NumericVector corrterm2(npers);
    
    corrterm1 = (THETA - mu)/sigma2;
    corrterm2 = 1/sigma2;
    
    l1l2M(_,1) = l1l2M(_,1) * (-1);
    l1l2M(_,2) = (l1l2M(_,0) - corrterm1) / (l1l2M(_,1)-corrterm2);
    
    
    // Korrekturterm um divergieren zu verhindern:
    for(int ai = 0; ai<npers; ai++)
      {
        if(std::abs(l1l2M(ai,2)) <= 5) 
        {
          continue;
        } else {
          l1l2M(ai,2) = l1l2M(ai,2)/std::abs(l1l2M(ai,2)) * 5;
        }
        
      }
      
    
    l1l2M(_,3) = THETA - l1l2M(_,2); 
    
  } else {
    
          l1l2M(_,1) = l1l2M(_,1) * (-1);
          l1l2M(_,2) = l1l2M(_,0)/l1l2M(_,1);
          
          
          for(int ai = 0; ai<npers; ai++)
            {
              if(std::abs(l1l2M(ai,2)) <= 5) 
              {
                continue;
              } else {
                l1l2M(ai,2) = l1l2M(ai,2)/std::abs(l1l2M(ai,2)) * 5;
              }
              
            }
            
          
          
          l1l2M(_,3) = THETA - l1l2M(_,2);
    
         }

return l1l2M;


}




// LIKELIHOOD - L1, L2 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// [[Rcpp::export]]
NumericMatrix L4pl_wle(IntegerMatrix awm, NumericMatrix DELTA, NumericVector ALPHA,
                NumericVector LOWA, NumericVector UPPA, NumericVector THETA) {

int npers = awm.nrow();
int nitem = awm.ncol();
int maxca = DELTA.nrow();
// ALPHA muss so lang sein wie nitem anzeigt
// THETA muss so lang sein wie npers anzeigt
// DELTA muss so viele Spalten haben wie nitem anzeigt
// maxca gibt an wieviele maximale kategorien


// 4 columns: 1st deriv of logL; 2nd deriv of logL; delta = 1st/2nd; theta - delta
NumericMatrix l1l2M(npers,6);

for(int it = 0; it < nitem; it++)
  {
    
  // Response vector of one item  
  IntegerVector respvec = awm(_,it);
  double alpha = ALPHA(it);
  NumericVector delta = DELTA(_,it);
  LogicalVector nas(maxca);
  double lowerA = LOWA(it);
  double upperA = UPPA(it);
  
  // find NA and kill them
  for (int fna = 0; fna < maxca; ++fna) {
    nas[fna] = NumericVector::is_na(delta[fna]);
  }
  // parameters without missing values. missing values should only be possible at the end of the matrix
  NumericVector delta1 = delta[!nas];
  
  //int kmax = delta1.size();
  
  for(int pe = 0; pe < npers; pe++)
    {
     
    int resp = respvec(pe);
    double theta = THETA(pe);
    //NumericVector ergP(3);
    
    
    // NA handling 

    // if the i,j obs is NA, add nothing
    if(IntegerVector::is_na(resp))
    { // in case of missing value as response
    continue;
    } else 
      {
        
     NumericVector ergP = P_4pl4wle(delta1, alpha, theta, lowerA, upperA);  
      
      // l1
      double Qj = 1 - ergP(0);
      l1l2M(pe,0) += (resp - ergP(0))/(ergP(0)*Qj) * ergP(1);
      l1l2M(pe,1) +=  ergP(2);
      
      // oberer korrekturterm
      l1l2M(pe,2) += ergP(1) * ergP(3) / (ergP(0)*Qj); // ersetzen?!
      // unterer korrekturterm
      l1l2M(pe,3) += ergP(5);
      // I' - w wont need this
      //l1l2M(pe,4) += ergP(6);
      }

    }
    
  }


//std::cout << "J' 1 = " << l1l2M(0,3) <<  std::endl ;
//std::cout << "J' 2 = " << l1l2M(1,3) <<  std::endl ;
//
// I * J' - I' * J
//NumericVector corru = (l1l2M(_,1)*l1l2M(_,3) - l1l2M(_,4)*l1l2M(_,2)) / (2 * pow(l1l2M(_,1),2));

NumericVector corru = l1l2M(_,3) / (2 * l1l2M(_,1)*l1l2M(_,1));

l1l2M(_,4) = (l1l2M(_,0) + l1l2M(_,2) / (2*l1l2M(_,1))) / (l1l2M(_,1) + corru);

// Korrekturterm um divergieren zu verhindern:
  for(int ai = 0; ai<npers; ai++)
    {
      if(std::abs(l1l2M(ai,4)) <= 5) 
        {
          continue;
        } else {
                l1l2M(ai,4) = l1l2M(ai,4)/std::abs(l1l2M(ai,4)) * 5;
               }
      
    }
    


l1l2M(_,5) = THETA + l1l2M(_,4);
//l1l2M(_,1) = l1l2M(_,1) * (-1); // erst jetzt * (-1)

return l1l2M;


}





// [[Rcpp::export]]
NumericMatrix L4pl_robust(IntegerMatrix awm, NumericMatrix DELTA, NumericVector ALPHA, 
                          NumericVector LOWA, NumericVector UPPA,
                          NumericVector THETA, double H) {

int npers = awm.nrow();
int nitem = awm.ncol();
int maxca = DELTA.nrow();
// ALPHA muss so lang sein wie nitem anzeigt
// THETA muss so lang sein wie npers anzeigt
// DELTA muss so viele Spalten haben wie nitem anzeigt
// maxca gibt an wieviele maximale kategorien


// 4 columns: 1st deriv of logL; 2nd deriv of logL; delta = 1st/2nd; theta - delta
NumericMatrix l1l2M(npers,4);

for(int it = 0; it < nitem; it++)
  {
    
  // Response vector of one item  
  IntegerVector respvec = awm(_,it);
  double alpha = ALPHA(it);
  NumericVector delta = DELTA(_,it);
  LogicalVector nas(maxca);
  double lowerA = LOWA(it);
  double upperA = UPPA(it);
  
  // find NA and kill them
  for (int fna = 0; fna < maxca; ++fna) {
    nas[fna] = NumericVector::is_na(delta[fna]);
  }
  // parameters without missing values. missing values should only be possible at the end of the matrix
  NumericVector delta1 = delta[!nas];
  
  //int kmax = delta1.size();
  
  for(int pe = 0; pe < npers; pe++)
    {
     
    int resp = respvec(pe);
    double theta = THETA(pe);
    NumericVector ergP(3);
    
    
    // NA handling 

    // if the i,j obs is NA, add nothing
    if(IntegerVector::is_na(resp))
    { // in case of missing value as response
    continue;
    } else 
      {
        
     ergP = P_4pl(delta1, alpha, theta, lowerA, upperA);
     // compute huber weight
     //NumericVector hub = r_huber_4pl(delta=delta1,alpha, theta, lowerA, upperA, H);
     double hub = r_huber_4pl(delta=delta1,alpha, theta, lowerA, upperA, H); 
      // l1
      double Qj = 1 - ergP(0);
      // huber weighted first deriv and Inf
      
      l1l2M(pe,0) += hub*(resp - ergP(0))/(ergP(0)*Qj) * ergP(1);
      //l1l2M(pe,0) += hub(1)*(resp - ergP(0))/(hub(1)*(ergP(0)*Qj)) * ergP(1);
      //l1l2M(pe,1) +=  ergP(2) * hub(1);
      l1l2M(pe,1) +=  ergP(2);

      }

    }
    
  }
  
l1l2M(_,1) = l1l2M(_,1) * (-1);
//l1l2M(_,1) = l1l2M(_,1);
l1l2M(_,2) = l1l2M(_,0)/l1l2M(_,1);

// Korrekturterm um divergieren zu verhindern:
for(int ai = 0; ai<npers; ai++)
  {
    if(std::abs(l1l2M(ai,2)) <= 5) 
    {
      continue;
    } else {
      l1l2M(ai,2) = l1l2M(ai,2)/std::abs(l1l2M(ai,2)) * 5;
    }
    
  }

l1l2M(_,3) = THETA - l1l2M(_,2);

//std::cout << "l1 = " << l1l2M(6,0) <<  std::endl ;

return l1l2M;


}



// NR - Algorithm --->>>  MLE + WLE + MAP + robust <<<<---- +++++++++++++++++++++++++++++++++++++++++++++++

// [[Rcpp::export]]
List NR_4PL(IntegerMatrix awm, NumericMatrix DELTA, NumericVector ALPHA, 
NumericVector LOWA, NumericVector UPPA, NumericVector THETA, String wm, 
int maxsteps, double exac, NumericVector mu, NumericVector sigma2, double H) {

int npers = awm.nrow();

NumericMatrix resPP(npers,2);
int howlong;


if(wm == "wle")
  {
    
  for(int newr = 0; newr < maxsteps; newr++)
    {
    NumericMatrix reso = L4pl_wle(awm,DELTA,ALPHA,LOWA,UPPA,THETA);
    THETA = reso(_,5);
    
    NumericVector diffs = reso(_,4);
    LogicalVector bxy = is_na(THETA);
    NumericVector diffs1 = diffs[!bxy];

    if( (is_true(all(abs(diffs1) < exac))) | (newr == (maxsteps-1)))
      {
        resPP(_,0) = THETA;
        //resPP(_,1) = pow(1/reso(_,1),0.5); // thank you solaris
        resPP(_,1) = 1/reso(_,1);
        howlong = newr + 1;
        break;
      }
  
    }
    
  } else if(wm == "mle")
    {
   bool map = FALSE;
   

    for(int newr = 0; newr < maxsteps; newr++)
      {
      NumericMatrix reso = L4pl(awm,DELTA,ALPHA,LOWA,UPPA,THETA,map,mu,sigma2);
      THETA = reso(_,3);
      
      NumericVector diffs = reso(_,2);
      LogicalVector bxy = is_na(THETA);
      NumericVector diffs1 = diffs[!bxy];
      
      if( (is_true(all(abs(diffs1) < exac))) | (newr == (maxsteps-1)))
        {
          resPP(_,0) = THETA;
          resPP(_,1) = 1/reso(_,1)*(-1);
          howlong = newr + 1;
          break;
        }
    
      }  

    } else if(wm == "map")
      {  
       bool map = TRUE; 
          for(int newr = 0; newr < maxsteps; newr++)
            {
            NumericMatrix reso = L4pl(awm,DELTA,ALPHA,LOWA,UPPA,THETA,map,mu,sigma2);
            THETA = reso(_,3);
            
            NumericVector diffs = reso(_,2);
            LogicalVector bxy = is_na(THETA);
            NumericVector diffs1 = diffs[!bxy];
            
            if( (is_true(all(abs(diffs1) < exac))) | (newr == (maxsteps-1)))
              {
                resPP(_,0) = THETA;
                resPP(_,1) = 1/reso(_,1)*(-1);
                howlong = newr + 1;
                break;
              }
          
            }   
        
      } else if(wm == "robust")
        {
            
            for(int newr = 0; newr < maxsteps; newr++)
              {
              NumericMatrix reso = L4pl_robust(awm,DELTA,ALPHA,LOWA,UPPA,THETA,H);
              THETA = reso(_,3);
              
              NumericVector diffs = reso(_,2);
              LogicalVector bxy = is_na(THETA);
              NumericVector diffs1 = diffs[!bxy];
              
              if( (is_true(all(abs(diffs1) < exac))) | (newr == (maxsteps-1)))
                {
                  resPP(_,0) = THETA;
                  resPP(_,1) = 1/reso(_,1)*(-1);
                  howlong = newr + 1;
                  break;
                }
            
              }    
          
        }

// ----
return List::create(_["resPP"] = resPP, _["nsteps"] = howlong);

}




// GPCM Stuff ***************************************************************************************




// P FUNCTION +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// [[Rcpp::export]]
double P_gpcm(NumericVector delta, double alpha, double theta, int resp) {

int nthres = delta.size();
double nenner = 0;
double zae = 0;

for(int aus = 0; aus < nthres; aus++)
  {
  
  double imex = 0;
    for(int i = 0; i < aus+1; i++)
        {
           
        imex += alpha * (theta - delta(i));
        }
  
  nenner += exp(imex);
  }
      
 
      
 for(int cat = 0; cat < resp+1; cat++)
  {
  zae += alpha * (theta - delta(cat));
  }

double P = exp(zae) / nenner;

return P;
}



// [[Rcpp::export]]
double r_huber_gpcm(NumericVector delta, double alpha,
                           double theta, double H) {

int nthres = delta.size();
//double zae = 0;
int nthresm1 = nthres - 1;
//double deltaw0 = delta(-0);
double resid = 0;
double HU = 0;
// compute residuum

// dont start with the first - because this is always zero!
for(int ru = 1; ru < nthres; ru++)
  {
  resid += alpha * (theta - delta(ru))/nthresm1;
  }

// compute Huber weight

if(resid < 0)
  {
  resid = resid * (-1);  
  }


if(resid <= H)
  {
  HU = 1;  
  } else 
    {
    HU = H/resid;
    }

//std::cout << "hu = " << HU <<  std::endl ;

return HU;
}






// P FUNCTION L1, L2 +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// [[Rcpp::export]]
NumericMatrix L12gpcm(IntegerMatrix awm, NumericMatrix DELTA, NumericVector ALPHA, 
                      NumericVector THETA, NumericVector mu, NumericVector sigma2, bool map) {
// awm = antwortmatrix
// 
int npers = awm.nrow();
int nitem = awm.ncol();
int maxca = DELTA.nrow();
// ALPHA muss so lang sein wie nitem anzeigt
// THETA muss so lang sein wie npers anzeigt
// DELTA muss so viele Spalten haben wie nitem anzeigt
// maxca gibt an wieviele maximale kategorien

// l1 fuer die personen
//NumericVector l1vec(npers);
//NumericVector l2vec(npers);
// 4 columns: 1st deriv of logL; 2nd deriv of logL; delta = 1st/2nd; theta - delta
NumericMatrix l1l2M(npers,4);

for(int it = 0; it < nitem; it++)
  {
    
  // Response vector of one item  
  IntegerVector respvec = awm(_,it);
  double alpha = ALPHA(it);
  NumericVector delta = DELTA(_,it);
  LogicalVector nas(maxca);
  
  // find NA and kill them
  for (int fna = 0; fna < maxca; ++fna) {
    nas[fna] = NumericVector::is_na(delta[fna]);
  }
  // parameters without missing values. missing values should only be possible at the end of the matrix
  NumericVector delta1 = delta[!nas];
  
  int kmax = delta1.size();
  
  
  
  for(int pe = 0; pe < npers; pe++)
    {
     
    int resp = respvec(pe);
    double theta = THETA(pe);
     
    // NA handling 

    // if the i,j obs is NA, add nothing
    if(IntegerVector::is_na(resp))
    {
    continue;
    // jumps to next iteration step
//    l1l2M(pe,0) += 0;
//    l1l2M(pe,1) += 0;

    } else 
      {
        
      // rattert die ks durch
      
      double rs = 0;
      double rs2 = 0;
      double ls2 = 0;
      double ergP = 0;
      
      for(int ks = 0; ks < kmax; ks++)
        {
      ergP = P_gpcm(delta1, alpha, theta, ks);
      
      rs += ks * alpha * ergP;
      // second derivates right and left side
      rs2 += ks * alpha * ergP;
      //ls2 += pow(ks,2) * pow(alpha,2) * ergP;
      ls2 += ks*ks * alpha*alpha * ergP;
        }
      rs2 = rs2*rs2;
      
      // write first and second derivs in 2 column matrix - for each person
      l1l2M(pe,0) += resp * alpha - rs;
      l1l2M(pe,1) += ls2 - rs2;
      }
//    double showme = l1l2M(pe,0); // weg
//    std::cout << "l1 = " << showme <<  std::endl ; // weg
    }
    
  }


if(map)
  {
  NumericVector corrterm1(npers);
  NumericVector corrterm2(npers);
  
  corrterm1 = (THETA - mu)/sigma2;
  corrterm2 = 1/sigma2;
  
  l1l2M(_,1) = l1l2M(_,1) * (-1);
  l1l2M(_,2) = (l1l2M(_,0) - corrterm1) / (l1l2M(_,1)-corrterm2);
  
  // Korrekturterm um divergieren zu verhindern:
  for(int ai = 0; ai<npers; ai++)
    {
      if(std::abs(l1l2M(ai,2)) <= 5) 
      {
        continue;
      } else {
        l1l2M(ai,2) = l1l2M(ai,2)/std::abs(l1l2M(ai,2)) * 5;
      }
      
    }
  
  
  l1l2M(_,3) = THETA - l1l2M(_,2); 
  
  } else 
    {
      
    l1l2M(_,1) = l1l2M(_,1) * (-1);
    l1l2M(_,2) = l1l2M(_,0)/l1l2M(_,1);
    
    // Korrekturterm um divergieren zu verhindern:
    for(int ai = 0; ai<npers; ai++)
      {
        if(std::abs(l1l2M(ai,2)) <= 5) 
        {
          continue;
        } else {
          l1l2M(ai,2) = l1l2M(ai,2)/std::abs(l1l2M(ai,2)) * 5;
        }
        
      }
    
    
    l1l2M(_,3) = THETA - l1l2M(_,2);
      
    }


return l1l2M;

}



// CORRECTION TERM FUNCTION +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// [[Rcpp::export]]
NumericVector Pcorr1_gpcm(NumericVector delta, double alpha, double theta) {

int nthres = delta.size();
// double which will contain the sum of weighted probs
double overac = 0;
double overac2 = 0;
double overac3 = 0;

NumericVector  ps(nthres);
NumericVector fds(nthres);
NumericVector sds(nthres);
NumericVector tds(nthres);


// Probs for each category  
for(int i = 0; i < nthres; i++)
  {
  ps(i) = P_gpcm(delta, alpha, theta, i);
  }
  

for(int i = 0; i < nthres; i++)
  {
  //overac += i * P_gpcm(delta, alpha, theta, i); 
  overac += i * ps(i);
  }

// first derivs of P for each category

for(int i = 0; i < nthres; i++)
  {
  fds(i) = alpha * P_gpcm(delta, alpha, theta, i) * (i - overac);
  }


// second derivates of P for each category
for(int i = 0; i < nthres; i++)
  {
  //overac2 += i * alpha * P_gpcm(delta, alpha, theta, i) * (i - overac);
  overac2 += i * fds(i);
  }
  
for(int i = 0; i < nthres; i++)
  {
  //sds(i) = alpha * fds(i) * (i - overac) + alpha * P_gpcm(delta, alpha, theta, i) * overac2;
  sds(i) = alpha * fds(i) * (i - overac) - alpha * P_gpcm(delta, alpha, theta, i) * overac2;
  }
  
// information
  NumericVector corrts(3); //corrterm1 und INF
  
  double INF = alpha * overac2;
  //double INF = sum(-sds * ps);
  //double INF = sum(pow(fds,2)/ps);

  //NumericVector corrts(2);
  
  // compute correction term  
  //double corrts = sum(fds * sds / ps)/(2*INF);
  corrts(0) = sum(fds * sds / ps);
  corrts(1) = INF;
  
// *******************************
// compute the first derivates of the correction term
// *******************************

// P - third deriv

for(int i = 0; i < nthres; i++)
  {
  overac3 += i * sds(i);
  }


for(int i = 0; i < nthres; i++)
  {
  //tds(i) = alpha * sds(i) * (i - overac) + alpha*fds(i)*overac2 + alpha*fds(i) * overac2 + alpha*ps(i) * overac3;
  tds(i) = alpha * sds(i) * (i - overac) - alpha*fds(i)*overac2 - alpha*fds(i) * overac2 + alpha*ps(i) * overac3;
  }

  
  double INF2 = alpha * overac3;
 
corrts(2) = sum((INF * fds * sds * fds - ps * (fds * (INF * tds - sds * INF2) + INF * sds *sds))/(ps*ps));


return corrts;


}



// P FUNCTION L1, L2 --->>>  WLE  <<<<---- +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// [[Rcpp::export]]
NumericMatrix L12gpcm_wle(IntegerMatrix awm, NumericMatrix DELTA, NumericVector ALPHA,
                          NumericVector THETA) {
// awm = antwortmatrix
// 
int npers = awm.nrow();
int nitem = awm.ncol();
int maxca = DELTA.nrow();
// ALPHA muss so lang sein wie nitem anzeigt
// THETA muss so lang sein wie npers anzeigt
// DELTA muss so viele Spalten haben wie nitem anzeigt
// maxca gibt an wieviele maximale kategorien

// l1 fuer die personen
//NumericVector l1vec(npers);
//NumericVector l2vec(npers);
// 
NumericMatrix l1l2M(npers,6);
NumericVector ct(npers); //weg

for(int it = 0; it < nitem; it++)
  {
    
  // Response vector of one item  
  IntegerVector respvec = awm(_,it);
  double alpha = ALPHA(it);
  NumericVector delta = DELTA(_,it);
  LogicalVector nas(maxca);
  
  // find NA and kill them
  for (int fna = 0; fna < maxca; fna++) {
    nas[fna] = NumericVector::is_na(delta[fna]);
  }
  // parameters without missing values. missing values should only be possible at the end of the matrix
  NumericVector delta1 = delta[!nas];
  
  int kmax = delta1.size();
  
  
  
  for(int pe = 0; pe < npers; pe++)
    {
     
    int resp = respvec(pe);
    double theta = THETA(pe);
     
    // NA handling 

    // if the i,j obs is NA, add nothing
    if(IntegerVector::is_na(resp))
    {
    continue;
    // jumps to next iteration step
//    l1l2M(pe,0) += 0;
//    l1l2M(pe,1) += 0;

    } else 
      {
        
      // rattert die ks durch
      
      double rs = 0;
      double rs2 = 0;
      double ls2 = 0;
      double ergP = 0;
      
      for(int ks = 0; ks < kmax; ks++)
        {
      ergP = P_gpcm(delta1, alpha, theta, ks);
      
      rs += ks * alpha * ergP;
      // second derivates right and left side
      rs2 += ks * alpha * ergP;
      //ls2 += pow(ks,2) * pow(alpha,2) * ergP;
      ls2 += ks*ks * alpha*alpha * ergP;
        }
      rs2 = rs2*rs2;
      
      // Pcorr2_gpcm(NumericVector delta, double alpha, double theta) 
      
      
      // ********** VERSION 1

      NumericVector ct = Pcorr1_gpcm(delta1,alpha,theta);
      l1l2M(pe,1) += ct(1); // summierte information
      l1l2M(pe,2) += ct(0); // summiert korrekturterm - korrekturterm oben
      //l1l2M(pe,3) += ct(1); // summierte information
      l1l2M(pe,3) += ct(2); // unterer korrekturterm 
      // *********************
      

      // ********** VERSION 2

//      double eps = 0.0001;
//      NumericVector ct = Pcorr1_gpcm(delta1,alpha,theta);
//      l1l2M(pe,2) += ct(0); // summiert korrekturterm - korrekturterm oben
//      l1l2M(pe,3) += ct(1); // summierte information
//      
//      //corrterms(0) = Pcorr1_gpcm(delta1,alpha,theta);
//      NumericVector crpep = Pcorr1_gpcm(delta1,alpha,theta+eps); //approx.
//      l1l2M(pe,4) += (crpep(0) - ct(0))/eps; // correkturterm unten
      
      // *********************

      l1l2M(pe,0) += (resp * alpha - rs);
      //l1l2M(pe,1) += (-1)*(ls2 - rs2); 
      }

//    std::cout << "l1 = " << showme <<  std::endl ; // weg
    }
    
  }


l1l2M(_,3) = l1l2M(_,3)/(2*l1l2M(_,1)*l1l2M(_,1));

//l1l2M(_,5) = (l1l2M(_,0) + l1l2M(_,2)/(2*l1l2M(_,3))) / (l1l2M(_,1) + l1l2M(_,4)); // here
//l1l2M(_,5) = (l1l2M(_,0) + l1l2M(_,2)/(2*l1l2M(_,3))) / (l1l2M(_,3)*(-1) + l1l2M(_,4)); // here

l1l2M(_,4) = (l1l2M(_,0) + l1l2M(_,2)/(2*l1l2M(_,1))) / (l1l2M(_,1)*(-1) + l1l2M(_,3));

// Korrekturterm um divergieren zu verhindern:
for(int ai = 0; ai<npers; ai++)
  {
    if(std::abs(l1l2M(ai,4)) <= 5) 
    {
      continue;
    } else {
      l1l2M(ai,4) = l1l2M(ai,4)/std::abs(l1l2M(ai,4)) * 5;
    }
    
  }


l1l2M(_,5) = THETA - l1l2M(_,4);



return l1l2M;

}




// [[Rcpp::export]]
NumericMatrix L12gpcm_robust(IntegerMatrix awm, NumericMatrix DELTA, NumericVector ALPHA, 
                      NumericVector THETA, double H) {
// awm = antwortmatrix
// 
int npers = awm.nrow();
int nitem = awm.ncol();
int maxca = DELTA.nrow();
// ALPHA muss so lang sein wie nitem anzeigt
// THETA muss so lang sein wie npers anzeigt
// DELTA muss so viele Spalten haben wie nitem anzeigt
// maxca gibt an wieviele maximale kategorien

// l1 fuer die personen
//NumericVector l1vec(npers);
//NumericVector l2vec(npers);
// 4 columns: 1st deriv of logL; 2nd deriv of logL; delta = 1st/2nd; theta - delta
NumericMatrix l1l2M(npers,4);

for(int it = 0; it < nitem; it++)
  {
    
  // Response vector of one item  
  IntegerVector respvec = awm(_,it);
  double alpha = ALPHA(it);
  NumericVector delta = DELTA(_,it);
  LogicalVector nas(maxca);
  
  // find NA and kill them
  for (int fna = 0; fna < maxca; ++fna) {
    nas[fna] = NumericVector::is_na(delta[fna]);
  }
  // parameters without missing values. missing values should only be possible at the end of the matrix
  NumericVector delta1 = delta[!nas];
  
  int kmax = delta1.size();
  
  
  
  for(int pe = 0; pe < npers; pe++)
    {
     
    int resp = respvec(pe);
    double theta = THETA(pe);
     
    // NA handling 

    // if the i,j obs is NA, add nothing
    if(IntegerVector::is_na(resp))
    {
    continue;
    // jumps to next iteration step
//    l1l2M(pe,0) += 0;
//    l1l2M(pe,1) += 0;

    } else 
      {
        
      // rattert die ks durch
      
      double rs = 0;
      double rs2 = 0;
      double ls2 = 0;
      double ergP = 0;
      
      
      double hu = r_huber_gpcm(delta1,alpha,theta,H);
      
      for(int ks = 0; ks < kmax; ks++)
        {
      ergP = P_gpcm(delta1, alpha, theta, ks);
      
      rs += ks * alpha * ergP;
      // second derivates right and left side
      rs2 += ks * alpha * ergP;
      //ls2 += pow(ks,2) * pow(alpha,2) * ergP;
      ls2 += ks*ks * alpha * alpha * ergP;
        }
      rs2 = rs2* rs2;
      
      // write first and second derivs in 2 column matrix - for each person
      l1l2M(pe,0) += (resp * alpha - rs)*hu;
      l1l2M(pe,1) += ls2 - rs2;
      
      //std::cout << "oben = " << l1l2M(pe,0) <<  std::endl ;
      }
//    double showme = l1l2M(pe,0); // weg
//    std::cout << "l1 = " << showme <<  std::endl ; // weg
    }
    
  }

   
    l1l2M(_,1) = l1l2M(_,1) * (-1);
    l1l2M(_,2) = l1l2M(_,0)/l1l2M(_,1);
    
    // Korrekturterm um divergieren zu verhindern:
    for(int ai = 0; ai<npers; ai++)
      {
        if(std::abs(l1l2M(ai,2)) <= 5) 
        {
          continue;
        } else {
          l1l2M(ai,2) = l1l2M(ai,2)/std::abs(l1l2M(ai,2)) * 5;
        }
        
      }
      
    l1l2M(_,3) = THETA - l1l2M(_,2);
      


return l1l2M;

}



// NR - Algorithm --->>>  MLE + WLE  <<<<---- +++++++++++++++++++++++++++++++++++++++++++++++

// [[Rcpp::export]]
List NR_GPCM(IntegerMatrix awm, NumericMatrix DELTA, NumericVector ALPHA, NumericVector THETA,
             String wm, int maxsteps, double exac, NumericVector mu, NumericVector sigma2, double H) {

int npers = awm.nrow();

NumericMatrix resPP(npers,2);
int howlong;


if(wm == "wle")
  {
    
  for(int newr = 0; newr < maxsteps; newr++)
    {
    NumericMatrix reso = L12gpcm_wle(awm,DELTA,ALPHA,THETA);
    THETA = reso(_,5);
    
    if( (is_true(all(abs(reso(_,4)) < exac))) | (newr == (maxsteps-1)))
      {
        resPP(_,0) = THETA;
        resPP(_,1) = 1/reso(_,1);
        howlong = newr;
        break;
      }
  
    }
    
  } else if(wm == "mle")
    {
    bool map = FALSE;  
    for(int newr = 0; newr < maxsteps; newr++)
      {
      NumericMatrix reso = L12gpcm(awm,DELTA,ALPHA,THETA,mu,sigma2,map);
      THETA = reso(_,3);
      
      if( (is_true(all(abs(reso(_,2)) < exac))) | (newr == (maxsteps-1)))
        {
          resPP(_,0) = THETA;
          resPP(_,1) = 1/reso(_,1)*(-1);
          howlong = newr;
          break;
        }
    
      }  

    } else if(wm == "map")
      {  
        bool map = TRUE;
          for(int newr = 0; newr < maxsteps; newr++)
            {
            NumericMatrix reso = L12gpcm(awm,DELTA,ALPHA,THETA,mu,sigma2,map);
            THETA = reso(_,3);
            
            if( (is_true(all(abs(reso(_,2)) < exac))) | (newr == (maxsteps-1)))
              {
                resPP(_,0) = THETA;
                resPP(_,1) = 1/reso(_,1)*(-1);
                howlong = newr;
                break;
              }
          
            }   
        
      } else if(wm == "robust")
          {  
              for(int newr = 0; newr < maxsteps; newr++)
                {
                NumericMatrix reso = L12gpcm_robust(awm,DELTA,ALPHA,THETA,H);
                THETA = reso(_,3);
                //std::cout << "theta = " << resPP(1,0) <<  std::endl ;

                if( (is_true(all(abs(reso(_,2)) < exac))) | (newr == (maxsteps-1)))
                  {
                    resPP(_,0) = THETA;
                    resPP(_,1) = 1/reso(_,1)*(-1);
                    howlong = newr;
                    break;
                  }
              
                }   
            
          }

// ----
return List::create(_["resPP"] = resPP, _["nsteps"] = howlong);

}




// GPCM - 4PL mixed - MLE**************************************************


// [[Rcpp::export]]
NumericMatrix Lgpcm4pl_mle(IntegerMatrix awm, NumericMatrix DELTA, NumericVector ALPHA,
                           NumericVector LOWA, NumericVector UPPA, NumericVector THETA, 
                           CharacterVector model, NumericVector mu, NumericVector sigma2, bool map) {
// awm = antwortmatrix
// 
int npers = awm.nrow();
int nitem = awm.ncol();
int maxca = DELTA.nrow();
// ALPHA muss so lang sein wie nitem anzeigt
// THETA muss so lang sein wie npers anzeigt
// DELTA muss so viele Spalten haben wie nitem anzeigt
// maxca gibt an wieviele maximale kategorien

// l1 fuer die personen
//NumericVector l1vec(npers);
//NumericVector l2vec(npers);
// 4 columns: 1st deriv of logL; 2nd deriv of logL; delta = 1st/2nd; theta - delta
NumericMatrix l1l2M(npers,4);

for(int it = 0; it < nitem; it++)
  {
    
  String modit = model(it);
  //std::cout << "modit = " << modit <<  std::endl ;
  IntegerVector respvec = awm(_,it);
  double alpha = ALPHA(it);
  NumericVector delta = DELTA(_,it);
  LogicalVector nas(maxca);
  
  if(modit == "4PL")  
    {

  double lowerA = LOWA(it);
  double upperA = UPPA(it);
  
  // find NA and kill them
  for (int fna = 0; fna < maxca; fna++) {
    nas[fna] = NumericVector::is_na(delta[fna]);
  }
  // parameters without missing values. missing values should only be possible at the end of the matrix
  NumericVector delta1 = delta[!nas];
  
  //int kmax = delta1.size();
  
  for(int pe = 0; pe < npers; pe++)
    {
     
    int resp = respvec(pe);
    double theta = THETA(pe);
    NumericVector ergP(3);
    
    
    // NA handling 

    // if the i,j obs is NA, add nothing
    if(IntegerVector::is_na(resp))
    { // in case of missing value as response
    continue;
    } else 
      {
        
     ergP = P_4pl(delta1, alpha, theta, lowerA, upperA);  
      
      // l1
      double Qj = 1 - ergP(0);
      l1l2M(pe,0) += (resp - ergP(0))/(ergP(0)*Qj) * ergP(1);
      l1l2M(pe,1) +=  ergP(2);

      }

//    std::cout << "l1 = " << showme <<  std::endl ; // weg
    }
      
    
      
    } else if(modit == "GPCM")
      {
        
        // find NA and kill them
        for (int fna = 0; fna < maxca; fna++) {
          nas[fna] = NumericVector::is_na(delta[fna]);
        }
        // parameters without missing values. missing values should only be possible at the end of the matrix
        NumericVector delta1 = delta[!nas];
        
        int kmax = delta1.size();
        
      
        for(int pe = 0; pe < npers; pe++)
          {
           
          int resp = respvec(pe);
          double theta = THETA(pe);
           
          // NA handling 
      
          // if the i,j obs is NA, add nothing
          if(IntegerVector::is_na(resp))
          {
          continue;
          // jumps to next iteration step
      //    l1l2M(pe,0) += 0;
      //    l1l2M(pe,1) += 0;
      
          } else 
            {
              
            // rattert die ks durch
            
            double rs = 0;
            double rs2 = 0;
            double ls2 = 0;
            double ergP = 0;
            
            for(int ks = 0; ks < kmax; ks++)
              {
            ergP = P_gpcm(delta1, alpha, theta, ks);
            
            rs += ks * alpha * ergP;
            // second derivates right and left side
            rs2 += ks * alpha * ergP;
            //ls2 += pow(ks,2) * pow(alpha,2) * ergP;
            ls2 += ks*ks * alpha * alpha * ergP;
              }
            rs2 = rs2 * rs2;
            
            // write first and second derivs in 2 column matrix - for each person
            l1l2M(pe,0) += resp * alpha - rs;
            l1l2M(pe,1) += ls2 - rs2;
            }
      //    double showme = l1l2M(pe,0); // weg
      //    std::cout << "l1 = " << showme <<  std::endl ; // weg
          }
              
              
      }
    

  }



if(map)
  {
  NumericVector corrterm1(npers);
  NumericVector corrterm2(npers);
  
  corrterm1 = (THETA - mu)/sigma2;
  corrterm2 = 1/sigma2;
  
  l1l2M(_,1) = l1l2M(_,1) * (-1);
  l1l2M(_,2) = (l1l2M(_,0) - corrterm1) / (l1l2M(_,1)-corrterm2);
  l1l2M(_,3) = THETA - l1l2M(_,2); 
  
  } else 
    {
      
    l1l2M(_,1) = l1l2M(_,1) * (-1);
    l1l2M(_,2) = l1l2M(_,0)/l1l2M(_,1);
    l1l2M(_,3) = THETA - l1l2M(_,2);
      
    }



return l1l2M;

}



// GPCM - 4PL mixed - WLE**************************************************



// [[Rcpp::export]]
NumericMatrix Lgpcm4pl_wle(IntegerMatrix awm, NumericMatrix DELTA, NumericVector ALPHA,
                           NumericVector LOWA, NumericVector UPPA, NumericVector THETA, 
                           CharacterVector model)
{
// awm = antwortmatrix
// 
int npers = awm.nrow();
int nitem = awm.ncol();
int maxca = DELTA.nrow();
// ALPHA muss so lang sein wie nitem anzeigt
// THETA muss so lang sein wie npers anzeigt
// DELTA muss so viele Spalten haben wie nitem anzeigt
// maxca gibt an wieviele maximale kategorien

// l1 fuer die personen
//NumericVector l1vec(npers);
//NumericVector l2vec(npers);
// 
NumericMatrix l1l2M(npers,6);
NumericVector ct(npers); //weg

for(int it = 0; it < nitem; it++)
  {
    
  String modit = model(it);  
  // Response vector of one item  
  IntegerVector respvec = awm(_,it);
  double alpha = ALPHA(it);
  NumericVector delta = DELTA(_,it);
  LogicalVector nas(maxca);
   
   
   if(modit == "4PL")  
      {

      double lowerA = LOWA(it);
      double upperA = UPPA(it);
      
      // find NA and kill them
      for (int fna = 0; fna < maxca; ++fna) {
        nas[fna] = NumericVector::is_na(delta[fna]);
      }
      // parameters without missing values. missing values should only be possible at the end of the matrix
      NumericVector delta1 = delta[!nas];
      
      //int kmax = delta1.size();
      
      for(int pe = 0; pe < npers; pe++)
        {
         
        int resp = respvec(pe);
        double theta = THETA(pe);
        //NumericVector ergP(3);
        
        
        // NA handling 
    
        // if the i,j obs is NA, add nothing
        if(IntegerVector::is_na(resp))
        { // in case of missing value as response
        continue;
        } else 
          {
            
         NumericVector ergP = P_4pl4wle(delta1, alpha, theta, lowerA, upperA);  
          
          // l1
          double Qj = 1 - ergP(0);
          l1l2M(pe,0) += (resp - ergP(0))/(ergP(0)*Qj) * ergP(1);
          l1l2M(pe,1) +=  ergP(2);
          
          // oberer korrekturterm
          l1l2M(pe,2) += ergP(1) * ergP(3) / (ergP(0)*Qj); // ersetzen?!
          // unterer korrekturterm
          l1l2M(pe,3) += ergP(5);
          // I' - w wont need this
          //l1l2M(pe,4) += ergP(6);
          }
    
        }  
            
     
        
        
      } else if(modit == "GPCM")
        {
          
               // find NA and kill them
      for (int fna = 0; fna < maxca; fna++) {
        nas[fna] = NumericVector::is_na(delta[fna]);
      }
      // parameters without missing values. missing values should only be possible at the end of the matrix
      NumericVector delta1 = delta[!nas];
      
      int kmax = delta1.size();
      
      
      
      for(int pe = 0; pe < npers; pe++)
        {
         
        int resp = respvec(pe);
        double theta = THETA(pe);
         
        // NA handling 
    
        // if the i,j obs is NA, add nothing
        if(IntegerVector::is_na(resp))
        {
        continue;
        // jumps to next iteration step
    //    l1l2M(pe,0) += 0;
    //    l1l2M(pe,1) += 0;
    
        } else 
          {
            
          // rattert die ks durch
          
          double rs = 0;
          double rs2 = 0;
          double ls2 = 0;
          double ergP = 0;
          
          for(int ks = 0; ks < kmax; ks++)
            {
          ergP = P_gpcm(delta1, alpha, theta, ks);
          
          rs += ks * alpha * ergP;
          // second derivates right and left side
          rs2 += ks * alpha * ergP;
          //ls2 += pow(ks,2) * pow(alpha,2) * ergP;
          ls2 += ks*ks * alpha * alpha * ergP;
            }
          rs2 = rs2*rs2;
          
          // Pcorr2_gpcm(NumericVector delta, double alpha, double theta) 
          
          NumericVector ct = Pcorr1_gpcm(delta1,alpha,theta);
          l1l2M(pe,1) += ct(1); // summierte information
          l1l2M(pe,2) += ct(0); // summiert korrekturterm - korrekturterm oben
          //l1l2M(pe,3) += ct(1); // summierte information
          l1l2M(pe,3) += ct(2); // unterer korrekturterm 
  
          
    
          l1l2M(pe,0) += (resp * alpha - rs);
          //l1l2M(pe,1) += (-1)*(ls2 - rs2); 
          }
    
    //    std::cout << "l1 = " << showme <<  std::endl ; // weg
        } 
                
      
        }           



  }


l1l2M(_,3) = l1l2M(_,3)/(2*l1l2M(_,1)*l1l2M(_,1));

l1l2M(_,4) = (l1l2M(_,0) + l1l2M(_,2)/(2*l1l2M(_,1))) / (l1l2M(_,1)*(-1) + l1l2M(_,3));
l1l2M(_,5) = THETA - l1l2M(_,4);

return l1l2M;

}



// GPCM - 4PL mixed - ROBUST*************************************************

// [[Rcpp::export]]
NumericMatrix Lgpcm4pl_robust(IntegerMatrix awm, NumericMatrix DELTA, NumericVector ALPHA,
                           NumericVector LOWA, NumericVector UPPA, NumericVector THETA, 
                           CharacterVector model, double H) {
// awm = antwortmatrix
// 
int npers = awm.nrow();
int nitem = awm.ncol();
int maxca = DELTA.nrow();
// ALPHA muss so lang sein wie nitem anzeigt
// THETA muss so lang sein wie npers anzeigt
// DELTA muss so viele Spalten haben wie nitem anzeigt
// maxca gibt an wieviele maximale kategorien

// l1 fuer die personen
//NumericVector l1vec(npers);
//NumericVector l2vec(npers);
// 4 columns: 1st deriv of logL; 2nd deriv of logL; delta = 1st/2nd; theta - delta
NumericMatrix l1l2M(npers,4);

for(int it = 0; it < nitem; it++)
  {
    
  String modit = model(it);
  //std::cout << "modit = " << modit <<  std::endl ;
  IntegerVector respvec = awm(_,it);
  double alpha = ALPHA(it);
  NumericVector delta = DELTA(_,it);
  LogicalVector nas(maxca);
  
  if(modit == "4PL")  
    {

  double lowerA = LOWA(it);
  double upperA = UPPA(it);
  
  // find NA and kill them
  for (int fna = 0; fna < maxca; fna++) {
    nas[fna] = NumericVector::is_na(delta[fna]);
  }
  // parameters without missing values. missing values should only be possible at the end of the matrix
  NumericVector delta1 = delta[!nas];
  
  //int kmax = delta1.size();
  
  for(int pe = 0; pe < npers; pe++)
    {
     
    int resp = respvec(pe);
    double theta = THETA(pe);
    NumericVector ergP(3);
    
    
    // NA handling 

    // if the i,j obs is NA, add nothing
    if(IntegerVector::is_na(resp))
    { // in case of missing value as response
    continue;
    } else 
      {
     double hub = r_huber_4pl(delta=delta1,alpha, theta, lowerA, upperA, H);   
     ergP = P_4pl(delta1, alpha, theta, lowerA, upperA);  
      
      // l1
      double Qj = 1 - ergP(0);
      l1l2M(pe,0) += hub*(resp - ergP(0))/(ergP(0)*Qj) * ergP(1);
      l1l2M(pe,1) +=  ergP(2);

      }

//    std::cout << "l1 = " << showme <<  std::endl ; // weg
    }
      
    
      
    } else if(modit == "GPCM")
      {
        
        // find NA and kill them
        for (int fna = 0; fna < maxca; fna++) {
          nas[fna] = NumericVector::is_na(delta[fna]);
        }
        // parameters without missing values. missing values should only be possible at the end of the matrix
        NumericVector delta1 = delta[!nas];
        
        int kmax = delta1.size();
        
      
        for(int pe = 0; pe < npers; pe++)
          {
           
          int resp = respvec(pe);
          double theta = THETA(pe);
           
          // NA handling 
      
          // if the i,j obs is NA, add nothing
          if(IntegerVector::is_na(resp))
          {
          continue;
          // jumps to next iteration step
      //    l1l2M(pe,0) += 0;
      //    l1l2M(pe,1) += 0;
      
          } else 
            {
              
            // rattert die ks durch
            
            double rs = 0;
            double rs2 = 0;
            double ls2 = 0;
            double ergP = 0;
            
            double hu = r_huber_gpcm(delta1,alpha,theta,H);
            
            for(int ks = 0; ks < kmax; ks++)
              {
            ergP = P_gpcm(delta1, alpha, theta, ks);
            
            rs += ks * alpha * ergP;
            // second derivates right and left side
            rs2 += ks * alpha * ergP;
            //ls2 += pow(ks,2) * pow(alpha,2) * ergP;
            ls2 += ks*ks * alpha * alpha * ergP;
              }
            rs2 = rs2 * rs2;
            
            // write first and second derivs in 2 column matrix - for each person
            l1l2M(pe,0) += (resp * alpha - rs) * hu;
            l1l2M(pe,1) += ls2 - rs2;
            }
      //    double showme = l1l2M(pe,0); // weg
      //    std::cout << "l1 = " << showme <<  std::endl ; // weg
          }
              
              
      }
    

  }

    
    l1l2M(_,1) = l1l2M(_,1) * (-1);
    l1l2M(_,2) = l1l2M(_,0)/l1l2M(_,1);
    l1l2M(_,3) = THETA - l1l2M(_,2);
      


return l1l2M;

}



// NR - Algorithm mixed --->>>  MLE + WLE + MAP <<<<---- +++++++++++++++++++++++++++++++++++++++++++++++

// [[Rcpp::export]]
List NR_mixed(IntegerMatrix awm, NumericMatrix DELTA, NumericVector ALPHA, 
              NumericVector LOWA, NumericVector UPPA, NumericVector THETA, CharacterVector model,
              String wm, int maxsteps, double exac, NumericVector mu, NumericVector sigma2, double H) {

int npers = awm.nrow();

NumericMatrix resPP(npers,2);
int howlong;

if(wm == "wle")
  {
    
  for(int newr = 0; newr < maxsteps; newr++)
    {
    NumericMatrix reso = Lgpcm4pl_wle(awm,DELTA,ALPHA,LOWA,UPPA,THETA,model);
    THETA = reso(_,5);
    
    if( (is_true(all(abs(reso(_,4)) < exac))) | (newr == (maxsteps-1)))
      {
        resPP(_,0) = THETA;
        resPP(_,1) = 1/reso(_,1);
        howlong = newr + 1;
        break;
      }
  
    }
    
  } else if(wm == "mle")
    {
      
    bool map = FALSE;
    for(int newr = 0; newr < maxsteps; newr++)
      {
      NumericMatrix reso = Lgpcm4pl_mle(awm,DELTA,ALPHA,LOWA,UPPA,THETA,model, mu, sigma2, map);
      THETA = reso(_,3);
      
      if( (is_true(all(abs(reso(_,2)) < exac))) | (newr == (maxsteps-1)))
        {
          resPP(_,0) = THETA;
          resPP(_,1) = 1/reso(_,1)*(-1);
          howlong = newr + 1;
          break;
        }
    
      }

    } else if(wm == "map")
      {  
        bool map = TRUE;
        
          for(int newr = 0; newr < maxsteps; newr++)
            {
            NumericMatrix reso = Lgpcm4pl_mle(awm,DELTA,ALPHA,LOWA,UPPA,THETA,model, mu, sigma2, map);
            THETA = reso(_,3);
            
            if( (is_true(all(abs(reso(_,2)) < exac))) | (newr == (maxsteps-1)))
              {
                resPP(_,0) = THETA;
                resPP(_,1) = 1/reso(_,1)*(-1);
                howlong = newr + 1;
                break;
              }
          
            }   
        
      } else if(wm == "robust")
    {
      

    for(int newr = 0; newr < maxsteps; newr++)
      {
      NumericMatrix reso = Lgpcm4pl_robust(awm,DELTA,ALPHA,LOWA,UPPA,THETA,model, H);
      THETA = reso(_,3);
      
      if( (is_true(all(abs(reso(_,2)) < exac))) | (newr == (maxsteps-1)))
        {
          resPP(_,0) = THETA;
          resPP(_,1) = 1/reso(_,1)*(-1);
          howlong = newr + 1;
          break;
        }
    
      }

    }

// ----
return List::create(_["resPP"] = resPP, _["nsteps"] = howlong);

}



// estimate likelihood for EAP estimation

// [[Rcpp::export]]
NumericVector Likgpcm(IntegerVector awv, NumericMatrix DELTA, NumericVector ALPHA, 
                      NumericVector nodes) 
{

//int npers = awm.nrow();
int nitem = awv.size();
int maxca = DELTA.nrow();
int nnodes = nodes.size();

NumericVector Likpernode(nnodes);

for(int no = 0; no < nnodes; no++)
  {
  double nod1 = nodes(no);  
  
  double singleLIK = 1;
   
  for(int it = 0; it < nitem; it++)  
    {
    NumericVector delta = DELTA(_,it); 
    double alpha = ALPHA(it);
    int resp = awv(it);
    LogicalVector nas(maxca);
    
      // find NA and kill them
  for (int fna = 0; fna < maxca; ++fna) 
        {
          nas[fna] = NumericVector::is_na(delta[fna]);
        }
  // parameters without missing values. missing values should only be possible at the end of the matrix
    NumericVector delta1 = delta[!nas];
    
    if(IntegerVector::is_na(resp))
      {
      continue;
      } else 
        {
        //singleLIK =* P_gpcm(delta=delta1,alpha=alpha,theta=nod1,resp=resp);
        singleLIK *= P_gpcm(delta1,alpha,nod1,resp); 
        }

    }
  
    Likpernode(no) = singleLIK;
    
}


return Likpernode;

}


// --------------------------- simulation ------------------------------------


//' Simulate data for 1/2/3/4-pl model
//' 
//' This function returns a dichotomous matrix of simulated responses under given item and person parameters.
//' 
//' 
//' 
//'@param beta A numeric vector which contains the difficulty parameters for each item.
//'@param alpha A numeric vector, which contains the slope parameters for each item.
//'@param lowerA A numeric vector, which contains the lower asymptote parameters (kind of guessing parameter) for each item.
//'@param upperA numeric vector, which contains the upper asymptote parameters for each item.
//'@param theta A numeric vector which contains person parameters.
//'
//'@seealso \link{sim_gpcm}, \link{PP_4pl}
//'@example ./R/.examples_sim.R
//'
//'@author Manuel Reif
//'
//'@export
// [[Rcpp::export]]
IntegerMatrix sim_4pl(NumericVector beta, NumericVector alpha,
                      NumericVector lowerA, NumericVector upperA, NumericVector theta) 
    {
    
    RNGScope scope;
    
    int dvs = beta.size();
    int ntheta = theta.size();
    IntegerMatrix simres(ntheta,dvs);
    
    for(int pers = 0; pers < ntheta; pers++)
      {
      double th = theta(pers);
      NumericVector Pl = lowerA + (upperA - lowerA) * exp(alpha*(th - beta))/(1+exp(alpha*(th - beta)));
      NumericVector ru = runif(dvs);
      simres(pers,_) = ifelse(Pl > ru,1,0);
      }
    
    return simres;
    }



















