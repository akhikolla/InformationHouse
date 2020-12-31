#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector CapplyB(NumericVector us, NumericMatrix X, NumericVector mu){
  int n = us.size();
  int p = mu.size();
  int ia, ib;
  NumericVector CapplyBout(p*p);

  for(ia = 0; ia < p*p; ia++){
    for(ib = 0; ib < n; ib++){
      CapplyBout(ia) += us(ib) * (X(ib,ia/p) - mu(ia/p)) * (X(ib,ia%p) - mu(ia%p));
    }
  }
  return CapplyBout;
}

// [[Rcpp::export]]
List myfitConGraphC(NumericMatrix amat, NumericMatrix S, double n, double tol = 0.000001) {
  int nrow = amat.nrow();
  int ia, ib, ic, id, sizeA, dummyInt = 0;
  NumericMatrix W(nrow, nrow), W0(nrow,nrow), W11(nrow-1,nrow-1);
  NumericVector w12(nrow-1), s12(nrow-1), paj(nrow-1), beta(nrow-1), w(nrow-1);
  bool dummyBool;
  double alpha, dummyDouble, dummyDouble2;

  for (ia = 0; ia < nrow; ia++){
    for (ib = 0; ib < nrow; ib++){
      W(ia,ib) = S(ia,ib);
      W0(ia,ib) = S(ia,ib);
    }
  }

  int nIterate = 0;
  bool Converge = false;


  while (Converge == false){
    nIterate += 1;
    for (ia = 0; ia < nrow; ia++){
      // Initialise W11
      for (ib = 0; ib < ia; ib++){
        for(ic = 0; ic < ia; ic++){
          W11(ib,ic) = W(ib,ic);
        }
      }
      for (ib = 0; ib < ia; ib++){
        for (ic = ia+1; ic < nrow; ic++){
          W11(ib,ic-1) = W(ib,ic);
        }
      }
      for (ib = ia+1; ib < nrow; ib++){
        for (ic = 0; ic < ia; ic++){
          W11(ib-1,ic) = W(ib,ic);
        }
      }
      for(ib = ia+1; ib < nrow; ib++){
        for(ic = ia+1; ic < nrow; ic++){
          W11(ib-1,ic-1) = W(ib,ic);
        }
      }
      // Initialise w12, s12
      for(ib = 0; ib < ia; ib++){
        w12(ib) = W(ib,ia);
        s12(ib) = W(ib,ia);
      }
      for(ib = ia+1; ib < nrow; ib++){
        w12(ib-1) = W(ib,ia);
        s12(ib-1) = W(ib,ia);
      }
      //Initialise paj
      for(ib = 0; ib < ia; ib++){
        paj(ib) = amat(ia,ib);
      }
      for(ib = ia+1; ib < nrow; ib++){
        paj(ib-1) = amat(ia,ib);
      }

      //Initialise beta

      for(ib = 0; ib < nrow-1; ib++){
        beta(ib) = 0;
      }


      //Initialise w
      dummyBool = false;
      for(ib = 0; ib < nrow-1; ib++){
        if(paj(ib) != 0){
          dummyBool = true;
          break;
        }
      }

      if(dummyBool == false){
        for(ib = 0; ib < nrow-1; ib++){
          w(ib) = 0;
        }
      }else{
        //create coeff matrix
        sizeA = 0;
        for(ib = 0; ib < nrow-1; ib++){
          if(paj(ib) == 1){
            sizeA +=1;
          }
        }
        NumericVector indicies(sizeA);
        dummyInt = 0;
        for(ib = 0; ib < nrow-1; ib++){
          if(paj(ib) == 1){
            dummyInt += 1;
            indicies(dummyInt-1) = ib;
            //indicies holds the index values where paj == 1.
          }
        }
        NumericMatrix A(sizeA, sizeA);
        for(ib = 0; ib < sizeA; ib++){
          for(ic = 0; ic < sizeA; ic++){
            A(ib,ic) = W11(indicies(ib),indicies(ic));
            //A is coefficient matrix of system to solve
          }
        }
        NumericVector b(sizeA);
        for(ib = 0; ib < sizeA; ib++){
          b(ib) = s12(indicies(ib));
          //b is right hand matrix of system to solve
        }

        //Gausian elimination. Get A in upper triangular form.

        for(ib = 0; ib < sizeA-1; ib++){
          if(A(ib,ib) == 0){
            stop("Zero found in covariance matrix diagonal.");
          }
          for(ic = ib+1; ic < sizeA; ic++){
            alpha = A(ic,ib) / A(ib,ib);
            for(id = 0; id < sizeA; id++){
              A(ic,id) -=alpha * A(ib,id);
            }
            b(ic) -= alpha*b(ib);
          }
        }

        //Back substitution.
        dummyInt = sizeA;
        NumericVector sol(sizeA);
        for(ib = sizeA-1; ib >= 0; ib--){
          dummyInt -= 1;
          dummyDouble = b(ib);
          for(ic = ib+1; ic < sizeA; ic++){
            dummyDouble -= A(ib,ic)*sol(ic);
          }
          sol(ib) = dummyDouble/A(ib,ib);
        }

        //Adjust beta
        for(ib = 0; ib < sizeA; ib++){
          beta(indicies(ib)) = sol(ib);
        }

        //Matrix multiply W11 with beta
        for(ib = 0; ib < nrow-1; ib++){
          w(ib) = 0;
          for(ic = 0; ic < nrow-1; ic++){
            w(ib) += W11(ib,ic)*beta(ic);
          }
        }

        for(ib = 0; ib < ia; ib++){
          W(ib,ia) = w(ib);
          W(ia,ib) = w(ib);
        }
        for(ib = ia+1; ib < nrow; ib++){
          W(ib,ia) = w(ib-1);
          W(ia,ib) = w(ib-1);
        }
      }
    }
    //Check for convergence
    Converge = true;
    for(ia = 0; ia < nrow; ia++){
      dummyDouble = 0;
      dummyDouble2 = 0;
      for(ib = 0; ib < nrow; ib++){
        dummyDouble2 = W0(ia,ib) - W(ia,ib);
        if(dummyDouble2 < 0){
          dummyDouble2 *= -1;
        }
        dummyDouble += dummyDouble2;
      }
      if(dummyDouble > tol){
        Converge = false;
        break;
      }
    }
    if(Converge == false){
      for(ia = 0; ia < nrow; ia++){
        for(ib = 0; ib < nrow; ib++){
          W0(ia,ib) = W(ia,ib);
        }
      }
    }
  }

  return List::create(Named("Shat") = W,
                      Named("iter") = nIterate);
}
