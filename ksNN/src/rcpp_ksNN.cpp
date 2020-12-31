//' @useDynLib ksNN, .registration = TRUE
//' @importFrom Rcpp evalCpp

#include <Rcpp.h>
using namespace Rcpp;

//' This function calculates the prediction value of k* nearest neighbors algorithm.
//' @param Label vectors of the known labels of the samples.
//' @param Distance vectors of the distance between the target sample we want to predict and the other samples.
//' @param L_C parameter of k* nearest neighbors algorithm.
//' @return the prediction value(pred) and the weight of the samples(alpha).
//' @note This algorithm is based on Anava and Levy(2017).
//' @export
//' @examples
//'   library(ksNN)
//'   set.seed(1)
//'
//'   #make the nonlinear regression problem
//'   X<-runif(100)
//'   Y<-X^6-3*X^3+5*X^2+2
//'
//'   suffle<-order(rnorm(length(X)))
//'   X<-X[suffle]
//'   Y<-Y[suffle]
//'
//'   test_X<-X[1]
//'   test_Y<-Y[1]
//'
//'   train_X<-X[-1]
//'   train_Y<-Y[-1]
//'
//'   Label<-train_Y
//'   Distance<-sqrt((test_X-train_X)^2)
//'
//'   pred_ksNN<-rcpp_ksNN(Label,Distance,L_C=1)
//'
//'   #the predicted value with k*NN
//'   pred_ksNN$pred
//'
//'   #the 'true' value
//'   test_Y
// [[Rcpp::export]]
List rcpp_ksNN(NumericVector Label, NumericVector Distance, double L_C=1.0){

	int dis_length = Distance.length();
  
  NumericVector sorted_distance = L_C * clone(Distance).sort();
  NumericVector sorted_label(Label.length());
  
  Function order("order");
  NumericVector order_num = order(Distance);
  for (int i = 0; i < max(order_num); i++){
    sorted_label[i] = Label[order_num[i] - 1];
  }
  
	double beta = 0;
	double beta2 = 0;
	double tmp = 0;

	//calculate Lambda
	NumericVector Lambda = NumericVector::create(sorted_distance[0] + 1);
	
	for (int k = 0; k < dis_length; k++) {
  		if(Lambda.length() > k){
  		  if(Lambda[k] > sorted_distance[k]){
      		beta = beta + sorted_distance[k];
      		beta2 = beta2 + pow(sorted_distance[k],2);
      		int k2 = k + 1;
      		tmp = (1.0/k2)*(beta + sqrt(k2 + (pow(beta,2) - k2*beta2 )));
      		Lambda.push_back(tmp);
  		  }
  		}
	}
	
	NumericVector alpha;
	double sum_lamda = 0;
	double tail_lambda = tail(Lambda, 1)[0];
	 
	for (int i = 0; i < dis_length; i++){
	  if(sorted_distance[i] < tail_lambda){
	    double tmp = tail_lambda - sorted_distance[i];
	    sum_lamda = sum_lamda + tmp;
	    alpha.push_back(tmp);
	  }
	}

	alpha = alpha/sum_lamda;
	double pred = sum(alpha * head(sorted_label, alpha.length()));
	
	List L = List::create(Named("pred") = pred, Named("alpha") = alpha);
	return L;
}
