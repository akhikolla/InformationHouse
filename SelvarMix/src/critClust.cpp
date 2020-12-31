#include "critClust.hpp"
//****************************************************************************//
//constructeur ***************************************************************//
CritClust::CritClust(){}

CritClust::CritClust(int k, S4 m, NumericMatrix data, string crit, IntegerVector knownlabels, bool DA)
{
    this->crit = crit;
    this->m = m;
    this->k = k;
    this->data = data;
    this->knownlabels = knownlabels;
    this->DA = DA;
}




List CritClust::ClustBestModel(vector<int> numExp)
{
  //.... les fonctions qui proviennent de R
  Environment Rmixmod("package:Rmixmod");
  Environment base("package:base");
  Function dataframe = base["data.frame"];
  Function RmixmodLearn = Rmixmod["mixmodLearn"];
  Function RmixmodCluster = Rmixmod["mixmodCluster"];
  //Function RmixmodStrategy = Rmixmod["mixmodStrategy"];
  
  
  if(DA == false){
    NumericMatrix dataAux(data.nrow(), numExp.size());
    for(int j = 0;  j < (int)numExp.size(); ++j)
    dataAux(_, j)  = data(_,numExp[j]-1);
    
    //S4 mixmodstrategy = RmixmodStrategy(Named("nbTry") = 2, 
    //Named("nbTryInInit") = 100, 
    //Named("nbIterationInInit") = 20);   
    
    S4 xem = RmixmodCluster(Named("data") = dataframe(dataAux),
    Named("nbCluster") = k,
    Named("models") = m,
    //Named("strategy") = mixmodstrategy,
    Named("criterion") = crit);
    S4 bestResult = xem.slot("bestResult");
    return List::create(Named("criterionValue") = -as<double>(bestResult.slot("criterionValue")),
    Named("criterion") = bestResult.slot("criterion"),
    Named("nbcluster") = bestResult.slot("nbCluster"),
    Named("model") = bestResult.slot("model"),
    Named("parameters") = bestResult.slot("parameters"),
    Named("proba") = bestResult.slot("proba"),
    Named("partition") = bestResult.slot("partition"),
    Named("error") = bestResult.slot("error"));
    
  }
  else
  {
    NumericMatrix dataAux(data.nrow(), numExp.size());
    for(int j = 0;  j < (int)numExp.size(); ++j)
    dataAux(_, j)  = data(_,numExp[j]-1);
    
    StringVector bic(2);
    bic[0] ="BIC";
    bic[1] ="CV";
    S4 xem = RmixmodLearn(Named("data") = dataframe(dataAux), 
    Named("knownLabels") = knownlabels,
    Named("models") = m,
    Named("criterion") = bic);
    S4 bestResult = xem.slot("bestResult");
    NumericVector critvalues(bestResult.slot("criterionValue"));
    return List::create(Named("criterionValue") = -critvalues[0],
                        Named("model") = bestResult.slot("model"),
                        Named("proba") = 0, 
                        Named("criterion") = "BIC",
                        Named("parameters") = bestResult.slot("parameters"),
                        Named("nbcluster") = bestResult.slot("nbCluster"),
                        Named("partition") = bestResult.slot("partition"),
                        Named("error") = bestResult.slot("error"));
  };
  
};




