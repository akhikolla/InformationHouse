#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <Rmath.h>
#include <Rdefines.h>
#include "select.hpp"
#include "selectRegGen.hpp"

//[[Rcpp::export]]
List rcppSelectS(NumericMatrix X, std::vector<int> Order, const int nbCluster, S4 CovForm, const int packSize, std::string Crit, IntegerVector knownlabels, IntegerVector DA){
    //cout<< " ... nbCluster =  "<<  nbCluster << " ... Crit =  " << Crit <<endl;
    Vect v(X);
    SelectReg sReg(v);
    vector<int> varSelectClust;
    CritClust b(nbCluster, CovForm, X, Crit, knownlabels, as<bool>(DA));
    Select s(v,b,sReg, packSize);
    return wrap(s.selectS(Order));
}


//[[Rcpp::export]]
IntegerVector rcppSelectW(NumericMatrix X, std::vector<int> Order, std::vector<int> OtherVar,const int packSize){
    Vect v(X);
    SelectReg sReg(v);
    Select s(v,sReg, packSize);
    return wrap(s.selectW(Order, OtherVar));
}


//[[Rcpp::export]]
IntegerVector rcppSelectR(NumericMatrix X, std::vector<int> S, std::vector<int> U, std::string regmodel){
    int nummodel = 0;
    if(regmodel == "LI")
        nummodel = 1;
    else
        if(regmodel == "LB")
            nummodel = 2;
        else
            nummodel = 3;
    
    Vect v(X);
    int InitialProjectsNb = v.experiments.size();
    SelectRegGen sRegGen(v);
    vector<int> varReg=sRegGen.selectReggen(S, U, nummodel,InitialProjectsNb);
    
    return wrap(varReg);
}

//[[Rcpp::export]]
List rcppCrit(NumericMatrix X, List MyList, std::vector<std::string> rgm, std::vector<std::string> idm){
    typedef vector<string> stdsvec;
    typedef vector<int> stdivec;
    typedef vector<double> stddvec;
    Vect v(X);
    int InitialProjectsNb = v.experiments.size();
    SelectReg sReg(v);
    SelectRegGen sRegGen(v);
    int rhat = 0, lhat = 0, initsave = 0;
    mat reg;
    stdivec varSelectClust, varIndep, varNonIndep, varReg, SFinal, RFinal, UFinal, WFinal, Empty, regmodel, indepmodel;
    long double critClustFinal, BicRegFinal, crit, Lmax;
    stddvec BicIndepFinal;
    BicIndepFinal.clear(); BicRegFinal=0.0; crit=0.0; Lmax = 0.0;
    varSelectClust = as<stdivec>(MyList["S"]);
    varNonIndep = as<stdivec>(MyList["U"]);
    varIndep = as<stdivec>(MyList["W"]);
    critClustFinal = as<double>(MyList["criterionValue"]);
    
    for(int p = 0; p < (int)rgm.size(); ++p)
    {
        if(rgm[p] == "LI")
            regmodel.push_back(1);
        if(rgm[p] == "LB")
            regmodel.push_back(2);
        if(rgm[p] == "LC")
            regmodel.push_back(3);
    }
    
    for(int p = 0; p < (int)idm.size(); ++p)
    {
        if(idm[p] == "LI")
            indepmodel.push_back(1);
        if(idm[p] == "LB")
            indepmodel.push_back(2);
    }
    if (varIndep.size()==0)                 //aucune variable W
        if (varNonIndep.size()==0)          //aucune variable dans U
        {
            crit = critClustFinal;
            if ((initsave==0) || ((initsave==1) & (crit>Lmax)))
            {
                initsave=1;
                SFinal=varSelectClust; RFinal=Empty; UFinal=Empty; WFinal=Empty;
                rhat=0; lhat=0;
                Lmax = crit;
            }
        }//
        else                                // else Card(U)!=0  && Card(W) = 0
            if (varNonIndep.size()==1)     // if Card(U)=1 et && Card(W) = 0
            {
                varReg=sReg.selectReg(varSelectClust,varNonIndep,InitialProjectsNb);  //r=1
                List mylist = v.bicReggen(varNonIndep,varReg,regmodel[0]);
                BicRegFinal= mylist["bicvalue"];
                crit = critClustFinal + BicRegFinal;
                if ((initsave==0) || ((initsave==1) & (crit>Lmax)))
                {
                    initsave=1;
                    SFinal=varSelectClust; RFinal=varReg; UFinal=varNonIndep; WFinal=Empty;
                    rhat=1; 
                    lhat=0;
                    reg = as<mat>(mylist("B"));
                    Lmax=crit;
                }
            } // end if Card(U)=1
            else                              // else Card(U)>1
            {
                for (int p=0; p < (int)regmodel.size();++p)
                {
                    varReg=sRegGen.selectReggen(varSelectClust,varNonIndep,regmodel[p],InitialProjectsNb);
                   List mylist =v.bicReggen(varNonIndep,varReg,regmodel[p]); 
                    BicRegFinal= mylist["bicvalue"];
                    crit = critClustFinal + BicRegFinal;
                    if ((initsave==0) || ((initsave==1) & (crit>Lmax)))
                    {
                        initsave=1;
                        SFinal=varSelectClust; RFinal=varReg; UFinal=varNonIndep; WFinal=Empty;
                        reg = as<mat>(mylist("B"));
                        rhat=regmodel[p]; lhat=0;
                        Lmax=crit;
                    }
                }// end for regmodel
            }
            else    // else Card(W)!=0
            {
                //calculation of BicIndepFinal
                for (int l=0; l < (int)indepmodel.size();++l)
                {   List mylist = v.bicReggen(varIndep,Empty,indepmodel[l]);
                    BicIndepFinal.push_back(mylist["bicvalue"]);
                }  
                if (varNonIndep.size()==0)         //The redundant variable set is empty
                {
                    for (int l=0; l < (int)indepmodel.size();++l)
                    {
                        crit = critClustFinal + BicIndepFinal[l];
                        if ((initsave==0) || ((initsave==1) & (crit>Lmax)))
                        {
                            initsave=1;
                            SFinal=varSelectClust; RFinal=Empty; UFinal=Empty; WFinal=varIndep;
                            rhat=0; lhat=indepmodel[l];
                            Lmax=crit;
                            
                        }
                    }
                }
                else                                 //The redundant variable set U is non-empty
                {
                    if (varNonIndep.size()==1)        //Card(U)=1
                    {                                                                    //card de U est 1
                        varReg=sReg.selectReg(varSelectClust,varNonIndep,InitialProjectsNb); //r=1
                        List mylist = v.bicReggen(varNonIndep,varReg,regmodel[0]); 
                        BicRegFinal= mylist["bicvalue"];
                        for (int l=0; l < (int)indepmodel.size();++l)
                        {
                            crit=critClustFinal + BicRegFinal + BicIndepFinal[l];
                            if ((initsave==0) || ((initsave==1) & (crit>Lmax)))
                            {
                                initsave=1;
                                SFinal=varSelectClust; RFinal=varReg; UFinal=varNonIndep; WFinal=varIndep;
                                reg = as<mat>(mylist("B"));
                                rhat=1; lhat=indepmodel[l];
                                Lmax=crit;
                                
                            }
                        }
                    }
                    else                              //Card(U)>1
                    {
                        for (int p=0; p < (int)regmodel.size();++p)
                        {
                            varReg=sRegGen.selectReggen(varSelectClust,varNonIndep,regmodel[p],InitialProjectsNb);
                            List mylist = v.bicReggen(varNonIndep,varReg,regmodel[p]); 
                            BicRegFinal = mylist["bicvalue"];
                            for (int l=0; l < (int)indepmodel.size();++l)
                            {
                                crit = critClustFinal + BicRegFinal+ BicIndepFinal[l]; 
                                if ((initsave==0) || ((initsave==1) & (crit>Lmax)))
                                {
                                    initsave=1;                            
                                    SFinal=varSelectClust; RFinal=varReg; UFinal=varNonIndep; WFinal=varIndep; 
                                    rhat=regmodel[p]; lhat=indepmodel[l];
                                    reg = as<mat>(mylist("B"));
                                    Lmax=crit;
                                } 
                            }
                        }
                    }//end else card(U)>1                                                                                      
                }//end else card(U)!=0                                                                                                 
            }//end else card(W)!=0  
    
    string rhats = "", lhats = ""; 
    if(rhat == 1)
        rhats = "LI";
    if(rhat == 2)
        rhats = "LB";
    if(rhat == 3)
        rhats = "LC";
    if(lhat == 1)
        lhats = "LI";
    if(lhat == 2)
        lhats = "LB";
    
    //cout << "rmodel " << rhats << "  imodel  " << lhats << endl; 
    return List::create(Named("S") = wrap(SFinal), 
                        Named("R") = wrap(RFinal), 
                        Named("U") = wrap(UFinal), 
                        Named("W") = wrap(WFinal),
                        Named("criterionValue") = Lmax, 
                        Named("criterion") = MyList["criterion"],
                        Named("nbcluster") = MyList["nbcluster"],
                        Named("model") = MyList["model"],
                        Named("rmodel") = rhats, 
                        Named("imodel") = lhats,
                        Named("parameters") = MyList["parameters"],
                        Named("proba") = MyList["proba"],
                        Named("partition") = MyList["partition"],
                        Named("regparameters")= wrap(reg));  
}
