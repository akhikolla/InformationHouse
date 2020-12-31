#include "selectRegGen.hpp"

SelectRegGen::SelectRegGen(){}
 
SelectRegGen::SelectRegGen(Vect v)
{ 
  this->v = v;
}

//****************************************************************************//
//***********    Exclusion step for regression  ******************************//
//****************************************************************************//
void SelectRegGen::exclusion_reggen(vector<int>& varSelectReg, vector<int>& varNonSig,vector<int>& jE, vector<int>& jI, int& stopreg, int& nummodel, int& InitialProjectsNb)
{
  //calculation of bicRegTotal
  List mylist = (this->v).bicReggen(varNonSig,varSelectReg,nummodel);  
  double bicRegTotal = as<double>(mylist["bicvalue"]);
  //initialization of bicDiffReg and aux
  double bicDiffReg = 0.0;
  vector<int> aux;
  aux.push_back(varSelectReg[0]); 
  vector<int> numExpAux = (this->v).enlever_var(varSelectReg,aux); 
  //Initialization of jEmin with the first element of varSelectReg
  vector<int> jEmin;
  jEmin.push_back(varSelectReg[0]);
  //calculation of bicDiffReg
  mylist = (this->v).bicReggen(varNonSig,numExpAux,nummodel);  
  bicDiffReg = bicRegTotal - as<double>(mylist["bicvalue"]);   
  
  aux.clear(); numExpAux.clear();
  //temporary variable for determining the minimal bicDiffReg
  double bicDiffReg_aux = 0.0;
  for (int j=1; j < (int)varSelectReg.size();++j)
     {
        aux.push_back(varSelectReg[j]);       
        numExpAux = (this->v).enlever_var(varSelectReg,aux);    
        List mylist = (this->v).bicReggen(varNonSig,numExpAux,nummodel);
        bicDiffReg_aux = bicRegTotal -  as<double>(mylist["bicvalue"]);
    
        if (bicDiffReg_aux<=bicDiffReg)
          {
             bicDiffReg = bicDiffReg_aux;
             jEmin.clear();
             jEmin.push_back(varSelectReg[j]);
          }
          
        aux.clear(); numExpAux.clear();
     }
 
  if (bicDiffReg<=0)
    {       
       varSelectReg = (this->v).enlever_var(varSelectReg,jEmin);  
       jE.clear();
       jE.push_back(jEmin[0]);    
       if (jE==jI)
         stopreg = 1; 
       else     
         stopreg = 0; 
    }
  else
    {
       jE.clear();
       if (jI.empty())
         stopreg = 1; 
       else
         stopreg = 0;   
    }
}//end SelectReg::exclusion_reg


//****************************************************************************//
//***********    Inclusion step for regression  ******************************//
//****************************************************************************//
void SelectRegGen::inclusion_reggen(vector<int> varSelect, vector<int>& varSelectReg, vector<int>& varNonSig,vector<int>& jE, vector<int>& jI, int& stopreg, int& nummodel, int& InitialProjectsNb)
{
  //calculation of bicRegTotal
  List mylist = (this->v).bicReggen(varNonSig,varSelectReg,nummodel);
  double bicRegTotal = as<double>(mylist["bicvalue"]);
  vector<int> varSelectRegBis = (this->v).enlever_var(varSelect,varSelectReg);   
  //Initialization of bicDiffReg and aux
  double bicDiffReg = 0.0;
  vector<int> aux;
  aux.push_back(varSelectRegBis[0]);   
  vector<int> numExpAux = (this->v).ajouter_var(varSelectReg,aux);  
  //Initialization of jImax with the first element of varSelectRegBis
  vector<int> jImax;
  jImax.push_back(varSelectRegBis[0]);
  //calculation of bicDiffReg
  mylist = (this->v).bicReggen(varNonSig,numExpAux,nummodel);
  bicDiffReg = -bicRegTotal + as<double>(mylist["bicvalue"]);   
 
  aux.clear(); numExpAux.clear();
  //temporary variable for determining the maximal  bicDiffReg
  double bicDiffReg_aux = 0.0;
  for (int j=1; j < (int)varSelectRegBis.size();++j)
     {       
        aux.push_back(varSelectRegBis[j]);
        numExpAux = (this->v).ajouter_var(varSelectReg,aux);
        mylist = (this->v).bicReggen(varNonSig,numExpAux,nummodel);
        bicDiffReg_aux = -bicRegTotal + as<double>(mylist["bicvalue"]); 
        if (bicDiffReg_aux>bicDiffReg)
          {
             bicDiffReg = bicDiffReg_aux;
             jImax.clear();
             jImax.push_back(varSelectRegBis[j]);
          }
     aux.clear();
   }

   if (bicDiffReg>0)
     {
       if (jImax==jE)
         stopreg = 1;
       else
         {        
           varSelectReg = (this->v).ajouter_var(varSelectReg,jImax);   
           jI.clear();             
           jI.push_back(jImax[0]);      
           stopreg = 0;    
         }
     }  
   else
     {
        stopreg = 0;  
        jI.clear();
     }
}//end SelectReg::inclusion_reg


//****************************************************************************//
//***************  Variable selection for regression  ************************//
//****************************************************************************//
vector<int> SelectRegGen::selectReggen(vector<int> varSelect, vector<int>& varNonSig, int nummodel, int& InitialProjectsNb)
{
  //Initialization
  vector<int> varSelectReg;
  varSelectReg = varSelect;
  vector<int> jI,jE;         
  int stopreg = 0;            
  //Variable selection
  while (stopreg==0 && !varSelectReg.empty())
       {
         //Exclusion step
         SelectRegGen::exclusion_reggen(varSelectReg,varNonSig,jE,jI,stopreg,nummodel,InitialProjectsNb); 
         //Inclusion step
         if (stopreg==0)
            SelectRegGen::inclusion_reggen(varSelect,varSelectReg,varNonSig,jE,jI,stopreg,nummodel,InitialProjectsNb); 
       }
return varSelectReg;
}//end SelectRegGen::selectReggen



