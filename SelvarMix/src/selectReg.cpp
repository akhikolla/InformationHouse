#include "selectReg.hpp"

//****************************************************************************//
SelectReg::SelectReg(){}
 
SelectReg::SelectReg(Vect v)
{ 
  this->v = v;
}

//****************************************************************************//
//***********    Exclusion step for regression  ******************************//
//****************************************************************************//
void SelectReg::exclusion_reg(vector<int>& varSelectReg, vector<int>& varNonSig, vector<int>& jE, vector<int>& jI, int& stop, int& InitialProjectsNb)
{
 
  //calculation of bicRegTotal
  const int numeromodeleaux=1;
  List mylist = v.bicReggen(varNonSig, varSelectReg, numeromodeleaux);
  double bicRegTotal = as<double>(mylist["bicvalue"]);
   
  //initialization of bicDiffReg and of aux
  double bicDiffReg = 0.0;
  vector<int> aux;
  aux.push_back(varSelectReg[0]); 
  vector<int> numProjets_aux = (this->v).enlever_var(varSelectReg,aux);
  //initialization of jEmin with the first element of varSelectReg
  vector<int> jEmin;
  jEmin.push_back(varSelectReg[0]);
  //calculation of bicDiffReg
  mylist = v.bicReggen(varNonSig, numProjets_aux, numeromodeleaux);
  bicDiffReg = bicRegTotal - as<double>(mylist["bicvalue"]);
  aux.clear(); numProjets_aux.clear();
  
  //temporary variable for determining the minimal  bicDiffReg
  double bicDiffReg_aux = 0.0;
  for (int j=1; j < (int)varSelectReg.size();++j)
     {
       aux.push_back(varSelectReg[j]);       
       numProjets_aux = (this->v).enlever_var(varSelectReg,aux);
       List mylist = v.bicReggen(varNonSig, numProjets_aux, numeromodeleaux);
       bicDiffReg_aux = bicRegTotal - as<double>(mylist["bicvalue"]);

       //the minimal bicDiffReg
       if (bicDiffReg_aux<=bicDiffReg)
         {
            bicDiffReg = bicDiffReg_aux;
            jEmin.clear();
            jEmin.push_back(varSelectReg[j]);
         }
    
       aux.clear(); numProjets_aux.clear();
     }
  //exclusion step decision
  if (bicDiffReg<=0)
    {       
      varSelectReg = (this->v).enlever_var(varSelectReg,jEmin);  
      jE.clear();
      jE.push_back(jEmin[0]);    
      if (jE==jI)
         stop = 1; 
      else     
         stop = 0;    
    }
  else
    {
      jE.clear();
      if (jI.empty())
         stop = 1; 
      else
         stop = 0;     
    }//end else

}//end SelectReg::exclusion_reg


//****************************************************************************//
//***********    Inclusion step for regression  ******************************//
//****************************************************************************//
void SelectReg::inclusion_reg(vector<int> varSelect, vector<int>& varSelectReg, vector<int>& varNonSig,vector<int>& jE,vector<int>& jI,int& stop, int& InitialProjectsNb)
{
  const int numeromodeleaux=1; 
  List mylist = v.bicReggen(varNonSig, varSelectReg, numeromodeleaux);
  double bicRegTotal = as<double>(mylist["bicvalue"]);   
  
  vector<int> varSelectRegBis = (this->v).enlever_var(varSelect,varSelectReg);   
  //initialization of bicDiffReg and of aux
  double bicDiffReg = 0.0;
  vector<int> aux;
  aux.push_back(varSelectRegBis[0]);   
  vector<int> numProjets_aux = (this->v).ajouter_var(varSelectReg,aux);  
 
  //initialization of jImax with the first valeur of varSelectRegBis
  vector<int> jImax;
  jImax.push_back(varSelectRegBis[0]);
  //calculation of bicDiffReg
  mylist = v.bicReggen(varNonSig, numProjets_aux, numeromodeleaux);
  bicDiffReg = -bicRegTotal + as<double>(mylist["bicvalue"]); 

  aux.clear(); numProjets_aux.clear();
  
  //temporary variable for determining the maximal bicDiffReg
  double bicDiffReg_aux = 0.0;
  for (int j=1; j < (int)varSelectRegBis.size();++j)
     {       
         aux.push_back(varSelectRegBis[j]);
         numProjets_aux = (this->v).ajouter_var(varSelectReg,aux);
         mylist = v.bicReggen(varNonSig, numProjets_aux, numeromodeleaux);
	       bicDiffReg_aux = -bicRegTotal + as<double>(mylist["bicvalue"]);
         //determination of the maximal bicDiffReg
         if (bicDiffReg_aux>bicDiffReg)
           {
             bicDiffReg = bicDiffReg_aux;
             jImax.clear();
             jImax.push_back(varSelectRegBis[j]);
           }
         aux.clear();
     }
   //inclusion step decision
   if (bicDiffReg>0)
     {
       if (jImax==jE)
         stop = 1;
       else
       {        
         varSelectReg =  (this->v).ajouter_var(varSelectReg,jImax);   
         jI.clear();             
         jI.push_back(jImax[0]);      
         stop = 0;    
       } 
     }  
   else
     {
        stop = 0;  
        jI.clear();
     } 
}//end SelectReg::inclusion_reg


//****************************************************************************//
//*******************s?lection dans la r?gression ****************************//
//****************************************************************************//
vector<int> SelectReg::selectReg(vector<int> varSelect,vector<int>& varNonSig, int& InitialProjectsNb)
{
  //Initialization
  vector<int> varSelectReg;
  varSelectReg = varSelect;
  vector<int> jI,jE;         
  int stop = 0;           
  
  //Variable selection for the regression 
  while (stop==0 && !varSelectReg.empty())
       {
         //Exclusion step
         SelectReg::exclusion_reg(varSelectReg,varNonSig,jE,jI,stop,InitialProjectsNb); 
         //Inclusion step
         if (stop==0)
            SelectReg::inclusion_reg(varSelect,varSelectReg,varNonSig,jE,jI,stop,InitialProjectsNb); 
       }
       
return varSelectReg;
}//end SelectReg::selectReg



