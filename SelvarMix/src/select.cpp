#include "select.hpp"

Select::Select(){}

Select::Select(Vect v, CritClust b, SelectReg sReg, int packSize)
{
    this->v = v;
    this->b = b;
    this->sReg = sReg;
    this->packSize = packSize;
}

Select::Select(Vect v, SelectReg sReg, int packSize)
{
    this->v = v;
    this->sReg = sReg;
    this->packSize = packSize;
}

List Select::selectS(vector<int> Order)
{
    int  InitialProjectsNb = this->v.experiments.size();
    //string m;
    double CritValue;
    const int numeromodeleaux = 1;
    int firstIndex, lastIndex, idx;
    firstIndex = 1;
    lastIndex = firstIndex + packSize;
    int ClustVar = 1000;
    vector<int> aux, varSelectReg, varSelectClust_aux, varSelectClust;
    double critClustaux = 0.0, critDiffClust = 0.0;
    List Mylist, Mylistaux;
    
    varSelectClust.clear();
    varSelectClust.push_back(Order[0]);
    Mylist = b.ClustBestModel(varSelectClust);
    CritValue = as<double>(Mylist["criterionValue"]);
    string s_target = as<string>(Mylist["error"]);
    //cout << "mixmod error  = "<< s_target << endl;
    while((ClustVar > 0) && (firstIndex < (int)Order.size()) && (s_target=="No error"))
    {
        ClustVar = 0;
        for(idx = firstIndex; idx < lastIndex; ++idx)
        {
             if(idx < (int)Order.size())
             {
                aux.clear(); varSelectReg.clear(); varSelectClust_aux.clear();
                aux.push_back(Order[idx]);
                varSelectReg = (this->sReg).selectReg(varSelectClust,aux,InitialProjectsNb);
                varSelectClust_aux = (this->v).ajouter_var(varSelectClust,aux);
                //cout << "idx ... " << idx << "...." << endl; 
                Mylistaux = b.ClustBestModel(varSelectClust_aux);
                s_target = as<string>(Mylistaux["error"]);
                //cout << "mixmod error  = "<< s_target << endl;
                List mylist = v.bicReggen(aux, varSelectReg, numeromodeleaux);
                //cout << "idx ... " << idx << " OK " << endl;
                critClustaux = as<double>(Mylistaux["criterionValue"]);
                critDiffClust = critClustaux - CritValue - as<double>(mylist["bicvalue"]);
                if(critDiffClust > 0 && (s_target=="No error"))
                {
                    varSelectClust =  varSelectClust_aux;
                    Mylist = Mylistaux;
                    CritValue = as<double>(Mylist["criterionValue"]);
                    ClustVar++;
                    
                }
            }
        }
        firstIndex = lastIndex;
        lastIndex += packSize;
        
    }
    return List::create(Named("S") = wrap(varSelectClust),
                        Named("model") = Mylist["model"],
                        Named("criterionValue") = Mylist["criterionValue"],
                        Named("criterion") = Mylist["criterion"],
                        Named("nbcluster") = Mylist["nbcluster"],
                        Named("parameters") = Mylist["parameters"],
                        Named("proba") = Mylist["proba"],
                        Named("partition") = Mylist["partition"]);
}

vector<int> Select::selectW(vector<int> Order, vector<int>OtherVar)
{
    vector<int> varIndep;
    int InitialProjectsNb = this->v.experiments.size();
    int firstIndex, lastIndex, idx;
    lastIndex = Order.size();
    firstIndex = lastIndex - packSize;
    int ClustVar = 1000;
    
    //vector<int> OtherVar, varIndep_Role(Order.size()), varRegAux, aux;
    vector<int> varIndep_Role(Order.size()), varRegAux, aux;
    for(int l = 0; l < (int)Order.size(); ++l)
        varIndep_Role[l] = 0;
    
    varIndep.clear();
    while((ClustVar > 0) && (firstIndex >= 0))
    {
        ClustVar =0;
        for(idx = (lastIndex - 1); idx >= firstIndex; idx--)
        {
            if(idx >= 0)
            {
                varRegAux.clear(); aux.clear();
                //OtherVar.clear();
                //for(int l = 0; l < idx; ++l)
                //    if((varIndep_Role[l] == 0) && (l != idx))
                //        OtherVar.push_back(Order[l]);
                
                aux.push_back(Order[idx]);
                varRegAux=sReg.selectReg(OtherVar,aux,InitialProjectsNb);              	    
                if(varRegAux.empty())
                {
                    varIndep.push_back(Order[idx]);
                    ClustVar++;
                    varIndep_Role[idx] = 1; 
                } 
                
            } 
        } 
        lastIndex = firstIndex; 
        firstIndex -= packSize; 
    }
    return(varIndep);
};




