//#include <Rcpp.h>
#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
int whichLowerEQThanX(NumericVector x, double y) {
  // Rcpp supports STL-style iterators
  int elem=0;
  int i=0;
  int stop=0;
  while ((i<=x.size())&&(stop==0)){
    if (x[i]>y){
      stop=1; 
      elem=i;
    }
    i++;
  }
  
  return elem;
}

// [[Rcpp::export]]
NumericVector COMP_Q_VECT(NumericVector x,NumericVector p, NumericVector vp) {
  
  std::sort(vp.begin(), vp.end());
  int vals;
  
  double p_val;
  NumericVector qua(vp.size());
  
  for (vals=0;vals<vp.size();vals++){
    p_val=vp[vals];
    double value=0;
    if (p_val<=0) value=x[0];
    if (p_val>=1) value=x[x.size()-1];  
    if ((p_val>0) && (p_val<1)){
      int pos1;
      double ini,fin;
      pos1=whichLowerEQThanX(p,p_val);
      ini=x[pos1-1];
      fin=x[pos1];
      value=ini+(fin-ini)*(p_val-p[pos1-1])/(p[pos1]-p[pos1-1]);
      
    }
    qua[vals]=value;
  }
  return qua;
}

// [[Rcpp::export]]
NumericVector concatenate_and_sort(NumericVector a,NumericVector b){
  // int DIG=14;
  std::vector<double> ttmp1=as<std::vector<double> >(a);
  std::vector<double> ttmp2=as<std::vector<double> >(b);
  
  ttmp1.insert(ttmp1.end(),ttmp2.begin(),ttmp2.end());
  NumericVector ttmp3;
  ttmp3=wrap(ttmp1);
  // NumericVector ttmp4;
  ttmp3=sort_unique(ttmp3);
  return ttmp3;
}
NumericVector concatenate_and_sort_not_unique(NumericVector a,NumericVector b){
  // int DIG=14;
  std::vector<double> ttmp1=as<std::vector<double> >(a);
  std::vector<double> ttmp2=as<std::vector<double> >(b);
  
  ttmp1.insert(ttmp1.end(),ttmp2.begin(),ttmp2.end());
  NumericVector ttmp3;
  ttmp3=wrap(ttmp1);
  // NumericVector ttmp4;
  std::sort(ttmp3.begin(),ttmp3.end());
  return ttmp3;
}

// [[Rcpp::export]]
List REGISTER2(S4 a, S4 b) {
  S4 x;
  S4 y;
  x=clone(a);
  y=clone(b);
  int DIG=14;
  int i;
  NumericVector p1=x.slot("p");
  NumericVector p2=y.slot("p");
  
  NumericVector tmp, tmp_p;
  int doit=0;
  if (p1.size()!=p2.size()){
    doit=1;
  }else{
    if (sum(abs(p1-p2))>1e-20) {doit=1;}
  }
  // Rcout<<"p1 "<<p1.size()<<"\n";
  // Rcout<<"p2 "<<p2.size()<<"\n";
  // Rcout<<"sab "<<p2.size()<<"\n";
  // Rcout<<" "<<doit<<"\n";
  // Rcout<<"doit "<<doit<<"\n";
      // Rcout<<"doit "<<doit<<"\n";
  if (doit>0){
  ///////////////////////////////////////////////////////////
  // check if there are empty bins in x
  tmp_p=x.slot("p");
  tmp_p=diff(tmp_p);
  
  if (min(tmp_p)<powf(0.1,DIG)){
    // Rcout<<"\n  zero weighted bin!! \n";
    tmp_p[tmp_p<powf(0.1,DIG)]=powf(0.1,DIG);
    tmp_p=tmp_p/(sum(as<NumericVector>(tmp_p)));
    NumericVector tt=cumsum(tmp_p);
    NumericVector tmp_p2(tt.size()+1);
    tmp_p2[0]=0;
    for (i=1;i<(tmp_p2.size()-1);i++)
      tmp_p2[i]=tt[i-1];
    tmp_p2[tmp_p2.size()-1]=1;
    x.slot("p")=tmp_p2;
  }
  //////////////////////////////////////////////////////////
  // check if there are empty bins in y
  tmp_p=y.slot("p");
  tmp_p=diff(tmp_p);
  
  if (min(tmp_p)<powf(0.1,DIG)){
    //Rcout<<"\n  zero weighted bin!! \n";
    tmp_p[tmp_p<powf(0.1,DIG)]=powf(0.1,DIG);
    tmp_p=tmp_p/(sum(as<NumericVector>(tmp_p)));
    NumericVector tt=cumsum(tmp_p);
    NumericVector tmp_p2(tt.size()+1);
    tmp_p2[0]=0;
    for (i=1;i<(tmp_p2.size()-1);i++)
      tmp_p2[i]=tt[i-1];
    tmp_p2[tmp_p2.size()-1]=1;
    y.slot("p")=tmp_p2;
  }
  
  ///////////////////////////////////////////////////////
  //create the commoncdf
  NumericVector ttmp3;
  
  ttmp3=concatenate_and_sort(x.slot("p"),y.slot("p"));
  tmp_p=y.slot("p");
  NumericVector ciccio;
  NumericVector t1,t2;
  t1=x.slot("p");
  t2=y.slot("p");
  
  ciccio=setdiff(ttmp3,t1);
  
  x.slot("p")=ttmp3;
  
  if (ciccio.size()>0){
    NumericVector tres,fintrex;
    tres=COMP_Q_VECT(x.slot("x"),t1,ciccio);
    fintrex=concatenate_and_sort_not_unique(x.slot("x"),tres);
    x.slot("x")=fintrex;
  }
  
  y.slot("p")=ttmp3;
  ciccio=setdiff(ttmp3,t2);
  
  if (ciccio.size()>0){
    NumericVector tres,fintrex;
    tres=COMP_Q_VECT(y.slot("x"),t2,ciccio);
    fintrex=concatenate_and_sort_not_unique(y.slot("x"),tres);
    y.slot("x")=fintrex;
  }
  }
  List My(3);
  My[0]=(as<S4>(x));
  My[1]=(as<S4>(y));
  NumericMatrix M(((as<NumericVector>(x.slot("x")))).size(),3);
  M(_,0)=(as<NumericVector>(x.slot("x")));
  M(_,1)=(as<NumericVector>(y.slot("x")));
  M(_,2)=(as<NumericVector>(x.slot("p")));
  My[2]=M;
  return My;
}

// [[Rcpp::export]]
List PREPARE_A_VEC_MAT(S4 MAT){
  //int DIG=14;
  List My(2);
  NumericVector a;
  List b=MAT.slot("M");
  int i,j;
  std::vector<double> tmp1;
  
  
  for (i=0;i<b.size();i++){
    S4 o=b[i];
    NumericVector t=o.slot("p");
    std::vector<double> tmp2=as<std::vector<double> >(t);
    tmp1.insert(tmp1.end(),tmp2.begin(),tmp2.end());
  }
  NumericVector tmp3;
  tmp3=wrap(tmp1);
  tmp3=sort_unique(tmp3);
  
  NumericMatrix M(b.size(),tmp3.size());
  //rows are individuals, columns are quantiles
  for (i=0;i<b.size();i++){
    S4 o=b[i];
    NumericVector x=o.slot("x");
    NumericVector t2=o.slot("p");
    NumericVector xf(clone(x));
    NumericVector t2f(clone(t2));  
    //std::vector<double> tmp2=as<std::vector<double> >(t);
    //tmp1.insert(tmp1.end(),tmp2.begin(),tmp2.end());
    NumericVector ciccio=setdiff(tmp3,t2f);
    
    if (ciccio.size()>0){
      NumericVector tres,fintrex;
      tres=COMP_Q_VECT(xf,t2f,ciccio);
      fintrex=concatenate_and_sort_not_unique(xf,tres);
      xf=fintrex;
    }
    //Rcout<<"\n STOASSEGNANDO \n";
    for (j=0;j<tmp3.size();j++){
      M(i,j)=xf(j);
    }
    
    
  }
  
  //Rcout<<"\n"<<tmp3<<"\n";
  My(0)=M;
  My(1)=tmp3;
  return My;
}
// [[Rcpp::export]]
S4 MEDIA_V(S4 MAT, NumericVector wei){ //compute the average distribution
  List TMP;
  TMP=PREPARE_A_VEC_MAT(MAT);
  NumericVector p;
  //assuming wei>0 and of size n
  //MAT has n rows and a column  
  p=TMP(1);
  //rows are individuals, columns are quantiles
  NumericMatrix A=TMP(0);
  int ind=A.nrow();
  int qua=A.ncol();
  NumericVector x(qua);
  double SumWei=0;
  int i,j;
  for (i=0;i<ind;i++){
    SumWei=SumWei+wei(i);
    for (j=0;j<qua;j++){
      //   Rcout<<"\n A(i,j) "<<A(i,j)<< " wei i " <<wei(i)<<" x(j) "<<x(j)<<"\n";
      x(j)=x(j)+A(i,j)*wei(i);
      
    }
  }
  
  NumericMatrix resu(qua,2);
  for (j=0;j<qua;j++){
    resu(j,0)=x(j)/SumWei;
    resu(j,1)=p(j);
  }
  
  S4 o("distributionH");
  NumericVector xx=resu(_,0);
  NumericVector pp=resu(_,1);
  
  o.slot("x")=xx;
  o.slot("p")=pp;
  NumericVector rr=diff(xx)*0.5;
  NumericVector cc(rr.size());
  NumericVector pp2=diff(pp);
  double m=0;
  double s=0;
  
  for (i=0;i<(xx.size()-1);i++){
    cc(i)=(xx[i]+xx[i+1])/2.0;
    m=m+cc[i]*pp2[i];
    s=s+pp2[i]*(cc[i]*cc[i]+rr[i]*rr[i]/3.0);
  }
  s=sqrt(s-m*m);
  
  o.slot("m")=m;
  o.slot("s")=s;
  return o;
}
// [[Rcpp::export]]
NumericMatrix SSQ_RCPP(S4 MAT,NumericVector wei){
  List resu;
  ListMatrix MM=MAT.slot("M");
  float Sumwei=sum(wei);
  arma::mat MIXP;
  MIXP.zeros(MM.ncol(),MM.ncol());
  int i;
  for (i=0;i<MM.nrow();i++){
    
    S4 TMPL("MatH");
    TMPL.slot("M")=MM(i,_);
    resu=PREPARE_A_VEC_MAT(TMPL);
    arma::mat M=as<arma::mat>(resu[0]);
    arma::vec p=resu[1];
    arma::mat Mc;
    arma::mat Mr;
    
    
    Mc=(M(span(0,(M.n_rows-1)),span(1,(M.n_cols-1)))+M(span(0,(M.n_rows-1)),span(0,(M.n_cols-2))))/2.0;
    Mr=(M(span(0,(M.n_rows-1)),span(1,(M.n_cols-1)))-M(span(0,(M.n_rows-1)),span(0,(M.n_cols-2))))/2.0;
    arma::mat W=diagmat(diff(p));
    
    
    MIXP=MIXP+Mc*W*Mc.t()*wei(i)+Mr*W*Mr.t()*wei(i)/3.0;
  }
  //now compute the second part
  
  S4 MEANS("MatH");
  ListMatrix TML(1,MM.ncol());
  int j;
  for (j=0;j<MM.ncol();j++){
    S4 TMPL("MatH");
    TMPL.slot("M")=MM(_,j);
    S4 tmp("distributionH");
    tmp=MEDIA_V(TMPL, wei);
    TML(0,j)=tmp;
  }
  MEANS.slot("M")=TML;
  List resu2;
  resu2=PREPARE_A_VEC_MAT(MEANS);
  
  ////////////////////////////////////////////////////////////////////////
  arma::mat M=as<arma::mat>(resu2[0]);
  arma::vec p=resu2[1];
  arma::mat Mc;
  arma::mat Mr;
  arma::mat W=diagmat(diff(p));
  Mc=(M(span(0,(M.n_rows-1)),span(1,(M.n_cols-1)))+M(span(0,(M.n_rows-1)),span(0,(M.n_cols-2))))/2.0;
  Mr=(M(span(0,(M.n_rows-1)),span(1,(M.n_cols-1)))-M(span(0,(M.n_rows-1)),span(0,(M.n_cols-2))))/2.0;
  arma::mat MIXP2;
  MIXP=MIXP-Mc*W*Mc.t()*Sumwei-Mr*W*Mr.t()*Sumwei/3.0;
  
  ////////////////////////////////////////////////////////////////////////
  NumericMatrix fin_resu;
  return fin_resu=wrap(MIXP);
}
// [[Rcpp::export]]
NumericMatrix COV_RCPP(S4 MAT,NumericVector wei){
  NumericMatrix cov;
  cov=SSQ_RCPP(MAT,wei);
  cov=cov/sum(wei);
  return cov;
}
// [[Rcpp::export]]
NumericMatrix CORR_RCPP(S4 MAT,NumericVector wei){
  arma::mat corr=as<arma::mat>(COV_RCPP(MAT,wei));
  arma::vec dia=sqrt(corr.diag());
  corr=corr/(dia*dia.t());
  NumericMatrix fin_resu;
  return fin_resu=wrap(corr);
}
// [[Rcpp::export]]
NumericVector M_STD_H(S4 o){ //compute mean and std of a distrbution 
  S4 o2;
  o2=clone(o);
  NumericVector xx=o2.slot("x");
  NumericVector pp=o2.slot("p");
  
  NumericVector rr=diff(xx);
  rr=rr/2.0;
  NumericVector cc(rr.size());
  NumericVector pp2=diff(pp);
  double m=0;
  int i;
  double s=0;
  
  for (i=0;i<(xx.size()-1);i++){
    cc(i)=(xx[i]+xx[i+1])/2.0;
    m=m+cc[i]*pp2[i];
    s=s+pp2[i]*(cc[i]*cc[i]+rr[i]*rr[i]/3.0);
  }
  
  s=sqrt(s-m*m);
  NumericVector resu(2);
  resu[0]=m;
  resu[1]=s;
  return resu;    
}
// [[Rcpp::export]]
double c_MM_L2_SQ_WASS_D(NumericMatrix MM){
  //int rows=MM.rows();
  double dist=0;
  NumericVector c1,r1,c2,r2,p;
  r1=diff(MM(_,0))/2.0;
  r2=diff(MM(_,1))/2.0;
  p=diff(MM(_,2));
  NumericVector tmp1=MM(_,0);
  c1=(tmp1[Range(0,(tmp1.size()-2))]+tmp1[Range(1,(tmp1.size()-1))])/2.0;
  NumericVector tmp2=MM(_,1);
  c2=(tmp2[Range(0,(tmp2.size()-2))]+tmp2[Range(1,(tmp2.size()-1))])/2.0;
  // NumericVector tmp3;
  // tmp3=(pow((c1-c2),2)+1/3.0*pow((r1-r2),2))*p;
  
  dist=sum((pow((c1-c2),2)+1/3.0*pow((r1-r2),2))*p);
  
  return dist;
}
// [[Rcpp::export]]
List c_DISTA_M(List MM, S4 Protos){// computes distances wrt prototypes
  
  ListMatrix prot_dis=(Protos.slot("M"));
  int npro=prot_dis.nrow();
  int vars=MM.size();
  int v,pr,ind;
  double dist;
  NumericMatrix GG=MM[0];
  NumericMatrix resu((GG.ncol()-1),npro);
  for(v=0;v<vars;v++){
    NumericMatrix M=MM[v];
    NumericMatrix TMP(M.nrow(),3);
    for(ind=0;ind<(M.ncol()-1);ind++){
      for(pr=0;pr<npro;pr++){
        S4 prok=prot_dis(pr,v);
        TMP(_,0)=M(_,ind);
        TMP(_,1)=as<NumericVector>(prok.slot("x"));
        TMP(_,2)=M(_,(M.ncol()-1));
        dist=c_MM_L2_SQ_WASS_D(TMP);
        resu(ind,pr)=resu(ind,pr)+dist;
               
      }
    }
    
  }
  //
  NumericVector wm(GG.ncol()-1);
  for(ind=0;ind<(GG.ncol()-1);ind++){
    wm[ind]=which_min(resu(ind,_))+1;
  }
  List L = List::create(Named("dist") = resu , _["wm"] = wm);
  return L;
}
// [[Rcpp::export]]
NumericMatrix c_DISTA_M2(List MM, S4 Protos){// computes distances wrt prototypes
  
  ListMatrix prot_dis=(Protos.slot("M"));
  int npro=prot_dis.nrow();
  int vars=MM.size();
  int v,pr,ind;
  double dist;
  
  NumericMatrix GG=MM[1];
  NumericMatrix resu((GG.ncol()-1),npro);
  for(v=0;v<vars;v++){
    NumericMatrix M=MM[v];
    NumericMatrix TMP(M.nrow(),3);
    for(ind=0;ind<(M.ncol()-1);ind++){
      for(pr=0;pr<npro;pr++){
        S4 prok=prot_dis(pr,v);
        TMP(_,0)=M(_,ind);
        TMP(_,1)=as<NumericVector>(prok.slot("x"));
        TMP(_,2)=M(_,(M.ncol()-1));
        
        dist=c_MM_L2_SQ_WASS_D(TMP);
        resu(ind,pr)=resu(ind,pr)+dist;
        //       
      }
    }
    
  }
  //
  return resu;
}
// [[Rcpp::export]]
List c_ComputeFASTSSQ(NumericMatrix subMM){
  
  int  ind=subMM.ncol()-1;
  int rr=subMM.nrow();
  double SSQ=0;
  int indiv,nqua;
  NumericVector c1,r1,cM,rM,mp,p;
  mp=subMM(_,(subMM.ncol()-1));
  p=diff(mp);
  NumericVector Mqua(rr);
  for(nqua=0;nqua<rr;nqua++){
    NumericVector tmp;
    tmp=subMM(nqua,_);
    Mqua[nqua]=mean(tmp[Range(0,(tmp.size()-2))]);
    
  }
  
  rM=diff(Mqua)/2.0;
  cM=(Mqua[Range(0,(Mqua.size()-2))]+Mqua[Range(1,(Mqua.size()-1))])/2.0;
  for(indiv=0;indiv<ind;indiv++){
    
    NumericVector tmp=subMM(_,indiv);
    r1=diff(tmp)/2.0;
    c1=(tmp[Range(0,(tmp.size()-2))]+tmp[Range(1,(tmp.size()-1))])/2.0;
    SSQ=SSQ+sum((pow((c1-cM),2)+ pow((r1-rM),2)/3.0)*p);
  }
  
  List L = List::create(Named("SSQ") = SSQ , _["mx"] = Mqua,_["mp"]=mp);
  
  return L;
}
// [[Rcpp::export]]
List c_Compute_M_from_MM(NumericMatrix subMM){
  
  int rr=subMM.nrow();
  int nqua;
  NumericVector cM,rM,mp,p;
  mp=subMM(_,(subMM.ncol()-1));
  p=diff(mp);
  NumericVector Mqua(rr);
  for(nqua=0;nqua<rr;nqua++){
    NumericVector tmp;
    tmp=subMM(nqua,_);
    Mqua[nqua]=mean(tmp[Range(0,(tmp.size()-2))]);
    
  }
  List L = List::create(Named("mx") = Mqua,_["mp"]=mp);
  
  return L;
}
// [[Rcpp::export]]
NumericMatrix c_Fast_D_Mat(List MM){
  int vars=MM.size();
  int v,ind,ind2,Nind;
  Nind=(as<NumericMatrix>(MM[0])).ncol()-1;
  NumericMatrix M=MM[0];
  NumericMatrix D(Nind,Nind);
  for(v=0;v<vars;v++){
    NumericMatrix M=MM[v];
    NumericMatrix TMP(M.nrow(),3);
    for(ind=0;ind<(M.ncol()-2);ind++){
      for(ind2=(ind+1);ind2<(M.ncol()-1);ind2++){
        TMP(_,0)=M(_,ind);
        TMP(_,1)=M(_,ind2);
        TMP(_,2)=M(_,(M.ncol()-1));
        D(ind,ind2)=D(ind,ind2)+c_MM_L2_SQ_WASS_D(TMP);
        D(ind2,ind)=D(ind,ind2);       
      }
    }
  }
  return D;
}
// [[Rcpp::export]]
List c_Prepare(S4 x,bool simplify,int qua,bool standardize){
  int ind=((as<ListMatrix>(x.slot("M"))).nrow());
  int vars=((as<ListMatrix>(x.slot("M"))).ncol());
  int i,j;
  
  List MM(vars);
  S4 y=clone(x);
  ListMatrix M=y.slot("M");
  if(simplify){
    NumericVector p(qua+1);
    for(i=0;i<(qua+1);i++){
      p[i]=i/double(qua);
    }
    // Rcout<<"\n p -->"<<p<<"\n";
    for (j=0;j<vars;j++){
      NumericMatrix TMP((qua+1),(ind+1));
      TMP(_,ind)=p;
      for(i=0;i<ind;i++){
        //  COMP_Q_VECT(NumericVector x,NumericVector p, NumericVector vp)
        S4 o=M(i,j);
        NumericVector dom,dom2;
        dom=COMP_Q_VECT((as<NumericVector>(o.slot("x"))),
                        (as<NumericVector>(o.slot("p"))),
                        p);
        
        
        o.slot("x")=dom;
        o.slot("p")=p;
        NumericVector MOM;
        MOM=M_STD_H(o);
        //transformation with invariance with respect mean and std
        
        dom2=(dom-MOM[0])*(((double)o.slot("s"))/MOM[1])+MOM[0]; 
        //
        
        
        //o.slot("x")=dom;
        o.slot("m")=MOM[0];
        o.slot("s")=MOM[1];
        M(i,j)=o;
        TMP(_,i)=dom;
      }
      
      MM[j]=TMP;
    }
  }else{
    for(j=0;j<vars;j++){
      
      List TMP_res;
      S4 T_MAT("MatH");
      T_MAT.slot("M")=M(_,j);
      TMP_res=PREPARE_A_VEC_MAT(T_MAT);
      NumericMatrix AA=as<NumericMatrix>(TMP_res(0));//Ã¨ la matrice
      NumericVector BB=as<NumericVector>(TMP_res(1));//sono i pesi
      NumericMatrix CC((AA.nrow()+1),AA.ncol());
      int kk;
      for(kk=0;kk<AA.nrow();kk++){
        CC(kk,_)=AA(kk,_);
      }
      CC(AA.nrow(),_)=BB;
      
      MM[j]=transpose(CC);
      for(i=0;i<ind;i++){
        S4 o=M(i,j);
        NumericVector dom;
        dom=CC(i,_);
        
        o.slot("x")=dom;
        o.slot("p")=BB;
        NumericVector MOM;
        MOM=M_STD_H(o);
        M(i,j)=o;
      }
      
    }
  }
  
  
  // standardize data if required
  if(standardize){
    NumericVector STAND(vars);
    NumericVector Mc(vars);
    for(j=0;j<vars;j++){
      NumericVector wei(ind,1.0);
      S4 T_MAT("MatH");
      ListMatrix GG(ind,1);
      for (i=0;i<ind;i++){
        GG(i,0)=M(i,j);
        
      }
      T_MAT.slot("M")=GG;
      NumericMatrix TTMMPP;
      TTMMPP=COV_RCPP(T_MAT,wei);
      
      STAND[j]=(sqrt(TTMMPP(0,0)));
      S4 tmp_mean("distributionH");
      tmp_mean=MEDIA_V(T_MAT,wei);
      Mc[j]=tmp_mean.slot("m");
      //Rcout<<"\n M "<<Mc[j]<<"\n STAND "<<STAND[j];
      NumericMatrix MTMP;
      MTMP=as<NumericMatrix>(MM[j]);
      for(i=0;i<ind;i++){
        S4 o=M(i,j);
        o.slot("x")=((as<NumericVector>(o.slot("x")))-Mc[j])/STAND[j];
        o.slot("m")=((as<double>(o.slot("m")))-Mc[j])/STAND[j];
        o.slot("s")=(as<double>(o.slot("s")))/STAND[j];
        M(i,j)=o;
        MTMP(_,i)=(MTMP(_,i)-Mc[j])/STAND[j];
        
      }
      MM[j]=MTMP;
    }
    
  }
  
  List resu = List::create(Named("MM") = MM , _["x"] = y);
  return resu;
}  
// [[Rcpp::export]]
List c_Prepare2(S4 x,bool simplify,int qua,bool standardize){
  int ind=((as<ListMatrix>(x.slot("M"))).nrow());
  int vars=((as<ListMatrix>(x.slot("M"))).ncol());
  int i,j;
  
  List MM(vars);
  //  S4 y=clone(x);
  ListMatrix M=x.slot("M");
  if(simplify){
    NumericVector p(qua+1);
    for(i=0;i<(qua+1);i++){
      p[i]=i/double(qua);
    }
    // Rcout<<"\n p -->"<<p<<"\n";
    for (j=0;j<vars;j++){
      NumericMatrix TMP((qua+1),(ind+1));
      TMP(_,ind)=p;
      for(i=0;i<ind;i++){
        //  COMP_Q_VECT(NumericVector x,NumericVector p, NumericVector vp)
        S4 o=M(i,j);
        NumericVector dom,dom2;
        dom=COMP_Q_VECT((as<NumericVector>(o.slot("x"))),
                        (as<NumericVector>(o.slot("p"))),
                        p);
        
        
        //o.slot("x")=dom;
        //o.slot("p")=p;
        //NumericVector MOM;
        //MOM=M_STD_H(o);
        //transformation with invariance with respect mean and std
        
        //dom2=(dom-MOM[0])*(((double)o.slot("s"))/MOM[1])+MOM[0]; 
        //
        
        
        //o.slot("x")=dom;
        //o.slot("m")=MOM[0];
        //o.slot("s")=MOM[1];
        //M(i,j)=o;
        TMP(_,i)=dom;
      }
      
      MM[j]=TMP;
    }
  }else{
    for(j=0;j<vars;j++){
      
      List TMP_res;
      S4 T_MAT("MatH");
      T_MAT.slot("M")=M(_,j);
      TMP_res=PREPARE_A_VEC_MAT(T_MAT);
      NumericMatrix AA=as<NumericMatrix>(TMP_res(0));//Ã¨ la matrice
      NumericVector BB=as<NumericVector>(TMP_res(1));//sono i pesi
      NumericMatrix CC((AA.nrow()+1),AA.ncol());
      int kk;
      for(kk=0;kk<AA.nrow();kk++){
        CC(kk,_)=AA(kk,_);
      }
      CC(AA.nrow(),_)=BB;
      
      MM[j]=transpose(CC);
      for(i=0;i<ind;i++){
        S4 o=M(i,j);
        NumericVector dom;
        dom=CC(i,_);
        
        // o.slot("x")=dom;
        //  o.slot("p")=BB;
        //  NumericVector MOM;
        //  MOM=M_STD_H(o);
        M(i,j)=o;
      }
      
    }
  }
  
  
  // standardize data if required
  if(standardize){
    NumericVector STAND(vars);
    NumericVector Mc(vars);
    for(j=0;j<vars;j++){
      NumericVector wei(ind,1.0);
      S4 T_MAT("MatH");
      ListMatrix GG(ind,1);
      for (i=0;i<ind;i++){
        GG(i,0)=M(i,j);
        
      }
      T_MAT.slot("M")=GG;
      NumericMatrix TTMMPP;
      TTMMPP=COV_RCPP(T_MAT,wei);
      
      STAND[j]=(sqrt(TTMMPP(0,0)));
      S4 tmp_mean("distributionH");
      tmp_mean=MEDIA_V(T_MAT,wei);
      Mc[j]=tmp_mean.slot("m");
      //Rcout<<"\n M "<<Mc[j]<<"\n STAND "<<STAND[j];
      NumericMatrix MTMP;
      MTMP=as<NumericMatrix>(MM[j]);
      for(i=0;i<ind;i++){
        //S4 o=M(i,j);
        //o.slot("x")=((as<NumericVector>(o.slot("x")))-Mc[j])/STAND[j];
        //o.slot("m")=((as<double>(o.slot("m")))-Mc[j])/STAND[j];
        //o.slot("s")=(as<double>(o.slot("s")))/STAND[j];
        //M(i,j)=o;
        MTMP(_,i)=(MTMP(_,i)-Mc[j])/STAND[j];
        
      }
      MM[j]=MTMP;
    }
    
  }
  
  List resu = List::create(Named("MM") = MM );
  return resu;
}
// [[Rcpp::export]]
List c_ComputeFastSSQ(NumericMatrix subMM){
  
  int ind=(subMM.ncol()-1);
  int rr=subMM.nrow();
  double SSQ=0;
  int i,j;
  NumericVector Mea(rr);
  for(i=0;i<rr;i++){
    NumericVector TMP;
    TMP=subMM(rr,_);
    double M=0;
    for(j=0;j<ind;j++){
      M=M+TMP[j];
    }
    Mea=M/double(ind);
  }
  NumericVector cum(rr);
  NumericVector p((rr-1));
  cum=subMM(_,ind);
  p=diff(cum);
  NumericVector cm((rr-1));
  NumericVector rm((rr-1));
  rm=diff(Mea)/2.0;
  for(i=0;i<(rr-1);i++){
    cm[i]=(Mea[i]+Mea[i+1])/2.0;
  }
  cm=0;
  for(i=0;i<ind;i++){
    NumericVector ci((rr-1));
    NumericVector ri((rr-1));
    NumericVector tmp(rr);
    tmp=subMM(_,i);
    ri=diff(tmp)/2.0;
    for(j=0;j<(rr-1);j++){
      ci[j]=(tmp[j]+tmp[j+1])/2.0;
    }
    SSQ=SSQ+sum(p*(pow((ci-cm),2)+pow((ri-rm),2)/3.0));
    
  }
  // for (indiv in 1:ind){
  //   SSQ=SSQ+sum(p1*((ci-cm)^2)+ p1/3*((ri-rm)^2))
  // }
  // return(list(SSQ=SSQ,mx=Mea, mp=cum)
  List resu = List::create(Named("SSQ") = SSQ , _["mx"] = Mea,_["mp"] = cum);
  return resu;
}
// [[Rcpp::export]]
double c_MOM_D(S4 o){//computes the integral of squares of a distrib
  int k;
  NumericVector dom,cum, w, tmp, radiiQ;
  dom=as<NumericVector>(o.slot("x"));
  cum=as<NumericVector>(o.slot("p"));
  w=diff(cum);
  radiiQ=pow(diff(dom)/2.0,2);
  int rr=dom.size();
  NumericVector cenQ(rr-1);
  for (k=0;k<(rr-1);k++){
    cenQ[k]=pow(((dom[k]+dom[k+1])*0.5),2);
  }
  tmp=w*(cenQ+radiiQ/3.0);
  double res;
  res=std::accumulate(tmp.begin(),
                      tmp.end(), 0.0);
  return res;
}
// [[Rcpp::export]]
NumericMatrix c_MOM_MAT(S4 x){ //computes the integral of squares of a MatH 
  int ind=((as<ListMatrix>(x.slot("M"))).nrow());
  int vars=((as<ListMatrix>(x.slot("M"))).ncol());
  int i,j;
  NumericMatrix res(ind,vars);
  
  ListMatrix M=x.slot("M");
  for(i=0;i<ind;i++){
    for (j=0;j<vars;j++){
      S4 o=M(i,j);
      res(i,j)=c_MOM_D(o);
    }
  }
  return res;
}
// [[Rcpp::export]]
List c_STEP_2_ADA_KMEANS(List MM,S4 proto, int schema, NumericMatrix memb){
  int clu,var,ind;
  int k=((as<ListMatrix>(proto.slot("M"))).nrow());
  int vars=((as<ListMatrix>(proto.slot("M"))).ncol());
  int inds=memb.nrow();
  ListMatrix proto_dis=(proto.slot("M"));// the matrix of protos
  
  List diINDtoPROT_m(k);
  List diINDtoPROT_v(k);
  NumericMatrix D1(vars,k);
  NumericMatrix D2(vars,k);
  for (clu=0;clu<k;clu++){
    NumericMatrix D_Im(inds,vars);
    NumericMatrix D_Iv(inds,vars);
    
    for(var=0;var<vars;var++){
      NumericMatrix M=MM[var];
      S4 prokv=proto_dis(clu,var);
      for(ind=0;ind<inds;ind++){
        NumericMatrix TMP(M.nrow(),3);
        TMP(_,0)=M(_,ind);
        TMP(_,1)=as<NumericVector>(prokv.slot("x"));
        TMP(_,2)=M(_,(M.ncol()-1));
        double tmpD;
        tmpD=c_MM_L2_SQ_WASS_D(TMP);
        if ((schema==1)||(schema==3)){
          D1(var,clu)=D1(var,clu)+(tmpD*memb(ind,clu));
          D_Im(ind,var)=tmpD;
        }
        if ((schema==2)||(schema==4)||(schema==5)||(schema==6)){
          NumericVector c;
          NumericVector tmp=TMP(_,0);
          NumericVector p;
          double m1,dv,dm;
          p=diff(TMP(_,2));
          c=(tmp[Range(0,(tmp.size()-2))]+tmp[Range(1,(tmp.size()-1))])/2.0;
          m1=sum(p*c);
          dm=pow((m1-as<double>(prokv.slot("m"))),2);
          dv=tmpD-dm;
          D1(var,clu)=D1(var,clu)+(dm*memb(ind,clu));
          D2(var,clu)=D2(var,clu)+(dv*memb(ind,clu));
          D_Im(ind,var)=dm;
          D_Iv(ind,var)=dv;
        }
      }
    }
    diINDtoPROT_m[clu]=D_Im;
    diINDtoPROT_v[clu]=D_Iv;
  }
  // distances[0]=D1;
  //distances[1]=D2;
  List resu = List::create(Named("distances1") = D1 ,_["distances2"]=D2,
                           _["dIpro_m"] = diINDtoPROT_m,
                           _["dIpro_v"] = diINDtoPROT_v);
  return resu;
}
// [[Rcpp::export]]
double Provec(NumericVector x){
  double res=1;
  int i;
  for (i=0;i<x.size();i++){
    res=res*x[i];
  }
  return res;
}
// [[Rcpp::export]]
List c_cen_rad(NumericVector x){
  int l=x.size();
  NumericVector cen(l-1,0.0);
  NumericVector rad(l-1,0.0);
  cen=(x[Range(0,l-2)]+x[Range(1,l-1)])*0.5;
  rad=diff(x)*0.5;
  List resu;
  resu["cen"]=cen;
  resu["rad"]=rad;
  return resu;
  
}
// [[Rcpp::export]]
NumericMatrix c_STEP_2_2_WEIGHTS_ADA_KMEANS(NumericMatrix D1,
                                            NumericMatrix D2,
                                            int PROSUM, int schema,
                                            double theta){
  int vars=D1.nrow();
  int k=D1.ncol();
  int var,clu;
  NumericMatrix lambdas((vars*2),k);
  double num,denom;
  for (var=0;var<vars;var++){
    for (clu=0;clu<k;clu++){
      if (PROSUM==1){//product
        if (schema==1){
          if (clu<1){
            //            num=(prod(apply(distances[,,1],MARGIN=c(1),sum)))^(1/vars)
            //            denom=max(sum(distances[variables,,1]),1e-10)
            
            num=pow(Provec(as<NumericVector>(rowSums(D1))),(1/double(vars)));
            denom=sum(D1(var,_));
            if (denom<1e-10){denom=1e-10;}
            lambdas(var*2,clu)=num/denom;
            lambdas((var*2+1),clu)=num/denom;}
          else{
            lambdas(var*2,clu)=lambdas(var*2,0);
            lambdas((var*2+1),clu)=lambdas((var*2+1),0);
          }
          
        }
        if (schema==2){
          if (clu<1){
            num=pow(Provec(as<NumericVector>(rowSums(D1))),(1/double(vars)));
            denom=sum(D1(var,_));
            if (denom<1e-10){denom=1e-10;}
            lambdas(var*2,clu)=num/denom;
            num=pow(Provec(as<NumericVector>(rowSums(D2))),(1/double(vars)));
            denom=sum(D2(var,_));
            if (denom<1e-10){denom=1e-10;}
            lambdas((var*2+1),clu)=num/denom;}
          else{
            lambdas(var*2,clu)=lambdas(var*2,0);
            lambdas((var*2+1),clu)=lambdas((var*2+1),0);
          }
          
        }
        if (schema==3){
          num=pow(Provec(D1(_,clu)),(1/double(vars)));
          denom=D1(var,clu);
          if (denom<1e-10){denom=1e-10;}
          lambdas(var*2,clu)=num/denom;
          lambdas((var*2+1),clu)=num/denom;
          
        }
        if (schema==4){
          num=pow(Provec(D1(_,clu)),(1/double(vars)));
          if (num<1e-10){num=1e-10;}
          denom=D1(var,clu);
          if (denom<1e-10){denom=1e-10;}
          lambdas(var*2,clu)=num/denom;
          num=pow(Provec(D2(_,clu)),(1/double(vars)));
          if (num<1e-10){num=1e-10;}
          denom=D2(var,clu);
          if (denom<1e-10){denom=1e-10;}
          lambdas((var*2+1),clu)=num/denom;
          
          
        }
        if (schema==5){// to be done
          lambdas(var*2,clu)=1;
          lambdas((var*2+1),clu)=1;
        }
        if (schema==6){// to be done
          lambdas(var*2,clu)=1;
          lambdas((var*2+1),clu)=1;
        }
      }
      if (PROSUM==2){//sum
        if (schema==1){
          if (clu<1){
            NumericVector num2(vars,sum(D1(var,_)));
            NumericVector den2;
            
            den2=rowSums(D1);//apply(distances[,,1],1,sum);
            lambdas(var*2,clu)=1/sum(pow((num2/den2),(1/(theta-1))));
            lambdas((var*2+1),clu)=lambdas(var*2,clu);
          }else{
            lambdas(var*2,clu)=lambdas(var*2,0);
            lambdas((var*2+1),clu)=lambdas((var*2+1),0);
          }
        }
        if (schema==2){
          if (clu<1){
            NumericVector num2(vars,sum(D1(var,_)));
            NumericVector den2;
            den2=rowSums(D1);// den=apply(distances[,,1],1,sum)
            lambdas(var*2,clu)=1/sum(pow((num2/den2),(1/(theta-1))));
            
            NumericVector num3(vars,sum(D2(var,_)));
            NumericVector den3;
            den3=rowSums(D2);// den=apply(distances[,,2],1,sum)
            lambdas((var*2+1),clu)=1/sum(pow((num3/den3),(1/(theta-1))));
          }
          else{
            lambdas(var*2,clu)=lambdas(var*2,0);
            lambdas((var*2+1),clu)=lambdas((var*2+1),0);
          }
          
        }
        if (schema==3){
          NumericVector num2(vars,D1(var,clu));
          NumericVector den2;
          
          den2=D1(_,clu);//den=distances[,cluster,1]
          lambdas(var*2,clu)=1/sum(pow((num2/den2),(1/(theta-1))));
          lambdas((var*2+1),clu)=lambdas(var*2,clu);
          
        }
        if (schema==4){
          NumericVector num2(vars,D1(var,clu));
          NumericVector den2;
          
          den2=D1(_,clu);//den=distances[,cluster,1]
          lambdas(var*2,clu)=1/sum(pow((num2/den2),(1/(theta-1))));
          
          NumericVector num3(vars,D2(var,clu));
          NumericVector den3;
          
          den3=D2(_,clu);//den=distances[,cluster,1]
          lambdas((var*2+1),clu)=1/sum(pow((num3/den3),(1/(theta-1))));
        }
        if (schema==5){//to be done
          lambdas(var*2,clu)=1/double(vars);
          lambdas((var*2+1),clu)=1/double(vars);
        }
        if (schema==6){
          //to be done
          lambdas(var*2,clu)=1/double(vars);
          lambdas((var*2+1),clu)=1/double(vars);
        }
      }
    }
  }
  
  return lambdas;
}
// [[Rcpp::export]]
List c_STEP_3_AFFECT_ADA_KMEANS(NumericMatrix lambdas, 
                                List dIpro_m, 
                                List dIpro_v, int ind, int k, int vars){
  NumericMatrix DiToClu(ind,k);
  int indiv,cluster,variable;
  NumericVector IDX(ind);
  double SSQ=0;
  
  for (indiv=0;indiv<ind;indiv++){
    for (cluster=0;cluster<k;cluster++){
      for (variable=0;variable<vars;variable++){
        DiToClu(indiv,cluster)=DiToClu(indiv,cluster)+
          (lambdas((variable*2),cluster)*(as<NumericMatrix>(dIpro_m[cluster]))(indiv,variable))+
          (lambdas((variable*2+1),cluster)*(as<NumericMatrix>(dIpro_v[cluster]))(indiv,variable));
      }
    }
    IDX[indiv]=which_min(DiToClu(indiv,_))+1;
    SSQ=SSQ+(min(DiToClu(indiv,_)));    
  }
   List resu = List::create(Named("DiToClu") = DiToClu , _["IDX"] = IDX,_["SSQ"]=SSQ);
   return resu;
}
/////////////////
// [[Rcpp::export]]
S4 c_COMP_SEVERAL_MEANS(S4 y, NumericMatrix w){//one weight for each histogram and each variable
  S4 x=clone(y);
  int ind=((as<ListMatrix>(x.slot("M"))).nrow());
  int vars=((as<ListMatrix>(x.slot("M"))).ncol());
  int j;
  ListMatrix cont=x.slot("M");
  CharacterVector ch = colnames(cont);
  S4 MEANS("MatH");
  ListMatrix mMEANS(1,vars);
  rownames(mMEANS)=CharacterVector::create("Average");
  colnames(mMEANS)=colnames(cont);
  for (j=0;j<vars;j++){
    //extract variable
    S4 part("MatH");
    ListMatrix tmp(ind,1);
    tmp(_,0)=cont(_,j);
    rownames(tmp)=rownames(cont);
    
    colnames(tmp)=CharacterVector::create(ch[j]);
    part.slot("M")=tmp;
    //end extract
    //extract weights
    NumericVector WEI;
    if (w.ncol()>1){
      WEI=w(_,j);
    }else{
      WEI=w(_,0);
    }
    //end extract weights
    S4 o("distributionH");
    o=MEDIA_V(part,WEI);
    mMEANS(0,j)=o;
  }
  MEANS.slot("M")=mMEANS;
  
  return MEANS;
}
// [[Rcpp::export]]
List c_COMP_SEVERAL_MEANS_M_D(S4 y, NumericMatrix wM, NumericMatrix wDis){
  S4 x=clone(y);
  //one weight for each histogram and each variable
  int ind=((as<ListMatrix>(x.slot("M"))).nrow());
  int vars=((as<ListMatrix>(x.slot("M"))).ncol());
  int i,j;
  ListMatrix cont=x.slot("M");
  CharacterVector ch = colnames(cont);
  S4 MEANS("MatH");
  ListMatrix mMEANS(1,vars);
  NumericMatrix centersM(1,vars);
  NumericMatrix MATMEA(ind,vars);
  rownames(mMEANS)=CharacterVector::create("Average");
  colnames(mMEANS)=colnames(cont);
  for (j=0;j<vars;j++){
    //extract variable
    S4 centrM("MatH");
    ListMatrix tmp(ind,1);
    tmp(_,0)=cont(_,j);
    rownames(tmp)=rownames(cont);
    
    colnames(tmp)=CharacterVector::create(ch[j]);
    centrM.slot("M")=tmp;
    
    NumericVector Means(ind);
    //center data
    for (i=0;i<ind;i++){
      S4 obj=cont(i,j);
      double m=obj.slot("m");
      Means[i]=m;
      MATMEA(i,j)=m;
      NumericVector dom=obj.slot("x");
      dom=dom-m;
      obj.slot("x")=dom;
      obj.slot("m")=0;
      cont(i,j)=obj;
    }
    //end center data
    //end extract
    //extract weights
    NumericVector WEI;
    if (wM.ncol()>1){
      WEI=wM(_,j);
    }else{
      WEI=wM(_,0);
    }
    NumericVector WEI2;
    if (wDis.ncol()>1){
      WEI2=wDis(_,j);
    }else{
      WEI2=wDis(_,0);
    }
    
    //end extract weights
    S4 o("distributionH");
    o=MEDIA_V(centrM,WEI2);
    mMEANS(0,j)=o;
    double tmpM;
    tmpM=sum(WEI*Means)/sum(WEI);
    centersM(0,j)=tmpM;
  }
  x.slot("M")=cont;
  MEANS.slot("M")=mMEANS;
  List resu; 
  resu["Cmeans"]=MEANS;//the MatH of centered mean histos
  resu["centr"]=centersM;//the centers of the mean histos
  resu["Cmath"]=x;//the centred matrix MatH
  resu["MatMEA"]=MATMEA;// the matrix of means
  return resu;
}

// [[Rcpp::export]]
List c_Wass_Q_dist_DET(S4 o1, S4 o2){
  List obj;
  
  NumericVector p1=o1.slot("p");
  NumericVector p2=o2.slot("p");
  NumericVector dom1,dom2,w,cen1,cen2,rad1,rad2,p;
   if ((p1.size()==p2.size()) && (sum(abs(p1-p2))<1e-20)){
     dom1=o1.slot("x");
     dom2=o2.slot("x");
     p=p1;
   }else{
  obj=as<List>(REGISTER2(o1,o2));
  
  //  Rcout<<"\n ok ";
  NumericMatrix M=obj[2];
  //S4 on1=obj[0];
  //S4 on2=obj[1];
  dom1=M(_,0);
  dom2=M(_,1);
  p=M(_,2);
   }
  //dom1=on1.slot("x");
  //dom2=on2.slot("x");
  //w=diff(as<NumericVector>(dom2.slot("p")));
  w=diff(p);
  rad1=diff(dom1)*0.5;
  rad2=diff(dom2)*0.5;
  cen1=(dom1[Range(0,(dom1.size()-2))]+dom1[Range(1,(dom1.size()-1))])*0.5;
  cen2=(dom2[Range(0,(dom1.size()-2))]+dom2[Range(1,(dom1.size()-1))])*0.5;
  
  double m1,m2,s1,s2;
  m1=o1.slot("m");
  m2=o2.slot("m");
  s1=o1.slot("s");
  s2=o2.slot("s");
  
  
  double dm2,ds2,ditot,rho;
  
  dm2=pow((m1-m2),2);
  ds2=pow((s1-s2),2);
  ditot=sum(w*((cen1-cen2)*(cen1-cen2)+ 1/3.0*(rad1-rad2)*(rad1-rad2)));
  rho=(sum(w*((cen1*cen2)+ 1/3.0*(rad1*rad2)))-m1*m2)/(s1*s2);
  
  List resu;resu["D"]=ditot;resu["Dm"]=dm2;resu["Dv"]=(ditot-dm2);resu["Ds"]=ds2;
  resu["Dsh"]=(ditot-dm2-ds2);resu["rho"]=rho;
  return resu;
}

// [[Rcpp::export]]
List c_WH_ADPT_KMEANS_TOTALSSQ(S4 x, NumericMatrix memb, double m,
                               NumericMatrix lambdas,S4 proto){//theta how to use it? 
  int ind=((as<ListMatrix>(x.slot("M"))).nrow());
  int vars=((as<ListMatrix>(proto.slot("M"))).ncol());
  ListMatrix proto_cont=proto.slot("M");
  ListMatrix x_cont=x.slot("M");
  int k=memb.ncol();
  int i,j,clu;
  // #Compute general weights for the global PROTOTYPE
  NumericVector mu(k);
  NumericMatrix MEMBM(ind,k);
  for (i=0;i<memb.nrow();i++){
    for (j=0;j<memb.ncol();j++){
      MEMBM(i,j)=pow(memb(i,j),m);
    }
  }
  mu=colSums(MEMBM);
  NumericMatrix mus(k,vars);
  for(clu=0;clu<k;clu++){
    NumericVector tmp(k, (mu[clu]));
    mus(clu,_)=tmp;  
  }
  
  S4 protoGEN("MatH");
  protoGEN=c_COMP_SEVERAL_MEANS(proto,mus);
  //#Compute general weights for the global PROTOTYPE2
  ListMatrix protoGEN_cont=protoGEN.slot("M");
  
  NumericMatrix W_m(ind,vars);
  NumericMatrix W_d(ind,vars);
  for (i=0;i<ind;i++){
    for (j=0;j<vars;j++){
      for (clu=0;clu<k;clu++){
        W_m(i,j)=W_m(i,j)+(MEMBM(i,clu)*lambdas(j*2,clu));
        W_d(i,j)=W_d(i,j)+(MEMBM(i,clu)*lambdas((j*2+1),clu));
      }
    }
  }
  
  List resu;
  resu=c_COMP_SEVERAL_MEANS_M_D(x, W_m, W_d);
  
  S4 MEANS=resu["Cmeans"];//the MatH of centered mean histos
  ListMatrix MEANS_cont=MEANS.slot("M");
  NumericMatrix centersM=resu["centr"];//the centers of the mean histos
  S4 cen_x=resu["Cmath"];//the centred matrix MatH
  ListMatrix cen_x_cont=cen_x.slot("M");
  NumericMatrix MATMEA=resu["MatMEA"];// the matrix of means
  
  //#Compute BETWEEN
  NumericMatrix BSQ_clu_m(vars,k);
  NumericMatrix BSQ_clu_d(vars,k);
  for (j=0;j<vars;j++){
    double GENPROTm=centersM(0,j);
    S4 o=MEANS_cont(0,j);
    for(clu=0;clu<k;clu++){
      S4 locprot=proto_cont(clu,j);
      double locm=locprot.slot("m");
      BSQ_clu_m(j,clu)=sum(memb(_,clu))*lambdas((j*2),clu)*(GENPROTm-locm)*(GENPROTm-locm);
      List tmpd= c_Wass_Q_dist_DET(o, locprot);
      double tmpdv=tmpd["Dv"];
      BSQ_clu_d(j,clu)=sum(memb(_,clu))*lambdas((j*2+1),clu)*tmpdv;
      
    }
  }
  
  // #Compute the total fuzzy sum of SQUARES
  NumericMatrix TSQ_clu_m(vars,k);
  NumericMatrix TSQ_clu_d(vars,k);
  for (i=0;i<ind;i++){
    for (clu=0;clu<k;clu++){
      for (j=0;j<vars;j++){
        S4 o1=x_cont(i,j);
        S4 prog=protoGEN_cont(0,j);
        List tmpd= c_Wass_Q_dist_DET(o1, prog);
        double dM=tmpd["Dm"];
        double dV=tmpd["Dv"];
        TSQ_clu_m(j,clu)=TSQ_clu_m(j,clu)+(pow(memb(i,clu),m)*lambdas((j*2),clu)*dM);
        TSQ_clu_d(j,clu)=TSQ_clu_d(j,clu)+(pow(memb(i,clu),m)*lambdas((j*2+1),clu)*dV);
        
        
      } 
    }  
  }
  
  NumericMatrix TSQ2_clu_m(vars,k);
  NumericMatrix TSQ2_clu_d(vars,k);
  for (i=0;i<ind;i++){
    for (clu=0;clu<k;clu++){
      for (j=0;j<vars;j++){
        S4 o1=x_cont(i,j);
        S4 prog2=MEANS_cont(0,j);
        List tmpd= c_Wass_Q_dist_DET(o1, prog2);
        double dM=pow((MATMEA(i,j)-centersM(0,j)),2);
        double dV=tmpd["Dv"];
        TSQ2_clu_m(j,clu)=TSQ2_clu_m(j,clu)+(pow(memb(i,clu),m)*lambdas((j*2),clu)*dM);
        TSQ2_clu_d(j,clu)=TSQ2_clu_d(j,clu)+(pow(memb(i,clu),m)*lambdas((j*2+1),clu)*dV);
        
        
      } 
    }  
  }
  List resu2; resu2["gen"]=protoGEN;  resu2["BSQ_m"]=BSQ_clu_m;resu2["BSQ_v"]=BSQ_clu_d;
  //resu2["SEC"]=resu;
  //resu2["Wm"]=W_m;resu2["Wd"]=W_d;
  resu2["TSQ_m"]=TSQ_clu_m;resu2["TSQ_v"]=TSQ_clu_d;
  resu2["TSQ2_m"]=TSQ2_clu_m;resu2["TSQ2_v"]=TSQ2_clu_d;
  return resu2;
}

// functions for fuzzycmeans
// [[Rcpp::export]]
double c_WH_ADPT_FCMEANS_SSQ(S4 x,NumericMatrix memb,double m,
                             NumericMatrix lambdas, S4 proto, double theta){
  int ind=((as<ListMatrix>(x.slot("M"))).nrow());
  int vars=((as<ListMatrix>(proto.slot("M"))).ncol());
  ListMatrix proto_cont=proto.slot("M");
  ListMatrix x_cont=x.slot("M");
  int k=memb.ncol();
  int i,j,clu;
  double SSQ=0;
  for (i=0;i<ind;i++){
    for (clu=0;clu<k;clu++){
      for (j=0;j<vars;j++){
        S4 o1=x_cont(i,j);
        S4 pro_clu=proto_cont(clu,j);
        List tmpd= c_Wass_Q_dist_DET(o1, pro_clu);
        double dM=tmpd["Dm"];
        double dV=tmpd["Dv"];
        SSQ=SSQ+(pow(memb(i,clu),m)*(pow(lambdas((j*2),clu),theta)*dM+
          pow(lambdas((j*2+1),clu),theta)*dV));          
      }
    }
  }
  return SSQ;
}
// [[Rcpp::export]]
NumericMatrix c_ADA_F_WHEIGHT(List distances,int k, int vars, int ind, int schema, int weight_sys,
                              double theta, double m){
  NumericMatrix lambdas((2*vars),k);
  NumericMatrix DiM=distances[0];
  NumericMatrix DiV=distances[1];
  int j,clu;
  double sr;
  for (clu=0;clu<k;clu++){
    sr=sum(DiM(_,clu));
    if (sr<1e-20){DiM(_,clu)=DiM(_,clu)+1e-20;}
    sr=sum(DiV(_,clu));
    if (sr<1e-20){DiV(_,clu)=DiV(_,clu)+1e-20;}
  }
  
  
  double num,denom;
  
  
  for (j=0;j<vars;j++){
    for (clu=0;clu<k;clu++){
      //product
      if (weight_sys==1){
        if (schema==1){//one weigth for one variable 1W GS Global-sum
          if (clu==0){
            NumericVector tD_M=rowSums(DiM);
            num=pow((Provec(tD_M)),(1/double(vars)));
            denom=sum(DiM(j,_));
            if (denom<1e-10){denom=1e-10;}
            lambdas(j*2,clu)=num/denom;
            lambdas((j*2+1),clu)=num/denom;
          }else{
            lambdas(j*2,clu)=lambdas(j*2,0);
            lambdas((j*2+1),clu)=lambdas((j*2+1),0);
          }
        }
        if (schema==2){//two weigths for the two components for each variable
          if (clu==0){
            NumericVector tD_M=rowSums(DiM);
            NumericVector tD_V=rowSums(DiV);
            num=pow((Provec(tD_M)),(1/double(vars)));
            denom=sum(DiM(j,_));
            if (denom<1e-10){denom=1e-10;}
            lambdas(j*2,clu)=num/denom;
            num=pow((Provec(tD_V)),(1/double(vars)));
            denom=sum(DiV(j,_));
            if (denom<1e-10){denom=1e-10;}
            lambdas((j*2+1),clu)=num/denom;
          }else{
            lambdas(j*2,clu)=lambdas(j*2,0);
            lambdas((j*2+1),clu)=lambdas((j*2+1),0);
          }
        }
        if (schema==3){//a weigth for one variable and for each cluster
          num=pow((Provec(DiM(_,clu))),(1/double(vars)));
          denom=DiM(j,clu);
          if (denom<1e-10){denom=1e-10;}
          lambdas(j*2,clu)=num/denom;
          lambdas((j*2+1),clu)=num/denom;
        }
        if (schema==4){//two weigths for the two components for each variable and each cluster
          num=pow((Provec(DiM(_,clu))),(1/double(vars)));
          denom=DiM(j,clu);
          if (denom<1e-10){denom=1e-10;}
          lambdas(j*2,clu)=num/denom;

          num=pow((Provec(DiV(_,clu))),(1/double(vars)));
          denom=DiV(j,clu);
          if (denom<1e-10){denom=1e-10;}
          lambdas((j*2+1),clu)=num/denom;
        }
        if (schema==5){//two weigths for the two components for each variable
          // multiplied all equal one
          if (clu==0){
            NumericVector tD_M=rowSums(DiM);
            NumericVector tD_V=rowSums(DiV);
            double num1,num2,denom1,denom2;

            num1=pow((Provec(tD_M)),(1/double(2*vars)));
            denom1=sum(DiM(j,_));
            if (denom1<1e-10){denom1=1e-10;}
            num2=pow((Provec(tD_V)),(1/double(2*vars)));
            denom2=sum(DiV(j,_));
            if (denom2<1e-10){denom2=1e-10;}
            lambdas(j*2,clu)=(num1*num2)/denom1;
            lambdas((j*2+1),clu)=(num1*num2)/denom2;
          }else{
            lambdas(j*2,clu)=lambdas(j*2,0);
            lambdas((j*2+1),clu)=lambdas((j*2+1),0);            }
        }
        if (schema==6){//two weigths for the two components for each variable and each cluster
          // multiplied all equal one
          num=pow((Provec(DiM(_,clu))*Provec(DiV(_,clu))),1/double(2*vars));
          denom=DiM(j,clu);
          if (denom<1e-10){denom=1e-10;}
          lambdas(j*2,clu)=num/denom;
          denom=DiV(j,clu);
          if (denom<1e-10){denom=1e-10;}
          lambdas((j*2+1),clu)=num/denom;
        }
      }
      //sum
      if (weight_sys==2){
        if (schema==1){//one weigth for one variable 1W GS Global-sum
          if (clu==0){
            double tmp1=sum(DiM(j,_));
            NumericVector tmpNUM(vars,tmp1);
            //             num=rep(sum(distances[variables,,1]),vars)
            NumericVector tmpDEN=rowSums(DiM);
            lambdas(j*2,clu)=1/(
              sum(
                pow((tmpNUM/tmpDEN),(1/(theta-1)))
              )
            );

            lambdas((j*2+1),clu)=lambdas(j*2,clu);
          }else{
            lambdas(j*2,clu)=lambdas(j*2,0);
            lambdas((j*2+1),clu)=lambdas((j*2+1),0);
          }


        }
        if (schema==2){//two weigths for the two components for each variable 2W GS Global-sum
          if (clu==0){
            double tmp1=sum(DiM(j,_));
            NumericVector tmpNUM(vars,tmp1);
            NumericVector tmpDEN=rowSums(DiM);
            lambdas(j*2,clu)=1/(
              sum(
                pow((tmpNUM/tmpDEN),(1/(theta-1)))
              )
            );

            double tmp2=sum(DiV(j,_));
            NumericVector tmpNUM2(vars,tmp2);
            NumericVector tmpDEN2=rowSums(DiV);
            lambdas((j*2+1),clu)=1/(
              sum(
                pow((tmpNUM2/tmpDEN2),(1/(theta-1)))
              )
            );

          }else{
            lambdas(j*2,clu)=lambdas(j*2,0);
            lambdas((j*2+1),clu)=lambdas((j*2+1),0);
          }
         
        }
        if (schema==3){//a weigth for one variable and for each cluster 1W LS Local-sum
          double tmp1=DiM(j,clu);
          NumericVector tmpNUM(vars,tmp1);
          //           num=rep(distances[variables,cluster,1],vars)
          NumericVector tmpDEN=(DiM(_,clu));
          //           den=distances[,cluster,1]
          lambdas(j*2,clu)=1/(
            sum(
              pow((tmpNUM/tmpDEN),(1/(theta-1)))
            )
          );
          lambdas((j*2+1),clu)=lambdas(j*2,clu);
          //         lambdas[(variables*2-1),cluster]=1/sum((num/den)^(1/(theta-1)))
          //           lambdas[(variables*2),cluster]=lambdas[(variables*2-1),cluster]
          
          
        }
        if (schema==4){//two weigths for the two components for each variable and each cluster 2W LS Local-sum
          double tmp1=DiM(j,clu);
          NumericVector tmpNUM(vars,tmp1);
          //           num=rep(distances[variables,cluster,1],vars)
          NumericVector tmpDEN=(DiM(_,clu));
          //           den=distances[,cluster,1]
          lambdas(j*2,clu)=1/(
            sum(
              pow((tmpNUM/tmpDEN),(1/(theta-1)))
            )
          );
          
          //           num=rep(distances[variables,cluster,1],vars)
          //           den=distances[,cluster,1]
          //
          //         lambdas[(variables*2-1),cluster]=1/sum((num/den)^(1/(theta-1)))
          double tmp2=DiV(j,clu);
          NumericVector tmpNUM2(vars,tmp2);
          //           num=rep(distances[variables,cluster,2],vars)
          NumericVector tmpDEN2=(DiV(_,clu));
          //           den=distances[,cluster,2]
          lambdas((j*2+1),clu)=1/(
            sum(
              pow((tmpNUM2/tmpDEN2),(1/(theta-1)))
            )
          );
        }
        if (schema==5){//two weigths for the two components for each variable 2W GS Global-sum
          if (clu==0){
            NumericVector unonum(vars,sum(DiM(j,_)));
            //             unonum=sum(distances[variables,,1])
            NumericVector unoden=rowSums(DiM);
            //             unoden=apply(distances[,,1],MARGIN=c(1),sum)
            double uno=sum(pow((unonum/unoden),(1/(theta-1))));
            //             uno=sum((unonum/unoden)^(1/(theta-1)))
            NumericVector dueden=rowSums(DiV);
            //             dueden=apply(distances[,,2],MARGIN=c(1),sum)
            double due=sum(pow((unonum/dueden),(1/(theta-1))));
            //             due=sum((unonum/dueden)^(1/(theta-1)))
            lambdas(j*2,clu)=1/(uno+due);
            //             lambdas[(variables*2-1),cluster]=1/(uno+due)
            
            NumericVector unonum2(vars,sum(DiV(j,_)));
            //             unonum=sum(distances[variables,,2])
            double tre=sum(pow((unonum2/unoden),(1/(theta-1))));
            //             tre=sum((unonum/unoden)^(1/(theta-1)))
            double quattro=sum(pow((unonum2/dueden),(1/(theta-1))));
            //             quattro=sum((unonum/dueden)^(1/(theta-1)))
            lambdas((j*2+1),clu)=1/(tre+quattro);
            //             lambdas[(variables*2),cluster]=1/(tre+quattro)
          }else{
            lambdas(j*2,clu)=lambdas(j*2,0);
            lambdas((j*2+1),clu)=lambdas((j*2+1),0);            }
          
        }
        if (schema==6){//two weigths for the two components for each variable
          //and each cluster 2W LS Local-sum
          NumericVector unonum(vars,DiM(j,clu));
          //           unonum=distances[variables,cluster,1]
          NumericVector unoden=DiM(_,clu);
          //         unoden=distances[,cluster,1]
          double uno=sum(pow((unonum/unoden),(1/(theta-1))));
          //         uno=sum((unonum/unoden)^(1/(theta-1)))
          NumericVector dueden=DiV(_,clu);
          //           dueden=distances[,cluster,2]
          double due=sum(pow((unonum/dueden),(1/(theta-1))));
          //         due=sum((unonum/dueden)^(1/(theta-1)))
          lambdas(j*2,clu)=1/(uno+due);
          //           lambdas[(variables*2-1),cluster]=1/(uno+due)
          NumericVector unonum2(vars,sum(DiV(j,clu)));
          //           unonum=distances[variables,cluster,2]
          double tre=sum(pow((unonum2/unoden),(1/(theta-1))));
          //         tre=sum((unonum/unoden)^(1/(theta-1)))
          double quattro=sum(pow((unonum2/dueden),(1/(theta-1))));
          //           quattro=sum((unonum/dueden)^(1/(theta-1)))
          lambdas((j*2+1),clu)=1/(tre+quattro);
          //           lambdas[(variables*2),cluster]=1/(tre+quattro)
          
        }
      }
    }
  }
  return lambdas;
}




// [[Rcpp::export]]
NumericMatrix c_MEMB_comp(int ind, int vars, int k, 
                          NumericMatrix lambdas,
                          List diINDtoPROT_M,
                          List diINDtoPROT_V,
                          double m,double  theta){
  NumericMatrix memb(ind,k);
  NumericMatrix TMP(ind,k);
  
  int i,j,clu;
  
  
  for (clu=0;clu<k;clu++){
    
    NumericMatrix DM_k=diINDtoPROT_M[clu];
    NumericMatrix DV_k=diINDtoPROT_V[clu];
    
    for (i=0;i<ind;i++){
      
      for (j=0;j<vars;j++){
        TMP(i,clu)=TMP(i,clu)+
          (pow(lambdas((j*2),clu),theta)*DM_k(i,j)+
          pow(lambdas((j*2+1),clu),theta)*DV_k(i,j));
      }
    }
  }
  
  for (i=0;i<ind;i++){
    for (clu=0;clu<k;clu++){
      if(TMP(i,clu)>0){
        NumericVector tmpi=TMP(i,clu)/TMP(i,_);
        memb(i,clu)=1/sum(pow(tmpi,(1/(m-1))));
        
        //         memb[indiv,cluster]=1/sum(
        //(d1[cluster]/d1)^(1/(m-1))
        //  )}        
      }else{
        memb(i,clu)=1;
      }
    }
  }
  // List resu; resu["memb"]=memb;resu["dis"]=TMP;
  // return resu;
  return memb;
}


// [[Rcpp::export]]
List c_ComputeFastSSQ_Fuzzy_1V(NumericMatrix M, NumericVector memb, double m){
  int ind=(M.ncol()-1);
  int rr=M.nrow();
  double SSQ=0;
  NumericVector W=pow(memb,m);
  NumericVector W2=W/sum(W);
  NumericVector domM(rr);
  int r,i;
  //compute the prototype
  for (r=0;r<rr;r++){
    NumericVector tmp=M(r,_);
    domM[r]=sum(tmp[Range(0,(ind-1))]*W2);
   }
  NumericVector p=diff(M(_,ind));
  NumericVector rm=diff(domM)*0.5;
  NumericVector cm=(domM[Range(0,rr-2)]+domM[Range(1,rr-1)])*0.5;
  double mm=sum(cm*p);
  double ss=sqrt((sum(cm*cm*p)+sum(p*rm*rm/3.0))-(mm*mm));
  
  S4 prot("distributionH");
  prot.slot("x")=domM;
  prot.slot("p")=M(_,ind);
  prot.slot("m")=mm;
  prot.slot("s")=ss;
  //end compute the prototype
  
  //compute SSQ
  for (i=0;i<ind;i++){
    NumericVector dom=M(_,i);
    NumericVector ri=diff(dom)*0.5;
    //     ri=(subMM[2:rr,indiv]-subMM[1:(rr-1),indiv])/2
    NumericVector ci=(dom[Range(0,rr-2)]+dom[Range(1,rr-1)])*0.5;
    //     ci=(subMM[2:rr,indiv]+subMM[1:(rr-1),indiv])/2
   
    SSQ=SSQ+((sum(p*((ci-cm)*(ci-cm)+(ri-rm)*(ri-rm)/3.0)))*W[i]);
    //     SSQ=SSQ+(sum(p1*((ci-cm)^2)+ p1/3*((ri-rm)^2)))*W[indiv]
    
  }
  //end compute SSQ
  
  List resu;resu["SSQ"]=SSQ;resu["prot"]=prot;//resu['W']=W;resu['W2']=W2;
  return resu;
}

// ComputeFastSSQ_Fuzzy=function(subMM,memb,m){
//   ind=ncol(subMM)-1
//   rr=nrow(subMM)
//   SSQ=0
//   W=memb^m
//   W2=W/sum(W)
//   m1=subMM[,1:ind]%*%as.matrix(W2)
//   
//   m1=as.vector(m1)
//   p1=subMM[2:rr,ind+1]-subMM[1:(rr-1),ind+1]
//   cm=(m1[2:rr]+m1[1:(rr-1)])/2
//   rm=(m1[2:rr]-m1[1:(rr-1)])/2
//   for (indiv in 1:ind){
//     ci=(subMM[2:rr,indiv]+subMM[1:(rr-1),indiv])/2
//     ri=(subMM[2:rr,indiv]-subMM[1:(rr-1),indiv])/2
//     SSQ=SSQ+(sum(p1*((ci-cm)^2)+ p1/3*((ri-rm)^2)))*W[indiv]
//   }
//   return(list(SSQ=SSQ,mx=m1, mp=subMM[,ind+1]))
// }


// [[Rcpp::export]]
List c_ComputeFastSSQ_Fuzzy(List MM,NumericMatrix memb, double m, SEXP nr, SEXP nc){
  int k=memb.ncol();
  int vars=MM.size();
  NumericMatrix WSSQ(k,vars);
  int j,clu;
  S4 protos("MatH");
  ListMatrix protos_d(k,vars);
  rownames(protos_d)=nr;
  colnames(protos_d)=nc;
  //AGGIUNGERE I NOMI RIGA E I NOMI COLONNA 
  for(clu=0;clu<k;clu++){
    for (j=0;j<vars;j++){
      NumericMatrix M=MM[j];
      List resup;
      resup=c_ComputeFastSSQ_Fuzzy_1V(M,memb(_,clu),m);
      double ttp=resup["SSQ"];
      WSSQ(clu,j)=ttp;
      protos_d(clu,j)=resup["prot"];
    }
  }
  protos.slot("M")=protos_d;
  List resu;resu["SSQ"]=WSSQ;resu["proto"]=protos;
  return resu;
}



// [[Rcpp::export]]
List c_ComputeFast_L2_SQ_WASS_DMAT(List MM,S4 prot){
  
  NumericMatrix M0=MM[0];
  int ind=(M0.ncol()-1);

  ListMatrix prot_con=prot.slot("M");
  int k=prot_con.nrow();
  int vars=prot_con.ncol();
  ListOf<List> DistaF(vars);
  NumericMatrix CompDista(ind,k);
  
  
  int i,j,clu;
  for (j=0;j<vars;j++){
    NumericMatrix DistM(ind,k);
    NumericMatrix DistV(ind,k);
    NumericMatrix Dista(ind,k);
    NumericMatrix M=MM[j];
    NumericVector cum=M(_,ind);
    NumericVector w=diff(cum);
    int rr=M.nrow();
    
    
    for (clu=0;clu<k;clu++){
      S4 pro=prot_con(clu,j);
      NumericVector domP=pro.slot("x");
      NumericVector cm=(domP[Range(0,rr-2)]+domP[Range(1,rr-1)])*0.5;
      NumericVector rm=diff(domP)*0.5;
      double mm=pro.slot("m");
      
      for (i=0;i<ind;i++){
        NumericVector domi=M(_,i);
        NumericVector ci=(domi[Range(0,rr-2)]+domi[Range(1,rr-1)])*0.5;
        NumericVector ri=diff(domi)*0.5;
        double mi=sum(ci*w);
        double tmpd=sum(w*((cm-ci)*(cm-ci)+(rm-ri)*(rm-ri)*1/3.0));
        Dista(i,clu)=tmpd;
        DistM(i,clu)=(mm-mi)*(mm-mi);
        DistV(i,clu)=tmpd-DistM(i,clu);
        CompDista(i,clu)=CompDista(i,clu)+tmpd;
      }
      
    }
    List D;
    D["D"]=Dista;
    D["DM"]=DistM;
    D["DV"]=DistV;
    DistaF[j]=D;
    
  }
  
  //Distances from a point and a center
  // DET with details for each variable  
  List resu; resu["Dist"]=CompDista; resu["DET"]=DistaF;
  return resu;
}

// [[Rcpp::export]]
List c_DISTA_ADA(S4 x, List MM, S4 proto,
                 int vars, int ind,int k, 
                 NumericMatrix memb,double m, int schema){
  NumericMatrix distances_M(vars,k);
  NumericMatrix distances_V(vars,k);
  List diINDtoPROT_M(k);
  List diINDtoPROT_V(k);
  
  
  
  int clu,j,i;
  for (clu=0;clu<k;clu++){
    NumericMatrix tmp(ind,vars);
    NumericMatrix tmp2(ind,vars);
    diINDtoPROT_M[clu]=tmp;
    diINDtoPROT_V[clu]=tmp2;
  }
  
  
  List resup;
  resup=c_ComputeFast_L2_SQ_WASS_DMAT(MM,proto);
  List resudet=resup["DET"];
  
  for (j=0;j<vars;j++){
    List DD=resudet[j];
    
    NumericMatrix tmpD=DD["D"];
    NumericMatrix tmpDM=DD["DM"];
    NumericMatrix tmpDV=DD["DV"];
    for (clu=0;clu<k;clu++){

      for (i=0;i<ind;i++){
        double tmp_d=tmpD(i,clu);
        double tmp_d_m=tmpDM(i,clu);
        double tmp_d_v=tmpDV(i,clu);
        if (schema==1 || schema==3){
          distances_M(j,clu)=distances_M(j,clu)+tmp_d*(pow(memb(i,clu),m));
          (as<NumericMatrix>(diINDtoPROT_M[clu]))(i,j)=tmp_d;
        }
        if (schema==2 || schema==4|| schema==5||schema==6){
          distances_M(j,clu)=distances_M(j,clu)+tmp_d_m*(pow(memb(i,clu),m));
          distances_V(j,clu)=distances_V(j,clu)+tmp_d_v*(pow(memb(i,clu),m));
          (as<NumericMatrix>(diINDtoPROT_M[clu]))(i,j)=tmp_d_m;
          (as<NumericMatrix>(diINDtoPROT_V[clu]))(i,j)=tmp_d_v;
        }

      }
    }
  }
  
  //   return(res=list(distances=distances,diINDtoPROT=diINDtoPROT))
  List resu;
  resu["dista_M"]=distances_M;
  resu["dista_V"]=distances_V;
  resu["dista_to_pro_M"]=diINDtoPROT_M;
  resu["dista_to_pro_V"]=diINDtoPROT_V;
  return resu;
}


// [[Rcpp::export]]
NumericMatrix c_ComputeFast_L2_SQ_WASS_DMAT_simple(List MM,S4 prot){
  NumericMatrix M0=MM[0];
  int ind=(M0.ncol()-1);
  ListMatrix prot_con=prot.slot("M");
  int k=prot_con.nrow();
  int vars=prot_con.ncol();
  NumericMatrix CompDista(ind,k);
  int i,j,clu;
  for (j=0;j<vars;j++){
    NumericMatrix M=MM[j];
    NumericVector cum=M(_,ind);
    NumericVector w=diff(cum);
    int rr=M.nrow();
    for (clu=0;clu<k;clu++){
      S4 pro=prot_con(clu,j);
      NumericVector domP=pro.slot("x");
      NumericVector cm=(domP[Range(0,rr-2)]+domP[Range(1,rr-1)])*0.5;
      NumericVector rm=diff(domP)*0.5;
      
      for (i=0;i<ind;i++){
        NumericVector domi=M(_,i);
        NumericVector ci=(domi[Range(0,rr-2)]+domi[Range(1,rr-1)])*0.5;
        NumericVector ri=diff(domi)*0.5;
        double tmpd=sum(w*((cm-ci)*(cm-ci)+(rm-ri)*(rm-ri)*1/3.0));
        CompDista(i,clu)=CompDista(i,clu)+tmpd;
      }
      
    }
  }
  return CompDista;
}

// [[Rcpp::export]]
List c_ComputeFast_Fuzzy_TOT_SSQ(List MM,
                                 S4 prot,
                                 NumericMatrix memb, 
                                 double m){
  NumericMatrix M0=MM[0];
  int ind=(M0.ncol()-1);
  ListMatrix prot_con=prot.slot("M");
  int vars=MM.size();
  int k=memb.ncol();
  
  // #computeGenProt
  NumericVector mus(k,0.0);
  
  int i,j,clu;
  for (i=0;i<ind;i++){
    for (clu=0;clu<k;clu++){
      mus[clu]=mus[clu]+(pow(memb(i,clu),m));
    }
  }
  mus=mus/sum(mus);
 
  S4 GP("MatH");
  ListMatrix GP_cont(1,vars);
  colnames(GP_cont)=colnames(prot_con);
  rownames(GP_cont)=CharacterVector::create("Average");
  // I NOMI
  for(j=0;j<vars;j++){
    S4 o1=prot_con(0,j);
    NumericVector dm=(o1.slot("x"));
    NumericVector p=(o1.slot("p"));
    dm=dm*mus[0];
    for(clu=1;clu<k;clu++){
      S4 o=prot_con(clu,j);
      NumericVector cum=o.slot("x");
      dm=dm+(cum*mus[clu]);

    }

    S4 mm("distributionH");
      mm.slot("x")=dm;
      mm.slot("p")=p;
    double m1,s1;
//    NumericVector ccm;
//    ccm=(dm[Range(0,(dm.size()-2))]+dm[Range(1,(dm.size()-1))])*0.5;
    NumericVector TMPR=M_STD_H(mm);
    m1=TMPR[0];
    s1=TMPR[1];
//    m1=sum(ccm*diff(p));
//    s1=sqrt(sum(ccm*ccm*diff(p)+pow((diff(dm)*0.5),2)*diff(p))-m1*m1);
   
    mm.slot("m")=m1;
    mm.slot("s")=s1;
    
    GP_cont(0,j)=mm;
  }
  // #compute TOTALSSQ
  NumericMatrix SSQ_M(vars,k);
  NumericMatrix SSQ_V(vars,k);
//  NumericMatrix SSQ_T(vars,k);
  //   SSQ_det=array(0,dim = c(2,var,clu), dimnames=list(c('Mean','Variability'),colnames(proto@M),rownames(proto@M)))
  for (j=0;j<vars;j++){
    NumericMatrix M=MM[j];
    NumericVector cum=M(_,ind);
    NumericVector w=diff(cum);
    S4 o=GP_cont(0,j);
    NumericVector dm=o.slot("x");
    NumericVector cm=(dm[Range(0,(dm.size()-2))]+dm[Range(1,(dm.size()-1))])*0.5;
    NumericVector rm=diff(dm)*0.5;
//    double mea_m=o.slot("m");

    for(i=0;i<ind;i++){
      NumericVector domi=M(_,i);
      NumericVector ci=(domi[Range(0,domi.size()-2)]+domi[Range(1,domi.size()-1)])*0.5;
      NumericVector ri=diff(domi)*0.5;
      double tmpd;//=sum(w*(pow((cm-ci),2)+pow((rm-ri),2)*1/3.0));
      tmpd=sum(w*((cm-ci)*(cm-ci)+ 1/3.0*(rm-ri)*(rm-ri)));
      double tmpd_m=pow((sum(cm*w)-sum(ci*w)),2);
      double tmpd_v=tmpd-tmpd_m;
      for (clu=0;clu<k;clu++){
//        SSQ_T(j,clu)=SSQ_T(j,clu)+(pow(memb(i,clu),m)*tmpd);
        SSQ_M(j,clu)=SSQ_M(j,clu)+(pow(memb(i,clu),m)*tmpd_m);
        SSQ_V(j,clu)=SSQ_V(j,clu)+(pow(memb(i,clu),m)*tmpd_v);
     }
      
    }
  }
 
  GP.slot("M")=GP_cont;
  List resu; resu["SSQ"]=sum(SSQ_M)+sum(SSQ_V);
//  resu["SSQ_T"]=SSQ_T;
  resu["SSQ_M"]=SSQ_M;
  resu["SSQ_V"]=SSQ_V;resu["ProtoGEN"]=GP;
  return resu;
}

// [[Rcpp::export]]
double c_WH_ADPT_FCMEANS_SSQ_FAST(List MM, S4 x,NumericMatrix memb,double m,
                             NumericMatrix lambdas, S4 proto, double theta){
  NumericMatrix M0=MM[0];
  int ind=(M0.ncol()-1);
  ListMatrix prot_con=proto.slot("M");
  int k=prot_con.nrow();
  int vars=prot_con.ncol();
  
  ListMatrix x_cont=x.slot("M");
  
  int i,j,clu;
  double SSQ=0;
  for (j=0;j<vars;j++){
    NumericMatrix M=MM[j];
    NumericVector cum=M(_,ind);
    NumericVector w=diff(cum);
    
    for (clu=0;clu<k;clu++){
      S4 pro_clu=prot_con(clu,j);
      NumericVector dom_p=pro_clu.slot("x");
      List tmpr=c_cen_rad(dom_p);
      NumericVector c_p=tmpr["cen"];
      NumericVector r_p=tmpr["rad"];

      
      
    for (i=0;i<ind;i++){
      S4 o1=x_cont(i,j);
      NumericVector dom_o=o1.slot("x");
      List tmpo=c_cen_rad(dom_o);
      NumericVector c_o=tmpo["cen"];
      NumericVector r_o=tmpo["rad"];
      double Tdist=sum(w*((c_p-c_o)*(c_p-c_o)+(r_p-r_o)*(r_p-r_o)*1/3.0));
  
  
  
        double dM=pow((double (o1.slot("m"))- double (pro_clu.slot("m"))),2.0);
        double dV=Tdist-dM;
        SSQ=SSQ+(pow(memb(i,clu),m)*(pow(lambdas((j*2),clu),theta)*dM+
          pow(lambdas((j*2+1),clu),theta)*dV));          
      }
    }
  }
  return SSQ;
}

// [[Rcpp::export]]
double c_WH_ADPT_FCMEANS_SSQ_FAST_NEW(List DM,List DV, NumericMatrix memb,double m,
                                  NumericMatrix lambdas, double theta){
  int ind=memb.nrow();
  int k=memb.ncol();
  int vars=lambdas.nrow()/2;//NO!!!
  

  int i,j,clu;
  double SSQ=0;
     
    for (clu=0;clu<k;clu++){
      NumericMatrix DM_k=DM[clu];
      NumericMatrix DV_k=DV[clu];
      for (i=0;i<ind;i++){
        for (j=0;j<vars;j++){
          
        double dM=DM_k(i,j);
        double dV=DV_k(i,j);
        SSQ=SSQ+(pow(memb(i,clu),m)*(pow(lambdas((j*2),clu),theta)*dM+
          pow(lambdas((j*2+1),clu),theta)*dV));          
      }
    }
  }
  return SSQ;
}

// [[Rcpp::export]]
double c_dotpW(S4 o1,S4 o2){
  double res=0;
  
  NumericVector p1=o1.slot("p");
  NumericVector p2=o2.slot("p");
  NumericVector dom1,dom2,w,cen1,cen2,rad1,rad2,p;
  if ((p1.size()==p2.size()) && (sum(abs(p1-p2))<1e-20)){
    dom1=o1.slot("x");
    dom2=o2.slot("x");
    p=p1;
  }else{
  List obj;
  obj=as<List>(REGISTER2(o1,o2));
    
    //  Rcout<<"\n ok ";
    NumericMatrix M=obj[2];
    
    //S4 on1=obj[0];
    //S4 on2=obj[1];
    dom1=M(_,0);
    dom2=M(_,1);
    p=M(_,2);
  }
    //dom1=on1.slot("x");
    //dom2=on2.slot("x");
    //w=diff(as<NumericVector>(dom2.slot("p")));
    w=diff(p);
    rad1=diff(dom1)*0.5;
    rad2=diff(dom2)*0.5;
    cen1=(dom1[Range(0,(dom1.size()-2))]+dom1[Range(1,(dom1.size()-1))])*0.5;
    cen2=(dom2[Range(0,(dom1.size()-2))]+dom2[Range(1,(dom1.size()-1))])*0.5;
    
    res=sum(w*((cen1*cen2)+ 1/3.0*(rad1*rad2)));
  
  return res;
  
}

// [[Rcpp::export]]
S4 c_PROTO_KOHONEN(S4 proto,int k, int ind, 
                   List MM,int vars,
                   NumericMatrix KT, NumericVector IDX){
  S4 x=clone(proto);
  ListMatrix x_cont=x.slot("M");
  
  int i,v,clu;
  for (clu=0;clu<k;clu++){
    NumericVector kern(ind);
    for (i=0;i<ind;i++){
      kern[i]=KT((IDX[i]-1),clu);
      if (kern[i]<1e-30) kern[i]=0;
    }
    double sk=sum(kern);
    if (sk>1e-30){
      kern=kern/sum(kern);
    for (v=0;v<vars;v++){
      NumericMatrix MM_tmp=MM[v];
      NumericVector p=MM_tmp(_,ind);
      NumericVector dm(MM_tmp.nrow());
      for (int cc=0;cc<ind;cc++){
        dm=dm+kern[cc]*MM_tmp(_,cc);
      }
      S4 tmp("distributionH");
      tmp.slot("x")=dm;
      tmp.slot("p")=p;
      double m1,s1;
      NumericVector TMPR=M_STD_H(tmp);
      m1=TMPR[0];
      s1=TMPR[1];
      tmp.slot("m")=m1;
      tmp.slot("s")=s1;
      x_cont(clu,v)=tmp;
      
    }
    }
  }
  x.slot("M")=x_cont;
  return x;
}

// [[Rcpp::export]]
List c_Wass_Q_dist_2P(S4 o1, S4 o2){
  List obj;
  
  double m1,m2,s1,s2;
  m1=o1.slot("m");
  m2=o2.slot("m");
  s1=o1.slot("s");
  s2=o2.slot("s");
  
  
  double dm2,ditot,cp;
  
  dm2=(m1-m2)*(m1-m2);
  cp=c_dotpW(o1,o2);
  
  
  ditot=s1*s1+m1*m1+s2*s2+m2*m2-2*cp;

  List resu;resu["D"]=ditot;resu["Dm"]=dm2;resu["Dv"]=(ditot-dm2);

  return resu;
}
