const int Imax=10000;
const double Tol=1e-5;
struct Param{  
  int i0;   // central site index
  const std::vector<std::vector<short> > &ai;
  const std::vector<int> &frq;
  const std::vector<bool> &qj;
  const std::vector<short> &L;
  double lambda;
  double lambdah;
  const std::vector<double> &f1;
  const std::vector<std::vector<double> > &f2;
  double &lzp;
  bool naive;
  bool lzhalf;
};

void f12(int i0, const std::vector<std::vector<short> > &si, 
         const std::vector<int> &frq, std::vector<double> &f1, 
         std::vector<std::vector<double> > &f2, 
         const std::vector<short> &L, bool naive, double pcount);

double pan3(std::vector<double> &peff, int nsnp, int i0, 
          const std::vector<short> &L, const std::vector<short> &ci, 
          std::vector<double> h1, const std::vector<std::vector<double> > &J1,
          bool naive, bool lzhalf);

double lnl_psl(const gsl_vector *v, void *params);

void dlnl_psl(const gsl_vector *v, void *params, gsl_vector *df);

void ln_dln_psl(const gsl_vector *x, void *params, double *f, gsl_vector *df);

double lpr_psl(int i0, const std::vector<std::vector<short> > &si,
               const std::vector<int> &frq, const std::vector<bool> &qj, 
    const std::vector<short> &L, double lambda, double lambdah, std::vector<double> &h, 
    std::vector<std::vector<double> > &J, int nprint,
    unsigned int imax, double tol, int verbose, double &lzp, 
    bool naive, bool &failed, bool lzhalf);

double pan2(int nsnp, int i0, int L, const std::vector<short> &ci, 
    const std::vector<double> &h1, const std::vector<std::vector<double> > &J1,
    double &lzp, bool naive);

void invC(const std::vector<std::vector<short> > &ai, 
          const std::vector<int> &frq, const std::vector<short> &L, 
          double &E, double &lnz, std::vector<std::vector<double> > &h,
            std::vector<std::vector<std::vector<double> > > &J, double eps,
            double priorCount);
