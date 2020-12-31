typedef double (* PerfFun)(NumericVector &, NumericVector &, NumericVector &);
enum Measure {AUC, TPR_AT_FPR, FPR_AT_TPR, P_AUC_OVER_FPR, P_AUC_OVER_TPR};

List roc_analysis(NumericVector pred, IntegerVector true_class);
PerfFun pick_measure(Measure measure);
NumericVector get_steps(int n_steps);
