class Bootstrapped_paired_ROC {
private:
  ROC roc [2];
  Sampler_base *sampler; //random number generator, see sampler
public:
  Bootstrapped_paired_ROC(NumericVector pred1, NumericVector pred2, IntegerVector true_class);
  ~Bootstrapped_paired_ROC();
  ROC & get_roc(int number);
  void bootstrap(); //performs a bootstrap sample on the ROC curve
};
