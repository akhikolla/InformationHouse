class Bootstrapped_ROC : public ROC{
private:
  Sampler_base *sampler; //random number generator, see sampler
public:
  Bootstrapped_ROC(NumericVector pred, IntegerVector true_class);
  ~Bootstrapped_ROC();
  void bootstrap(); //performs a bootstrap sample on the ROC curve
};
