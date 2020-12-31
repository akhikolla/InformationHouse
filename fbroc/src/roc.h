//Class models a ROC curve
class ROC{
private:
  void find_thresholds(NumericVector pred, IntegerVector true_class); //Finds ROC thresholds
  void build_pred(NumericVector pred, IntegerVector true_class); //Separates predictions by class
  IntegerVector build_index(NumericVector pred); //Builds initial link between samples and thresholds
  void reset_delta(); //Helper function, resets delta to 0 for next bootstrap iteration
  void get_positives_delta(); //Helper function calculates 
  void get_positives(); //Helper functions, count total positives at each thresholds 
  void get_rate(); //Calculate TPR and FPR at each threshold
protected:
  NumericVector pred_pos; //Predictions for positive samples
  NumericVector pred_neg; //Predictions for negative samples
   //Thresholds of the roc curve. Uses fact, that bootstrap thresholds are subset of original
  NumericVector thresholds;
  IntegerVector index_pos; //Link sample-threshold. Updates while bootstrapping.
  IntegerVector index_neg;  //Link sample-threshold. Updates while bootstrapping.
  IntegerVector original_index_pos; //Link sample-threshold for original data
  IntegerVector original_index_neg; //Link sample-threshold for original data
  IntegerVector delta_pos; // Contaings change in number of true positives at each threshold
  IntegerVector delta_neg; // Contaings change in number of false positives at each threshold
  IntegerVector true_positives; // Contains number of true positives at each threshold
  IntegerVector false_positives; // Contains number of false positives at each threshold
  NumericVector tpr; // true positive rate
  NumericVector fpr; // false positive rate
  int n; //total number of observations
  int n_thresholds; //number of unique thresholds
  int n_pos; //number of positive samples
  int n_neg; //number of negative samples
public:
  ROC(NumericVector pred, IntegerVector true_class); //constructor
  ROC(); // empty constructor
  void shuffle(IntegerVector &shuffle_pos, IntegerVector &shuffle_neg); //Shuffle (bootstrap)
  void strat_shuffle(IntegerVector &shuffle_pos, IntegerVector &shuffle_neg); //Stratified shuffle
  NumericVector & get_tpr(); // return TPR
  NumericVector & get_fpr(); // return FPR
  NumericVector & get_thresholds(); // return thresholds
  NumericVector get_tpr_at_fpr(NumericVector &steps) const; // return TPR at FPR
  NumericVector get_fpr_at_tpr(NumericVector &steps) const; // return FPR at TPR
  // return TPR at FPR if TPR and FPR are pre-cached in a matrix
  static NumericVector get_tpr_at_fpr(NumericVector &tpr_in, NumericVector &fpr_in, NumericVector &steps);
  static NumericVector get_fpr_at_tpr(NumericVector &tpr_in, NumericVector &fpr_in, NumericVector &steps);
  int get_n_thres() const; // get number of thresholds
};

