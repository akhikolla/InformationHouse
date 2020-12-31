class Sampler_base { //abstract base class
protected:
  int n_pos;
  int n_neg;
  int n;
  IntegerVector shuffled_pos_index;
  IntegerVector shuffled_neg_index;
public:
  virtual void generate() = 0; //sample with replacement
  IntegerVector get_shuffled_index(bool which_class) const; // get index 
  virtual ~Sampler_base();
};

class Sampler_Stratified : public Sampler_base { //class for stratified sampling
public:
  virtual void generate();
  Sampler_Stratified(IntegerVector true_class);
  virtual ~Sampler_Stratified();
};
