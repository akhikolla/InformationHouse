extern "C" {
  #include "cl_min.h"
  #include <attributes.h>
  #include "utils.h"
}


  
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export(name=".cwb_makeall")]]
int cwb_makeall(SEXP x, SEXP registry_dir, SEXP p_attribute){
  
  /* char *progname = "RcppCWB"; */
  
  char *registry_directory = strdup(Rcpp::as<std::string>(registry_dir).c_str());
  char *attr_name = strdup(Rcpp::as<std::string>(p_attribute).c_str());
  char * corpus_id = strdup(Rcpp::as<std::string>(x).c_str());
  int validate = 1;

  ComponentID cid = CompLast;

  corpus = cl_new_corpus(registry_directory, corpus_id);
  
  Rprintf("=== Makeall: processing corpus %s ===\n", corpus_id);
  Rprintf("Registry directory: %s\n", corpus->registry_dir);
  
  Attribute *attribute = cl_new_attribute(corpus, attr_name, ATT_POS);
  do_attribute(attribute, cid, validate);

  Rprintf("========================================\n");
  return 0;
}


// [[Rcpp::export(name=".cwb_huffcode")]]
int cwb_huffcode(SEXP x, SEXP registry_dir, SEXP p_attribute) {
  
  char *registry_directory = strdup(Rcpp::as<std::string>(registry_dir).c_str());
  char *attr_name = strdup(Rcpp::as<std::string>(p_attribute).c_str());
  char * corpus_id = strdup(Rcpp::as<std::string>(x).c_str());
  
  char *output_fn = NULL;
  Attribute *attr;
  
  HCD hc;
  
  int i_want_to_believe = 0;        /* skip error checks? */
  /* int all_attributes = 0; */
  
  /* protocol = stdout;   */             /* 'delayed' init (see top of file) */
  
  
  if ((corpus = cl_new_corpus(registry_directory, corpus_id)) == NULL) {
    Rprintf("Corpus %s not found in registry %s . Aborted.\n", 
            corpus_id,
            (registry_directory ? registry_directory
               : central_corpus_directory()));
    return 1;
  }
  
  if ((attr = cl_new_attribute(corpus, attr_name, ATT_POS)) == NULL) {
    Rprintf("Attribute %s.%s doesn't exist. Aborted.\n",  corpus_id, attr_name);
    return 1;
  }
  
  compute_code_lengths(attr, &hc, output_fn);
  if (! i_want_to_believe) decode_check_huff(attr, corpus_id, output_fn);
  
  cl_delete_corpus(corpus);
  
  return 0;
}


// [[Rcpp::export(name=".cwb_compress_rdx")]]
int cwb_compress_rdx(SEXP x, SEXP registry_dir, SEXP p_attribute) {
  
  char *registry_directory = strdup(Rcpp::as<std::string>(registry_dir).c_str());
  char *attr_name = strdup(Rcpp::as<std::string>(p_attribute).c_str());
  char *corpus_id = strdup(Rcpp::as<std::string>(x).c_str());
  
  Attribute *attr;
  
  char *output_fn = NULL;
  
  #ifdef _WIN32
    int i_want_to_believe = 1;        /* skip error on Windows, for the time being */
  #else
    int i_want_to_believe = 0;        /* do not skip error checks on macOS and Linux */
  #endif
  
  
  
  
  /* debug_output = stderr; */        /* 'delayed' init (see top of file) */
  int debug = 0;
  
  if ((corpus = cl_new_corpus(registry_directory, corpus_id)) == NULL) {
    Rprintf("Corpus %s not found in registry %s . Aborted.\n", 
            corpus_id,
            (registry_directory ? registry_directory
               : central_corpus_directory()));
    cleanup(1);
  }
  
  if ((attr = find_attribute(corpus, attr_name, ATT_POS, NULL)) == NULL) {
    Rprintf("Attribute %s.%s doesn't exist. Aborted.\n", corpus_id, attr_name);
    cleanup(1);
  }
  
  compress_reversed_index(attr, output_fn, corpus_id, debug);
  if (! i_want_to_believe) decompress_check_reversed_index(attr, output_fn, corpus_id, debug);
  
  cleanup(0);
  return 0;                        /* to keep gcc from complaining */
}
