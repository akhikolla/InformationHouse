extern "C" {
  #include <stdio.h>
  #include <stdlib.h>
  #include <unistd.h>
  #include <string.h>
  #include "cl_min.h"
  #include <pcre.h>
/*  #include "globals.h" */
/*  #include "context_descriptor.h" */
/*  #include "_options.h" */
/*  #include "cqp.h" */
/*  #include "corpmanag.h" */
  #include "server.h"
  #include "globalvars.h"
  #include "env.h"
}

#include <Rcpp.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>


using namespace Rcpp;

int cqp_initialization_status = 0;

Attribute* make_s_attribute(SEXP corpus, SEXP s_attribute, SEXP registry){
  
  char* reg_dir = strdup(Rcpp::as<std::string>(registry).c_str());
  char* s_attr = strdup(Rcpp::as<std::string>(s_attribute).c_str());
  char* corpus_pointer  = strdup(Rcpp::as<std::string>(corpus).c_str());
  
  Corpus *corpus_obj = cl_new_corpus(reg_dir, corpus_pointer);
  Attribute* att = cl_new_attribute(corpus_obj, s_attr, ATT_STRUC);
  
  return att;
}


Attribute* make_p_attribute(SEXP corpus, SEXP p_attribute, SEXP registry){
  
  char* reg_dir = strdup(Rcpp::as<std::string>(registry).c_str());
  char* p_attr = strdup(Rcpp::as<std::string>(p_attribute).c_str());
  char* corpus_pointer  = strdup(Rcpp::as<std::string>(corpus).c_str());
  
  Corpus *corpus_obj = cl_new_corpus(reg_dir, corpus_pointer);
  Attribute* att = cl_new_attribute(corpus_obj, p_attr, ATT_POS);
  
  return att;
}

/* these are the wrappers for the functions of the corpus library (CL) */


// [[Rcpp::export(name=".cl_attribute_size")]]
int _cl_attribute_size(SEXP corpus, SEXP attribute, SEXP attribute_type, SEXP registry) {
  int size;
  std::string atype = Rcpp::as<std::string>(attribute_type);
  if (atype == "p"){
    Attribute* att = make_p_attribute(corpus, attribute, registry);
    size = cl_max_cpos(att);
  } else {
    Attribute* att = make_s_attribute(corpus, attribute, registry);
    size = cl_max_struc(att);
  }
  return(size);
}


// [[Rcpp::export(name=".cl_lexicon_size")]]
int _cl_lexicon_size(SEXP corpus, SEXP p_attribute, SEXP registry){
  Attribute* att = make_p_attribute(corpus, p_attribute, registry);
  int size = cl_max_id(att);
  return( size );
}


// [[Rcpp::export(name=".cl_cpos2struc")]]
Rcpp::IntegerVector _cl_cpos2struc(SEXP corpus, SEXP s_attribute, Rcpp::IntegerVector cpos, SEXP registry){
  Attribute* att = make_s_attribute(corpus, s_attribute, registry);
  int i;
  int len = cpos.length();
  Rcpp::IntegerVector strucs(len);
  for (i = 0; i < len; i++){
    strucs(i) = cl_cpos2struc(att, cpos(i));
  }
  return( strucs );
}


// [[Rcpp::export(name=".cl_cpos2str")]]
Rcpp::StringVector _cl_cpos2str(SEXP corpus, SEXP p_attribute, SEXP registry, Rcpp::IntegerVector cpos){
  Attribute* att = make_p_attribute(corpus, p_attribute, registry);
  int i;
  int len;
  len = cpos.length();
  Rcpp::StringVector result(len);
  for (i = 0; i < len; i++){
    result(i) = cl_cpos2str(att, cpos(i));
  }
  return(result);
}


// [[Rcpp::export(name=".cl_cpos2id")]]
Rcpp::IntegerVector _cl_cpos2id(SEXP corpus, SEXP p_attribute, SEXP registry, Rcpp::IntegerVector cpos){
  Attribute* att = make_p_attribute(corpus, p_attribute, registry);
  int i;
  int len = cpos.length();
  Rcpp::IntegerVector ids(len);
  for (i = 0; i < len; i++){
    ids(i) = cl_cpos2id(att, cpos(i));
  }
  return( ids );
}

// [[Rcpp::export(name=".cl_struc2cpos")]]
Rcpp::IntegerVector _cl_struc2cpos(SEXP corpus, SEXP s_attribute, SEXP registry, int struc){
  Attribute* att = make_s_attribute(corpus, s_attribute, registry);
  Rcpp::IntegerVector cpos(2);
  int lb, rb;
  cl_struc2cpos(att, struc, &lb, &rb);
  cpos(0) = lb;
  cpos(1) = rb;
  return( cpos );
}


// [[Rcpp::export(name=".cl_id2str")]]
Rcpp::StringVector _cl_id2str(SEXP corpus, SEXP p_attribute, SEXP registry, Rcpp::IntegerVector id){
  /* potentially cpos > max cpos causing a crash */
  int len = id.length();
  Attribute* att = make_p_attribute(corpus, p_attribute, registry);
  Rcpp::StringVector result(len);
  int i;
  for (i = 0; i < len; i++){
    result(i) = cl_id2str(att, id(i));
  }
  return ( result );
}


// [[Rcpp::export(name=".cl_struc2str")]]
Rcpp::StringVector _cl_struc2str(SEXP corpus, SEXP s_attribute, Rcpp::IntegerVector struc, SEXP registry){
  Attribute* att = make_s_attribute(corpus, s_attribute, registry);
  int len = struc.length();
  Rcpp::StringVector result(len);
  if ( cl_struc_values(att) ){
    int i;
    for (i = 0; i < len; i++){
      result(i) = cl_struc2str(att, struc(i));
    }
  }
  return ( result );
}


// [[Rcpp::export(name=".cl_regex2id")]]
Rcpp::IntegerVector _cl_regex2id(SEXP corpus, SEXP p_attribute, SEXP regex, SEXP registry){
  Attribute* att = make_p_attribute(corpus, p_attribute, registry);
  char *r = strdup(Rcpp::as<std::string>(regex).c_str());
  int *idlist;
  int len;
  int i;
  idlist = collect_matching_ids(att, r, 0, &len);
  Rcpp::IntegerVector result(len);
  for (i = 0; i < len; i++){
    result(i) = idlist[i];
  }
  return( result );
}


// [[Rcpp::export(name=".cl_str2id")]]
Rcpp::IntegerVector _cl_str2id(SEXP corpus, SEXP p_attribute, Rcpp::StringVector str, SEXP registry){
  Attribute* att = make_p_attribute(corpus, p_attribute, registry);
  int len = str.length();
  Rcpp::IntegerVector ids(len);
  int i;
  for (i = 0; i < len; i++){
    ids(i) = cl_str2id(att, str(i));
  }
  return( ids );
}

// [[Rcpp::export(name=".cl_id2freq")]]
Rcpp::IntegerVector _cl_id2freq(SEXP corpus, SEXP p_attribute, Rcpp::IntegerVector id, SEXP registry){
  Attribute* att = make_p_attribute(corpus, p_attribute, registry);
  int len = id.length();
  Rcpp::IntegerVector result(len);
  int i;
  for (i = 0; i < len; i++){
    result(i) = cl_id2freq(att, id(i));
  }
  return( result );
}


// [[Rcpp::export(name=".cl_id2cpos")]]
Rcpp::IntegerVector _cl_id2cpos(SEXP corpus, SEXP p_attribute, SEXP id, SEXP registry){
  Attribute* att = make_p_attribute(corpus, p_attribute, registry);
  int *cposlist;
  int len;
  int idx = Rcpp::as<int>(id);
  cposlist = cl_id2cpos(att, idx, &len);
  Rcpp::IntegerVector cpos(len);
  int i;
  for (i = 0; i < len; i++){
    cpos(i) = cposlist[i];
  }
  return( cpos );
}


// [[Rcpp::export(name=".cl_cpos2lbound")]]
int _cl_cpos2lbound(SEXP corpus, SEXP s_attribute, SEXP cpos, SEXP registry){
  Attribute* att = make_s_attribute(corpus, s_attribute, registry);
  int cpos_int = Rcpp::as<int>(cpos);
  int struc = cl_cpos2struc(att, cpos_int);
  int lb, rb;
  cl_struc2cpos(att, struc, &lb, &rb);
  return( lb );
}

// [[Rcpp::export(name=".cl_cpos2rbound")]]
int _cl_cpos2rbound(SEXP corpus, SEXP s_attribute, SEXP cpos, SEXP registry){
  Attribute* att = make_s_attribute(corpus, s_attribute, registry);
  int cpos_int = Rcpp::as<int>(cpos);
  int struc = cl_cpos2struc(att, cpos_int);
  int lb, rb;
  cl_struc2cpos(att, struc, &lb, &rb);
  return( rb );
}


// [[Rcpp::export(name=".cl_delete_corpus")]]
int _cl_delete_corpus(SEXP corpus, SEXP registry){
  
  Corpus * c;
  static char *canonical_name = NULL;
  
  char* registry_dir = strdup(Rcpp::as<std::string>(registry).c_str());
  char* registry_name  = strdup(Rcpp::as<std::string>(corpus).c_str());
  
  /* code copied from cl_new_corpus in corpus.c */
  cl_free(canonical_name);
  canonical_name = cl_strdup(registry_name);
  cl_id_tolower(canonical_name);
  if (!cl_id_validate(canonical_name)) {
    Rprintf("cl_new_corpus: <%s> is not a valid corpus name\n", registry_name);
  }
  
  Rprintf("Corpus to delete (ID): %s\n", registry_name);
  c = find_corpus(registry_dir, canonical_name); 
  Rprintf("Corpus name: %s\n", c->name);
  Rprintf("Number of loads before reset: %d\n", c->nr_of_loads);
  c->nr_of_loads = 1;
  Rprintf("Number of loads resetted: %d\n", c->nr_of_loads);
  cl_delete_corpus(c);
  
  return( 0 );
}



// [[Rcpp::export(name=".cl_charset_name")]]
Rcpp::StringVector _cl_charset_name(SEXP corpus, SEXP registry){
  
  char* corpus_pointer  = strdup(Rcpp::as<std::string>(corpus).c_str());
  char* reg_dir = strdup(Rcpp::as<std::string>(registry).c_str());
  Corpus *corpus_obj = cl_new_corpus(reg_dir, corpus_pointer);
  
  Rcpp::StringVector result(1);
  
  result(0) = cl_charset_name(cl_corpus_charset(corpus_obj));
  
  return( result );
}




// [[Rcpp::export(name=".init_cqp")]]
void init_cqp() {
	int		ac = 1;
	char *		av[1];
	av[0] = (char *)"RcppCWB";
	which_app = cqp;
	silent = 1; 
	paging = 0;
	autoshow = 0;
	auto_save = 0;
	server_log = 0;
	enable_macros = 0;

	initialize_cqp(ac, av);
	cqp_initialization_status = 1;
	make_attribute_hash(16384);
}


// [[Rcpp::export(name=".cqp_get_registry")]]
Rcpp::StringVector cqp_get_registry(){
  Rcpp::StringVector result(1);
  result(0) = cl_standard_registry();
  return result;
}

// [[Rcpp::export(name=".cqp_get_status")]]
int cqp_get_status(){
  return cqp_initialization_status;
}


// [[Rcpp::export(name=".cqp_set_registry")]]
SEXP cqp_set_registry(SEXP registry_dir){
  char * registry_new;
  registry_new = (char*)CHAR(STRING_ELT(registry_dir,0));
  registry = cl_strdup(registry_new);
  int		ac = 1;
  char *		av[1];
  av[0] = (char *)"RcppCWB";
  set_current_corpus(NULL, 0); /* required to avoid crash! */ 
  
  initialize_cqp(ac, av);
  make_attribute_hash(16384);
  
  SEXP result = R_NilValue;
  return result;
}


// [[Rcpp::export(name=".cqp_list_corpora")]]
Rcpp::StringVector cqp_list_corpora(){
  
  CorpusList *	cl;
  int	i = 0, n = 0;
  
  /* First count corpora */
  for (cl = FirstCorpusFromList(); cl != NULL; cl = NextCorpusFromList(cl)) {
    if (cl->type == SYSTEM) n++;
  }
  Rcpp::StringVector result(n);

  /* Then build list of names */
  for (cl = FirstCorpusFromList(); cl != NULL; cl = NextCorpusFromList(cl)) {
    if (cl->type == SYSTEM) {
      result(i) = cl->name;
      i++;
    }
  }
  return result;
  
}


// [[Rcpp::export(name=".cqp_query")]]
SEXP cqp_query(SEXP corpus, SEXP subcorpus, SEXP query){
  
  char * mother = (char*)CHAR(STRING_ELT(corpus,0));
  char * child = (char*)CHAR(STRING_ELT(subcorpus,0));
  char * q = (char*)CHAR(STRING_ELT(query,0));
  char * cqp_query;
  CorpusList *cl;
  
  /* is this necessary */
  char	*c, *sc;
  if (!split_subcorpus_spec(mother, &c, &sc)) {
    Rprintf("ERROR (function: split_subcorpus_spec)");
  }
  
  cl = cqi_find_corpus(mother);
  set_current_corpus(cl, 0);

  int len = strlen(child) + strlen(q) + 10;
  cqp_query = (char *) cl_malloc(len);
  
  sprintf(cqp_query, "%s = %s", child, q);

  if (!cqi_activate_corpus(mother)){
    Rprintf("activation failed");
  }
  if (!check_subcorpus_name(child)){
    Rprintf("checking subcorpus name failed \n");
  }
  
  if (!cqp_parse_string(cqp_query)){
    Rprintf("ERROR: Cannot parse the CQP query.\n");
  } else {
    char *			full_child;
    CorpusList *	childcl;
    
    full_child = combine_subcorpus_spec(mother, child); /* c is the 'physical' part of the mother corpus */

    childcl = cqi_find_corpus(full_child);
    if ((childcl) == NULL) {
      Rprintf("subcorpus not found\n");
    } 
  }

  
  SEXP result = R_NilValue;
  return result;
}


// [[Rcpp::export(name=".cqp_subcorpus_size")]]
int cqp_subcorpus_size(SEXP scorpus)
{
  int result;
  char * subcorpus;
  CorpusList * cl;
  
  subcorpus = (char*)CHAR(STRING_ELT(scorpus,0));
  cl = cqi_find_corpus(subcorpus);
  
  if (cl == NULL) {
    result = 0;
  } else {
    result = cl->size;
  }
  return result;
}

// [[Rcpp::export(name=".cqp_list_subcorpora")]]
Rcpp::StringVector cqp_list_subcorpora(SEXP inCorpus)
{
  char * corpus;
  CorpusList *cl, *mother;
  int i = 0, n = 0;
  Rcpp::StringVector result;

  corpus = (char*)CHAR(STRING_ELT(inCorpus,0));
  
  mother = cqi_find_corpus(corpus);
  if (!check_corpus_name(corpus) || mother == NULL) {
    /* Rcpp::StringVector result(1); */
  } else {
    /* First count subcorpora */
    for (cl = FirstCorpusFromList(); cl != NULL; cl = NextCorpusFromList(cl)) {
      if (cl->type == SUB && cl->corpus == mother->corpus){
        n++;
      }
    }
    Rprintf("number of subcorpora: %d\n", n);
    Rcpp::StringVector result(n);
    
    /* Then build list of names */
    for (cl = FirstCorpusFromList(); cl != NULL; cl = NextCorpusFromList(cl)) {
      if (cl->type == SUB && cl->corpus == mother->corpus) {
        result(i) = cl->name;
        Rprintf("subcorpus name: %s\n", cl->name);
        /* printf("added to result: %s", result(i)); */
        i++;
      }
    }
  }
  return result;
}

// [[Rcpp::export(name=".cqp_dump_subcorpus")]]
Rcpp::IntegerMatrix cqp_dump_subcorpus(SEXP inSubcorpus)
{
  char * subcorpus;
  CorpusList * cl;
  int i;
  int nrows = cqp_subcorpus_size(inSubcorpus);

  subcorpus = (char*)CHAR(STRING_ELT(inSubcorpus,0));
  cl = cqi_find_corpus(subcorpus);
  if (cl == NULL) {
    Rprintf("subcorpus not found\n");
  }
  
  Rcpp::IntegerMatrix result(nrows,2);

  for (i = 0; i < nrows; i++) {
    result(i,0) = cl->range[i].start;
    result(i,1) = cl->range[i].end;
  }
  
  return result;
}

// [[Rcpp::export(name=".cqp_drop_subcorpus")]]
SEXP cqp_drop_subcorpus(SEXP inSubcorpus)
{
  SEXP			result = R_NilValue;
  char *			subcorpus;
  char 			*c, *sc;
  CorpusList *	cl;
  
  PROTECT(inSubcorpus);
  
  subcorpus = (char*)CHAR(STRING_ELT(inSubcorpus,0));
  
  /* Make sure it is a subcorpus, not a root corpus */
  if (!split_subcorpus_spec(subcorpus, &c, &sc)) {
    UNPROTECT(1);
    /* rcqp_error_code(cqi_errno); */
  } else if (sc == NULL) {
    free(c);
    UNPROTECT(1);
    /* error("can't drop a root corpus."); */
  } else {
    free(c); free(sc);
    cl = cqi_find_corpus(subcorpus);
    if (cl == NULL) {
      UNPROTECT(1);
      /* rcqp_error_code(cqi_errno); */
    } else {
      dropcorpus(cl);
    }
  }
  
  UNPROTECT(1);
  
  return result;
}



// [[Rcpp::export(name=".check_corpus")]]
int check_corpus(SEXP corpus){
  
  char * c;
  CorpusList * cl;
  
  c = (char*)CHAR(STRING_ELT(corpus,0));
  cl = findcorpus(c, SYSTEM, 0);
  
  if (cl == NULL || !access_corpus(cl)) {
    return 0;
  } else {
    return 1;
  }
}


