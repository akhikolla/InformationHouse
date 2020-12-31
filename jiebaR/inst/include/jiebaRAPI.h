#ifndef __JIEBARAPI_H__
#define __JIEBARAPI_H__

// #include <R.h>
// #include <Rinternals.h>

#include <stddef.h>
#include <R_ext/Rdynload.h>

#ifdef HAVE_VISIBILITY_ATTRIBUTE
    # define attribute_hidden __attribute__ ((visibility ("hidden")))
#else
    # define attribute_hidden
#endif

#ifdef __cplusplus
extern "C" {
#endif

/**
 *  **********  all text input should be UTF-8 encoding string**************
 *  For windows platform, the string output of C++ to be printed by R should have Encoding() setted first
 *
 */

  /**
   *   removed in jiebaR v0.8 cppjieba v4.4.1
   */
SEXP attribute_hidden jiebaR_mp_ptr(SEXP dictSEXP, SEXP userSEXP,SEXP stopSEXP){
    static SEXP(*f)(SEXP,SEXP,SEXP) = NULL;
    if (!f) {
        f = (SEXP(*)(SEXP,SEXP,SEXP)) R_GetCCallable("jiebaR", "jiebaR_mp_ptr");
    }
    return f(dictSEXP,userSEXP,stopSEXP);
}
  /**
   *   removed in jiebaR v0.8 cppjieba v4.4.1
   */
SEXP attribute_hidden jiebaR_mp_cut(SEXP xSEXP, SEXP cutterSEXP){
    static SEXP(*f)(SEXP,SEXP) = NULL;
    if (!f) {
        f = (SEXP(*)(SEXP,SEXP)) R_GetCCallable("jiebaR", "jiebaR_mp_cut");
    }
    return f(xSEXP,cutterSEXP);
}
  /**
   *   removed in jiebaR v0.8 cppjieba v4.4.1
   */
SEXP attribute_hidden jiebaR_mix_ptr(SEXP dictSEXP, SEXP modelSEXP, SEXP userSEXP, SEXP stopSEXP){
    static SEXP(*f)(SEXP,SEXP,SEXP,SEXP) = NULL;
    if (!f) {
        f = (SEXP(*)(SEXP,SEXP,SEXP,SEXP)) R_GetCCallable("jiebaR", "jiebaR_mix_ptr");
    }
    return f(dictSEXP, modelSEXP, userSEXP, stopSEXP);
}
  /**
   *   removed in jiebaR v0.8 cppjieba v4.4.1
   */
SEXP attribute_hidden jiebaR_mix_cut(SEXP xSEXP, SEXP cutterSEXP){
    static SEXP(*f)(SEXP,SEXP) = NULL;
    if (!f) {
        f = (SEXP(*)(SEXP,SEXP)) R_GetCCallable("jiebaR", "jiebaR_mix_cut");
    }
    return f(xSEXP,cutterSEXP);
}
  /**
   *   removed in jiebaR v0.8 cppjieba v4.4.1
   */
SEXP attribute_hidden jiebaR_query_ptr(SEXP dictSEXP, SEXP modelSEXP, SEXP nSEXP, SEXP stopSEXP){
    static SEXP(*f)(SEXP,SEXP,SEXP,SEXP) = NULL;
    if (!f) {
        f = (SEXP(*)(SEXP,SEXP,SEXP,SEXP)) R_GetCCallable("jiebaR", "jiebaR_query_ptr");
    }
    return f(dictSEXP, modelSEXP, nSEXP, stopSEXP);
}
  /**
   *   removed in jiebaR v0.8 cppjieba v4.4.1
   */
SEXP attribute_hidden jiebaR_query_cut(SEXP xSEXP, SEXP cutterSEXP){
    static SEXP(*f)(SEXP,SEXP) = NULL;
    if (!f) {
        f = (SEXP(*)(SEXP,SEXP)) R_GetCCallable("jiebaR", "jiebaR_query_cut");
    }
    return f(xSEXP,cutterSEXP);
}
  /**
   *   removed in jiebaR v0.8 cppjieba v4.4.1
   */
SEXP attribute_hidden jiebaR_hmm_ptr(SEXP modelSEXP, SEXP stopSEXP){
    static SEXP(*f)(SEXP,SEXP) = NULL;
    if (!f) {
        f = (SEXP(*)(SEXP,SEXP)) R_GetCCallable("jiebaR", "jiebaR_hmm_ptr");
    }
    return f(modelSEXP, stopSEXP);
}
  /**
   *   removed in jiebaR v0.8 cppjieba v4.4.1
   */
SEXP attribute_hidden jiebaR_hmm_cut(SEXP xSEXP, SEXP cutterSEXP){
    static SEXP(*f)(SEXP,SEXP) = NULL;
    if (!f) {
        f = (SEXP(*)(SEXP,SEXP)) R_GetCCallable("jiebaR", "jiebaR_hmm_cut");
    }
    return f(xSEXP,cutterSEXP);
}
  /**
   *   removed in jiebaR v0.8 cppjieba v4.4.1
   */
SEXP attribute_hidden jiebaR_tag_ptr(SEXP dictSEXP, SEXP modelSEXP, SEXP userSEXP, SEXP stopSEXP){
    static SEXP(*f)(SEXP,SEXP,SEXP,SEXP) = NULL;
    if (!f) {
        f = (SEXP(*)(SEXP,SEXP,SEXP,SEXP)) R_GetCCallable("jiebaR", "jiebaR_tag_ptr");
    }
    return f(dictSEXP, modelSEXP, userSEXP, stopSEXP);
}
  /**
   *   removed in jiebaR v0.8 cppjieba v4.4.1
   */
SEXP attribute_hidden jiebaR_tag_tag(SEXP xSEXP, SEXP cutterSEXP){
    static SEXP(*f)(SEXP,SEXP) = NULL;
    if (!f) {
        f = (SEXP(*)(SEXP,SEXP)) R_GetCCallable("jiebaR", "jiebaR_tag_tag");
    }
    return f(xSEXP,cutterSEXP);
}
  /**
   *   removed in jiebaR v0.8 cppjieba v4.4.1
   */
SEXP attribute_hidden jiebaR_tag_file(SEXP xSEXP, SEXP cutterSEXP){
    static SEXP(*f)(SEXP,SEXP) = NULL;
    if (!f) {
        f = (SEXP(*)(SEXP,SEXP)) R_GetCCallable("jiebaR", "jiebaR_tag_file");
    }
    return f(xSEXP,cutterSEXP);
}

/**
 * [jiebaR_key_ptr description]
 * @param  nSEXP     numbers of keywords
 * @param  dictSEXP  a file path
 * @param  modelSEXP a file path
 * @param  idfSEXP   a file path
 * @param  stopSEXP  a file path
 * @param  userSEXP  a file path
 * @return           a key ptr
 */
SEXP attribute_hidden jiebaR_key_ptr(SEXP nSEXP, SEXP dictSEXP, SEXP modelSEXP, SEXP idfSEXP, SEXP stopSEXP, SEXP userSEXP){
    static SEXP(*f)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP) = NULL;
    if (!f) {
        f = (SEXP(*)(SEXP,SEXP,SEXP,SEXP,SEXP,SEXP)) R_GetCCallable("jiebaR", "jiebaR_key_ptr");
    }
    return f(nSEXP, dictSEXP, modelSEXP, idfSEXP, stopSEXP,userSEXP);
}

/**
 * [jiebaR_key_tag description]
 * @param  xSEXP      a string in R, only the first element will be used
 * @param  cutterSEXP a key ptr, warning no checking for ptr type in runtime
 * @return            a vector with attributte
 */
SEXP attribute_hidden jiebaR_key_tag(SEXP xSEXP, SEXP cutterSEXP){
    static SEXP(*f)(SEXP,SEXP) = NULL;
    if (!f) {
        f = (SEXP(*)(SEXP,SEXP)) R_GetCCallable("jiebaR", "jiebaR_key_tag");
    }
    return f(xSEXP,cutterSEXP);
}

/**
 * [jiebaR_key_cut description]
 * @param  xSEXP      a string in R, only the first element will be used
 * @param  cutterSEXP a key ptr, warning no checking for ptr type in runtime
 * @return            a vector only with keywords
 */
SEXP attribute_hidden jiebaR_key_cut(SEXP xSEXP, SEXP cutterSEXP){
    static SEXP(*f)(SEXP,SEXP) = NULL;
    if (!f) {
        f = (SEXP(*)(SEXP,SEXP)) R_GetCCallable("jiebaR", "jiebaR_key_cut");
    }
    return f(xSEXP,cutterSEXP);
}

/**
 * [jiebaR_key_keys description]
 * @param  xSEXP      a string in R, only the first element will be used
 * @param  cutterSEXP a key ptr, warning no checking for ptr type in runtime
 * @return            a vector only with keywords
 */
SEXP attribute_hidden jiebaR_key_keys(SEXP xSEXP, SEXP cutterSEXP){
    static SEXP(*f)(SEXP,SEXP) = NULL;
    if (!f) {
        f = (SEXP(*)(SEXP,SEXP)) R_GetCCallable("jiebaR", "jiebaR_key_keys");
    }
    return f(xSEXP,cutterSEXP);
}

/**
 * [jiebaR_sim_ptr description]
 * @param  dictSEXP  a file path
 * @param  modelSEXP a file path
 * @param  idfSEXP   a file path
 * @param  stopSEXP  a file path
 * @param  userSEXP  a file path
 * @return           a sim ptr
 */
SEXP attribute_hidden jiebaR_sim_ptr(SEXP dictSEXP, SEXP modelSEXP, SEXP idfSEXP, SEXP stopSEXP,SEXP userSEXP){
    static SEXP(*f)(SEXP,SEXP,SEXP,SEXP,SEXP) = NULL;
    if (!f) {
        f = (SEXP(*)(SEXP,SEXP,SEXP,SEXP,SEXP)) R_GetCCallable("jiebaR", "jiebaR_sim_ptr");
    }
    return f(dictSEXP, modelSEXP, idfSEXP, stopSEXP,userSEXP);
}

/**
 * [jiebaR_sim_sim description]
 * @param  codeSEXP   a string in R, only the first element will be used
 * @param  topnSEXP   numbers of keywords
 * @param  cutterSEXP a sim ptr
 * @return            a R list of result
 */
SEXP attribute_hidden jiebaR_sim_sim(SEXP codeSEXP, SEXP topnSEXP, SEXP cutterSEXP){
    static SEXP(*f)(SEXP,SEXP,SEXP) = NULL;
    if (!f) {
        f = (SEXP(*)(SEXP,SEXP,SEXP)) R_GetCCallable("jiebaR", "jiebaR_sim_sim");
    }
    return f(codeSEXP, topnSEXP, cutterSEXP);
}

/**
 * [jiebaR_sim_distance description]
 * @param  lhsSEXP    a string in R, only the first element will be used
 * @param  rhsSEXP    a string in R, only the first element will be used
 * @param  topnSEXP   numbers of keywords
 * @param  cutterSEXP a sim ptr
 * @return            a R list of result
 */
SEXP attribute_hidden jiebaR_sim_distance(SEXP lhsSEXP, SEXP rhsSEXP, SEXP topnSEXP, SEXP cutterSEXP){
    static SEXP(*f)(SEXP,SEXP,SEXP,SEXP) = NULL;
    if (!f) {
        f = (SEXP(*)(SEXP,SEXP,SEXP,SEXP)) R_GetCCallable("jiebaR", "jiebaR_sim_distance");
    }
    return f(lhsSEXP, rhsSEXP, topnSEXP, cutterSEXP);
}


// v4 new

/**
  * [jiebaR_mix_ptr description]
  * @param  dictSEXP  a file path
  * @param  modelSEXP a file path
  * @param  userSEXP  a file path
  * @param  stopSEXP  a file path
  * @return           a jieba ptr
  */
SEXP attribute_hidden jiebaR_jiebaclass_ptr(SEXP dict, SEXP model, SEXP user,SEXP stop){
    static SEXP(*f)(SEXP,SEXP,SEXP,SEXP) = NULL;
    if (!f) {
      f = (SEXP(*)(SEXP,SEXP,SEXP,SEXP)) R_GetCCallable("jiebaR", "jiebaR_jiebaclass_ptr");
    }
    return f(dict, model, user,stop);
}


  SEXP attribute_hidden jiebaR_jiebaclass_ptr_v2(SEXP dict, SEXP model, SEXP user,SEXP stop,SEXP uw){
    static SEXP(*f)(SEXP,SEXP,SEXP,SEXP,SEXP) = NULL;
    if (!f) {
      f = (SEXP(*)(SEXP,SEXP,SEXP,SEXP,SEXP)) R_GetCCallable("jiebaR", "jiebaR_jiebaclass_ptr_v2");
    }
    return f(dict, model, user,stop, uw);
  }
  /**
   * [description]
   * @param  xSEXP      a string in R, only the first element will be used
   * @param  cutterSEXP a jieba ptr, warning: no checking for ptr type in runtime
   * @return            a string vector
   */
SEXP attribute_hidden jiebaR_jiebaclass_mix_cut(SEXP x, SEXP cutter){
  static SEXP(*f)(SEXP,SEXP) = NULL;
  if (!f) {
    f = (SEXP(*)(SEXP,SEXP)) R_GetCCallable("jiebaR", "jiebaR_jiebaclass_mix_cut");
  }
  return f(x, cutter);
}
  /**
   * [description]
   * @param  xSEXP      a string in R, only the first element will be used
   * @param  cutterSEXP a jieba ptr, warning: no checking for ptr type in runtime
   * @return            a string vector
   */
  SEXP attribute_hidden jiebaR_jiebaclass_mp_cut(SEXP x, SEXP cutter){
    static SEXP(*f)(SEXP,SEXP) = NULL;
    if (!f) {
      f = (SEXP(*)(SEXP,SEXP)) R_GetCCallable("jiebaR", "jiebaR_jiebaclass_mp_cut");
    }
    return f(x, cutter);
  }
  /**
   * [description]
   * @param  xSEXP      a string in R, only the first element will be used
   * @param  cutterSEXP a jieba ptr, warning: no checking for ptr type in runtime
   * @return            a string vector
   */
  SEXP attribute_hidden jiebaR_jiebaclass_hmm_cut(SEXP x, SEXP cutter){
    static SEXP(*f)(SEXP,SEXP) = NULL;
    if (!f) {
      f = (SEXP(*)(SEXP,SEXP)) R_GetCCallable("jiebaR", "jiebaR_jiebaclass_hmm_cut");
    }
    return f(x, cutter);
  }
  /**
   * [description]
   * @param  xSEXP      a string in R, only the first element will be used
   * @param  cutterSEXP a jieba ptr, warning: no checking for ptr type in runtime
   * @return            a string vector
   */
  SEXP attribute_hidden jiebaR_jiebaclass_query_cut(SEXP x, SEXP cutter){
    static SEXP(*f)(SEXP,SEXP) = NULL;
    if (!f) {
      f = (SEXP(*)(SEXP,SEXP)) R_GetCCallable("jiebaR", "jiebaR_jiebaclass_query_cut");
    }
    return f(x, cutter);
  }
  /**
   * [description]
   * @param  xSEXP      a string in R, only the first element will be used
   * @param  cutterSEXP a jieba ptr, warning: no checking for ptr type in runtime
   * @return            a string vector
   */
  SEXP attribute_hidden jiebaR_jiebaclass_full_cut(SEXP x, SEXP cutter){
    static SEXP(*f)(SEXP,SEXP) = NULL;
    if (!f) {
      f = (SEXP(*)(SEXP,SEXP)) R_GetCCallable("jiebaR", "jiebaR_jiebaclass_full_cut");
    }
    return f(x, cutter);
  }

  /**
   * [description]
   * @param  xSEXP      a string in R, only the first element will be used
   * @param  cutterSEXP a jieba ptr, warning: no checking for ptr type in runtime
   * @return            a string vector with names
   */
  SEXP attribute_hidden jiebaR_jiebaclass_tag_tag(SEXP x, SEXP cutter){
    static SEXP(*f)(SEXP,SEXP) = NULL;
    if (!f) {
      f = (SEXP(*)(SEXP,SEXP)) R_GetCCallable("jiebaR", "jiebaR_jiebaclass_tag_tag");
    }
    return f(x, cutter);
  }
  /**
   * [description]
   * @param  xSEXP      a string in R, only the first element will be used
   * @param  cutterSEXP a jieba ptr, warning: no checking for ptr type in runtime
   * @return            a string vector
   */
  SEXP attribute_hidden jiebaR_jiebaclass_tag_file(SEXP x, SEXP cutter){
    static SEXP(*f)(SEXP,SEXP) = NULL;
    if (!f) {
      f = (SEXP(*)(SEXP,SEXP)) R_GetCCallable("jiebaR", "jiebaR_jiebaclass_tag_file");
    }
    return f(x, cutter);
  }

  /**
   * [description]
   * @param  x      a string
   * @param tag     a string
   * @param  cutterSEXP a jieba ptr, warning: no checking for ptr type in runtime
   * @return        bool
   */
  SEXP attribute_hidden jiebaR_add_user_word(SEXP x,SEXP tag, SEXP cutter){
    static SEXP(*f)(SEXP,SEXP,SEXP) = NULL;
    if (!f) {
      f = (SEXP(*)(SEXP,SEXP,SEXP)) R_GetCCallable("jiebaR", "jiebaR_add_user_word");
    }
    return f(x, tag, cutter);
  }
  /**
   * [description]
   * @param  x      a string
   * @return        a string
   */
  SEXP attribute_hidden jiebaR_u64tobin(SEXP x){
    static SEXP(*f)(SEXP) = NULL;
    if (!f) {
      f = (SEXP(*)(SEXP)) R_GetCCallable("jiebaR", "jiebaR_u64tobin");
    }
    return f(x);
  }



  /**
   * [jiebaR_filecoding description]
   * @param  fileSEXP a file path
   * @return          a string about encoding
   */
  SEXP attribute_hidden jiebaR_filecoding(SEXP fileSEXP){
    static SEXP(*f)(SEXP) = NULL;
    if (!f) {
      f = (SEXP(*)(SEXP)) R_GetCCallable("jiebaR", "jiebaR_filecoding");
    }
    return f(fileSEXP);
  }



#ifdef __cplusplus
}
#endif

#endif /* __JIEBARAPI_H__ */
