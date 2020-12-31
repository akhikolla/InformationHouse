#' Get Attribute Size (of Positional/Structural Attribute).
#' 
#' Use \code{cl_attribute_size} to get the total number of values of a
#' positional attribute (param \code{attribute_type} = "p"), or structural
#' attribute (param \code{attribute_type} = "s"). Note that indices are
#' zero-based, i.e. the maximum position of a positional / structural
#' attribute is attribute size minus 1 (see examples).
#' @rdname cl_attribute_size
#' @param corpus name of a CWB corpus (upper case)
#' @param attribute name of a p- or s-attribute
#' @param attribute_type either "p" or "s", for structural/positional attribute
#' @param registry path to the registry directory, defaults to the value of the
#'   environment variable CORPUS_REGISTRY
#' @examples 
#' registry <- if (!check_pkg_registry_files()) use_tmp_registry() else get_pkg_registry()
#' 
#' Sys.setenv(CORPUS_REGISTRY = registry)
#' token_no <- cl_attribute_size("REUTERS", attribute = "word", attribute_type = "p")
#' corpus_positions <- seq.int(from = 0, to = token_no - 1)
#' cl_cpos2id("REUTERS", "word", cpos = corpus_positions)
#' 
#' places_no <- cl_attribute_size("REUTERS", attribute = "places", attribute_type = "s")
#' strucs <- seq.int(from = 0, to = places_no - 1)
#' cl_struc2str("REUTERS", "places", struc = strucs)
cl_attribute_size <- function(corpus, attribute, attribute_type, registry = Sys.getenv("CORPUS_REGISTRY")){
  .cl_attribute_size(corpus = corpus, attribute = attribute, attribute_type = attribute_type, registry = registry)
}

#' Get Lexicon Size.
#' 
#' Get the total number of unique tokens/ids of a positional attribute. Note
#' that token ids are zero-based, i.e. when iterating through tokens, start at
#' 0, the maximum will be \code{cl_lexicon_size()} minus 1.
#' 
#' @param corpus name of a CWB corpus (upper case)
#' @param p_attribute name of positional attribute
#' @param registry path to the registry directory, defaults to the value of the
#'   environment variable CORPUS_REGISTRY
#' @rdname cl_lexicon_size
#' @examples 
#' registry <- if (!check_pkg_registry_files()) use_tmp_registry() else get_pkg_registry()
#' Sys.setenv(CORPUS_REGISTRY = registry)
#' lexicon_size <- cl_lexicon_size("REUTERS", p_attribute = "word")
#' token_ids <- seq.int(from = 0, to = lexicon_size - 1)
#' cl_id2str("REUTERS", p_attribute = "word", id = token_ids)
cl_lexicon_size <- function(corpus, p_attribute, registry = Sys.getenv("CORPUS_REGISTRY")){
  .cl_lexicon_size(corpus = corpus, p_attribute = p_attribute, registry = registry)
}


#' @title Using Structural Attributes.
#' 
#' @description Structural attributes store the metadata of texts in a CWB
#'   corpus and/or any kind of annotation of a region of text. The fundamental
#'   unit are so-called strucs, i.e. indices of regions identified by a left and
#'   a right corpus position. The corpus library (CL) offers a set of functions
#'   to make the translations between corpus positions (cpos) and strucs
#'   (struc).
#' 
#' @param corpus name of a CWB corpus (upper case)
#' @param s_attribute name of structural attribute (character vector)
#' @param cpos corpus positions (integer vector)
#' @param struc a struc identifying a region
#' @param registry path to the registry directory, defaults to the value of the
#'   environment variable CORPUS_REGISTRY
#' @rdname s_attributes
#' @name CL: s_attributes
#' @examples
#' registry <- if (!check_pkg_registry_files()) use_tmp_registry() else get_pkg_registry()
#' 
#' # get metadata for matches of token
#' # scenario: id of the texts with occurrence of 'oil'
#' token_to_get <- "oil"
#' token_id <- cl_str2id("REUTERS", p_attribute = "word", str = "oil")
#' token_cpos <- cl_id2cpos("REUTERS", p_attribute = "word", id = token_id)
#' strucs <- cl_cpos2struc("REUTERS", s_attribute = "id", cpos = token_cpos)
#' strucs_unique <- unique(strucs)
#' text_ids <- cl_struc2str("REUTERS", s_attribute = "id", struc = strucs_unique)
#' 
#' # get the full text of the first text with match for 'oil'
#' left_cpos <- cl_cpos2lbound("REUTERS", s_attribute = "id", cpos = min(token_cpos))
#' right_cpos <- cl_cpos2rbound("REUTERS", s_attribute = "id", cpos = min(token_cpos))
#' txt <- cl_cpos2str("REUTERS", p_attribute = "word", cpos = left_cpos:right_cpos)
#' fulltext <- paste(txt, collapse = " ")
#' 
#' # alternativ approach to achieve same result
#' first_struc_match_oil <- cl_cpos2struc("REUTERS", s_attribute = "id", cpos = min(token_cpos))
#' cpos_struc <- cl_struc2cpos("REUTERS", s_attribute = "id", struc = first_struc_match_oil)
#' txt <- cl_cpos2str("REUTERS", p_attribute = "word", cpos = cpos_struc[1]:cpos_struc[2])
#' fulltext <- paste(txt, collapse = " ")
cl_cpos2struc <- function(corpus, s_attribute, cpos, registry = Sys.getenv("CORPUS_REGISTRY")){
  check_registry(registry)
  check_corpus(corpus, registry)
  check_s_attribute(corpus = corpus, registry = registry, s_attribute = s_attribute)
  check_cpos(corpus = corpus, p_attribute = "word", cpos = cpos, registry = registry)
  .cl_cpos2struc(corpus = corpus, s_attribute = s_attribute, cpos = cpos, registry = registry)
}

#' @rdname s_attributes
cl_struc2cpos <- function(corpus, s_attribute, registry = Sys.getenv("CORPUS_REGISTRY"), struc){
  check_registry(registry)
  check_corpus(corpus, registry)
  check_s_attribute(corpus = corpus, registry = registry, s_attribute = s_attribute)
  check_strucs(corpus = corpus, s_attribute = s_attribute, strucs = struc, registry = registry)
  .cl_struc2cpos(corpus = corpus, s_attribute = s_attribute, registry = registry, struc = struc)
}

#' @rdname s_attributes
cl_struc2str <- function(corpus, s_attribute, struc, registry = Sys.getenv("CORPUS_REGISTRY")){
  check_registry(registry)
  check_corpus(corpus, registry)
  check_s_attribute(corpus = corpus, registry = registry, s_attribute = s_attribute)
  check_strucs(corpus = corpus, s_attribute = s_attribute, strucs = struc, registry = registry)
  .cl_struc2str(corpus = corpus, s_attribute = s_attribute, struc = struc, registry = registry)
}

#' @rdname s_attributes
cl_cpos2lbound <- function(corpus, s_attribute, cpos, registry = Sys.getenv("CORPUS_REGISTRY")){
  check_registry(registry)
  check_corpus(corpus, registry)
  check_s_attribute(corpus = corpus, registry = registry, s_attribute = s_attribute)
  check_cpos(corpus = corpus, p_attribute = "word", cpos = cpos, registry = registry)
  .cl_cpos2lbound(corpus = corpus, s_attribute = s_attribute, cpos = cpos, registry = registry)
}

#' @rdname s_attributes
cl_cpos2rbound <- function(corpus, s_attribute, cpos, registry = Sys.getenv("CORPUS_REGISTRY")){
  check_registry(registry)
  check_corpus(corpus, registry)
  check_s_attribute(corpus = corpus, registry = registry, s_attribute = s_attribute)
  check_cpos(corpus = corpus, p_attribute = "word", cpos = cpos, registry = registry)
  .cl_cpos2rbound(corpus = corpus, s_attribute = s_attribute, cpos = cpos, registry = registry)
}






#' @title Using Positional Attributes.
#' 
#' @description CWB indexed corpora store the text of a corpus as numbers: Every token
#' in the token stream of the corpus is identified by a unique corpus
#' position. The string value of every token is identified by a unique integer
#' id. The corpus library (CL) offers a set of functions to make the transitions
#' between corpus positions, token ids, and the character string of tokens.
#' 
#' @param corpus name of a CWB corpus (upper case)
#' @param registry path to the registry directory, defaults to the value of the
#'   environment variable CORPUS_REGISTRY
#' @param p_attribute a p-attribute (positional attribute)
#' @param cpos corpus positions (integer vector)
#' @param id id of a token
#' @param regex a regular expression
#' @param str a character string
#' @rdname p_attributes
#' @name CL: p_attributes
#' @examples 
#' # registry directory and cpos_total will be needed in examples
#' registry <- if (!check_pkg_registry_files()) use_tmp_registry() else get_pkg_registry()
#' Sys.setenv(CORPUS_REGISTRY = registry)
#' cpos_total <- cl_attribute_size(
#'   corpus = "REUTERS", attribute = "word",
#'   attribute_type = "p", registry = registry
#'   )
#' 
#' # decode the token stream of the corpus (the quick way)
#' token_stream_str <- cl_cpos2str(
#'   corpus = "REUTERS", p_attribute = "word",
#'   cpos = seq.int(from = 0, to = cpos_total - 1),
#'   registry = registry
#'   )
#'   
#' # decode the token stream (cpos2id first, then id2str)
#' token_stream_ids <- cl_cpos2id(
#'   corpus = "REUTERS", p_attribute = "word",
#'   cpos = seq.int(from = 0, to = cpos_total - 1),
#'   registry = registry
#'   )
#' token_stream_str <- cl_id2str(
#'   corpus = "REUTERS", p_attribute = "word",
#'   id = token_stream_ids, registry = registry
#' )
#' 
#' # get corpus positions of a token
#' token_to_get <- "oil"
#' id_oil <- cl_str2id(
#'   corpus = "REUTERS", p_attribute = "word",
#'   str = token_to_get
#'   )
#' cpos_oil <- cl_id2cpos <- cl_id2cpos(
#'   corpus = "REUTERS", p_attribute = "word",
#'   id = id_oil
#' )
#' 
#' # get frequency of token
#' oil_freq <- cl_id2freq(
#'   corpus = "REUTERS", p_attribute = "word", id = id_oil
#' )
#' length(cpos_oil) # needs to be the same as oil_freq
#' 
#' # use regular expressions 
#' ids <- cl_regex2id(
#'   corpus = "REUTERS", p_attribute = "word",
#'   regex = "M.*"
#' )
#' m_words <- cl_id2str(
#'   corpus = "REUTERS", p_attribute = "word",
#'   id = ids
#' )
#' 
cl_cpos2str <- function(corpus, p_attribute, registry = Sys.getenv("CORPUS_REGISTRY"), cpos){
  check_registry(registry)
  check_corpus(corpus, registry)
  .cl_cpos2str(corpus = corpus, p_attribute = p_attribute, registry = registry, cpos = cpos)
}

#' @rdname p_attributes
cl_cpos2id <- function(corpus, p_attribute, registry = Sys.getenv("CORPUS_REGISTRY"), cpos){
  check_registry(registry)
  check_corpus(corpus, registry)
  .cl_cpos2id(corpus = corpus, p_attribute = p_attribute, registry = registry, cpos = cpos)
}

#' @rdname p_attributes
cl_id2str <- function(corpus, p_attribute, registry = Sys.getenv("CORPUS_REGISTRY"), id){
  check_registry(registry)
  check_corpus(corpus, registry)
  check_id(corpus = corpus, p_attribute = p_attribute, id = id, registry = registry)
  .cl_id2str(corpus = corpus, p_attribute = p_attribute, registry = registry, id = id)
}

#' @rdname p_attributes
cl_regex2id <- function(corpus, p_attribute, regex, registry = Sys.getenv("CORPUS_REGISTRY")){
  check_registry(registry)
  check_corpus(corpus, registry)
  .cl_regex2id(corpus = corpus, p_attribute = p_attribute, regex = regex, registry = registry)
}

#' @rdname p_attributes
cl_str2id <- function(corpus, p_attribute, str, registry = Sys.getenv("CORPUS_REGISTRY")){
  check_registry(registry)
  check_corpus(corpus, registry)
  .cl_str2id(corpus = corpus, p_attribute = p_attribute, str = str, registry = registry)
}

#' @rdname p_attributes
cl_id2freq <- function(corpus, p_attribute, id, registry = Sys.getenv("CORPUS_REGISTRY")){
  check_registry(registry)
  check_corpus(corpus, registry)
  check_p_attribute(p_attribute = p_attribute, corpus = corpus, registry = registry)
  check_id(corpus = corpus, p_attribute = p_attribute, id = id, registry = registry)
  .cl_id2freq(corpus = corpus, p_attribute = p_attribute, id = id, registry = registry)
}


#' @rdname p_attributes
cl_id2cpos <- function(corpus, p_attribute, id, registry = Sys.getenv("CORPUS_REGISTRY")){
  check_registry(registry)
  check_corpus(corpus, registry)
  check_p_attribute(p_attribute = p_attribute, corpus = corpus, registry = registry)
  check_id(corpus = corpus, p_attribute = p_attribute, id = id, registry = registry)
  .cl_id2cpos(corpus = corpus, p_attribute = p_attribute, id = id, registry = registry)
}

#' Drop loaded corpus.
#' 
#' Remove a corpus from the list of loaded corpora of the corpus library (CL).
#' 
#' The corpus library (CL) internally maintains a list of corpora including
#' information on positional and structural attributes so that the registry file
#' needs not be parsed again and again. However, when an attribute has been
#' added to the corpus, it will not yet be visible, because it is not part of
#' the data that has been loaded. The \code{cl_delete_corpus} function exposes a
#' CL function named identically, to force reloading the corpus (after it has
#' been deleted), which will include parsing an updated registry file.
#' 
#' @param corpus name of a CWB corpus (upper case) 
#' @param registry path to the registry directory, defaults to the value of the
#'   environment variable CORPUS_REGISTRY
#' @export cl_delete_corpus
cl_delete_corpus <- function(corpus, registry = Sys.getenv("CORPUS_REGISTRY")){
  .cl_delete_corpus(corpus = corpus, registry = registry)
}

#' Get charset of a corpus.
#' 
#' The encoding of a corpus is declared in the registry file (corpus property
#' "charset"). Once a corpus is loaded, this information is available without
#' parsing the registry file again and again. The \code{cl_charset_name} offers
#' a quick access to this information.
#' 
#' @param corpus Name of a CWB corpus (upper case).
#' @param registry Path to the registry directory, defaults to the value of the
#'   environment variable CORPUS_REGISTRY
#' @export cl_charset_name
#' @examples
#' cl_charset_name(
#'   corpus = "REUTERS",
#'   registry = system.file(package = "RcppCWB", "extdata", "cwb", "registry")
#' )
cl_charset_name <- function(corpus, registry = Sys.getenv("CORPUS_REGISTRY")){
  .cl_charset_name(corpus = corpus, registry = registry)
}

