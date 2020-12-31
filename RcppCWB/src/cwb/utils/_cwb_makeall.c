void Rprintf(const char *, ...); /* alternative to include R_ext/Print.h */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../cl/cl.h"
#include "../cl/globals.h"
#include "../cl/attributes.h"
#include "../cl/fileutils.h"

#include <netinet/in.h>

Corpus *corpus;


/* COPIED AND SLIGHTLY MODIFIED FROM utils-makeall.c 
 * Modifications:
 * - printf -> Rprintf;
 * - fprintf -> Rprintf (without stderr)
 * - exit -> return
 * - fflush(stdout); commented out (twice)
 * - do_attribute and make_component return integer value;
 */

/**
 * Checks whether a component has already been created.
 *
 * @param attr  The attribute of the component to check.
 * @param cid   The component ID of the component to check.
 * @return      RUE iff the component has already been created.
 */
int
  component_ok(Attribute *attr, ComponentID cid)
  {
    ComponentState state;
    
    state = component_state(attr, cid);
    if ((state == ComponentLoaded) || (state == ComponentUnloaded)) {
      return 1;
    }
    else {
      if (state != ComponentDefined) {
        Rprintf("Internal Error: Illegal state %d/component ID %d ???\n", state, cid);
        return 1;
      }
      return 0;
    }
  }


/**
* Creates a component for the specified attribute.
*
* This function will create the component if it doesn't already exist;
* it aborts on error.
*
* @see create_component
*
* @param attr  The attribute of the component to create.
* @param cid   The component ID of the component to create.
*/
int
  make_component(Attribute *attr, ComponentID cid)
  {
    int state;
    
    if (! component_ok(attr, cid)) {
      
      Rprintf(" + creating %s ... ", cid_name(cid));
      /* fflush(stdout); */
      (void) create_component(attr, cid);
      
      state = component_state(attr, cid);
      if (!(state == ComponentLoaded || state == ComponentUnloaded)) {
        Rprintf("FAILED\n");
        Rprintf("ERROR. Aborted.\n");
        return 1;
      }
      
      Rprintf("OK\n");
    }
    return 0;
  }



/**
* Validates the REVCORP component of the given attribute.
*
* This function validates a REVCORP (i.e. an uncompressed index).
* It assumes that a lexicon, frequencies and (compressed or
* uncompressed) token stream are available for CL access for the
* given attribute.
*
* @param attr  The attribute whose REVCORP should be checked.
* @return      True for all OK, false for a problem.
*/
int
  validate_revcorp(Attribute *attr)
  {
    
    Component *revcorp = ensure_component(attr, CompRevCorpus, 0);
    int *ptab;                        /* table of index offsets for each lexicon entry */
int lexsize, corpsize;
int i, offset, cpos, id;

Rprintf(" ? validating %s ... ", cid_name(CompRevCorpus));
/* fflush(stdout); */

if (revcorp == NULL) {
  Rprintf("FAILED (no data)\n");
  return 0;
}
lexsize = cl_max_id(attr);
corpsize = cl_max_cpos(attr);
if ((lexsize <= 0) || (corpsize <= 0)) {
  Rprintf("FAILED (corpus access error)\n");
  return 0;
}
if (revcorp->size != corpsize) {
  Rprintf("FAILED (wrong size)\n");
  return 0;
}

/* init offsets by calculating REVIDX component from token frequencies */
ptab = (int *) cl_calloc(lexsize, sizeof(int));
offset = 0;
for (i = 0; i < lexsize; i++) {
  ptab[i] = offset;
  offset += cl_id2freq(attr, i);
}

/* now read token stream, check each token id against REVCORP, and increment its pointer */
for (cpos = 0; cpos < corpsize; cpos++) {
  id = cl_cpos2id(attr, cpos);
  if ((id < 0) || (id >= lexsize)) {
    Rprintf("FAILED (inconsistency in token stream)\n");
    cl_free(ptab);
    return 0;
  }
  if (ntohl(revcorp->data.data[ptab[id]]) != cpos) {
    Rprintf("FAILED\n");
    cl_free(ptab);
    return 0;
  }
  ptab[id]++;
}

/* validate frequencies by comparing final offsets against those calculated from token frequencies */
offset = 0;
for (i = 0; i < lexsize; i++) {
  offset += cl_id2freq(attr, i);
  if (ptab[i] != offset) {
    Rprintf("FAILED (token frequencies incorrect)\n");
    cl_free(ptab);
    return 0;
  }
}

cl_free(ptab);

Rprintf("OK\n");
return 1;
  }

/**
* Create a given component (or all components) for an attribute.
*
* @param attr      The attribute to work on.
* @param cid       If this is CompLast, all components will be created.
*                  Otherwise, it specifies the single component that will
*                  be created.
* @param validate  boolean - if true, validate_revcorp is called to check
*                  the resulting revcorp.
*/
int
  do_attribute(Attribute *attr, ComponentID cid, int validate)
  {
    assert(attr);
    
    if (cid == CompLast) {
      Rprintf("ATTRIBUTE %s\n", attr->any.name);
      /* automatically create all necessary components */
      
      /* check whether directory for data files exists (may be misspelt in registry) */
      if (! is_directory(attr->any.path)) {
        Rprintf("WARNING. I cannot find the data directory of the '%s' attribute.\n",
                attr->any.name);
        Rprintf("WARNING  Directory: %s/ \n", attr->any.path);
        Rprintf("WARNING  Perhaps you misspelt the directory name in the registry file?\n");
      }
      
      /* lexicon and lexicon offsets must have been created by encode */
      if (! (component_ok(attr, CompLexicon) && component_ok(attr, CompLexiconIdx))) {
        /* if none of the components exits, we assume that the attribute will be created later & skip it */
        if (!component_ok(attr, CompLexicon) && !component_ok(attr, CompLexiconIdx) &&
        !component_ok(attr, CompLexiconSrt) &&
        !component_ok(attr, CompCorpus) && !component_ok(attr, CompCorpusFreqs) &&
        !component_ok(attr, CompHuffSeq) && !component_ok(attr, CompHuffCodes) &&
        !component_ok(attr, CompHuffSync) &&
        !component_ok(attr, CompRevCorpus) && !component_ok(attr, CompRevCorpusIdx) &&
        !component_ok(attr, CompCompRF) && !component_ok(attr, CompCompRFX))
        {
          /* issue a warning message & return */
          Rprintf(" ! attribute not created yet (skipped)\n");
          if (strcmp(attr->any.name, "word") == 0) {
            Rprintf("WARNING. The 'word' attribute must be created before using CQP on this corpus!\n");
          }
          return 0;
        }
        else {
          Rprintf("ERROR. Lexicon is missing. You must use the 'encode' tool first!\n");
          return 1;
        }
      }
      else {
        /* may need to create "alphabetically" sorted lexicon */
        make_component(attr, CompLexiconSrt);
        Rprintf(" - lexicon      OK\n");
      }
      
      /* create token frequencies if necessary (must be able to do so if they aren't already there) */
      make_component(attr, CompCorpusFreqs);
      Rprintf(" - frequencies  OK\n");
      
      /* check if token sequence has been compressed, otherwise create CompCorpus (if necessary) */
      if (component_ok(attr, CompHuffSeq) && component_ok(attr, CompHuffCodes) && component_ok(attr, CompHuffSync)) {
        Rprintf(" - token stream OK (COMPRESSED)\n");
      }
      else {
        make_component(attr, CompCorpus);
        Rprintf(" - token stream OK\n");
      }
      
      /* same for index (check if compressed, otherwise create if not already there) */
      if (component_ok(attr, CompCompRF) && component_ok(attr, CompCompRFX)) {
        Rprintf(" - index        OK (COMPRESSED)\n");
      }
      else {
        make_component(attr, CompRevCorpusIdx);
        if (! component_ok(attr, CompRevCorpus)) { /* need this check to avoid validation of existing revcorp  */
      make_component(attr, CompRevCorpus);
          if (validate) {
            /* validate the index, i.e. the REVCORP component we just created */
            if (! validate_revcorp(attr)) {
              Rprintf("ERROR. Validation failed.\n");
              return 1;
            }
          }
        }
        Rprintf(" - index        OK\n");
      }
    }
    else {
      /* create requested component only */
      Rprintf("Processing component %s of ATTRIBUTE %s\n",
             cid_name(cid), attr->any.name);
      make_component(attr, cid);
      if (validate && (cid == CompRevCorpus)) { /* validates even if REVCORP already existed -> useful trick for validating later */
      if (! validate_revcorp(attr)) {
        Rprintf("ERROR. Validation failed.\n");
        return 1;
      }
      }
    }
    return 0;
    
  }


