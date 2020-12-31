/* 
 *  IMS Open Corpus Workbench (CWB)
 *  Copyright (C) 1993-2006 by IMS, University of Stuttgart
 *  Copyright (C) 2007-     by the respective contributers (see file AUTHORS)
 * 
 *  This program is free software; you can redistribute it and/or modify it
 *  under the terms of the GNU General Public License as published by the
 *  Free Software Foundation; either version 2, or (at your option) any later
 *  version.
 * 
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
 *  Public License for more details (in the file "COPYING", or available via
 *  WWW at http://www.gnu.org/copyleft/gpl.html).
 */

#include "../cl/globals.h"
#include <math.h>
/* apparently we don't need to include <search.h> for tfind() etc. and it doesn't */
/* exist in Gnuwin32, so we'll just leave it out */
/* #include <search.h> */

#include "../cl/attributes.h"
#include "../cl/cdaccess.h"
#include "../cl/macros.h"

#include "../cl/binsert.h"

#include "options.h"
#include "html-print.h"
#include "ascii-print.h"
#include "sgml-print.h"
#include "latex-print.h"

#include "corpmanag.h"
#include "groups.h"
#include "output.h"

#include "html-print.h"


#define GROUP_DEBUG 0

#define ANY_ID -2

/* ---------------------------------------------------------------------- */

/* support function for qsort() of grouping cells, controlled by global variables */
Group *compare_cells_group = NULL;

static int 
compare_cells(const void *p1, const void *p2)
{
  ID_Count_Mapping *cell1 = (ID_Count_Mapping *)p1;
  ID_Count_Mapping *cell2 = (ID_Count_Mapping *)p2;
  char *w1, *w2;
  int res, f1, f2;
  
  if (compare_cells_group->is_grouped) {
    /* grouped sort: order by s_freq, then source, then freq, then target */
    f1 = cell1->s_freq;
    f2 = cell2->s_freq;
    res = (f2 > f1) - (f2 < f1); /* corresponds to f2 <=> f1 in Perl */
    if (res != 0) return res;

    w1 = Group_id2str(compare_cells_group, cell1->s, 0);
    w2 = Group_id2str(compare_cells_group, cell2->s, 0);
    res = cl_strcmp(w1, w2);
    if (res != 0) return res;
    
    f1 = cell1->freq;
    f2 = cell2->freq;
    res = (f2 > f1) - (f2 < f1);
    if (res != 0) return res;

    w1 = Group_id2str(compare_cells_group, cell1->t, 1);
    w2 = Group_id2str(compare_cells_group, cell2->t, 1);
    res = cl_strcmp(w1, w2);
    return res;
  }
  else {
    /* ungrouped sort: order by freq, then source, then target */
    f1 = cell1->freq;
    f2 = cell2->freq;
    res = (f2 > f1) - (f2 < f1);
    if (res != 0) return res;

    w1 = Group_id2str(compare_cells_group, cell1->s, 0);
    w2 = Group_id2str(compare_cells_group, cell2->s, 0);
    res = cl_strcmp(w1, w2);
    if (res != 0) return res;

    w1 = Group_id2str(compare_cells_group, cell1->t, 1);
    w2 = Group_id2str(compare_cells_group, cell2->t, 1);        
    res = cl_strcmp(w1, w2);
    return res;
  }
} 


int
get_group_id(Group *group, int i, int target) {
  CorpusList *cl = group->my_corpus;
  int field_type = (target) ? group->target_field : group->source_field;
  int offset = (target) ? group->target_offset : group->source_offset;
  Attribute *attr = (target) ? group->target_attribute : group->source_attribute;
  int is_struc = (target) ? group->target_is_struc : group->source_is_struc;
  char *base = (target) ? group->target_base : group->source_base;
  int pos, id;
  
  switch (field_type) {
  case KeywordField:
    pos = cl->keywords[i];
    break;
  case TargetField:
    pos = cl->targets[i];
    break;
  case MatchField:
    pos = cl->range[i].start;
    break;
  case MatchEndField:
    pos = cl->range[i].end;
    break;
  case NoField:
    pos = ANY_ID;
    break;
  default:
    assert(0 && "Can't be");
    break;
  }
  if (pos >= 0)
    pos += offset;
  if (pos < 0) 
    id = -1;
  else {
    if (is_struc) {
      char *str = cl_cpos2struc2str(attr, pos);
      if (str) 
        id = str - base;
      else
        id = -1;
    }
    else {
      id = cl_cpos2id(attr, pos);
    }
  }
  return id;
}

char *
Group_id2str(Group *group, int id, int target) {
  Attribute *attr = (target) ? group->target_attribute : group->source_attribute;
  int is_struc = (target) ? group->target_is_struc : group->source_is_struc;
  char *base = (target) ? group->target_base : group->source_base;
  if (id == ANY_ID) 
    return "(all)";
  else if (id < 0)
    return "(none)";
  else if (is_struc) 
    return base+id;             /* keep fingers crossed ... */
  else 
    return cl_id2str(attr, id);
}

Group *
ComputeGroupInternally(Group *group)
{
  cl_ngram_hash pairs, groups; /* frequency counts for (source, target) pairs and groups (= source ID) */
  cl_ngram_hash_entry item;

  int i;
  int s_t[2]; /* pair of source and target IDs */
  size_t nr_nodes;
  int percentage, new_percentage; /* for ProgressBar */
  int size = group->my_corpus->size;

  /* ---------------------------------------------------------------------- */

  if (progress_bar)
    progress_bar_clear_line();
  percentage = -1;

  pairs = cl_new_ngram_hash(2, 0);
  groups = cl_new_ngram_hash(1, 0);
  EvaluationIsRunning = 1;

  for (i = 0; i < size; i++) {
    if (! EvaluationIsRunning)
      break;                    /* user abort (Ctrl-C) */

    if (progress_bar) {
      new_percentage = floor(0.5 + (100.0 * i) / size);
      if (new_percentage > percentage) {
        percentage = new_percentage;
        progress_bar_percentage(1, 2, percentage);
      }
    }

    s_t[0] = get_group_id(group, i, 0);       /* source ID */
    s_t[1] = get_group_id(group, i, 1);       /* target ID */
    cl_ngram_hash_add(pairs, s_t, 1);         /* count frequency of (source, target) pair */
    cl_ngram_hash_add(groups, s_t, 1);        /* frequency counts for groups (source) */  
  }

  if (EvaluationIsRunning) {

    if (progress_bar) 
      progress_bar_message(2, 2, " cutoff freq.");

    /* extract vector of pairs above the specified frequency threshold */
    group->count_cells = (ID_Count_Mapping *)cl_malloc(cl_ngram_hash_size(pairs) * sizeof(ID_Count_Mapping));
    nr_nodes = 0;
    cl_ngram_hash_iterator_reset(pairs);
    while ((item = cl_ngram_hash_iterator_next(pairs)) != NULL) {
      if (item->freq >= group->cutoff_frequency) {
        group->count_cells[nr_nodes].s = item->ngram[0];
        group->count_cells[nr_nodes].t = item->ngram[1];
        group->count_cells[nr_nodes].freq = item->freq;
        group->count_cells[nr_nodes].s_freq = cl_ngram_hash_freq(groups, item->ngram);
        nr_nodes++;
      }
    }
    
    /* free unused memory if frequency threshold has filtered out items */
    if (nr_nodes < cl_ngram_hash_size(pairs))
      group->count_cells = (ID_Count_Mapping *)cl_realloc(group->count_cells, (nr_nodes * sizeof(ID_Count_Mapping)));
    group->nr_cells = nr_nodes;

    if (progress_bar) 
      progress_bar_message(2, 2, " sorting rslt");

    /* now sort entries by decreasing frequency, breaking ties in cl_strcmp() order */
    compare_cells_group = group;
    qsort(group->count_cells, nr_nodes, sizeof(ID_Count_Mapping), compare_cells);
    
    if (progress_bar) {
      progress_bar_percentage(2, 2, 100); /* so total percentage runs up to 100% */
      progress_bar_clear_line();
    }
    
  }
  else {
    cqpmessage(Warning, "Group operation aborted by user.");
    if (which_app == cqp) install_signal_handler();
    free_group(&group);         /* sets return value to NULL to indicate failure */
  }
  EvaluationIsRunning = 0;
  
  cl_delete_ngram_hash(pairs); /* free memory */
  cl_delete_ngram_hash(groups);
  
  return group;
}

/*
Group *
ComputeGroupExternally(Group *group)
{
  int i;
  int size = group->my_corpus->size;
  int cutoff_freq = group->cutoff_frequency;

  char temporary_name[TEMP_FILENAME_BUFSIZE];
  FILE *fd;
  FILE *pipe;
  char sort_call[CL_MAX_LINE_LENGTH];

  if ((fd = open_temporary_file(temporary_name)) == NULL) {
    perror("Error while opening temporary file");
    cqpmessage(Warning, "Can't open temporary file");
    return group;
  }

  for (i = 0; i < size; i++) {
    Rprintf("%d %d\n", get_group_id(group, i, 0), get_group_id(group, i, 1));
  }
  fclose(fd);

  sprintf(sort_call, ExternalGroupingCommand, temporary_name);
  if (GROUP_DEBUG)
    Rprintf("Running grouping sort: \n\t%s\n",
            sort_call);
  if ((pipe = popen(sort_call, "r")) == NULL) {
    perror("Failure opening grouping pipe");
    cqpmessage(Warning, "Can't open grouping pipe:\n%s\n"
               "Disable external grouping by\n"
               "  set UseExternalGrouping off;", 
               sort_call);
  }
  else {
    int freq, p1, p2, tokens;
#define GROUP_REALLOC 1024

    while ((tokens = fscanf(pipe, "%d%d%d", &freq, &p1, &p2)) == 3) {
      if (freq > cutoff_freq) {
        if ((group->nr_cells % GROUP_REALLOC) == 0) {
          if (group->count_cells == NULL) {
            group->count_cells = 
              (ID_Count_Mapping *)cl_malloc(GROUP_REALLOC *
                                         sizeof(ID_Count_Mapping));
          }
          else {
            group->count_cells = 
              (ID_Count_Mapping *)cl_realloc(group->count_cells,
                                          (group->nr_cells + GROUP_REALLOC) *
                                          sizeof(ID_Count_Mapping));
          }
          assert(group->count_cells);
        }

        group->count_cells[group->nr_cells].s = p1;
        group->count_cells[group->nr_cells].t = p2;
        group->count_cells[group->nr_cells].freq = freq;
        group->count_cells[group->nr_cells].s_freq = 0;

        group->nr_cells = group->nr_cells + 1;
      }
    }

    if (tokens != EOF) {
      Rprintf("Warning: could not reach EOF of temporary file!\n");
    }

    pclose(pipe);
  }

  if (GROUP_DEBUG) {
    Rprintf("Keeping temporary file %s -- delete manually\n",
            temporary_name);
  }
  else if (unlink(temporary_name) != 0) {
    perror(temporary_name);
    Rprintf("Can't remove temporary file %s -- \n\tI will continue, "
            "but you should remove that file.\n", temporary_name);
  }
  
  return group;
}
*/

Group *compute_grouping(CorpusList *cl,
                        FieldType source_field,
                        int source_offset,
                        char *source_attr_name,
                        FieldType target_field,
                        int target_offset,
                        char *target_attr_name,
                        int cutoff_freq,
                        int is_grouped)
{
  Group *group;
  Attribute *source_attr, *target_attr;
  int source_is_struc = 0, target_is_struc = 0;
  char *source_base = NULL, *target_base = 0;

  if ((cl == NULL) || (cl->corpus == NULL)) {
    cqpmessage(Warning, "Grouping:\nCan't access corpus.");
    return NULL;
  }

  if ((cl->size == 0) || (cl->range == NULL)) {
    cqpmessage(Warning, "Corpus %s is empty, no grouping possible",
               cl->name);
    return NULL;
  }

  if ((source_attr_name == NULL) && (source_field == NoField)) {
    source_attr = NULL;
  }
  else {
    source_attr = find_attribute(cl->corpus, source_attr_name, ATT_POS, NULL);
    if (source_attr == NULL) {
      source_attr = find_attribute(cl->corpus, source_attr_name, ATT_STRUC, NULL);
      source_is_struc = 1;
    }
    if (source_attr == NULL) {
      cqpmessage(Error, "Can't find attribute ``%s'' for named query %s",
                 source_attr_name, cl->name);
      return NULL;
    }
    if (source_is_struc) {
      if (cl_struc_values(source_attr)) {
        source_base = cl_struc2str(source_attr, 0); /* should be beginning of the attribute's lexicon */
        assert(source_base && "Internal error. Please don't use s-attributes in group command.");
      }
      else {
        cqpmessage(Error, "No annotated values for s-attribute ``%s'' in named query %s",
                   source_attr_name, cl->name);
        return NULL;
      }
    }

    switch (source_field) {
    case KeywordField:
      if (cl->keywords == NULL) {
        cqpmessage(Error, "No keyword anchors defined for %s", cl->name);
        return NULL;
      }
      break;
      
    case TargetField:
      if (cl->targets == NULL) {
        cqpmessage(Error, "No target anchors defined for %s", cl->name);
        return NULL;
      }
      break;
      
    case MatchField:
    case MatchEndField:
      assert(cl->range && cl->size > 0);
      break;
      
    case NoField:
    default:
      cqpmessage(Error, "Illegal second anchor in group command");
      return NULL;
      break;
    }
  }

  target_attr = find_attribute(cl->corpus, target_attr_name, ATT_POS, NULL);
  if (target_attr == NULL) {
      target_attr = find_attribute(cl->corpus, target_attr_name, ATT_STRUC, NULL);
      target_is_struc = 1;
  }
  if (target_attr == NULL) {
    cqpmessage(Error, "Can't find attribute ``%s'' for named query %s",
               target_attr_name, cl->name);
    return NULL;
  }
  if (target_is_struc) {
    if (cl_struc_values(target_attr)) {
      target_base = cl_struc2str(target_attr, 0); /* should be beginning of the attribute's lexicon */
      assert(target_base && "Internal error. Please don't use s-attributes in group command.");
    }
    else {
      cqpmessage(Error, "No annotated values for s-attribute ``%s'' in named query %s",
                 target_attr_name, cl->name);
      return NULL;
    }
  }

  switch (target_field) {
  case KeywordField:
    if (cl->keywords == NULL) {
      cqpmessage(Error, "No keyword anchors defined for %s", cl->name);
      return NULL;
    }
    break;
    
  case TargetField:
    if (cl->targets == NULL) {
      cqpmessage(Error, "No target anchors defined for %s", cl->name);
      return NULL;
    }
    break;
    
  case MatchField:
  case MatchEndField:
    assert(cl->range && cl->size > 0);
    break;
    
  case NoField:
  default:
    cqpmessage(Error, "Illegal anchor in group command");
    return NULL;
    break;
  }

  /* set up Group object */
  group = (Group *) cl_malloc(sizeof(Group));
  group->my_corpus = cl;
  group->source_attribute = source_attr;
  group->source_offset = source_offset;
  group->source_is_struc = source_is_struc;
  group->source_base = source_base;
  group->source_field = source_field;
  group->target_attribute = target_attr;
  group->target_offset = target_offset;
  group->target_is_struc = target_is_struc;
  group->target_base = target_base;
  group->target_field = target_field;
  group->nr_cells = 0;
  group->count_cells = NULL;
  group->cutoff_frequency = cutoff_freq;
  group->is_grouped = is_grouped;


  /* if (UseExternalGrouping && !insecure && !(source_is_struc || target_is_struc || is_grouped))
    return ComputeGroupExternally(group);
  else
    return ComputeGroupInternally(group);
   */
  return group;
}

void
free_group(Group **group)
{
  cl_free((*group)->count_cells);
  (*group)->my_corpus = NULL;
  (*group)->nr_cells = 0;
  (*group)->source_attribute = NULL;
  (*group)->target_attribute = NULL;

  free(*group);
  *group = NULL;
}

void print_group(Group *group, int expand, struct Redir *rd)
{
  if (group && open_stream(rd, group->my_corpus->corpus->charset)) {

    switch (GlobalPrintMode) {
    
    case PrintSGML:
      sgml_print_group(group, expand, rd->stream);
      break;
      
    case PrintHTML:
      html_print_group(group, expand, rd->stream);
      break;
      
    case PrintLATEX:
      latex_print_group(group, expand, rd->stream);
      break;
      
    case PrintASCII:
      ascii_print_group(group, expand, rd->stream);
      break;
    
    default:
      cqpmessage(Error, "Unknown print mode");
      break;
    }

    close_stream(rd);
    
  }
}


