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



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "../cl/macros.h"

#include "symtab.h"
#include "output.h"


/*
 * some internal helper functions 
 */

/**
 * Frees the whole of a linked list of labels entries.
 */
void
free_labellist(LabelEntry l)
{
  LabelEntry tmp;

  while (l) {
    tmp = l;
    l = tmp->next;

    if (tmp->name)
      free(tmp->name);
    free(tmp);
  }
}


/*
 * The SYMBOL LOOKUP part: SymbolTable and LabelEntry 
 */

/** create new symbol table */
SymbolTable
new_symbol_table(void)
{
  SymbolTable st = cl_malloc(sizeof(struct _symbol_table));
  st->user = (LabelEntry) NULL;
  st->rdat = (LabelEntry) NULL;
  st->next_index = 0;
  return st;
}

/** delete symbol table (free all entries) */
void
delete_symbol_table(SymbolTable st)
{
  if (st != NULL) {
    free_labellist(st->user);
    free_labellist(st->rdat);
    free(st);
  }
}

/** Returns label entry, or NULL if undefined (flags are used _only_ to determine namespace) */
LabelEntry 
findlabel(SymbolTable st, char *s, int flags)
{
  LabelEntry l;

  if (flags & LAB_RDAT)
    l = st->rdat;
  else
    l = st->user;
  
  while (l != NULL) {
    if (strcmp(l->name, s) == 0)
      return(l);
    l = l->next;
  }

  return ((LabelEntry) NULL);
}

/**
 * Look up a label, add flags, and return label entry (NULL if undefined).
 * if create is set and label does not exist, it is added to the symbol table
 */
LabelEntry 
labellookup(SymbolTable st, char *s, int flags, int create)
{
  LabelEntry l;
  int user_namespace, this_label, is_special;
  
  if ((l = findlabel(st, s, flags)) != NULL) {
    l->flags |= flags;                /* add flags from this call */
    return l;
  }

  /* are we in user namespace? (otherwise namespace is rdat, so far) */
  user_namespace = (flags & LAB_RDAT) ? 0 : 1;

  /* check if the requested label is the 'this' label (label _ in user namespace) */
  this_label = (user_namespace && (strcmp(s, "_") == 0)) ? 1 : 0;

  /* add label to symbol table if create is set ('_' is always defined and created automatically) */
  if (create || this_label) {

    /* special labels: 'this' and all field names (used to keep track of target anchors in queries) */
    is_special = (this_label || (field_name_to_type(s) != NoField)) ? 1 : 0;
    if (is_special) flags |= LAB_SPECIAL;
    
    /* allocate new label entry and initialise data fields */
    l = cl_malloc(sizeof(struct _label_entry));
    l->name = cl_strdup(s);
    l->flags = flags;
    if (this_label)
      l->ref = -1;
    else
      l->ref = st->next_index++;

    /* insert into correct namespace (so far only user or rdat) */
    if (user_namespace) {
      l->next = st->user;
      st->user = l;
    }
    else {
      l->next = st->rdat;
      st->rdat = l;
    }

    return l;

  }
  else {
    return NULL;
  }
}

/**
 * Drops a label from the symbol table (NB: its reference index can't be re-used).
 * WARNING: this function is not implemented (doesn't seem to be useful at the moment)
 */
void 
droplabel(SymbolTable st, LabelEntry l)
{
  assert(0 && "droplabel(): function not implemented");
}

/**
 * Checks whether all used labels are defined (and vice versa).
 * [only non-special labels in the user namespace will be checked]
 */
int 
check_labels(SymbolTable st)
{
  LabelEntry l = st->user;
  int result = 1;

  while (l != NULL) {
    if (! (l->flags & LAB_SPECIAL)) {
      if (! (l->flags & LAB_USED)) {
        cqpmessage(Warning, "Label %s defined but not used", l->name);
        result = 0;
      }
      if (! (l->flags & LAB_DEFINED)) {
        cqpmessage(Warning, "Label %s used but not defined", l->name);
        result++;
      }
    }
    l = l->next;
  }
  return result;
}

/** print symbol table contents (for debugging purposes) */
void 
print_symbol_table(SymbolTable st)
{
  LabelEntry l;
  char *namespace;
  int i;

  Rprintf("Contents of SYMBOL TABLE:\n");
  
  for (i = 1; i <= 2; i++) {
    switch (i) {
    case 1:
      l = st->user;
      namespace = "USER";
      break;
    case 2:
      l = st->rdat;
      namespace = "RDAT";
      break;
    default:
      l = NULL;
      namespace = "???";
    }
    while (l != NULL) {
      Rprintf("\t%s\t%s(flags: %d)  ->  RefTab[%d]\n",
              namespace, l->name, l->flags, l->ref);
      l = l->next;
    }
  }
}

LabelEntry 
symbol_table_new_iterator(SymbolTable st, int flags)
{
  int user_namespace = (flags & LAB_RDAT) ? 0 : 1;
  LabelEntry lab;
  
  if (st != NULL) {
    lab = (user_namespace) ? st->user : st->rdat;
    if ((lab != NULL) && ((lab->flags & flags) != flags)) {
      lab = symbol_table_iterator(lab, flags); /* if first label in list doesn't match flags, find next matching label */
    }
    return lab;
  }
  else 
    return NULL;
}

LabelEntry 
symbol_table_iterator(LabelEntry prev, int flags)
{
  if (prev == NULL) 
    return NULL;
  else {
    prev = prev->next;
    while (prev != NULL) {
      if ((prev->flags & flags) == flags) 
        break;
      prev = prev->next;
    }
    return prev;
  }
}




/* 
 * The DATA ARRAY part: RefTab 
 */

/**
 * Create new reference table of required size for the given symbol table.
 *
 * NB If further labels are added to st, you must reallocate the reference table
 * to make room for the new reference indices
 */
RefTab 
new_reftab(SymbolTable st)
{
  RefTab rt;

  rt = (RefTab) cl_malloc(sizeof(struct _RefTab));
  rt->size = st->next_index;
  rt->data = (int *) cl_malloc(rt->size * sizeof(int));

  return rt;
}

/** Deletes (and frees) a reference table. */
void 
delete_reftab(RefTab rt)
{
  if (rt) {
    if (rt->data)
      free(rt->data);
    free(rt);
  }
}

/** Copies rt1 to rt2; doesn't allocate new reftab for efficiency reasons. */
void 
dup_reftab(RefTab rt1, RefTab rt2)
{
  assert(rt1);
  assert(rt2);
  if (rt1->size != rt2->size) {
    Rprintf("dup_reftab()<symtab.c>: Tried to dup() RefTab (%d entries) to RefTab of different size (%d entries)\n", rt1->size, rt2->size);
    exit(1);
  }
  (void) memcpy(rt2->data, rt1->data, rt1->size * sizeof(int));
}

/** resets all referenced corpus position to -1 -> undefine all references */
void 
reset_reftab(RefTab rt)
{
  int i;
  
  assert(rt);
  for (i = 0; i < rt->size; i++)
    rt->data[i] = -1;
}

/** set references (cpos value in get_reftab is returned for 'this' label (_), set to -1 if n/a) */
void 
set_reftab(RefTab rt, int index, int value)
{
  if (rt != NULL) {
    if ((index < 0) || (index >= rt->size)) {
      cqpmessage(Error, "RefTab index #%d not in range 0 .. %d", index, rt->size - 1);
      exit(1);
    }
    else {
      rt->data[index] = value;
    }
  }
}

/** read references (cpos value in get_reftab is returned for 'this' label (_), set to -1 if n/a) */
int 
get_reftab(RefTab rt, int index, int cpos)
{
  if (index == -1)              /* -1 == 'this' label returns <cpos> value */
    return cpos;
  else if (rt == NULL)          /* NULL is used for dummy reftabs */
    return -1;
  else if ((index < 0) || (index >= rt->size)) {
    Rprintf("get_reftab()<symtab.c>: RefTab index #%d not in range 0 .. %d", index, rt->size - 1);
    return -1;
  }
  else
    return (rt->data[index]);
}

/**
 * Prints the current label values (for debugging).
 *
 * @param st    The SymbolTable
 * @param rt
 * @param cpos  The corpus position
 */
void 
print_label_values(SymbolTable st, RefTab rt, int cpos)
{
  LabelEntry l;
  int i;

  Rprintf("Label values:\n");
  if ((st == NULL) || (rt == NULL) || (st->next_index != rt->size)) {
    Rprintf("ERROR\n");
    return;
  }

  for (i = 1; i <= 2; i++) {
    switch (i) {
    case 1:
      l = st->user;
      Rprintf("USER:\t");
      break;
    case 2:
      l = st->rdat;
      Rprintf("RDAT:\t");
      break;
    default:
      l = NULL;
      Rprintf("???");
    }
    while (l != NULL) {
      Rprintf("%s=%d  ", l->name, get_reftab(rt, l->ref, cpos));
      l = l->next;
    }
    Rprintf("\n");
  }
}

