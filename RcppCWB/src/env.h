#include "corpmanag.h" /* get ContextDescriptor,AttributeList */

#define MAXPATTERNS 5000
#define MAXENVIRONMENT 10

enum ctxtdir { leftright, left, right };
enum spacet { word, structure };

/**
 * Evaltree object
 *
 * (TODO specify here what on earth one of these is)
 */
typedef union e_tree *Evaltree;


typedef struct _label_entry {
  int        flags;
  char      *name;
  int        ref;             /**< array index the label refers to */
  struct _label_entry *next;
} *LabelEntry;

typedef struct ctxtsp {
  enum ctxtdir   direction;     /**< direction of context expansion (if valid).
 Might be left, right, or leftright*/
  enum spacet    type;          /**< kind of space (word or structure)         */
  Attribute     *attrib;        /**< attribute representing the structure.     */
  int            size;          /**< size of space in number of structures.    */
  int            size2;         /**< only for meet-context                     */
} Context;

typedef enum _avstype {
  Pattern, Tag, MatchAll, Anchor
} AVSType;


enum bnodetype { bnode,                 /**< boolean evaluation node            */
                 cnode,                 /**< constant node                      */
                 func,                  /**< function call                      */
                 sbound,                /**< structure boundary (open or close) */
                 pa_ref,                /**< reference to positional attribute  */
                 sa_ref,                /**< reference to structural attribute  */

                 string_leaf,           /**< string constant */
                 int_leaf,              /**< integer constant */
                 float_leaf,            /**< float constant */

                 id_list,               /**< list of IDs */
                 var_ref                /**< variable reference */
               };


enum b_ops { b_and,	 /**< boolean and operator           */
             b_or,	 /**< boolean or operator            */
             b_implies,  /**< boolean implication (->) operator */
             b_not,	 /**< boolean negation               */

             cmp_gt,     /**< compare: greater than          */
             cmp_lt,     /**< compare: less than             */
             cmp_get,    /**< compare: greater or equal than */
             cmp_let,    /**< compare: less or equal than    */
             cmp_eq,     /**< compare: equal                 */
             cmp_neq,    /**< compare: not equal             */

             cmp_ex	 /**< is value present? bool exprs   */
           };

enum wf_type { NORMAL, REGEXP, CID };



/**
 * Union of structures underlying the Constraint / Constrainttree objects.
 *
 * Each Constraint is a node in the Constrainttree, i.e. a single element of a compiled CQP query.
 */
typedef union c_tree {

  /** The type of this particular node.
   * Allows the type member of the other structures within the union to be accessed. */
  enum bnodetype type;

  /** "standard" operand node in the evaluation tree; type is "bnode" */
  struct {
    enum bnodetype type;                  /**< must be bnode                     */
    enum b_ops     op_id;                 /**< identifier of the bool operator   */
    union c_tree  *left,                  /**< points to the first operand       */
                  *right;                 /**< points to the second operand,
                                               if present                        */
  }                node;

  /** "constant" node in the evaluation tree */
  struct {
    enum bnodetype type;                  /**< must be cnode                     */
    int            val;                   /**< Value of the constant: 1 or 0 for true or false */
  }                constnode;

  /** function call (dynamic attribute), type is "func" */
  struct {
    enum bnodetype type;                  /**< must be func                  */
    int            predef;
    Attribute     *dynattr;
    struct _ActualParamList *args;        /**< arguments of the function     */
    int            nr_args;               /**< nr of arguments for this call */
  }                func;

  /** structure boundary */
  struct {
    enum bnodetype type;                  /**< must be sbound                */
    Attribute     *strucattr;             /**< the attribute which corresponds to the structure */
    Boolean        is_closing;            /**< True if closing tag, False for opening tag */
  }                sbound;

  /** reference to positional attribute */
  struct {
    enum bnodetype type;                  /**< must be pa_ref */
    LabelEntry     label;                 /**< may be empty (NULL) */
    Attribute     *attr;                  /**< the P-attribute we are referring to */
    int            del;                /**< delete label after using it ? */
  }                pa_ref;

  /**
   * reference to structural attribute.
   *
   * If label is empty, this checks if the current position is at start
   * or end of structural_attribute and returns INT value (this is kept for
   * backward compatibility regarding lbound() and rbound() builtins; the new
   * syntax is to use {s} and {/s}, which are represented as 'Tag' nodes.
   *
   * If label is non-empty, the referenced S-attribute must have values, and
   * the value of the enclosing region is returned as a string; in short,
   * values of attributes can be accessed through label references .
   */
  struct {
    enum bnodetype type;                  /**< must be sa_ref */
    LabelEntry     label;                 /**< may be empty (NULL) */
    Attribute     *attr;                  /**< the s-attribute we are referring to */
    int            del;                /**< delete label after using it ? */
  }                sa_ref;

  struct {
    enum bnodetype type;                  /**< must be var_ref */
    char          *varName;
  }                varref;

  struct {
    enum bnodetype type;                  /**< must be id_list */
    Attribute     *attr;
    LabelEntry     label;                 /**< may be empty (NULL) */
    int            negated;
    int            nr_items;
    int           *items;                 /**< an array of item IDs of size nr_items */
    int            del;                /**< delete label after using it ? */
  }                idlist;

  /** constant (string, int, float, ...) */
  struct {
    enum bnodetype type;                  /**< string_leaf, int_leaf, or float_leaf */

    int            canon;                 /**< canonicalization mode (i.e. flags)         */
    enum wf_type   pat_type;              /**< pattern type: normal wordform or reg. exp. */
    CL_Regex       rx;                    /**< compiled regular expression (using CL frontend) */

    /** Union containing the constant type. */
    union {
      char        *sconst;               /**< operand is a string constant.           */
      int          iconst;               /**< operand is a integer constant.          */
      int          cidconst;             /**< operand is {?? corpus position?? corpus lexicon id??} constant */
      double       fconst;               /**< operand is a float (well, double) constant */
    }              ctype;
  }                leaf;
} Constraint;

/**
 * The Constrainttree object.
 */
typedef Constraint *Constrainttree;





/**
 * The AVStructure object.
 *
 * A union of structures with the type member always accessible.
 */
typedef union _avs {

  /** What type of AV structure does this union represent? */
  AVSType type;

  /** a matchall item */
  struct {
    AVSType type;                /* set to MatchAll */
    LabelEntry label;
    Boolean is_target;
    Boolean lookahead;           /**< whether pattern is just a lookahead constraint */
  } matchall;

  /** a constraint tree */
  struct {
    AVSType type;                /* set to Pattern */
    LabelEntry label;
    Constrainttree constraint;
    Boolean is_target;
    Boolean lookahead;           /**< whether pattern is just a lookahead constraint */
  } con;

  /** a structure describing tag */
  struct {
    AVSType type;                /* set to Tag */
    int is_closing;
    Attribute *attr;
    char *constraint;            /**< constraint for annotated value of region (string or regexp); NULL = no constraint */
    int flags;                   /**< flags passed to regexp or string constraint (information purposes only) */
    CL_Regex rx;                 /**< if constraint is a regexp, this holds the compiled regexp; otherwise NULL */
    int negated;                 /**< whether constraint is negated (!=, not matches, not contains) */
    LabelEntry right_boundary;   /**< label in RDAT namespace: contains right boundary of constraining region (in StrictRegions mode) */
  } tag;

  /* an anchor point tag (used in subqueries) */
  struct {
    AVSType type;                /* set to Anchor */
    int is_closing;
    FieldType field;
  } anchor;
} AVStructure;

/** AVS is a pointer type for AVStructure */
typedef AVStructure *AVS;

/** Patternlist is an array of AVStructures */
typedef AVStructure Patternlist[MAXPATTERNS];




typedef struct _symbol_table {
  LabelEntry  user;                /**< user namespace */
  LabelEntry  rdat;                /**< namespace for LAB_RDAT labels */
  int next_index;                  /**< next free reference table index */
} *SymbolTable;


typedef struct dfa {
  int Max_States;         /**< max number of states of the current dfa;
 state no. 0 is the initial state.             */
int Max_Input;          /**< max number of input chars of the current dfa. */
int **TransTable;       /**< state transition table of the current dfa.    */
Boolean *Final;         /**< set of final states.                          */
int E_State;            /**< Error State -- it is introduced in order to
 *   make the dfa complete, so the state transition
 *   is a total mapping. The value of this variable
 *   is Max_States.
 */
} DFA;





/**
 * The EvalEnvironment object: environment variables for the evaluation of
 * a corpus query.
 */
typedef struct evalenv {
  
  CorpusList *query_corpus;         /**< the search corpus for this query part */

int rp;                           /**< index of current range (in subqueries) */

SymbolTable labels;               /**< symbol table for labels */

int MaxPatIndex;                  /**< the current number of patterns */
Patternlist patternlist;          /**< global variable which holds the pattern list */

Constrainttree gconstraint;       /**< the "global constraint" */

Evaltree evaltree;                /**< the evaluation tree (with regular exprs) */

DFA  dfa;                         /**< the regex DFA for the current query */

int has_target_indicator;         /**< is there a target mark ('@') in the query? */
LabelEntry target_label;          /**< targets are implemented as a special label "target" now */

LabelEntry match_label;           /**< special "match" and "matchend"-Labels for access
 to start & end of match within query */
LabelEntry matchend_label;

Context search_context;           /**< the search context (within...) */

Attribute *aligned;               /**< the attribute holding the alignment info */

int negated;                      /**< 1 iff we should negate alignment constr */

enum _matching_strategy matching_strategy; /**< copied from global option unless overwritten by (?...) directive */

} EvalEnvironment;

/**
 * EEPs are Eval Environment pointers.
 */
typedef EvalEnvironment *EEP;

/** A global array of EvalEnvironment structures */
EvalEnvironment Environment[MAXENVIRONMENT];

EEP CurEnv, evalenv;
