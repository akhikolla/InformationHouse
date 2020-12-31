enum _which_app { undef, cqp, cqpcl, cqpserver} which_app;
int insecure;                     /**< Boolean: != 0 means we should not allow pipes etc. (For example, in CGI.) */
int inhibit_activation;           /**< Boolean: inhibit corpus activations in parser */
int parseonly;                    /**< if true, queries are only parsed, not evaluated. */
int verbose_parser;               /**< if true, absolutely all messages from the parser get printed (inc Message-level). */
int show_symtab;                  /**< Doesn't seem to be used anywhere; should show_environment use it? if not, remove? TODO  */
int show_gconstraints;            /**< if true, the tree of global contraints is printed when an EvalEnvironment is displayed */
int show_evaltree;                /**< if true, the evaluation tree is printed when an EvalEnvironment is displayed */
int show_patlist;                 /**< if true, the pattern list is printed when an EvalEnvironment is displayed */
int show_compdfa;                 /**< if true, the complete DFA is printed when an EvalEnvironment is displayed */
int show_dfa;                     /**< if true, the regex2dfa module will print out the states of the DFA after it is parsed. */
int symtab_debug;                 /**< if this AND debug_simulation are true, print extra messages relating to eval */
int parser_debug;
int tree_debug;                   /**< if true, extra messages are embedded when an evaluation tree is pretty-printed */
int eval_debug;                   /**< if true, assorted debug messages related to query evaluation are printed */
int search_debug;                 /**< if true, the evaltree of a pattern is pretty-printed before the DFA is created. */
int initial_matchlist_debug;      /**< if true, debug messages relating to the initial set of candidate matches are printed. */
int debug_simulation;             /**< if true, debug messages are printed when simulating an NFA. @see simulate */
int activate_cl_debug;            /**< if true, the CL's debug message setting is set to On. */
int server_log;                   /**< cqpserver option: logging (print log messages to standard output) */
int server_debug;                 /**< cqpserver option: debugging output (print debug messages to standard error) */
int snoop;                        /**< cqpserver option: monitor CQi network communication */
int private_server;               /**< cqpserver option: makes CQPserver accept a single connection only */
int server_port;                  /**< cqpserver option: CQPserver's listening port (if 0, listens on CQI_PORT) */
int localhost;                    /**< cqpserver option: accept local connections (loopback) only */
int server_quit;                  /**< cqpserver option: spawn server and return to caller (for CQI::Server.pm) */
int query_lock;                   /**< cqpserver option: safe mode for network/HTTP servers (allow query execution only) */
int query_lock_violation;         /**< cqpserver option: set for CQPserver's sake to detect attempted query lock violation */
int enable_macros;                /**< enable macros only at user request in case they introduce compatibility problems */
int macro_debug;                  /**< enable debugging of macros (and print macro hash stats on shutdown). */
int hard_boundary;                /**< Query option: use implicit 'within' clause (unless overridden by explicit spec) */
int hard_cut;                     /**< Query option: use hard cut value for all queries (cannot be changed) */
int subquery;                     /**< Query option: use auto-subquery mode (TODO rename to auto_subquery for clarity) */
char *def_unbr_attr;              /**< Query option: unbracketed attribute (attribute matched by "..." patterns) */
int query_optimize;               /**< Query option: use query optimisation (untested and expensive optimisations) */
enum _matching_strategy { traditional, shortest_match, standard_match, longest_match } matching_strategy;
char *matching_strategy_name;     /**< The matching strategy option: which is implemented as a vstring option with side-effect */
int strict_regions;               /**< boolean: expression between {s} ... {/s} tags is constrained to single {s} region  */
int use_readline;                 /**< UI option: use GNU Readline for input line editing if available */
int highlighting;                 /**< UI option: highlight match / fields in terminal output? (default = yes) */
int paging;                       /**< UI option: activate/deactivate paging of query results */
char *pager;                      /**< UI option: pager program to used for paged kwic display */
char *tested_pager;               /**< UI option: CQP tests if selected pager works & will fall back to "more" if it doesn't */
char *less_charset_variable;      /**< UI option: name of environment variable for controlling less charset (usually LESSCHARSET) */
int use_colour;                   /**< UI option: use colours for terminal output (experimental) */
int progress_bar;                 /**< UI option: show progress bar during query execution */
int pretty_print;                 /**< UI option: pretty-print most of CQP's output (turn off to simplify parsing of CQP output) */
int autoshow;                     /**< UI option: show query results after evaluation (otherwise, just print number of matches) */
int timing;                       /**< UI option: time queries (printed after execution) */
int show_tag_attributes;          /**< kwic option: show values of s-attributes as SGML tag attributes in kwic lines */
int show_targets;                 /**< kwic option: show numbers of target anchors in brackets */
char *printModeString;            /**< kwic option: string of current printmode */
char *printModeOptions;           /**< kwic option: some printing options */
int printNrMatches;               /**< kwic option: -> 'cat' prints number of matches in first line (do we need this?) */
char *printStructure;             /**< kwic option: show annotations of structures containing match */
char *left_delimiter;             /**< kwic option: the match start prefix (defaults to '<') */
char *right_delimiter;            /**< kwic option: the match end suffix   (defaults to '>') */
char *registry;                   /**< registry directory */
char *LOCAL_CORP_PATH;            /**< directory where subcorpora are stored (saved & loaded) */
int auto_save;                    /**< automatically save subcorpora */
int save_on_exit;                 /**< save unsaved subcorpora upon exit */
char *cqp_init_file;              /**< changed from 'init_file' because of clash with a # define in {term.h} */
char *macro_init_file;            /**< secondary init file for loading macro definitions (not read if macros are disabled) */
char *cqp_history_file;           /**< filename where CQP command history will be saved */
int write_history_file;           /**< Controls whether CQP command history is written to file */
int batchmode;                    /**< set by -f {file} option (don't read ~/.cqprc, then process input from {file}) */
int silent;                       /**< Disables some messages & warnings (used rather inconsistently). */
char *default_corpus;             /**< corpus specified with -D {corpus} */
char *query_string;               /**< query specified on command line (-E {string}, cqpcl only) */
int UseExternalSorting;           /**< (option which should not exist) use external sorting algorithm */
char *ExternalSortingCommand;     /**< (option which should not exist) external sort command to use */
int UseExternalGrouping;          /**< (option which should not exist) use external grouping algorithm */
char *ExternalGroupingCommand;    /**< (option which should not exist) external group command to use */
int user_level;                   /**< (option which should not exist) user level: 0 == normal, 1 == advanced, 2 == expert) */
int rangeoutput;                  /**< (option which should not exist) */
int child_process;
ContextDescriptor CD;
int handle_sigpipe;
char *progname;
char *licensee;
FILE *batchfd;
int initialize_cqp(int argc, char **argv);
int cqp_parse_file(FILE *fd, int exit_on_parse_errors);
int cqp_parse_string(char *s);
int setInterruptCallback(InterruptCheckProc f);
void CheckForInterrupts(void);
void install_signal_handler(void);
int eep;

#ifdef _WIN32
  /* nothing to be done */
#elif _WIN64
  /* nothing to be done */
#else
  CorpusList *current_corpus;
  CorpusList *corpuslist;
  CYCtype LastExpression;
  int exit_cqp;                   /**< 1 iff exit-command was issued while parsing */
  char *cqp_input_string;
  int cqp_input_string_position;
  int EvaluationIsRunning;
  int signal_handler_is_installed;
#endif
