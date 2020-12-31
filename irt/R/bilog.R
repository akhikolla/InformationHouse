
#' Prepare BILOG-MG data file
#'
#' @description Prepare \code{".dat"} data file for BILOG-MG.
#'
#' @param x Either a \\code{data.frame} or \code{matrix} object. Each row
#'   should represent a subject (examinee) and each column represent a
#'   variable (usually an item). The response values should be either 0, 1 or
#'   \code{NA}.
#' @param items A vector of column names or numbers of the \code{x} that
#'   represents the responses.
#' @param id_var The column name or number that contains individual subject IDs.
#'   If none is provided (i.e. \code{id_var = NULL}), the program will check
#'   whether the data provided has row names.
#' @param group_var The column name or number that contains group membership
#'   information if multi-group calibration is desired. Ideally, it grouping
#'   variable is represented by single digit integers. If other type of data
#'   provided, an integer value will automatically assigned to the variables.
#'   The default value is \code{NULL}, where no multi-group analysis is
#'   performed.
#' @param target_path The location where the BILOG-MG data file fill be saved.
#'   For example: \code{file.path(getwd(), "data", "mydata.dat")}.
#' @param create_np_key_file If \code{TRUE}, a file that contains
#'   the non-presented key will be created. This key is the same format as
#'   the corresponding response data file.
#' @param overwrite If TRUE and there is already a BILOG-MG data file in the
#'   target path with the same name, the file will be overwritten.
#'
#' @author Emre Gonulates
#'
#' @noRd
#'
#' @examples
#' \dontrun{
#' resp <- sim_resp(generate_ip(n = 15), rnorm(200), prop_missing = .2)
#' create_bilog_datafile(x = resp)
#' create_bilog_datafile(x = resp, items = paste0("Item-", 1:7))
#' }
#'
create_bilog_datafile <- function(
  x,
  items = NULL,
  id_var = NULL,
  group_var = NULL,
  target_path = file.path(getwd(), "bilog_data.dat"),
  create_np_key_file = FALSE,
  overwrite = FALSE) {

  if (!(inherits(x, c("matrix", "data.frame"))))
    stop("Invalid data file. Data file should inherit from either a matrix ",
         "or a data.frame", call. = FALSE)

  # Check path
  target_path <- normalizePath(target_path, mustWork = FALSE)
  target_dir <- dirname(target_path)
  if (!grepl("(.*)\\.dat$", target_path))
    stop("Invalid 'target_path' argument. Please provide a valid file name ",
         "with '.dat' extension as a target_path", call. = FALSE)
  # Make sure the directory exists
  if (!dir.exists(target_dir))
    dir.create(path = target_dir, recursive = TRUE)
  if (!dir.exists(target_dir))
    stop(paste0("The directory for BILOG-MG data file cannot be created at: \n",
                target_dir, "\nPlease create directory manually.",
                call. = FALSE))

  # If items are NULL, then use all of the items.
  if (is.null(items)) {
    # Check if the column names are NULL
    if (is.null(colnames(x)))
      colnames(x) <- paste0("Item-", 1:ncol(x))
    items <- colnames(x)
  }
  # Make sure items are valid
  if (!is.null(items) && (
    any(duplicated(items)) ||
    !(all(items %in% colnames(x)) || all(items %in% 1:ncol(x)) )
    ))
    stop("Invalid 'items' argument. The elements of 'items' argument should ",
         "be a character vector of column names of the 'x'. All elements ",
         "should be unique.", call. = FALSE)

  ### Prepare responses ###
  resp <- apply(x[, items], 2, as.character)
  # Check if the responses are valid:
  if (!all(sapply(unique(as.vector(resp)), function(x) is.na(x) ||
                  x %in% c("0", "1"))))
    stop("Invalid response data. The response values should be either 0, 1 ",
         "or missing (i.e. NA).", call. = FALSE)
  # Convert missing data to "."
  resp <- ifelse(is.na(resp), ".", resp)
  num_of_items <- ncol(resp)
  resp <- apply(resp, 1, paste0, collapse = "")

  ### ID ###
  # Create the id field, if it is NULL, automatically the row names of the data
  # will be the subject ids.
  id <- rownames(x)
  # Check whether id_var is valid
  if (!is.null(id_var) && (
    length(id_var) != 1 || !(id_var %in% 1:ncol(x) ||
                             (!is.null(colnames(x)) &&
                              id_var %in% colnames(x)))
    ))
    stop("Invalid 'id_var' argument. 'id_var' value should be a column ",
         "number or a column name.", call. = FALSE)
  if (!is.null(id_var)) id <- paste0(x[, id_var, drop = TRUE])
  if (is.null(id)) id <- paste0(1:nrow(x))
  # find the maximum number or characters in an ID
  num_id_char <- max(nchar(id))
  if (!num_id_char %in% 1:30)
    stop("Invalid IDs. Number of characters alloted to IDs cannot be more ",
         "than 30 characters.", call. = FALSE)
  id <- sprintf(paste0("%", num_id_char, "s"), id)

  ### Group ###
  # group: In BILOG-MG, "the group identifier has to be a single digit
  # (integer), starting with 1."
  group <- NULL
  groups <- NULL
  num_of_groups <- 0
  # Check whether group_var is valid
  if (!is.null(group_var) && (
    length(group_var) != 1 || !(group_var %in% 1:ncol(x) ||
                                (!is.null(colnames(x)) &&
                                 group_var %in% colnames(x)))
    ))
    stop("Invalid 'group_var' argument. 'group_var' value should be a column ",
         "number or a column name.", call. = FALSE)
  if (!is.null(group_var)) {
    group_names <- as.character(x[, group_var, drop = TRUE])
    unique_groups <- unique(group_names)
    if (any(is.na(unique_groups))) {
      stop("Grouping variable cannot have missing (NA) values.", call. = FALSE)
    }
    if (any(nchar(unique_groups) > 8)) {
      unique_groups <- substr(unique_groups, 1, 8)
      if (length(unique(unique_groups)) != length(unique_groups)) {
        stop("Invalid group values. Values that designate group membership ",
             "should be up to eight (8) characters and be unique.",
             call. = FALSE)
      } else
        message(paste0("Some of the group names are larger than eight ",
                       "characters. Group names are shortened to fit ",
                       "within eight characters:\n",
                       paste0("\"", unique_groups, "\"", collapse = ", ")))
    }
    num_of_groups <- length(unique_groups)
    # Make sure group variable is integers from 1:9
    if (!all(unique_groups %in% paste0(1:num_of_groups))) {
      group <- as.character(as.integer(factor(group_names,
                                              levels = unique_groups)))
    } else group <- group_names
    groups <- data.frame(name = unique_groups, code = unique(group))
  }

  ### Prepare the data file and formal statement ###
  # Add ID
  data_text <- paste0(id, " ")
  formal_statement <- paste0("(", num_id_char, "A1", ",", "1X")

  # Add group if there is
  if (!is.null(group)) {
    data_text <- paste0(data_text, group, " ")
    formal_statement <- paste0(formal_statement,  ",I1", ",1X")
  }

  # Create not-presented key file:
  if (create_np_key_file) {
    np_text <- paste0(
      # Add id
      paste0(rep(" ", num_id_char), collapse = ""), " ",
      # Add groups
      ifelse(num_of_groups == 0, "", "  "),
      # Add items
      paste0(rep(".", num_of_items), collapse = ""))

    np_key_file_path <- file.path(
      target_dir,
      paste0(gsub("\\.(.*)$", "", basename(target_path)), ".NFN"))
    if (overwrite || !file.exists(np_key_file_path))
      writeLines(text = np_text, con = np_key_file_path)
  } else
    np_key_file_path <- NULL

  # Add response data:
  data_text <- paste0(data_text, resp)
  formal_statement <- paste0(formal_statement,  ",", num_of_items, "A1)")

  if (overwrite || !file.exists(target_path)) writeLines(data_text, target_path)


  return(list(formal_statement = formal_statement,
              data_file_path = target_path,
              num_of_items = num_of_items,
              num_of_groups = num_of_groups,
              groups = groups,
              np_key_file_path = np_key_file_path,
              num_id_char = num_id_char))
}


#' Read the parameters of BILOG-MG Calibration
#'
#' @param par_file Path for a BILOG-MG parameter file with '.PAR' extension.
#' @param items The names of the items. This will be given assigned as the id's
#'   of the items. The default value is \code{NULL}, where the item id's
#'   assigned by BILOG-MG will be assigned as item parameters.
#' @param model The model of the items. The value is one of the
#'   following:
#'   One-parameter logistic model (\code{"1PL"}),
#'   Two-parameter logistic model (\code{"2PL"}),
#'   Three-parameter logistic model (\code{"3PL"}).
#'
#'   The default value is \code{"3PL"}.
#' @param D Scaling constant. The default value is \code{1}. If the item
#'   parameters needs to be converted to the commonly used normal scale where
#'   \code{D = 1.7} or \code{D = 1.702}, change this value accordingly.
#'
#' @noRd
#'
read_bilog_pars <- function(
  par_file,
  items = NULL,
  model = "3PL",
  D = 1
  ) {
  result <- list(ip = NULL,
                 failed_items = NULL)

  if (is.null(D)) D <- 1.7

  # Wait three seconds for the parameter file and
  counter <- 1
  while (!file.exists(par_file) && counter < 4) {
    message(paste0("Waiting for the parameter file... (", counter, ")\n"))
    Sys.sleep(1)
    counter <- counter + 1
  }
  if (!file.exists(par_file))
    stop(paste0("Parameter file cannot be found at: \n\"", par_file, "\"\n"),
         call. = FALSE)
  # Locate the parameter file and read it
  # In BILOG-MG's terms:
  # "a" = slope
  # "b" = threshold
  # "c" = lower asymptote
  pars <- utils::read.fwf(
    file = par_file,
    widths = c(8, 8, rep(10, 13), 4, 1, 1), skip = 4,  header = FALSE,
    col.names = c("id", "subtest_name", "intercept", "intercept_se",
                  "a", "a_se", "b", "b_se", "dispersion", "dispersion_se",
                  "c", "c_se", "drift", "drift_se", "unused", "item_location",
                  "answer_key", "dummy_values"
                                 ))

  if (model == "3PL") {
    ipdf <- pars[, c("id", "a", "b", "c")]
    ipdf_se <- pars[, c("a_se", "b_se", "c_se")]
  } else if (model == "2PL") {
    ipdf <- pars[, c("id", "a", "b")]
    ipdf_se <- pars[, c("a_se", "b_se")]
  } else if (model == "1PL") {
    ipdf <- pars[, c("id", "b"), drop = FALSE]
    ipdf_se <- pars[, "b_se", drop = FALSE]
  }
  if (!is.null(items)) ipdf$id <- items

  # Check whether all item parameter values are valid:
  if (any(sapply(ipdf[, -1], class) != "numeric") ||
      any(sapply(ipdf_se[, -1], class) != "numeric")) {
    invalid_items <- !apply(sapply(ipdf[, -1], function(x) grepl(
      "[-]?[0-9]+[.]?[0-9]*|[-]?[0-9]+[l]?|[-]?[0-9]+[.]?[0-9]*[ee][0-9]+", x)),
      1, all)
    invalid_items <- invalid_items | !apply(sapply(ipdf_se[, -1], function(x)
      grepl(
        "[-]?[0-9]+[.]?[0-9]*|[-]?[0-9]+[l]?|[-]?[0-9]+[.]?[0-9]*[ee][0-9]+",
        x)), 1, all)
    warning(paste0("\nEstimation failed for following item(s):\n\n",
                   paste0(utils::capture.output(print(ipdf[invalid_items, ])),
                          collapse = "\n"),
                   "\n\nThese items will be remove from the output."))
    result$failed_items <- ipdf[invalid_items, ]
    # Remove invalid items
    ipdf <- ipdf[!invalid_items, ]
    ipdf_se <- ipdf_se[!invalid_items, ]
  }

  # Convert the cleaned data frame to numeric.
  for (i in which(sapply(ipdf[, -1], class) != "numeric"))
    ipdf[, i+1] <- as.numeric(ipdf[, i+1])
  for (i in which(sapply(ipdf_se[, -1], class) != "numeric"))
    ipdf_se[, i+1] <- as.numeric(ipdf_se[, i+1])

  # if ("a" %in% colnames(ipdf)) ipdf$a <- ipdf$a / D
  ipdf$id <- gsub("^\\s+|\\s+$", "", ipdf$id)
  result$ip <- itempool(ipdf[, which(colnames(ipdf) %in% c("a", "b", "c"))],
                        id = ipdf$id, model = model, D = D,
                        se_parameters = ipdf_se)
  return(result)
}


#' This function reads the next CTT table from the row given
#' @param text BILOG-MG text grabbed
#' @param row The row number to start reading from 'text'
#'
#' @noRd
read_ctt_table <- function(text, row) {
  sub_text <- text[row:length(text)]
  sub_text <- sub_text[6:(which(grepl(paste0(
    "^ ", paste0(rep("-", 73), collapse = "")), sub_text))[2]-1)]

  pars <- utils::read.fwf(textConnection(sub_text),
                   widths = c(5, 11, 9, 10, 8, 9, 10, 9),
                   col.names = c("order", "name", "tried", "right", "pvalue",
                                 "logit", "pbis", "bis"),
                   header = FALSE,
                   colClasses = c("integer", "character", "numeric",
                                  "numeric", "numeric", "numeric", "numeric",
                                  "numeric")
                   )
  if (any(pars$pvalue > 1))
    pars$pvalue <- pars$pvalue/100
  pars$name <- gsub("^\\s+|\\s+$", "", pars$name)
  return(pars)
}


#' Read the CTT parameters of BILOG-MG Calibration
#'
#' @param ph1_file Path for a BILOG-MG ".PH1" file which holds CTT statistics.
#' @return A list of CTT parameters. If there are groups, then the CTT
#'   statistics for groups can be found in \code{$group$GROUP-NAME}.
#'   Overall statistics for the whole group is at \code{$overall}.
#'
#' @noRd
read_bilog_ctt <- function(ph1_file) {
  text <- readLines(ph1_file)
  result <- list(overall = NULL)
  # Check whether there is a group
  if (any(grepl("ITEM STATISTICS FOR MULTIPLE GROUPS", text))) {
    ### Multi-group ###
    result$group <- NULL
    # Read the CTT statistics for groups
    for (i in which(grepl("ITEM STATISTICS FOR GROUP:", text))) {
      group_name <- sub("( )*$", "", text[i])
      group_name <- sub("^(.*) ", "", group_name)
      result$group[[group_name]] <- read_ctt_table(text, i)
    }
    # Read the CTT statistics for the overall test
    i <- which(grepl("ITEM STATISTICS FOR MULTIPLE GROUPS", text))
    result$overall <- read_ctt_table(text, i)
  } else if (any(grepl("ITEM STATISTICS FOR SUBTEST", text))) {
    ### Single-Group ###
    i <- which(grepl("ITEM STATISTICS FOR SUBTEST", text))
    result$overall <- read_ctt_table(text, i)
  }
  if (is.null(result$group)) return(result$overall)
  return(result)
}



#' Read the group means of multigroup BILOG-MG Calibration
#'
#' @param ph3_file Path for a BILOG-MG ".PH3" file which holds group means.
#' @return Group means
#'
#' @noRd
read_bilog_group_means <- function(ph3_file, group_info) {
  text <- readLines(ph3_file, skipNul = TRUE)
  row <- which(grepl(pattern = "PRIOR MEANS AND STANDARD DEVIATIONS", text))
  result <- group_info
  if (length(row) == 1) {
    text <- text[row:length(text)]
    row <- which(grepl("^ ---------------------------$", text))[1:2]
    text <- text[(row[1]+1):(row[2]-1)]
    pars <- utils::read.fwf(textConnection(text), widths = c(3, 13, 10),
                            col.names = c("code", "mean", "sd"))
    if (is.null(result)) return(pars)
    result <- merge(x = result, pars, by = "code")
  }
  return(result)
}


#' Read the scale scores from BILOG-MG Calibration
#'
#' @param score_file Path for a BILOG-MG ".SCO" file which holds estimated
#'   scale scores.
#' @return A data frame consist of scores.
#'
#' @noRd
read_bilog_scores <- function(score_file) {
  text <- readLines(score_file)
  # Remove first two lines which does not contain any information
  text <- text[-c(1:2)]
  # Get the width of the first line that contains group number and id
  # widths <- c(3, nchar(text[1])-2, 6, 11, 3, 5, 10, 12, 12, 11, 10)
  n <- 1:(length(text)/2)
  # text <- paste(text[2*n-1], text[2*n])
  text <- text[2*n]

  scores <- utils::read.fwf(
    textConnection(text),
    widths = c(6, 11, 3, 5, 10, 12, 12, 11, 10),
    col.names = c("weight", "test", "tried", "right", "percent",
                  "ability", "se", "prob", "unknown1")
                            )
  return(scores[, c("tried", "right", "ability", "se", "prob")])
}


#' Run BILOG-MG in batch mode
#'
#' @description \code{est_bilog} runs BILOG-MG in batch mode or reads BILOG-MG
#'   output generated by BILOG-MG program. In the first case, this function
#'   requires BILOG-MG already installed on your computer under
#'   \code{bilog_exe_folder} directory.
#'
#'   In the latter case, where appropriate BILOG-MG files are present (i.e.
#'   \code{"<analysis_name>.PAR"}, \code{"<analysis_name>.PH1"},
#'   \code{"<analysis_name>.PH2"} and \code{"<analysis_name>.PH3"} files exist)
#'   and \code{overwrite = FALSE}, there is no need for BILOG-MG program. This
#'   function can read BILOG-MG output without BILOG-MG program.
#'
#' @param x Either a \code{data.frame} or \code{matrix} object. When the data
#'   is not necessary, i.e. user only wants to read the BILOG-MG output from
#'   the \code{target_dir}, then this can be set to \code{NULL}.
#' @param model The model of the items. The value is one of the
#'   following:
#'   \describe{
#'     \item{\code{"1PL"}}{One-parameter logistic model.}
#'     \item{\code{"2PL"}}{Two-parameter logistic model.}
#'     \item{\code{"3PL"}}{Three-parameter logistic model.}
#'     \item{\code{"CTT"}}{Return only Classical Test theory statistics such as
#'       p-values, point-biserial and  biserial correlations.}
#'   }
#'
#'   The default value is \code{"3PL"}.
#' @param target_dir The directory/folder where the BILOG-MG analysis and data
#'   files will be saved. The default value is the current working directory,
#'   i.e. \code{get_wd()}.
#' @param analysis_name A short file name that will be used for the data files
#'   created for the analysis.
#' @param items A vector of column names or numbers of the \code{x} that
#'   represents the responses. If, in the syntax file, no entry for item
#'   names are desired, then, simply write \code{items = "none"}.
#' @param id_var The column name or number that contains individual subject IDs.
#'   If none is provided (i.e. \code{id_var = NULL}), the program will check
#'   whether the data provided has row names.
#' @param group_var The column name or number that contains group membership
#'   information if multi-group calibration is desired. Ideally, it grouping
#'   variable is represented by single digit integers. If other type of data
#'   provided, an integer value will automatically assigned to the variables.
#'   The default value is \code{NULL}, where no multi-group analysis will be
#'   performed.
#' @param logistic A logical value. If \code{TRUE}, \code{LOGISTIC} keyword will
#'   be added to the BILOG-MG command file which means the calibration will
#'   assume the natural metric of the logistic response function in all
#'   calculations. If \code{FALSE}, the logit is multiplied by D = 1.7 to obtain
#'   the metric of the normal-ogive model. The default value is \code{TRUE}.
#' @param num_of_alternatives An integer specifying the maximum number of
#'   response alternatives in the raw data. \code{1/num_of_alternatives} is
#'   used by the analysis as automatic starting value for estimating the
#'   pseudo-guessing parameters.
#'
#'   The default value is \code{NULL}. In this case, for 3PL, 5 will be used
#'   and for 1PL and 2PL, 1000 will be used.
#'
#'   This value will be represented in BILOG-MG control file as:
#'   \code{NALT = num_of_alternatives}.
#'
#' @param criterion Convergence criterion for EM and Newton iterations. The
#'   default value is 0.01.
#' @param num_of_quadrature The number of quadrature points in MML estimation.
#'   The default value is 81. This value will be represented in BILOG-MG control
#'   file as: \code{NQPT = num_of_quadrature}. The BILOG-MG default value is
#'   20 if there are more than one group, 10 otherwise.
#' @param max_em_cycles An integer (0, 1, ...) representing the maximum number
#'   of EM cycles. This value will be represented in BILOG-MG control file as:
#'   \code{CYCLES = max_em_cycles}.
#'   The default value is 10.
#' @param newton An integer (0, 1, ...) representing the number of Gauss-Newton
#'   iterations following EM cycles. This value will be represented in BILOG-MG
#'   control file as: \code{NEWTON = newton}.
#' @param reference_group Represent which group's ability distribution will be
#'   set to mean = 0 and standard deviation = 1. For example, if the value is 1,
#'   then the group whose code is 1 will have ability distribution with mean 0
#'   and standard deviation 1. When groups are assumed to coming from a single
#'   population, set this value to 0.
#'
#'   The default value is \code{NULL}.
#'
#'   This value will be represented in BILOG-MG control file as:
#'   \code{REFERENCE = reference_group}.
#' @param scoring_options A string vector of keywords/options that will be added
#'   to the \code{SCORE} section in BILOG-MG syntax. Set the value of
#'   \code{scoring_options} to \code{NULL} if scoring of individual examinees is
#'   not necessary.
#'
#'   The default value is \code{c("METHOD=1", "NOPRINT")} where scale scores
#'   will be estimated using Maximum Likelihood estimation and the scoring
#'   process will not be printed to the R console (if
#'   \code{show_output_on_console = TRUE}).
#'
#'   The main option to be added to this vector is \code{"METHOD=n"}.
#'   Following options are available:
#'
#'   \describe{
#'     \item{"METHOD=1"}{Maximum Likelihood (ML)}
#'     \item{"METHOD=2"}{Expected a Posteriori (EAP)}
#'     \item{"METHOD=3"}{Maximum a Posteriori (MAP)}
#'   }
#'
#'   In addition to \code{"METHOD=n"} keyword, following keywords can be added:
#'
#'   \code{"NOPRINT"}: Suppresses the display of the scores on the R console.
#'
#'   \code{"FIT"}: likelihood ratio chi-square goodness-of-fit statistic for
#'     each response pattern will be computed.
#'
#'   \code{"NQPT=(list)"}, \code{"IDIST=n"}, \code{"PMN=(list)"},
#'   \code{"PSD=(list)"}, \code{"RSCTYPE=n"}, \code{"LOCATION=(list)"},
#'   \code{"SCALE=(list)"}, \code{"INFO=n"}, \code{"BIWEIGHT"},
#'   \code{"YCOMMON"}, \code{"POP"}, \code{"MOMENTS"},
#'   \code{"FILE"}, \code{"READF"}, \code{"REFERENCE=n"}, \code{"NFORMS=n"}
#'
#'   See BILOG-MG manual for more details about these keywords/options.
#' @param calib_options A string vector of keywords/options that will be added
#'   to the \code{CALIB} section in BILOG-MG syntax in addition to the keywords
#'   \code{NQPT}, \code{CYCLES}, \code{NEWTON}, \code{CRIT}, \code{REFERENCE}.
#'
#'   The default value is \code{c("NORMAL")}.
#'
#'   When \code{"NORMAL"} is included in \code{calib_options}, the prior
#'   distributions of ability in the population is assumed to have normal
#'   distribution.
#'
#'   When \code{"COMMON"} is included in \code{calib_options}, a common value
#'   for the lower asymptote for all items in the 3PL model will be estimated.
#'
#'   Following keywords/options can be added to \code{calib_options}:
#'
#'   \code{"PRINT=n"}, \code{"IDIST=n"}, \code{"PLOT=n"}, \code{"DIAGNOSIS=n"},
#'   \code{"REFERENCE=n"}, \code{"SELECT=(list)"}, \code{"RIDGE=(a,b,c)"},
#'   \code{"ACCEL=n"}, \code{"NSD=n"}, \code{"COMMON"}, \code{"EMPIRICAL"},
#'   \code{"NORMAL"}, \code{"FIXED"}, \code{"TPRIOR"}, \code{"SPRIOR"},
#'   \code{"GPRIOR"}, \code{"NOTPRIOR"}, \code{"NOSPRIOR"}, \code{"NOGPRIOR"},
#'   \code{"READPRIOR"}, \code{"NOFLOAT"}, \code{"FLOAT"}, \code{"NOADJUST"},
#'   \code{"GROUP-PLOT"}, \code{"RASCH"}, \code{"NFULL"}, \code{"CHI=(a,b)"}.
#'
#'   See BILOG-MG manual for more details about these keywords/options.
#'
#'   NOTE: Do not add any of the following keywords to \code{calib_options}
#'   since they will already be included:
#'
#'   \code{NQPT}, \code{CYCLES}, \code{NEWTON}, \code{CRIT}, \code{REFERENCE}
#'
#' @param overwrite If \code{TRUE} and there are already a BILOG-MG analysis
#'   files in the target path with the same name, these file will be
#'   overwritten.
#' @param show_output_on_console logical (not NA), indicates whether to capture
#'   the output of the command and show it on the R console. The default value
#'   is \code{TRUE}.
#' @param bilog_exe_folder The location of the \code{"blm1.exe"},
#'   \code{"blm2.exe"} and \code{"blm3.exe"} files. The default location is
#'   \code{file.path("C:/Program Files/BILOGMG")}.
#'
#' @return A list of following objects:
#'   \describe{
#'     \item{"ip"}{An \code{\link{Itempool-class}} object holding the item
#'       parameters. This element will not be created when
#'       \code{model = "CTT"}.}
#'     \item{"score"}{A data frame object that holds the number of item
#'       examinee has attempted (\code{tried}), the number of item examinee got
#'       right (\code{right}), the estimated scores of examinees
#'       (\code{ability}), the standard errors of ability estimates (\code{se}),
#'       and the probability of the response string (\code{prob}). This element
#'       will not be created when \code{model = "CTT"}.}
#'     \item{"ctt"}{The Classical Test Theory (CTT) stats such as p-value,
#'       biserial, point-biserial estimated by BILOG-MG. If there are groups,
#'       then the CTT statistics for groups can be found in
#'       \code{ctt$group$GROUP-NAME}. Overall statistics for the whole group is
#'       at \code{ctt$overall}.
#'       }
#'     \item{"failed_items"}{A data frame consist of items that cannot be
#'       estimated.}
#'     \item{"syntax"}{The syntax file.}
#'     \item{"converged"}{A logical value indicating whether a model has been
#'       converged or not. If the value is \code{TRUE}, model has been
#'       converged. This element will not be created when \code{model = "CTT"}.}
#'     \item{"neg_2_log_likelihood"}{-2 Log Likelihood value. This value is
#'       \code{NULL}, when model does not converge. This element will not be
#'       created when \code{model = "CTT"}.}
#'     \item{"input"}{A list object that stores the arguments that are passed
#'       to the function.}
#'   }
#'
#' @export
#'
#' @author Emre Gonulates
#'
#' @examples
#' \dontrun{
#' ### Example 1 ###
#' # Create responses to be used in BILOG-MG estimation
#' true_ip <- generate_ip(n = 30, model = "2PL")
#' resp <- sim_resp(true_ip, rnorm(4000))
#'
#' # The following line will run BILOG-MG, estimate 2PL model and put the
#' # analysis results under the target directory:
#' bilog_calib <- est_bilog(x = resp, model = "2PL",
#'                          target_dir = "C:/Temp/Analysis", overwrite = TRUE)
#' # Check whether the calibration converged
#' bilog_calib$converged
#'
#' # Get the estimated item pool
#' bilog_calib$ip
#'
#' # See the BILOG-MG syntax
#' cat(bilog_calib$syntax)
#'
#' # See the classical test theory statistics estimated by BILOG-MG:
#' bilog_calib$ctt
#'
#' # Get -2LogLikelihood for the model (mainly for model comparison purposes):
#' bilog_calib$neg_2_log_likelihood
#'
#'
#' ### Example 2 ###
#' # Get Expected-a-posteriori theta scores:
#' result <- est_bilog(x = resp, model = "2PL",
#'                     scoring_options = c("METHOD=2", "NOPRINT"),
#'                     target_dir = "C:/Temp/Analysis", overwrite = TRUE)
#'
#'
#' ### Example 3 ###
#' # Multi-group calibration
#' ip <- generate_ip(n = 35, model = "3PL", D = 1.7)
#' n_upper <- sample(1200:3000, 1)
#' n_lower <- sample(1900:2800, 1)
#' theta_upper <- rnorm(n_upper, 1.5, .25)
#' theta_lower <- rnorm(n_lower)
#' resp <- sim_resp(ip = ip, theta = c(theta_lower, theta_upper))
#' dt <- data.frame(level = c(rep("Lower", n_lower), rep("Upper", n_upper)), resp)
#'
#' mg_calib <- est_bilog(x = dt, model = "3PL",
#'                       group_var = "level",
#'                       reference_group = "Lower",
#'                       items = 2:ncol(dt), # Exclude the 'group' column
#'                       num_of_alternatives = 5,
#'                       # Use MAP ability estimation.
#'                       # "FIT": calculate GOF for response patterns
#'                       scoring_options = c("METHOD=3", "NOPRINT", "FIT"),
#'                       target_dir = "C:/Temp/Analysis", overwrite = TRUE,
#'                       show_output_on_console = FALSE)
#' # Estimated item pool
#' mg_calib$ip
#' # Print group means
#' mg_calib$group_info
#' # Check Convergence
#' mg_calib$converged
#' # Print estimated scores of first five examinees
#' head(mg_calib$score)
#'
#'
#' ### Example 4 ###
#' # When user wants to read BILOG-MG output saved in the directory "Analysis/"
#' # with file names "my_analysis", use the following syntax:
#' # (The following code does not require an installed BILOG-MG program.)
#' result <- est_bilog(target_dir = file.path("Analysis/"), model = "3PL",
#'                     analysis_name = "my_analysis", overwrite = FALSE)
#'
#' }
est_bilog <- function(
  x = NULL,
  model = "3PL",
  target_dir = getwd(),
  analysis_name = "bilog_calibration",
  items = NULL,
  id_var = NULL,
  group_var = NULL,
  logistic = TRUE,
  num_of_alternatives = NULL,
  criterion = 0.01,
  num_of_quadrature = 81,
  max_em_cycles = 100,
  newton = 20,
  reference_group = NULL,
  scoring_options = c("METHOD=1", "NOPRINT"),
  calib_options = c("NORMAL"),
  overwrite = FALSE,
  show_output_on_console = TRUE,
  bilog_exe_folder = file.path("C:/Program Files/BILOGMG")
  ) {

  # model = "3PL"
  # target_dir = getwd()
  # analysis_name = "bilog_calibration"
  # items = NULL
  # id_var = NULL
  # group_var = NULL
  # logistic = TRUE
  # num_of_alternatives = NULL
  # overwrite = TRUE
  # criterion = 0.01
  # num_of_quadrature = 31
  # max_em_cycles = 100
  # newton = 20
  # reference_group = NULL
  # scoring_options = c("METHOD=1", "NOPRINT")
  # calib_options = c("NORMAL")
  # bilog_exe_folder = file.path("C:/Program Files/BILOGMG")
  # show_output_on_console = TRUE

  result <- list(ip = NULL,
                 score = NULL,
                 ctt = NULL,
                 failed_items = NULL,
                 syntax = NULL
                 )

  ### ### ### ### ### GLOBAL ### ### ### ### ###
  if (!model %in% c("CTT", "1PL", "2PL", "3PL"))
    stop("Invalid 'model' argument. 'model' should be either '1PL', '2PL' ",
         "'3PL' or 'CTT'.", call. = FALSE)

  par_file <- file.path(target_dir, paste0(analysis_name, ".PAR"))
  ph1_file <- file.path(target_dir, paste0(analysis_name, ".PH1"))
  ph2_file <- file.path(target_dir, paste0(analysis_name, ".PH2"))
  ph3_file <- file.path(target_dir, paste0(analysis_name, ".PH3"))
  score_file <- file.path(target_dir, paste0(analysis_name, ".SCO"))

  D <- ifelse(logistic, 1, 1.7)

  target_path <- file.path(target_dir, paste0(analysis_name, ".blm"))

  if (overwrite || any(!file.exists(c(ph1_file, ph2_file, ph3_file,
                                      par_file)))) {
    # Check whether the filename is less than 128 characters:
    if (nchar(target_path) > 128)
      stop(paste0("\nThe following path for control file is ",
                  nchar(target_path), " character long, i.e. longer than the ",
                  "maximum length of 128 characters:\n\"", target_path,
                  "\"\nPlease choose a shorter 'analysis_name' or change the ",
                  "'target_dir'."), call. = FALSE)

    tab <- "       " # Determines the tab size, space added to lines following
                     # main analysis commands
    text <- paste0(analysis_name, "\n\n")

    # If overwrite is TRUE, delete important output files:
    if (overwrite) {
      extensions <- c("blm", "dat", "PH1", "PH2", "PH3", "COV", "NFN", "PAR",
                      "SCO")
      suppressWarnings(file.remove(
        file.path(target_dir, paste0(analysis_name, ".", extensions))))
    }
    # Create the data file:
    data_output <- create_bilog_datafile(
      x = x, items = items, id_var = id_var, group_var = group_var,
      target_path = file.path(target_dir, paste0(analysis_name, ".dat")),
      create_np_key_file = TRUE,
      overwrite = overwrite)

    # This functions wraps text within 80 characters. If it overflows, the
    # remaining text will be flowed the next line beginning column 1.
    wrap_text <- function(t, width = 80)
      gsub(paste0("(.{", width, "})"), "\\1\n", t)

    temp_text <- wrap_text(paste0(
      ">GLOBAL DFNAME='", data_output$data_file_path, "',\n"))
    # Add model
    temp_text <- paste0(
      temp_text, tab, "NPARM = ",
      gsub("(.*)PL", "\\1", ifelse(model == "CTT", "1PL", model)),
      ",\n", tab)
    # Add Logistic: when added the natural metric of the logistic response
    # function is assumed in all calculations.
    if (logistic) temp_text <- paste0(temp_text, "LOGISTIC,\n", tab)

    # Add SAVE statement
    # This indicates that the SAVE command will follow the GLOBAL command.
    temp_text <- paste0(temp_text, "SAVE;\n")

    text <- c(text, temp_text)

    ### ### ### ### ### SAVE ### ### ### ### ###
    temp_text <- paste0(">SAVE ")
    # Add PARM
    temp_text <- wrap_text(paste0(
      temp_text, "PARM='",
      normalizePath(file.path(target_dir, paste0(analysis_name, ".par")),
                    mustWork = FALSE), "',\n"))
    # Add SCORE
    temp_text <- paste0(
      temp_text, wrap_text(paste0(tab, "SCORE='",
      normalizePath(file.path(target_dir, paste0(analysis_name, ".sco")),
                    mustWork = FALSE), "',\n")))
    # Add COVARIANCE
    temp_text <- paste0(
      temp_text, wrap_text(paste0(tab, "COVARIANCE='",
      normalizePath(file.path(target_dir, paste0(analysis_name, ".cov")),
                    mustWork = FALSE), "';\n")))
    # # Add TSTAT
    # temp_text <- paste0(
    #   temp_text, "',\n", tab, "TSTAT = '",
    #   normalizePath(file.path(target_dir, paste0(analysis_name, ".tst")),
    #                 mustWork = FALSE))
    # # Add PDISTRIB
    # temp_text <- paste0(
    #   temp_text, "',\n", tab, "PDISTRIB = '",
    #   normalizePath(file.path(target_dir, paste0(analysis_name, ".dst")),
    #                 mustWork = FALSE))

    # temp_text <- paste0(temp_text, "';\n")

    text <- c(text, temp_text)

    ### ### ### ### ### LENGTH ### ### ### ### ###
    temp_text <- paste0(">LENGTH NITEMS = (", data_output$num_of_items, ");\n")

    text <- c(text, temp_text)

    ### ### ### ### ### INPUT ### ### ### ### ###
    # This section provides information that describes the raw data file.
    temp_text <- paste0(">INPUT NTOTAL = ", data_output$num_of_items)
    # Add number of alternatives
    if (!is.null(num_of_alternatives) && is.numeric(num_of_alternatives))
      temp_text <- paste0(temp_text, ",\n", tab, "NALT = ", num_of_alternatives)
    # Add NFNAME:
    temp_text <- paste0(
      temp_text, ",\n", wrap_text(paste0(tab, "NFNAME = '",
      normalizePath(data_output$np_key_file_path), "',\n")))
    # Add NIDCHAR
    temp_text <- paste0(temp_text, tab,
                        "NIDCHAR = ", data_output$num_id_char)

    # Add number of groups NGROUP:
    if (data_output$num_of_groups > 1)
      temp_text <- paste0(temp_text, ",\n", tab,
                          "NGROUP = ", data_output$num_of_groups)

    temp_text <- paste0(temp_text, ";\n")

    text <- c(text, temp_text)

    ### ### ### ### ### ITEMS ### ### ### ### ###
    # Setup the items to be printed.
    if (is.null(items)) {
      items <- colnames(x)
    } else if (length(items) == 1 && items == "none") {
      items <- NULL
    }
    if (is.null(items)) {
      text <- c(text, paste0(">ITEMS ;\n"))
    } else {
      temp_text <- wrap_text(paste0(
        ">ITEMS INAMES=(", paste0("'", items, "'", collapse = ","), ");\n"))
      text <- c(text, temp_text)
    }

    ### ### ### ### ### TEST ### ### ### ### ###
    # This section provides information that describes the raw data file.
    temp_text <- paste0(">TEST1 TNAME = 'TEST01'")
    temp_text <- paste0(temp_text, ",\n", tab, "INUMBER = (1(1)",
                        data_output$num_of_items, ")")
    temp_text <- paste0(temp_text, ";\n")

    text <- c(text, temp_text)

    ### ### ### ### ### GROUPS ### ### ### ### ###
    temp_text <- ""
    for (g in seq_len(data_output$num_of_groups)) {
      groups <- data_output$groups
      temp_text <- paste0(temp_text,
                          ">GROUP", g, " GNAME = '", groups$name[g], "', ",
                          "LENGTH = ", data_output$num_of_items, ",\n", tab,
                          "INUMBER = (1(1)", data_output$num_of_items, ")")
      temp_text <- paste0(temp_text, ";\n")

      # Check whether reference_group is NULL, if yes, assign a default
      # reference group:
      if (is.null(reference_group))
        reference_group <- groups$name[1]
    }
    text <- c(text, temp_text)

    ### Formal Statement ###
    text <- c(text, paste0(data_output$formal_statement, "\n"))

    ### ### ### ### ### CALIB ### ### ### ### ###
    temp_text <- paste0(">CALIB ")

    # Add NQPT
    temp_text <- paste0(temp_text, "NQPT = ", num_of_quadrature)

    # Add CYCLES
    temp_text <- paste0(temp_text, ",\n", tab, "CYCLES = ", max_em_cycles)

    # Add NEWTON
    temp_text <- paste0(temp_text, ",\n", tab, "NEWTON = ", newton)

    # Add CRIT
    temp_text <- paste0(temp_text, ",\n", tab, "CRIT = ", criterion)

    # Add calib_options
    if (length(calib_options) > 0)
      for (i in calib_options) temp_text <- paste0(temp_text, ",\n", tab, i)

    # Add REFERENCE
    if (!is.null(reference_group) && data_output$num_of_groups > 1) {
      if (length(reference_group) == 1 &&
          reference_group %in% data_output$groups$name) {
        temp_text <- paste0(
          temp_text, ",\n", tab, "REFERENCE = ",
          data_output$groups$code[data_output$groups$name == reference_group])
      } else {
        stop(paste0("'reference_group' argument should be one of the following ",
                    "group names:\n",
                    paste0("\"", data_output$groups$name, "\"", collapse = ",")),
             call. = FALSE)
      }
    }

    temp_text <- paste0(temp_text, ";\n")

    text <- c(text, temp_text)

    ### ### ### ### ### SCORE ### ### ### ### ###
    if (!is.null(scoring_options) && length(scoring_options) > 0) {
      temp_text <- paste0(">SCORE ")
      for (i in 1:length(scoring_options)) {
        temp_text <- paste0(temp_text, ifelse(i == 1, "", paste0(",\n", tab)),
                            scoring_options[i])
      }
      temp_text <- paste0(temp_text, ";\n")
      text <- c(text, temp_text)
    }

    ### ### ### ### ### RUN CALIBRATION ### ### ### ### ###
    result$syntax <- paste0(text, collapse = "\n")
    command <- paste0(
      "cd ", normalizePath(target_dir), " && \"",
      file.path(bilog_exe_folder, "blm1.exe"), "\" ", analysis_name ,
      " NUM=8900 CHAR=2000"
    )
    if (model != "CTT")
      command <- paste0(
        command, " && \"",
        file.path(bilog_exe_folder, "blm2.exe"), "\" ", analysis_name ,
        " NUM=8900 CHAR=2000 && \"",
        file.path(bilog_exe_folder, "blm3.exe"), "\" ", analysis_name,
        " NUM=8900 CHAR=2000"
      )

    if (overwrite || !file.exists(par_file)) {
      if (!all(file.exists(file.path(bilog_exe_folder,
                                     paste0("blm", 1:3, ".exe"))))) {
        stop(paste0("The required BILOG-MG executable files does not exist at ",
                    "\"", bilog_exe_folder, "\"."),
             call. = FALSE)
      }
      writeLines(text = text, con = target_path, sep = "")
      system("cmd.exe", input = command, wait = TRUE,
             show.output.on.console = show_output_on_console)
      if (show_output_on_console) cat("\n")
    }
  } else if (file.exists(target_path)) {
    result$syntax <- paste0(readLines(con = target_path, skipNul = TRUE),
                            collapse = "\n")
  }

  # Read PH2 file to determine whether the calibration converged, the method
  # below is not tested extensively and based on small number of experiences.
  if (file.exists(ph2_file)) {
    ph2_content <- readLines(con = ph2_file)
    result$converged <- ifelse(
      any(grepl(pattern = "BYTES OF NUMERICAL WORKSPACE USED",
                x = ph2_content)),
      TRUE, FALSE)
    # If the result converged get the -2 Log Likelihood
    if (result$converged) {
      result$neg_2_log_likelihood <- as.numeric(
        gsub("^ -2 LOG LIKELIHOOD = ", "",
             utils::tail(ph2_content[grepl(pattern = "^ -2 LOG LIKELIHOOD = ",
                                           ph2_content)],1)))
    } else result$neg_2_log_likelihood <- NULL
  }
  # Read item parameters
  if (model == "CTT") {
    result <- result[-which(names(result) %in%
                              c("ip", "converged", "failed_items", "score"))]
  } else {
    par_results <- read_bilog_pars(par_file = par_file, items = items,
                                   model = model, D = D)

    result$ip <- par_results$ip
    result$failed_items <- par_results$failed_items
  }

  result$ctt <- read_bilog_ctt(ph1_file)

  # Check if there is valid group means

  if (model != "CTT" && file.exists(ph3_file) && result$converged) {
    if (exists("data_output") && data_output$num_of_groups > 1) {
      result$group_info <- read_bilog_group_means(
        ph3_file, group_info = data_output$groups)
    # In case reading directly from the Bilog output and data is not available
    } else if (!is.null(result$ctt$group)) {
      result$group_info <- read_bilog_group_means(
        ph3_file, group_info = NULL)
    }
  }

  if (model != "CTT" && !is.null(scoring_options) &&
      (length(scoring_options) > 0) && result$converged) {
    result$score <- read_bilog_scores(score_file)
  }

  # set the input list
  result$input <- list(
    model = model,
    target_dir = target_dir,
    analysis_name = analysis_name,
    items = items,
    id_var = id_var,
    group_var = group_var,
    logistic = logistic,
    num_of_alternatives = num_of_alternatives,
    criterion = criterion,
    num_of_quadrature = num_of_quadrature,
    max_em_cycles = max_em_cycles,
    newton = newton,
    reference_group = reference_group,
    scoring_options = scoring_options,
    calib_options = calib_options,
    bilog_exe_folder = bilog_exe_folder
    )

  return(result)
}


