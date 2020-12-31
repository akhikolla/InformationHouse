context("Testing bilog")

#' Check to see whether the BILOG-MG is installed in this machine. If no, skip
#' the tests on this file.
#'
#' @noRd
skip_if_bilog_exe_not_found <- function() {
  bilog_exe_paths <- c(file.path("C:/Program Files/BILOGMG"))
  if (!any(dir.exists(bilog_exe_paths)))
    skip("Bilog is not installed on this compter.")
}



# library(testthat)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%% est_bilog %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


test_that("est_bilog", {
  skip_if_bilog_exe_not_found()

  target_dir <- "C:/Temp/testthat-bilog"

  n_item <- sample(30:36, 1)
  n_theta <- sample(1900:2200, 1)
  true_ip <- generate_ip(n = n_item, model = "2PL")
  resp <- sim_resp(ip = true_ip, theta = rnorm(n_theta), prop_missing = .1)

  # The following line will run BILOG-MG, estimate 2PL model and put the
  # analysis results under the target directory:
  c1 <- est_bilog(x = resp,
                  model = "2PL",
                  target_dir = target_dir, overwrite = TRUE,
                  show_output_on_console = FALSE)

  expect_true(c1$converged)
  expect_is(c1$ip, "Itempool")
  expect_equal(nrow(c1$score), n_theta)
  expect_equal(nrow(c1$ctt), n_item)

  # -------------------------------------------------------------------------- #
  # model = "CTT"
  c2 <- est_bilog(x = resp,
                  model = "CTT",
                  target_dir = target_dir,
                  overwrite = TRUE,
                  show_output_on_console = FALSE)
  expect_null(c2$ip)
  expect_null(c2$converged)
  expect_null(c2$score)
  expect_equal(nrow(c2$ctt), n_item)
  expect_equal(c2$input$model, "CTT")
  expect_equal(c2$input$target_dir, target_dir)
})


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%% create_bilog_datafile %%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# test_that("create_bilog_datafile", {

  # skip_if_bilog_exe_not_found()
  # expect_equal("Uncomment below", "Uncomment below")
  # skip("Skip 'create_bilog_datafile'.")
  # ip <- generate_ip(n = 30, model = "2PL")
  # resp <- sim_resp(ip, rnorm(4000), prop_missing = .2)
  #
  # # Use all of the data in the analysis
  # data_output <- create_bilog_datafile(x = resp)
  #
  # # Only use a subset of items in the analysis.
  # create_bilog_datafile(data = resp, items = paste0("Item-", 1:7))
  #
  # # Response data without column names
  # resp_matrix <- matrix(sample(0:1, 500, T), ncol=10)
  # create_bilog_datafile(data = resp_matrix)
  #
  # # 'items' can be column numbers
  # create_bilog_datafile(data = resp_matrix, items = 1:4)
  #
  #
  # ### Invalid items argument
  # # Duplicated items
  # expect_error(
  #   create_bilog_datafile(data = resp, items = c("Item0", "Item-2","Item-2")),
  #   "Invalid 'items' argument.")
  #
  # # Invalid column names given in 'items' argument
  # expect_error(create_bilog_datafile(data = resp, items = paste0("I-", 1:7)),
  #              "Invalid 'items' argument.")
  #
  # # cat(est_bilog(data = data), sep = "")
  #
  # D <- 1.7
  # ip <- generate_ip(n = 30, model = "2PL", D = D)
  # resp <- sim_resp(ip, rnorm(4000))
  # ip_est <- est_bilog(x = resp, model = "2PL", D = D,
  #                     target_dir = "C:/Users/EGonulates/Desktop/NF2")
  # cbind(a_est = ip_est$a, a = ip$a, a_new = ip_est$a/D)
  # cbind(b_est = ip_est$b, b = ip$b)
  # plot(x = ip$a, y = ip_est$a)
  # abline(0,1)
  #
  #
  # # Multi Group
  # require(irt)
  # D = 1.7
  # n_item <- sample(40:50, 1)
  # ip <- generate_ip(n = n_item, D= D)
  # n_upper <- sample(1200:3000, 1)
  # n_lower <- sample(1900:2800, 1)
  # theta_upper <- rnorm(n_upper, 1.5, .25)
  # theta_lower <- rnorm(n_lower)
  # resp <- sim_resp(ip = ip, theta = c(theta_lower, theta_upper))
  # dt <- data.frame(group = c(rep("Lower", n_lower), rep("Upper", n_upper)),
  #                  resp)
  # head(dt)
  #
  # x = dt
  # model = "3PL"
  # items = colnames(x)[-1]
  # group_var = "group"
  # target_dir = "C:/Users/EGonulates/Desktop/NF2"
  # scoring_method = 3
  # reference_group = "Upper"
  # overwrite = TRUE
  # analysis_name = "bilog_calibration"
  # id_var = NULL
  # num_of_alternatives = NULL
  # overwrite = TRUE
  # criterion = 0.01
  # common_guessing = FALSE
  # num_of_quadrature = 81
  # max_em_cycles = 100
  # newton = 20
  # normal_ability = TRUE
  # bilog_exe_folder = file.path("C:/Program Files/BILOGMG")
  #
  #
  # data_output <- create_bilog_datafile(
  #   x = x, items = items, id_var = id_var,
  #   group_var = group_var,
  #   target_path = file.path(target_dir, paste0(analysis_name, ".dat")),
  #   create_np_key_file = TRUE,
  #   overwrite = TRUE)
  #
  #
  # ip_est <- est_bilog(x = x, model = model, items = items,
  #                     group_var = group_var,
  #                     D = D,
  #                     target_dir = target_dir,
  #                     scoring_method = scoring_method,
  #                     reference_group = reference_group,
  #                     max_em_cycles = 500,
  #                     num_of_quadrature = 81,
  #                     overwrite = overwrite)
  # cbind(a_est = ip_est$ip$a, a = ip$a, a_new = ip_est$ip$a * D)
  #
  # ip_est$ctt$overall
  # ip_est$group_info
# })
