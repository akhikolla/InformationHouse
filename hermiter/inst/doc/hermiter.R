## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  fig.width=6, 
  fig.height=4,
  collapse = TRUE,
  comment = "#>"
)
set.seed(10)

## ----setup, message=FALSE, warning=FALSE--------------------------------------
library(hermiter)

## ----message=FALSE, warning=FALSE---------------------------------------------
library(magrittr)
library(ggplot2)
library(dplyr)
library(data.table)
library(DT)
library(mvtnorm)
library(patchwork)

## -----------------------------------------------------------------------------
hermite_est <- hermite_estimator(N=10, standardize=TRUE, 
                                 est_type = "univariate")

## -----------------------------------------------------------------------------
hermite_est <- hermite_estimator(N=10, standardize=TRUE, 
                                 est_type = "bivariate")

## -----------------------------------------------------------------------------
observations <- rlogis(n=1000)
hermite_est <- hermite_estimator(N=10, standardize=TRUE)
hermite_est <- update_batch(hermite_est,observations)

## -----------------------------------------------------------------------------
observations <- matrix(data = rnorm(2000),nrow = 1000, ncol=2)
hermite_est <- hermite_estimator(N=10, standardize=TRUE, 
                                 est_type = "bivariate")
hermite_est <- update_batch(hermite_est,observations)

## -----------------------------------------------------------------------------
observations <- rlogis(n=1000)
hermite_est <- hermite_estimator(N=10, standardize=TRUE)
hermite_est <- hermite_est %>% update_batch(observations)

## -----------------------------------------------------------------------------
observations <- matrix(data = rnorm(2000),nrow = 1000, ncol=2)
hermite_est <- hermite_estimator(N=10, standardize=TRUE, 
                                 est_type = "bivariate")
hermite_est <- hermite_est %>% update_batch(observations)

## -----------------------------------------------------------------------------
observations <- rlogis(n=1000)
hermite_est <- hermite_estimator(N=10, standardize=TRUE)
for (idx in seq_along(observations)) {
  hermite_est <- update_sequential(hermite_est,observations[idx])
}

## -----------------------------------------------------------------------------
observations <- matrix(data = rnorm(2000),nrow = 1000, ncol=2)
hermite_est <- hermite_estimator(N=10, standardize=TRUE, 
                                 est_type = "bivariate")
for (idx in seq_len(nrow(observations))) {
  hermite_est <- update_sequential(hermite_est,observations[idx,])
}

## -----------------------------------------------------------------------------
observations <- rlogis(n=1000)
hermite_est <- hermite_estimator(N=10, standardize=TRUE)
for (idx in seq_along(observations)) {
  hermite_est <- hermite_est %>% update_sequential(observations[idx])
}

## -----------------------------------------------------------------------------
observations <- matrix(data = rnorm(2000),nrow = 1000, ncol=2)
hermite_est <- hermite_estimator(N=10, standardize=TRUE, 
                                 est_type = "bivariate")
for (idx in seq_len(nrow(observations))) {
  hermite_est <- hermite_est %>% update_sequential(observations[idx,])
}

## -----------------------------------------------------------------------------
observations_1 <- rlogis(n=1000)
observations_2 <- rlogis(n=1000)
hermite_est_1 <- hermite_estimator(N=10, standardize=TRUE) %>% 
  update_batch(observations_1)
hermite_est_2 <- hermite_estimator(N=10, standardize=TRUE) %>% 
  update_batch(observations_2)
hermite_est_merged <- merge_hermite(list(hermite_est_1,hermite_est_2))

## -----------------------------------------------------------------------------
observations_1 <- matrix(data = rnorm(2000),nrow = 1000, ncol=2)
observations_2 <- matrix(data = rnorm(2000),nrow = 1000, ncol=2)
hermite_est_1 <- hermite_estimator(N=10, standardize=TRUE, 
                                 est_type = "bivariate") %>% 
  update_batch(observations_1)
hermite_est_2 <- hermite_estimator(N=10, standardize=TRUE, 
                                 est_type = "bivariate") %>% 
  update_batch(observations_2)
hermite_est_merged <- merge_hermite(list(hermite_est_1,hermite_est_2))

## -----------------------------------------------------------------------------
observations <- rlogis(n=2000)
hermite_est <- hermite_estimator(N=10, standardize=TRUE)
hermite_est <- update_batch(hermite_est, observations)

x <- seq(-15,15,0.1)
pdf_est <- dens(hermite_est,x)
cdf_est <- cum_prob(hermite_est,x)

p <- seq(0.05,1,0.05)
quantile_est <- quant(hermite_est,p)

## -----------------------------------------------------------------------------
observations <- rlogis(n=2000)
hermite_est <- hermite_estimator(N=10, standardize=TRUE)
hermite_est <- hermite_est %>% update_batch(observations)

x <- seq(-15,15,0.1)
pdf_est <- hermite_est %>% dens(x)
cdf_est <- hermite_est %>% cum_prob(x)

p <- seq(0.05,0.95,0.05)
quantile_est <- hermite_est %>% quant(p)

## -----------------------------------------------------------------------------
actual_pdf <- dlogis(x)
actual_cdf <- plogis(x)
df_pdf_cdf <- data.frame(x,pdf_est,cdf_est,actual_pdf,actual_cdf)

actual_quantiles <- qlogis(p)
df_quant <- data.frame(p,quantile_est,actual_quantiles)

## -----------------------------------------------------------------------------
ggplot(df_pdf_cdf,aes(x=x)) + geom_line(aes(y=pdf_est, colour="Estimated")) +
  geom_line(aes(y=actual_pdf, colour="Actual")) +
  scale_colour_manual("", 
                      breaks = c("Estimated", "Actual"),
                      values = c("blue", "black")) + ylab("Probability Density")

## -----------------------------------------------------------------------------
ggplot(df_pdf_cdf,aes(x=x)) + geom_line(aes(y=cdf_est, colour="Estimated")) +
  geom_line(aes(y=actual_cdf, colour="Actual")) +
  scale_colour_manual("", 
                      breaks = c("Estimated", "Actual"),
                      values = c("blue", "black")) +
  ylab("Cumulative Probability")

## -----------------------------------------------------------------------------
ggplot(df_quant,aes(x=actual_quantiles)) + geom_point(aes(y=quantile_est),
                                                      color="blue") +
  geom_abline(slope=1,intercept = 0) +xlab("Theoretical Quantiles") +
  ylab("Estimated Quantiles")

## -----------------------------------------------------------------------------
# Prepare bivariate normal data
sig_x <- 1
sig_y <- 1
num_obs <- 4000
rho <- 0.5
observations_mat <- mvtnorm::rmvnorm(n=num_obs,mean=rep(0,2),sigma = matrix(c(sig_x^2,rho*sig_x*sig_y,rho*sig_x*sig_y,sig_y^2), 
                                                      nrow=2,ncol=2, 
                                                      byrow = TRUE))

hermite_est <- hermite_estimator(N = 20, standardize = TRUE, 
                                 est_type = "bivariate") 
hermite_est <-  update_batch(hermite_est,observations_mat)

vals <- seq(-5,5,by=0.25)
x_grid <- as.matrix(expand.grid(X=vals, Y=vals))
pdf_est <- dens(hermite_est,x_grid)
cdf_est <- cum_prob(hermite_est,x_grid)

## -----------------------------------------------------------------------------
sig_x <- 1
sig_y <- 1
num_obs <- 4000
rho <- 0.5
observations_mat <- mvtnorm::rmvnorm(n=num_obs,mean=rep(0,2),sigma = matrix(c(sig_x^2,rho*sig_x*sig_y,rho*sig_x*sig_y,sig_y^2), nrow=2, ncol=2, 
                                                                byrow = TRUE))

hermite_est <- hermite_estimator(N = 20, standardize = TRUE, 
                                 est_type = "bivariate") 
hermite_est <-  hermite_est %>% update_batch(observations_mat)

vals <- seq(-5,5,by=0.25)
x_grid <- as.matrix(expand.grid(X=vals, Y=vals))
pdf_est <- hermite_est %>% dens(x_grid, clipped = TRUE)
cdf_est <- hermite_est %>% cum_prob(x_grid, clipped = TRUE)

spear_est <- hermite_est %>% spearmans()

## -----------------------------------------------------------------------------
actual_pdf <-mvtnorm::dmvnorm(x_grid,mean=rep(0,2),sigma = matrix(c(sig_x^2,rho*sig_x*sig_y,rho*sig_x*sig_y,sig_y^2), nrow=2,ncol=2, 
                                                                byrow = TRUE))

actual_cdf <- rep(NA,nrow(x_grid))
for (row_idx in seq_len(nrow(x_grid))) {
  actual_cdf[row_idx] <-  mvtnorm::pmvnorm(lower = c(-Inf,-Inf),upper=as.numeric(x_grid[row_idx,]),mean=rep(0,2),sigma = matrix(c(sig_x^2,rho*sig_x*sig_y,rho*sig_x*sig_y,sig_y^2), nrow=2,ncol=2, 
                                                                                                                                byrow = TRUE))
}

actual_spearmans <- cor(observations_mat,method = "spearman")[1,2]

df_pdf_cdf <- data.frame(x_grid,pdf_est,cdf_est,actual_pdf,actual_cdf)

## -----------------------------------------------------------------------------
p1 <- ggplot(df_pdf_cdf) + geom_tile(aes(X, Y, fill= actual_pdf)) +
  scale_fill_gradient2(low="blue", mid="cyan", high="purple",
                       midpoint=.1,    
                       breaks=seq(0,.2,by=.05), 
                       limits=c(0,.2))  

p2 <- ggplot(df_pdf_cdf) + geom_tile(aes(X, Y, fill= pdf_est)) +
  scale_fill_gradient2(low="blue", mid="cyan", high="purple",
                       midpoint=.1,
                       breaks=seq(0,.2,by=.05), 
                       limits=c(0,.2))

p1+ ggtitle("Actual PDF")+ theme(legend.title = element_blank()) + p2 +
  ggtitle("Estimated PDF") +theme(legend.title = element_blank()) +
  plot_layout(guides = 'collect')

## -----------------------------------------------------------------------------
p1 <- ggplot(df_pdf_cdf) + geom_tile(aes(X, Y, fill= actual_cdf)) +
  scale_fill_gradient2(low="blue", mid="cyan", high="purple", 
                       midpoint=0.5,    
                       breaks=seq(0,1,by=.2), 
                       limits=c(0,1)) 

p2 <- ggplot(df_pdf_cdf) + geom_tile(aes(X, Y, fill= cdf_est)) +
  scale_fill_gradient2(low="blue", mid="cyan", high="purple",
                       midpoint=0.5,
                       breaks=seq(0,1,by=.2), #breaks in the scale bar
                       limits=c(0,1))

p1+ ggtitle("Actual CDF") + theme(legend.title = element_blank()) + p2 +
  ggtitle("Estimated CDF") + theme(legend.title = element_blank())+
  plot_layout(guides = 'collect')

## -----------------------------------------------------------------------------
# Prepare Test Data
test_data <- data.frame()
for (i in seq_len(5)) {
  exponential_data <- rexp(n=1000)
  logistic_data <- rlogis(n=1000)
  logn_data <- rlnorm(n=1000)
  test_data <- rbind(test_data,data.frame(dist_name=rep("exp",
                length(exponential_data)),idx=i,observations=exponential_data))
  test_data <- rbind(test_data,data.frame(dist_name=rep("logis",
                      length(logistic_data)),idx=i,observations=logistic_data))
  test_data <- rbind(test_data,data.frame(dist_name=rep("lnorm",
                              length(logn_data)),idx=i,observations=logn_data))
}
setDT(test_data)

## -----------------------------------------------------------------------------
# Group observations by distribution and idx and create Hermite estimators
estimates <- test_data[,.(herm_est = list(hermite_estimator(N=10,
             standardize = TRUE) %>% update_batch(observations))),
             by=.(dist_name,idx)]
estimates

## -----------------------------------------------------------------------------
# Group observations by distribution and merge Hermite estimators
merged_estimates <- estimates[,.(herm_comb = list(merge_hermite(herm_est))),
                                by=.(dist_name)]
merged_estimates

## -----------------------------------------------------------------------------
# Estimate probability densities, cumulative probabilities and quantiles
dens_vals <- merged_estimates[,.(dens_est = list(dens(herm_comb[[1]],
                                            c(0.5,1,1.5,2)))),by=.(dist_name)]
cum_prob_vals <- merged_estimates[,.(cum_prob_est = list(cum_prob(herm_comb[[1]],c(0.5,1,1.5,2)))),by=.(dist_name)]
quantile_vals <- merged_estimates[,.(quantile_est = list(quant(herm_comb[[1]],c(0.25,0.5,0.75)))),by=.(dist_name)]

## -----------------------------------------------------------------------------
# Prepare Test Data
test_data <- data.frame()
for (i in seq_len(5)) {
  exponential_data <- rexp(n=1000)
  logistic_data <- rlogis(n=1000)
  logn_data <- rlnorm(n=1000)
  test_data <- rbind(test_data,data.frame(dist_name=rep("exp",
                length(exponential_data)),idx=i,observations=exponential_data))
  test_data <- rbind(test_data,data.frame(dist_name=rep("logis",
                      length(logistic_data)),idx=i,observations=logistic_data))
  test_data <- rbind(test_data,data.frame(dist_name=rep("lnorm",
                              length(logn_data)),idx=i,observations=logn_data))
}

## -----------------------------------------------------------------------------
# Group observations by distribution and idx and create Hermite estimators
estimates <- test_data %>% group_by(dist_name,idx) %>% summarise(herm_est = list(hermite_estimator(N=10,standardize = TRUE) %>% update_batch(observations)))
estimates

## -----------------------------------------------------------------------------
# Group observations by distribution and merge Hermite estimators
merged_estimates <- estimates %>% group_by(dist_name) %>% summarise(herm_comb
                                              = list(merge_hermite(herm_est)))
merged_estimates

## -----------------------------------------------------------------------------
# Estimate probability densities, cumulative probabilities and quantiles
dens_vals <- merged_estimates %>%
  rowwise() %>% mutate(dens_est = list(dens(herm_comb,c(0.5,1,1.5,2))))
cum_prob_vals <- merged_estimates %>%
  rowwise() %>% mutate(cum_prob_est = list(cum_prob(herm_comb,c(0.5,1,1.5,2))))
quantile_vals <- merged_estimates %>%
  rowwise() %>% mutate(quantile_est = list(quant(herm_comb,c(0.25,0.5,0.75))))

## -----------------------------------------------------------------------------
# Compute Mean Absolute Error
dens_vals <- dens_vals %>%
  rowwise() %>% mutate(dens_actual = list(do.call(paste0("d",dist_name),
                list(c(0.5,1,1.5,2))))) %>% mutate(mean_abs_error_density =
                                              mean(abs(dens_est-dens_actual)))
cum_prob_vals <- cum_prob_vals %>%
  rowwise() %>% mutate(cum_prob_actual = list(do.call(paste0("p",dist_name),
                list(c(0.5,1,1.5,2)))))%>% mutate(mean_abs_error_cum_prob = mean(abs(cum_prob_est-cum_prob_actual)))
quantile_vals <- quantile_vals %>%
  rowwise() %>% mutate(quantile_actual= list(do.call(paste0("q",dist_name),
            list(c(0.25,0.5,0.75)))))%>% mutate(mean_abs_error_quantiles = mean(abs(quantile_est-quantile_actual)))
mean_abs_error_summary <- data.frame(dist_name=dens_vals$dist_name, mean_abs_error_density=dens_vals$mean_abs_error_density, mean_abs_error_cum_prob=cum_prob_vals$mean_abs_error_cum_prob,
              mean_abs_error_quantiles=quantile_vals$mean_abs_error_quantiles)

## -----------------------------------------------------------------------------
datatable(mean_abs_error_summary) %>% formatRound(columns =c("mean_abs_error_density","mean_abs_error_cum_prob",
                                        "mean_abs_error_quantiles"),digits = 3)

## ----eval=FALSE---------------------------------------------------------------
#  # Not Run. Copy and paste into app.R and run.
#  library(shiny)
#  library(hermiter)
#  library(ggplot2)
#  library(magrittr)
#  
#  ui <- fluidPage(
#      titlePanel("Streaming Statistics Analysis Example: Exponential
#                 i.i.d. stream"),
#      sidebarLayout(
#          sidebarPanel(
#              sliderInput("percentile", "Percentile:",
#                          min = 0.01, max = 0.99,
#                          value = 0.5, step = 0.01)
#          ),
#          mainPanel(
#             plotOutput("plot"),
#             textOutput("quantile_text")
#          )
#      )
#  )
#  
#  server <- function(input, output) {
#      values <- reactiveValues(hermite_est =
#                                   hermite_estimator(N = 10, standardize = TRUE))
#      x <- seq(-15, 15, 0.1)
#      # Note that the stub below could be replaced with code that reads streaming
#      # data from various sources, Kafka etc.
#      read_stream_stub_micro_batch <- reactive({
#          invalidateLater(1000)
#          new_observation <- rexp(10)
#          return(new_observation)
#      })
#      updated_cdf_calc <- reactive({
#          micro_batch <- read_stream_stub_micro_batch()
#          for (idx in seq_along(micro_batch)) {
#              values[["hermite_est"]] <- isolate(values[["hermite_est"]]) %>%
#                  update_sequential(micro_batch[idx])
#          }
#          cdf_est <- isolate(values[["hermite_est"]]) %>%
#              cum_prob(x, clipped = TRUE)
#          df_cdf <- data.frame(x, cdf_est)
#          return(df_cdf)
#      })
#      updated_quantile_calc <- reactive({
#          values[["hermite_est"]]  %>% quant(input$percentile)
#      })
#      output$plot <- renderPlot({
#          ggplot(updated_cdf_calc(), aes(x = x)) + geom_line(aes(y = cdf_est)) +
#              ylab("Cumulative Probability")
#      }
#      )
#      output$quantile_text <- renderText({
#          return(paste(input$percentile * 100, "th Percentile:",
#                       round(updated_quantile_calc(), 2)))
#      })
#  }
#  shinyApp(ui = ui, server = server)

## -----------------------------------------------------------------------------
# Prepare Test Data
num_obs <-2000
test <- rchisq(num_obs,5)
test <- c(test,rlogis(num_obs))
test <- c(test,rnorm(num_obs))

## -----------------------------------------------------------------------------
# Calculate theoretical pdf, cdf and quantile values for comparison
x <- seq(-15,15,by=0.1)
actual_pdf_lognorm <- dchisq(x,5)
actual_pdf_logis <- dlogis(x)
actual_pdf_norm <- dnorm(x)
actual_cdf_lognorm <- pchisq(x,5)
actual_cdf_logis <- plogis(x)
actual_cdf_norm <- pnorm(x)
p <- seq(0.05,0.95,by=0.05)
actual_quantiles_lognorm <- qchisq(p,5)
actual_quantiles_logis <- qlogis(p)
actual_quantiles_norm <- qnorm(p)

## -----------------------------------------------------------------------------
# Construct Hermite Estimator 
h_est <- hermite_estimator(N=20,standardize = TRUE,exp_weight_lambda = 0.005)

## -----------------------------------------------------------------------------
# Loop through test data and update h_est to simulate observations arriving 
# sequentially
count <- 1
res <- data.frame()
res_q <- data.frame()
for (idx in seq_along(test)) {
  h_est <- h_est %>% update_sequential(test[idx])
  if (idx %% 100 == 0){
    if (floor(idx/num_obs)==0){
      actual_cdf_vals <- actual_cdf_lognorm
      actual_pdf_vals <-actual_pdf_lognorm
      actual_quantile_vals <- actual_quantiles_lognorm
    }
    if (floor(idx/num_obs)==1){
      actual_cdf_vals <- actual_cdf_logis
      actual_pdf_vals <-actual_pdf_logis
      actual_quantile_vals <- actual_quantiles_logis
    }
    if (floor(idx/num_obs)==2){
      actual_cdf_vals <- actual_cdf_norm
      actual_pdf_vals <- actual_pdf_norm
      actual_quantile_vals <- actual_quantiles_norm
    }
    idx_vals <- rep(count,length(x))
    cdf_est_vals <- h_est %>% cum_prob(x, clipped=TRUE)
    pdf_est_vals <- h_est %>% dens(x, clipped=TRUE)
    quantile_est_vals <- h_est %>% quant(p)
    res <- rbind(res,data.frame(idx_vals,x,cdf_est_vals,actual_cdf_vals,
                                pdf_est_vals,actual_pdf_vals))
    res_q <- rbind(res_q,data.frame(idx_vals=rep(count,length(p)),p,
                                    quantile_est_vals,actual_quantile_vals))
    count <- count +1
  }
}
res <- res %>% mutate(idx_vals=idx_vals*100)
res_q <- res_q %>% mutate(idx_vals=idx_vals*100)

## ----eval=FALSE---------------------------------------------------------------
#  # Visualize Results for PDF (Not run, requires gganimate, gifski and transformr
#  # packages)
#  p <- ggplot(res,aes(x=x)) + geom_line(aes(y=pdf_est_vals, colour="Estimated")) + geom_line(aes(y=actual_pdf_vals, colour="Actual")) +
#    scale_colour_manual("",
#                        breaks = c("Estimated", "Actual"),
#                        values = c("blue", "black")) + ylab("Probability Density") +transition_states(idx_vals,transition_length = 2,state_length = 1) +
#    ggtitle('Observation index {closest_state}')
#  anim_save("pdf.gif",p)

## ----eval=FALSE---------------------------------------------------------------
#  # Visualize Results for CDF (Not run, requires gganimate, gifski and transformr
#  # packages)
#  p <- ggplot(res,aes(x=x)) + geom_line(aes(y=cdf_est_vals, colour="Estimated")) + geom_line(aes(y=actual_cdf_vals, colour="Actual")) +
#    scale_colour_manual("",
#                        breaks = c("Estimated", "Actual"),
#                        values = c("blue", "black")) +
#    ylab("Cumulative Probability") +
#    transition_states(idx_vals, transition_length = 2,state_length = 1) +
#    ggtitle('Observation index {closest_state}')
#  anim_save("cdf.gif", p)

## ----eval=FALSE---------------------------------------------------------------
#  # Visualize Results for Quantiles (Not run, requires gganimate, gifski and
#  # transformr packages)
#  p <- ggplot(res_q,aes(x=actual_quantile_vals)) +
#    geom_point(aes(y=quantile_est_vals), color="blue") +
#    geom_abline(slope=1,intercept = 0) +xlab("Theoretical Quantiles") +
#    ylab("Estimated Quantiles") +
#    transition_states(idx_vals,transition_length = 2, state_length = 1) +
#    ggtitle('Observation index {closest_state}')
#  anim_save("quant.gif",p)

## ----eval=FALSE---------------------------------------------------------------
#  # Not Run. Copy and paste into app.R and run.
#  library(shiny)
#  library(hermiter)
#  library(ggplot2)
#  library(magrittr)
#  
#  ui <- fluidPage(
#    titlePanel("Bivariate Streaming Statistics Analysis Example"),
#    sidebarLayout(
#      sidebarPanel(
#        sliderInput("spearmans", "True Spearman's Correlation:",
#                    min = -0.9, max = 0.9,
#                    value = 0, step = 0.1)
#      ),
#      mainPanel(
#        plotOutput("plot"),
#        textOutput("spearman_text")
#      )
#    )
#  )
#  
#  server <- function(input, output) {
#    values <- reactiveValues(hermite_est =
#                               hermite_estimator(N = 10, standardize = TRUE,
#                                                 exp_weight_lambda = 0.01,
#                                                 est_type="bivariate"))
#    # Note that the stub below could be replaced with code that reads streaming
#    # data from various sources, Kafka etc.
#    read_stream_stub_micro_batch <- reactive({
#      invalidateLater(1000)
#      sig_x <- 1
#      sig_y <- 1
#      num_obs <- 100
#      rho <- 2 *sin(pi/6 * input$spearmans)
#      observations_mat <- mvtnorm::rmvnorm(n=num_obs,mean=rep(0,2),sigma = matrix(c(sig_x^2,rho*sig_x*sig_y,rho*sig_x*sig_y,sig_y^2), nrow=2,ncol=2,
#                                                                  byrow = TRUE))
#      return(observations_mat)
#    })
#    updated_spear_calc <- reactive({
#      micro_batch <- read_stream_stub_micro_batch()
#      for (idx in seq_len(nrow(micro_batch))) {
#        values[["hermite_est"]] <- isolate(values[["hermite_est"]]) %>%
#          update_sequential(micro_batch[idx,])
#      }
#      spear_est <- isolate(values[["hermite_est"]]) %>%
#        spearmans(clipped = TRUE)
#      return(spear_est)
#    })
#    output$plot <- renderPlot({
#      vals <- seq(-5,5,by=0.25)
#      x_grid <- as.matrix(expand.grid(X=vals, Y=vals))
#      rho <- 2 *sin(pi/6 * input$spearmans)
#      actual_pdf <-mvtnorm::dmvnorm(x_grid,mean=rep(0,2),sigma = matrix(c(sig_x^2,rho*sig_x*sig_y,rho*sig_x*sig_y,sig_y^2), nrow=2,ncol=2,
#                                                                  byrow = TRUE))
#      df_pdf <- data.frame(x_grid,actual_pdf)
#      p1 <- ggplot(df_pdf) + geom_tile(aes(X, Y, fill= actual_pdf)) +
#        scale_fill_gradient2(low="blue", mid="cyan", high="purple",
#                             midpoint=.2,
#                             breaks=seq(0,.4,by=.1),
#                             limits=c(0,.4)) +ggtitle(paste("True Bivariate
#                      Normal Density with matched Spearman's correlation")) +
#         theme(legend.title = element_blank())
#      p1
#    }
#    )
#    output$spearman_text <- renderText({
#      return(paste("Spearman's Correlation Estimate from Hermite Estimator:",
#                   round(updated_spear_calc(), 1)))
#    })
#  }
#  shinyApp(ui = ui, server = server)

## ----eval=FALSE---------------------------------------------------------------
#  citation("hermiter")

