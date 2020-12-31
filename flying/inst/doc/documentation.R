## ----X1, echo=FALSE, fig.height=5, fig.width=6---------------------------
#load("~/Documents/R projects/Flight project/flying/data/birds.rda")
data("birds", package = "flying")
prof.pow.ratio <- function(ws, wa) {
  # ws = wing span
  # wa = wing area
  X1 <- 8.4 / (ws ^ 2 / wa)

  return(X1)
}

boxplot(prof.pow.ratio(birds$Wing.span, birds$Wing.area),
        main = "Profile power ratio distribution among preset birds",
        col = "#58CCED")

## ----factor_table2, echo=FALSE, message=FALSE, warning=FALSE-------------
gen.table2 <- function() {
  x1Plusx2 <- c(0.00, 0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00,
                2.25, 2.50, 2.75, 3.00, 3.25, 3.50, 3.75, 4.00, 4.25,
                4.50, 4.75, 5.00)
  B <- c(1.360, 1.386, 1.452, 1.515, 1.574, 1.631, 1.684, 1.735, 1.784, 1.830,
         1.875, 1.918, 1.959, 1.999, 2.038, 2.075, 2.111, 2.146, 2.180, 2.213,
         2.246)
  C <- c(1.140, 1.458, 1.783, 2.115, 2.453, 2.795, 3.141, 3.490, 3.841, 4.195,
         4.550, 4.907, 5.266, 5.625, 5.986, 6.348, 6.711, 7.074, 7.438,
         7.803, 8.168)
  D <- c(1.000, 0.824, 0.706, 0.621, 0.556, 0.506, 0.465, 0.431, 0.402, 0.378,
         0.357, 0.339, 0.322, 0.308, 0.295, 0.283, 0.273, 0.263, 0.254, 0.246,
         0.238)
  table2 <- as.data.frame(cbind(x1Plusx2, B, C, D))

  return(table2)
}

table1 <- knitr::kable(gen.table2(),format = "latex", label = "Table 1", caption = "Factor Table")

kableExtra::kable_styling(table1, latex_options = "hold_position")

## ----data, echo = FALSE, message=FALSE, warning=FALSE, paged.print=FALSE----
data("birds", package = "flying")
table2 <- knitr::kable(birds, format = "latex", caption = "Preset birds from Flight program")
kableExtra::kable_styling(table2, latex_options = "hold_position")
#birds

## ------------------------------------------------------------------------
results_fixed_wing <- flying::flysim(data = birds) 

# range 
results_fixed_wing$range

## ------------------------------------------------------------------------
# constants used
results_fixed_wing$constants

## ------------------------------------------------------------------------
# constant speed
results_cmm_cs <- flying::migrate(data = birds, speed_control = "constant_speed")

results_cmm_cs$range

## ------------------------------------------------------------------------
# constant ratio between true air-speed and minimum power speed
results_cmm_cratio <- flying::migrate(data = birds, speed_control = "vvmp_constant")

results_cmm_cratio$range

## ----tables, echo=FALSE--------------------------------------------------
# results simulated from flight true air-speed constant
flight_cons_speed <- c(3090, 2670, 3702, 1115, 3821, 3519, 11418, 3383, 2922,
                            NA ,2248, 1867, 2031, 5186, 9880, 5308, 5498, 2606, 
                            6439, 5776, NA, NA, 2385, 2221, 2684, 5658, 5071, 6427)
flight_cons_ratio <- c(3007, 2586, 3566, 1076, 3679, 3417, 10647, 3267, 2789,
                       3016, 2146, 1798, 1975, 4946, 9499, 5120, 5378, 2541, 
                       6140, 5544, NA, 3374, 2275, 2119, 2571, 5381, 4871, 6153
                       )

results_table1 <- as.data.frame(results_cmm_cs$range)

results_table1$flight_results <- flight_cons_speed


results_table2 <-  as.data.frame(results_cmm_cratio$range)
results_table2$fight_results <- flight_cons_ratio  

colnames(results_table1) <-  c("package_cons_speed", "flight_cons_speed")

colnames(results_table2) <- c("package_cons_ratio","flight_cons_ratio")


## ----table1, echo=FALSE--------------------------------------------------
table3 <- knitr::kable(results_table1, format = "latex", label = "Table 2", caption = "Comparison table constant speed")
kableExtra::kable_styling(table3, latex_options = "hold_position")
#results_table1

## ----table2, echo=FALSE--------------------------------------------------
table4 <- knitr::kable(results_table2, format = "latex", label = "Table 3", caption = "Comparison table constant speed ratio")
kableExtra::kable_styling(table4, latex_options = "hold_position")
#results_table2

