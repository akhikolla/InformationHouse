## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----message=FALSE------------------------------------------------------------
library(magrittr)
library(dplyr)
library(tidyr)
library(purrr)
library(fable)
library(ProbReco)
set.seed(1983)
sim_hierarchy

## -----------------------------------------------------------------------------
#Length of window
N<-500 

#Make data windows

data_windows<-purrr::map(1:(nrow(sim_hierarchy)-N+1),
                  function(i){return(sim_hierarchy[i:(i+N-1),])})


## -----------------------------------------------------------------------------
data_windows[[1]]%>%head(3)
data_windows[[1]]%>%tail(3)
data_windows[[2]]%>%head(3)
data_windows[[2]]%>%tail(3)
data_windows[[500]]%>%head(3)
data_windows[[500]]%>%tail(3)

## -----------------------------------------------------------------------------

forecast_window <- function(data_w){
  data_w%>%
    tidyr::pivot_longer(cols = -Time,
               names_to = 'var',
               values_to = 'value')%>%
    as_tsibble(index = 'Time',key='var')%>%
    model(arma11=ARIMA(value~1+pdq(1,0,1)+PDQ(0,0,0)))%>%
    forecast(h=1)%>%
    dplyr::arrange(match(var,c("Tot","A","B","AA","AB","BA","BB")))->f
    mean<-map_dbl(f$.distribution,use_series,mean)
    sd<-map_dbl(f$.distribution,use_series,sd)
 return(list(mean=mean,sd=sd))
}




## -----------------------------------------------------------------------------
#Number of windows for training
W<-10
all_fc<-purrr::map(data_windows[1:W],forecast_window)

## -----------------------------------------------------------------------------
S<-matrix(c(1,1,1,1,
            1,1,0,0,
            0,0,1,1,
            1,0,0,0,
            0,1,0,0,
            0,0,1,0,
            0,0,0,1),7,4,byrow = TRUE)

## -----------------------------------------------------------------------------

obs_i<-function(i){
  sim_hierarchy%>%
  dplyr::filter(Time==i)%>%
  tidyr::pivot_longer(-Time,names_to = 'var')%>%
  dplyr::arrange(match(var,c("Tot","A","B","AA","AB","BA","BB")))%>%
  dplyr::pull(value)
}

all_y<-purrr::map((N+1):(N+W),obs_i)


## -----------------------------------------------------------------------------
make_genfunc<-function(input){
  f<-function(){
    fc_mean<-input$mean
    fc_sd<-input$sd
    out<-matrix(rnorm(350,mean=fc_mean,sd=fc_sd),7,50)
    return(out)
  }
  return(f)
}


## -----------------------------------------------------------------------------
all_prob<-purrr::map(all_fc,make_genfunc)

## ---- eval=FALSE--------------------------------------------------------------
#    G_bu<-as.vector(rbind(matrix(0,4,4),diag(rep(1,4))))
#    es_bu<-total_score(all_y,all_prob,S,G_bu)

## ----eval=FALSE---------------------------------------------------------------
#    opt<-scoreopt(all_y,all_prob,S)

