#' -----------------------------------------------------------------------------
#' VN555_CCCmix.R
#'
#' Averaging CCC using the bts base forecasts from both models (SA and ETS)
#'
#' Base forecasts: ETS lev
#'
#' Reconcile forecasts (+ nn):
#'       - mix-CCC (mCCC and bCCC)
#'
#' Input files: VN555_cslccd.RData VN555_cslccd_mean.RData VN555_base.RData
#' Output files: VN555_CCCmix.RData
#'
#' This code is written by Daniele Girolimetto
#' Department of Statistics, University of Padua (Italy)
#' -----------------------------------------------------------------------------
rm(list = ls(all = TRUE))
library(tidyverse)

load("./Reconciliation/ETSlev/VN555_cslccd.RData")
DFb <- DFcslev %>% filter(`R-comb` == "bCCCexod") %>% 
  mutate(K = factor(K,c(12,6,4,3,2,1), ordered = TRUE))
rm(list=setdiff(ls(), c("DFb")))
load("./Reconciliation/ETSlev/VN555_cslccd_mean.RData")
DFm <- DFcslcc %>% filter(`R-comb` == "mCCCexod") %>% 
  mutate(K = factor(K,c(12,6,4,3,2,1), ordered = TRUE))
rm(list=setdiff(ls(), c("DFb","DFm")))
load("./BaseForecasts/ETSlev/555/VN555_base.RData")
rm(list=setdiff(ls(), c("DFb","DFm", "name555", "DFbase")))

DFb$Series <- factor(DFb$Series, name555, ordered = TRUE)
DFm$Series <- factor(DFm$Series, name555, ordered = TRUE)
DFbase$Series <- factor(DFbase$Series, name555, ordered = TRUE)

test_length <- as.numeric(max(DFb$Replication))
DF <- NULL

for(j in 1:test_length){
  brecf <- DFb %>% 
    filter(Replication==j, FoReco == "direct") %>% 
    arrange(Series) %>%
    select(Series, Forecasts, `Forecast Horizon`, K) %>% 
    pivot_wider(names_from = Series, values_from = Forecasts) %>%
    arrange(K, `Forecast Horizon`) %>% 
    select(-`Forecast Horizon`, -K) %>% as.matrix() %>% t()
  
  mrecf <- DFm %>% 
    filter(Replication==j, FoReco == "direct") %>% 
    arrange(Series) %>%
    select(Series, Forecasts, `Forecast Horizon`, K) %>% 
    pivot_wider(names_from = Series, values_from = Forecasts) %>%
    arrange(K, `Forecast Horizon`) %>% 
    select(-`Forecast Horizon`, -K) %>% as.matrix() %>% t()
  
  Fltr <- DFbase  %>% 
    mutate(K = factor(K,c(12,6,4,3,2,1), ordered = TRUE)) %>% 
    filter(Replication==j) %>% 
    dplyr::select(-"Forecasts", -"R-method", -"R-comb") %>% 
    arrange(Series, K, `Forecast Horizon`)
  
  Recon_PointF <- (brecf+mrecf)/2
  nn_val <- all(Recon_PointF>=0)
  Df1 <- cbind(Fltr, "Forecasts" = as.vector(t(Recon_PointF)),
               "R-method" = "cs", "R-comb" = "mix-CCC", 
               nn = nn_val, 
               FoReco = "direct")
  Df1 <- Df1[names(DFb)]
  DF <- rbind(DF, Df1)
  
  if(!nn_val){
    brecfnn <- DFb  %>% 
      filter(Replication==j, nn == TRUE) %>% 
      arrange(Series) %>%
      select(Series, Forecasts, `Forecast Horizon`, K) %>% 
      pivot_wider(names_from = Series, values_from = Forecasts) %>%
      arrange(K, `Forecast Horizon`) %>% 
      select(-`Forecast Horizon`, -K) %>% as.matrix() %>% t()
    
    mrecfnn <- DFm %>% 
      filter(Replication==j, nn == TRUE) %>% 
      arrange(Series) %>%
      select(Series, Forecasts, `Forecast Horizon`, K) %>% 
      pivot_wider(names_from = Series, values_from = Forecasts) %>%
      arrange(K, `Forecast Horizon`) %>% 
      select(-`Forecast Horizon`, -K) %>% as.matrix() %>% t()
    
    Recon_PointFnn <- (brecfnn+mrecfnn)/2
    
    Df1 <- cbind(Fltr, "Forecasts" = as.vector(t(Recon_PointFnn)),
                 "R-method" = "cs", "R-comb" = "mix-CCC", 
                 nn = all(Recon_PointFnn>=0), 
                 FoReco = "osqp")
    Df1 <- Df1[names(DFb)]
    DF <- rbind(DF, Df1)
  }
  cat(j, " ")
}
DFmix <- DF
save(DFmix, 
     file="./Reconciliation/ETSlev/VN555_CCCmix.RData")
