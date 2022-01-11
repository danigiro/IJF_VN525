#' -----------------------------------------------------------------------------
#' VN525_csmix.R
#'
#' Averaging hts using the bts base forecasts from both models (SA and ETS)
#'
#' Base forecasts: ETS lev
#'
#' Reconcile forecasts (+ nn):
#'       - mix-wls (wls and m-wls)
#'       - mix-ols (ols and m-ols)
#'       - mix-shr (shr and m-shr)
#'
#' Input files: VN525_htsrec.RData VN525_htsrec_mean.RData
#'              VN525.RData
#' Output files: VN525_csmix.RData
#'
#' This code is written by Daniele Girolimetto
#' Department of Statistics, University of Padua (Italy)
#' -----------------------------------------------------------------------------
rm(list = ls(all = TRUE))
library(tidyverse)

load("./Reconciliation/ETSlev/VN525_htsrec.RData")
DFb <- DFhts %>% filter(`R-comb` %in% c("shr", "wls", "ols")) %>% 
  mutate(K = factor(K,c(12,6,4,3,2,1), ordered = TRUE))
rm(list=setdiff(ls(), c("DFb")))
load("./Reconciliation/ETSlev/VN525_htsrec_mean.RData")
DFm <- DFhts %>% filter(`R-comb` %in% c("m-shr", "m-wls", "m-ols")) %>% 
  mutate(K = factor(K,c(12,6,4,3,2,1), ordered = TRUE))
rm(list=setdiff(ls(), c("DFb","DFm")))
load("./VN525.RData")

DFb$Series <- factor(DFb$Series, colnames(VNdata), ordered = TRUE)
DFm$Series <- factor(DFm$Series, colnames(VNdata), ordered = TRUE)
test_length <- as.numeric(max(DFb$Replication))
DF <- NULL
ty <- unique(DFb$`R-comb`)
names(ty) <- paste0("m-", ty)
for(j in 1:test_length){
  for(i in 1:length(ty)){
    brecf <- DFb %>% 
      filter(Replication==j, `R-comb` == ty[i], FoReco == "direct") %>% 
      arrange(Series) %>%
      select(Series, Forecasts, `Forecast Horizon`, K) %>% 
      pivot_wider(names_from = Series, values_from = Forecasts) %>%
      arrange(K, `Forecast Horizon`) %>% 
      select(-`Forecast Horizon`, -K) %>% as.matrix() %>% t()
    
    mrecf <- DFm %>% 
      filter(Replication==j, `R-comb` == names(ty[i]), FoReco == "direct") %>% 
      arrange(Series) %>%
      select(Series, Forecasts, `Forecast Horizon`, K) %>% 
      pivot_wider(names_from = Series, values_from = Forecasts) %>%
      arrange(K, `Forecast Horizon`) %>% 
      select(-`Forecast Horizon`, -K) %>% as.matrix() %>% t()
    
    Fltr <- DFb %>% 
      filter(Replication==j, `R-comb` == ty[i], FoReco == "direct") %>% 
      dplyr::select(-"Forecasts", -"R-method", -"R-comb", -"FoReco", -"nn") %>% 
      arrange(Series, K, `Forecast Horizon`)
    
    Recon_PointF <- (brecf+mrecf)/2
    nn_val <- all(Recon_PointF>=0)
    Df1 <- cbind(Fltr, "Forecasts" = as.vector(t(Recon_PointF)),
                 "R-method" = "cs", "R-comb" = paste0("mix-", ty[i]), 
                 nn = nn_val, 
                 FoReco = "direct")
    Df1 <- Df1[names(DFb)]
    DF <- rbind(DF, Df1)
    
    if(!nn_val){
      brecfnn <- DFb  %>% 
        filter(Replication==j, `R-comb` == ty[i], nn == TRUE) %>% 
        arrange(Series) %>%
        select(Series, Forecasts, `Forecast Horizon`, K) %>% 
        pivot_wider(names_from = Series, values_from = Forecasts) %>%
        arrange(K, `Forecast Horizon`) %>% 
        select(-`Forecast Horizon`, -K) %>% as.matrix() %>% t()
      
      mrecfnn <- DFm %>% 
        filter(Replication==j, `R-comb` == names(ty[i]), nn == TRUE) %>% 
        arrange(Series) %>%
        select(Series, Forecasts, `Forecast Horizon`, K) %>% 
        pivot_wider(names_from = Series, values_from = Forecasts) %>%
        arrange(K, `Forecast Horizon`) %>% 
        select(-`Forecast Horizon`, -K) %>% as.matrix() %>% t()
      
      Recon_PointFnn <- (brecfnn+mrecfnn)/2
      
      Df1 <- cbind(Fltr, "Forecasts" = as.vector(t(Recon_PointFnn)),
                   "R-method" = "cs", "R-comb" = paste0("mix-", ty[i]), 
                   nn = all(Recon_PointFnn>=0), 
                   FoReco = "osqp")
      Df1 <- Df1[names(DFb)]
      DF <- rbind(DF, Df1)
    }
  }
  cat(j, " ")
}
DFmix <- DF
save(DFmix, 
     file="./Reconciliation/ETSlev/VN525_csmix.RData")
