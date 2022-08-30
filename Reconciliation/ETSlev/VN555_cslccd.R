#' -----------------------------------------------------------------------------
#' VN555_cslccd.R
#'
#' Creating an RData file of conditional levels reconciled forecasts for VN555 
#' (cross-sectional framework, balance version) - diagonal covariance matrix 
#' 
#' Base forecasts: ETS lev
#' 
#' Reconcile forecasts (+ nn):
#'       - HCCCexod (Hollyman + CCC + exogenous + diagonal cov)
#'       - bCCCexod (basef + CCC + exogenous + diagonal cov)
#'       - HCCCendod (Hollyman + CCC + endogenous + diagonal cov)
#'       - bCCCendod (basef + CCC + endogenous + diagonal cov)
#'
#' Input files: VN555_base.RData
#' Output files: VN525_cslccd.RData
#'
#' This code is written by Daniele Girolimetto
#' Department of Statistics, University of Padua (Italy)
#' -----------------------------------------------------------------------------
rm(list = ls(all = TRUE))
library(FoReco)
library(tidyverse)
library(progress)

load("./BaseForecasts/ETSlev/555/VN555_base.RData")
DF <- NULL
C <- C555
nl <- c(1,7,27,76,4,28,108)
test_length <- as.numeric(max(DFbase$Replication))

mvdf$Series <- factor(mvdf$Series, name555, ordered = TRUE)
DFbase$Series <- factor(DFbase$Series, name555, ordered = TRUE)
resmat_all <- Residuals_ETS
K_vec <- unique(DFbase$K)
time_cslev <- array(NA, dim = c(test_length, length(K_vec), 4),
                    dimnames = list(NULL, K_vec, c("HCCCexod", "bCCCexod", 
                                                   "HCCCendod", "bCCCendod")))

pb <- progress_bar$new(
  format = "lccrec Hollyman: [:bar] :percent eta: :eta", 
  show_after = 1,
  total = test_length, clear = FALSE, width= 80)

for (j in 1:test_length) { #test_length
  for(i in 1:length(K_vec)){
    resmat <- resmat_all[[j]][[i]]
    mse <- diag(crossprod(resmat)/NROW(resmat))
    
    basef <- DFbase %>% filter(Replication==j, K==K_vec[i]) %>% 
      arrange(Series) %>%
      select(Series, Forecasts, `Forecast Horizon`) %>% 
      pivot_wider(names_from = Series, values_from = Forecasts) %>%
      arrange(`Forecast Horizon`) %>% 
      select(-`Forecast Horizon`) %>% as.matrix()
    
    bnaive <- mvdf %>% filter(Replication==j, K==K_vec[i]) %>% 
      arrange(Series) %>%
      select(Series, Mean, `Forecast Horizon`) %>% 
      pivot_wider(names_from = Series, values_from = Mean) %>%
      arrange(`Forecast Horizon`) %>% 
      select(-`Forecast Horizon`) %>% as.matrix()
    
    fixv <- mvdf %>% filter(Replication==j, K==K_vec[i]) %>% 
      arrange(Series) %>%
      select(Series, MeanVar, `Forecast Horizon`) %>% 
      pivot_wider(names_from = Series, values_from = MeanVar) %>%
      arrange(`Forecast Horizon`) %>% 
      select(-`Forecast Horizon`) %>% as.matrix()
    fixv <- fixv[1,]
    
    seasv <- mvdf %>% filter(Replication==j, K==K_vec[i]) %>% 
      arrange(Series) %>%
      select(Series, SeasVar, `Forecast Horizon`) %>% 
      pivot_wider(names_from = Series, values_from = SeasVar) %>%
      arrange(`Forecast Horizon`) %>% 
      select(-`Forecast Horizon`) %>% as.matrix()
    
    # Allocate the matrix for the reconcile forecasts with origin j 
    Fltr <- DFbase %>% filter(Replication==j, K==K_vec[i]) %>% 
      arrange(Series) %>%
      dplyr::select(-"Forecasts", -"R-method", -"R-comb") %>% 
      arrange(Series, `Forecast Horizon`)
    
    ## Reconciliation ----
    # HCCCexod (Hollyman + CCC + exogenous + diagonal cov) ----
    Start <- Sys.time()
    objH <- suppressMessages(lccrec(basef = basef, C = C, nl = nl, 
                                    bnaive = bnaive[,-c(1:NROW(C))],
                                    CCC = TRUE, weights = fixv[-c(1:NROW(C))]))
    End <- Sys.time()
    time_cslev[j,i,1] <- as.numeric(difftime(End, Start, units = "secs"))
    
    Recon_PointF <- objH$recf
    
    Df1 <- cbind(Fltr, "Forecasts" = as.vector(Recon_PointF),
                 "R-method" = "cslcc", "R-comb" = "HCCCexod", 
                 nn = all(Recon_PointF>=0), 
                 FoReco = "direct")
    Df1 <- Df1[c(names(DFbase), "nn", "FoReco")]
    DF <- rbind(DF, Df1)
    
    # bCCCexod (basef + CCC + exogenous + diagonal cov) ----
    Start <- Sys.time()
    objb <- suppressMessages(lccrec(basef = basef, C = C, nl = nl, 
                                    CCC = TRUE, weights = fixv[-c(1:NROW(C))]))
    End <- Sys.time()
    time_cslev[j,i,2] <- as.numeric(difftime(End, Start, units = "secs"))
    
    Recon_PointF <- objb$recf
    
    Df1 <- cbind(Fltr, "Forecasts" = as.vector(Recon_PointF),
                 "R-method" = "cslcc", "R-comb" = "bCCCexod", 
                 nn = all(Recon_PointF>=0), 
                 FoReco = "direct")
    Df1 <- Df1[c(names(DFbase), "nn", "FoReco")]
    DF <- rbind(DF, Df1)
    
    # HCCCendod (Hollyman + CCC + endogenous + diagonal cov) ----
    Start <- Sys.time()
    objHe <- suppressMessages(lccrec(basef = basef, C = C, nl = nl, 
                                     bnaive = bnaive[,-c(1:NROW(C))], const = "endo",
                                     CCC = TRUE, weights = fixv))
    End <- Sys.time()
    time_cslev[j,i,3] <- as.numeric(difftime(End, Start, units = "secs"))
    
    Recon_PointF <- objHe$recf
    
    Df1 <- cbind(Fltr, "Forecasts" = as.vector(Recon_PointF),
                 "R-method" = "cslcc", "R-comb" = "HCCCendod", 
                 nn = all(Recon_PointF>=0), 
                 FoReco = "direct")
    Df1 <- Df1[c(names(DFbase), "nn", "FoReco")]
    DF <- rbind(DF, Df1)
    
    # bCCCendod (basef + CCC + endogenous + diagonal cov) ----
    Start <- Sys.time()
    objbe <- suppressMessages(lccrec(basef = basef, C = C, nl = nl, const = "endo",
                                     CCC = TRUE, weights = fixv))
    End <- Sys.time()
    time_cslev[j,i,4] <- as.numeric(difftime(End, Start, units = "secs"))
    
    Recon_PointF <- objbe$recf
    
    Df1 <- cbind(Fltr, "Forecasts" = as.vector(Recon_PointF),
                 "R-method" = "cslcc", "R-comb" = "bCCCendod", 
                 nn = all(Recon_PointF>=0), 
                 FoReco = "direct")
    Df1 <- Df1[c(names(DFbase), "nn", "FoReco")]
    DF <- rbind(DF, Df1)
  }
  pb$tick()
}

nn_data <- DF %>% filter(nn=="FALSE")%>%
  select(`R-comb`,K, Replication) %>% unique()

nn_data$time <-  NA
save.image(file = "./Reconciliation/ETSlev/cslcc_555_preosqp.RData")
problem <- c()
for(row_nn in 1:NROW(nn_data)){
  if(row_nn == 1) 
    cat(" A.O. | F.O. | Comb  | Iteration", sep="")
  
  j <- nn_data$Replication[row_nn]
  i <- which(K_vec %in% nn_data$K[row_nn])
  ty <- nn_data$`R-comb`[row_nn]
  
  resmat <- resmat_all[[j]][[i]]
  mse <- diag(crossprod(resmat)/NROW(resmat))
  
  basef <- DFbase %>% filter(Replication==j, K==K_vec[i]) %>% 
    arrange(Series) %>%
    select(Series, Forecasts, `Forecast Horizon`) %>% 
    pivot_wider(names_from = Series, values_from = Forecasts) %>%
    arrange(`Forecast Horizon`) %>% 
    select(-`Forecast Horizon`) %>% as.matrix()
  
  bnaive <- mvdf %>% filter(Replication==j, K==K_vec[i]) %>% 
    arrange(Series) %>%
    select(Series, Mean, `Forecast Horizon`) %>% 
    pivot_wider(names_from = Series, values_from = Mean) %>%
    arrange(`Forecast Horizon`) %>% 
    select(-`Forecast Horizon`) %>% as.matrix()
  
  fixv <- mvdf %>% filter(Replication==j, K==K_vec[i]) %>% 
    arrange(Series) %>%
    select(Series, MeanVar, `Forecast Horizon`) %>% 
    pivot_wider(names_from = Series, values_from = MeanVar) %>%
    arrange(`Forecast Horizon`) %>% 
    select(-`Forecast Horizon`) %>% as.matrix()
  fixv <- fixv[1,]
  
  # Allocate the matrix for the reconcile forecasts with origin j
  Fltr <- DFbase %>% filter(Replication==j, K==K_vec[i]) %>% 
    arrange(Series) %>%
    dplyr::select(-"Forecasts", -"R-method", -"R-comb") %>% 
    arrange(Series, `Forecast Horizon`)
  
  ## NN reconciliation ----
  if(ty == "HCCCexod"){
    # HCCCexod (Hollyman + CCC + exogenous + diagonal cov + nn) ----
    Start <- Sys.time()
    obj <- suppressMessages(lccrec(basef = basef, C = C, nl = nl, nn= TRUE,
                                   bnaive = bnaive[,-c(1:NROW(C))],
                                   CCC = TRUE, weights = fixv[-c(1:NROW(C))],
                                   settings = osqpSettings(verbose = FALSE,
                                                           max_iter = 100000L,
                                                           check_termination = 5,
                                                           eps_abs = 1e-5, 
                                                           eps_rel = 1e-10,
                                                           eps_dual_inf = 1e-07,
                                                           #delta = 1e-5,
                                                           polish_refine_iter = 10000, 
                                                           polish = TRUE)))
    End <- Sys.time()
    if(any(sapply(obj$info, function(x) any(x[,5]!=1 | x[,4]>1e-4)))){
      cat("\n")
      problem = c(problem, row_nn)
      print(obj$info)
      cat("\n")
    }
  }else if(ty == "bCCCexod"){
    # bCCCexod (basef + CCC + exogenous + diagonal cov + nn) ----
    Start <- Sys.time()
    obj <- suppressMessages(lccrec(basef = basef, C = C, nl = nl, nn= TRUE,
                                   CCC = TRUE, weights = fixv[-c(1:NROW(C))],
                                   settings = osqpSettings(verbose = FALSE,
                                                           max_iter = 100000L,
                                                           check_termination = 5,
                                                           eps_abs = 1e-5, 
                                                           eps_rel = 1e-10,
                                                           #eps_prim_inf = 1e-16,
                                                           eps_dual_inf = 1e-07,
                                                           #delta = 1e-5,
                                                           polish_refine_iter = 10000, 
                                                           polish = TRUE)))
    End <- Sys.time()
    if(any(sapply(obj$info, function(x) any(x[,5]!=1 | x[,4]>1e-4)))){
      cat("\n")
      problem = c(problem, row_nn)
      print(obj$info)
      cat("\n")
    }
  }else if(ty == "HCCCendod"){
    # HCCCendod (Hollyman + CCC + endogenous + diagonal cov + nn) ----
    Start <- Sys.time()
    obj <- suppressMessages(lccrec(basef = basef, C = C, nl = nl, nn= TRUE,
                                   bnaive = bnaive[,-c(1:NROW(C))], const = "endo",
                                   CCC = TRUE, weights = fixv,
                                   settings = osqpSettings(verbose = FALSE,
                                                           max_iter = 100000L,
                                                           check_termination = 5,
                                                           eps_abs = 1e-5, 
                                                           eps_rel = 1e-10,
                                                           eps_dual_inf = 1e-07,
                                                           #delta = 1e-5,
                                                           polish_refine_iter = 10000, 
                                                           polish = TRUE)))
    End <- Sys.time()
    if(any(sapply(obj$info, function(x) any(x[,5]!=1 | x[,4]>1e-4)))){
      cat("\n")
      problem = c(problem, row_nn)
      print(obj$info)
      cat("\n")
    }
  }else if(ty == "bCCCendod"){
    # bCCCendod (basef + CCC + endogenous + diagonal cov) ----
    Start <- Sys.time()
    obj <- suppressMessages(lccrec(basef = basef, C = C, nl = nl, nn= TRUE,
                                   CCC = TRUE, weights = fixv, const = "endo",
                                   settings = osqpSettings(verbose = FALSE,
                                                           max_iter = 100000L,
                                                           check_termination = 5,
                                                           eps_abs = 1e-5, 
                                                           eps_rel = 1e-10,
                                                           eps_dual_inf = 1e-07,
                                                           #delta = 1e-5,
                                                           polish_refine_iter = 10000, 
                                                           polish = TRUE)))
    End <- Sys.time()
    if(any(sapply(obj$info, function(x) any(x[,5]!=1 | x[,4]>1e-4)))){
      cat("\n")
      problem = c(problem, row_nn)
      print(obj$info)
      cat("\n")
    }
  }
  
  Recon_PointF <- obj$recf
  
  nn_data$time[row_nn] <- as.numeric(difftime(End, Start, units = "secs"))
  
  # Add rows to DF
  Df1 <- cbind(Fltr, "Forecasts" = as.vector(Recon_PointF),
               "R-method" = "cslcc", "R-comb" = ty, 
               nn = all(Recon_PointF>=0), 
               FoReco = "osqp")
  Df1 <- Df1[c(names(DFbase), "nn", "FoReco")]
  DF <- rbind(DF, Df1)
  cat("\n  ", formatC(i, width=2, flag="0"), "  |  ", 
      formatC(j, width=2, flag="0"), "  | ", 
      format(ty, width=5), " |  ", row_nn, "/", NROW(nn_data), sep="")
}
nn_data_cslev <- nn_data
DFcslev <- DF
save(DFcslev, time_cslev, nn_data_cslev, 
     file="./Reconciliation/ETSlev/VN555_cslccd.RData")