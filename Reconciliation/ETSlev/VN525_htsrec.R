#' -----------------------------------------------------------------------------
#' VN525_htsrec.R
#'
#' Creating an RData file of cross sectional reconciled forecasts for VN525
#' 
#' Base forecasts: ETS lev
#' 
#' Reconcile forecasts (+ nn):
#'       - bu (Bottom-up)
#'       - ols (Identity)
#'       - struc (Structural variances)
#'       - wls (Series variances) 
#'       - shr (Shrunk covariance matrix)
#'
#' Input files: DFbase.RData, Residuals.RData, VN525.RData
#' Output files: VN525_htsrec.RData
#'
#' This code is written by Daniele Girolimetto
#' Department of Statistics, University of Padua (Italy)
#' -----------------------------------------------------------------------------
rm(list = ls(all = TRUE))
libs <- c("tidyverse")
invisible(lapply(libs, library, character.only = TRUE))
rm(libs)

library(FoReco)
load("./BaseForecasts/ETSlev/DFbase.RData")
load("./BaseForecasts/ETSlev/Residuals.RData")
load("./VN525.RData")
DF <- NULL
resmat_all <- Residuals_ETS

DFbase$Series <- factor(DFbase$Series, colnames(VNdata), ordered = TRUE)
test_length <- as.numeric(max(DFbase$Replication))
ty_hts <- c("bu", "ols", "struc" ,"wls","shr")

K_vec <- unique(DFbase$K)
time_hts <- array(0, dim=c(test_length, length(K_vec), NROW(ty_hts)))
for (j in 1:test_length) { #test_length
  for(i in 1:length(K_vec)){
    resmat <- resmat_all[[j]][[i]]
    basef <- DFbase %>% filter(Replication==j, K==K_vec[i]) %>% 
      select(Series, Forecasts, `Forecast Horizon`) %>% 
      pivot_wider(names_from = Series, values_from = Forecasts) %>%
      arrange(`Forecast Horizon`) %>% 
      select(-`Forecast Horizon`) %>% as.matrix()

    Fltr <- DFbase %>% filter(Replication==j, K==K_vec[i]) %>% 
      dplyr::select(-"Forecasts", -"R-method", -"R-comb") %>% arrange(Series, `Forecast Horizon`)

    ## Reconciliation ----
    for(l in 1:NROW(ty_hts)){
      Start <- Sys.time()
      Recon_PointF <- htsrec(basef = basef, C = C, comb = ty_hts[l], res = resmat, 
                             type = "M", keep = "recf")
      End <- Sys.time()
      time_hts[j,i,l] <- as.numeric(difftime(End, Start, units = "secs"))
      
      # Add rows to DF
      Df1 <- cbind(Fltr, "Forecasts" = as.vector(Recon_PointF),
                   "R-method" = "hts", "R-comb" = ty_hts[l], 
                   nn = all(Recon_PointF>=0), 
                   FoReco = "direct")
      Df1 <- Df1[c(names(DFbase), "nn", "FoReco")]
      DF <- rbind(DF, Df1)
      cat(l)
    }
    cat("\nAggregate order number ", i, " (out of ", length(K_vec), ")\n", sep = "")
  }
  cat("Forecast origin number ", j, " (out of ", test_length, ", ",j/test_length*100,"%)\n", sep = "")
}

nn_data <- DF %>% filter(nn=="FALSE")%>%
  select(`R-comb`,K, Replication) %>% unique()

nn_data$time <-  NA
osqp_check <- list()
save.image(file = "./Reconciliation/ETSlev/hts_preosqp.RData")
for(row_nn in 1:NROW(nn_data)){
  if(row_nn == 1) 
    cat(" A.O. | F.O. | Comb  | Iteration", sep="")
  
  j <- nn_data$Replication[row_nn]
  i <- which(K_vec %in% nn_data$K[row_nn])
  ty <- nn_data$`R-comb`[row_nn]
  
  resmat <- resmat_all[[j]][[i]]
  basef <- DFbase %>% filter(Replication==j, K==K_vec[i]) %>% 
    select(Series, Forecasts, `Forecast Horizon`) %>% 
    pivot_wider(names_from = Series, values_from = Forecasts) %>%
    arrange(`Forecast Horizon`) %>% 
    select(-`Forecast Horizon`) %>% as.matrix()
  
  Fltr <- DFbase %>% filter(Replication==j, K==K_vec[i]) %>% 
    dplyr::select(-"Forecasts", -"R-method", -"R-comb") %>% arrange(Series, `Forecast Horizon`)
  
  ## NN reconciliation ----
  Start <- Sys.time()
  obj <- htsrec(basef = basef, C = C, comb = ty, res = resmat, nn=TRUE,
                type = "M", keep = "recf",
                settings = osqpSettings(verbose = FALSE,
                                        max_iter = 100000L,
                                        check_termination = 5,
                                        eps_abs = 1e-6, 
                                        eps_rel = 0,
                                        #eps_prim_inf = 1e-07,
                                        eps_dual_inf = 1e-07,
                                        #delta = 1e-5,
                                        polish_refine_iter = 500, polish = TRUE))
  End <- Sys.time()
  if(is.matrix(obj)){
    osqp_check[[row_nn]] <- 0
    Recon_PointF <- obj
  }else{
    osqp_check[[row_nn]] <- obj$info
    
    if(any(obj$info[,5]!=1 | obj$info[,4]>1e-6)){
      cat("\n")
      print(obj$info)
      cat("\n")
    }
    
    Recon_PointF <- obj$recf
  }
  nn_data$time[row_nn] <- as.numeric(difftime(End, Start, units = "secs"))
  
  
  # Add rows to DF
  Df1 <- cbind(Fltr, "Forecasts" = as.vector(Recon_PointF),
               "R-method" = "hts", "R-comb" = ty, 
               nn = all(Recon_PointF>=0), 
               FoReco = "osqp")
  Df1 <- Df1[c(names(DFbase), "nn", "FoReco")]
  DF <- rbind(DF, Df1)
  cat("\n  ", formatC(i, width=2, flag="0"), "  |  ", 
      formatC(j, width=2, flag="0"), "  | ", 
      format(ty, width=5), " |  ", row_nn, "/", NROW(nn_data), sep="")
}
osqp_check_hts <- osqp_check
nn_data_hts <- nn_data
DFhts <- DF
save(DFhts, time_hts, nn_data_hts, osqp_check_hts, 
     file="./Reconciliation/ETSlev/VN525_htsrec.RData")
