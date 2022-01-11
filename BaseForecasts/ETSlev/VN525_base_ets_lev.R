#' -----------------------------------------------------------------------------
#' VN525_base_ets_lev.R
#'
#' Creating an RData file of VN525 ETS levels base forecasts.
#' 
#' Base forecasts: ETS lev
#'
#' Input files: VN525.RData
#' Output files: ALL.RData, DFbase.RData and residuals.RData
#'
#' This code is written by Daniele Girolimetto
#' Department of Statistics, University of Padua (Italy)
#' -----------------------------------------------------------------------------
rm(list = ls(all = TRUE))
library(forecast)
library(tidyverse)
library(zoo)

load("./VN525.RData")
rm(C)
DF <- tibble("Date" = yearmon(NULL),
             "Series" = character(),
             "F-method" = character(),
             "F-model" = character(),
             "R-method" = character(),
             "R-comb" = character(),
             "Forecast Horizon" = integer(),
             "Forecasts" = double(),
             "Actual" = double(),
             "Training window_length" = integer(),
             "Replication" = integer(),
             "K" = integer())

# Seconds to d h m s
dhms <- function(t){
  paste(t %/% (60*60*24), "d ", 
        paste(formatC(t %/% (60*60) %% 24, width = 2, format = "d", flag = "0"), "h ",
              formatC(t %/% 60 %% 60, width = 2, format = "d", flag = "0"), "m ",
              formatC(t %% 60, width = 2, format = "d", flag = "0"), "s",
              sep = ""), 
        sep="")
}

k <- 1
H <- 12
fixed_length <- 96 # 1998.01 - 2005.12 (first)
end_traing <- NROW(VNdata)-fixed_length-1 # 10 years

# Test
# end_traing = 2
# VNdata = VNdata[,1:3]

month_all <- as.yearmon(time(VNdata[,1]))

fit_ETS_all <- list()
Residuals_ETS <- list()

time_start <- Sys.time()
for(j in 0:end_traing){
  time_start_par <- Sys.time()
  cat("--------------------------- Aggregation order k = 1 ----------------------------\n")
  fit_ETS <- list()
  month_test1 <- month_all[(j+fixed_length+1):(j+fixed_length+H)]
  train1 <- window(VNdata, start=c(1998,(j+1)), end=c(1998,j+fixed_length))
  test1 <- window(VNdata, start=c(1998,j+fixed_length+1), end=c(1998,j+fixed_length+H))
  
  Residuals_ETS1 <- matrix(NA, NROW(train1), NCOL(train1))
  for(i in 1:NCOL(VNdata)){
    ts_train <- na.omit(train1[,i])
    # ETS
    fit_ETS$k1[[i]] <- ets(ts_train)                                     # Model
    Forecast_ETS <- forecast(fit_ETS$k1[[i]], h = H)                          # Forecasts
    Residuals_ETS1[which(!is.na(train1[,i])),i] <- as.vector(ts_train - fitted(fit_ETS$k1[[i]]))    # Residuals
    
    # Add rows line to DF tibble
    for (h in 1:sum(!is.na(month_test1))) {
      DF <- DF %>% add_row("Date" = month_test1[h],
                           "Series" = colnames(VNdata)[i],
                           "F-method" = "ETS_lev",
                           "F-model" = fit_ETS$k1[[i]]$method,
                           "R-method" = "base",
                           "R-comb" = "none",
                           "Forecast Horizon" = h,
                           "Forecasts" = Forecast_ETS$mean[h],
                           "Actual" = as.numeric(test1[h,i]),
                           "Training window_length" = fixed_length,
                           "Replication" = j+1,
                           "K" = 1)
    }
    if(i%%25==0) cat(i, " ", sep="")
  }
  
  Residuals_ETS[[j+1]] <- list(k1 = Residuals_ETS1)
  time_end <- Sys.time()
  cat("\n")
  cat("Training window n.", (j+1), "(out of", 
      end_traing, ",",
      (j+1)/end_traing*100,"%)\n")
  
  parJ <- end_traing-(j+1)
  par_time <- as.numeric(difftime(time_end, time_start_par, units = "secs"))
  tot_time <- as.numeric(difftime(time_end, time_start, units = "secs"))
  cat("Elapsed time: ", dhms(tot_time),"\nEstimated time remaining: ", 
      dhms(mean(par_time)*parJ),"\n", sep="")
}

save.image(file = "./BaseForecasts/ETSlev/ALL.RData")
DFbase <- DF
save(DFbase, file = "./BaseForecasts/ETSlev/DFbase.RData")
save(Residuals_ETS, file = "./BaseForecasts/ETSlev/residuals.RData")
rm(list = ls(all = TRUE))
