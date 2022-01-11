#' -----------------------------------------------------------------------------
#' VN525_features.R
#'
#' Creating an RData file with some usefull tools for the level conditional 
#' forecast reconciliation.
#' 
#' Input files: VN525.RData
#' Output files: VN525_features.RData
#'
#' This code is written by Daniele Girolimetto
#' Department of Statistics, University of Padua (Italy)
#' -----------------------------------------------------------------------------
rm(list = ls(all = TRUE))
library(forecast)
library(tidyverse)
library(zoo)
library(progress)

load("./VN525.RData")
rm(C)
mvdf <- tibble("Date" = yearmon(NULL),
               "Series" = character(),
               "Forecast Horizon" = integer(),
               "Mean" = double(),
               "SeasVar" = double(),
               "MeanVar" = double(),
               "fc_scale" = double(),
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
end_traing <- NROW(VNdata)-1-fixed_length # 10 years

# Test
# end_traing = 2
# VNdata = VNdata[,1:3]

month_all <- as.yearmon(time(VNdata[,1]))

time_start <- Sys.time()
residuals_mean <- list()
for(j in 0:end_traing){
  time_start_par <- Sys.time()
  month_test1 <- month_all[(j+fixed_length+1):(j+fixed_length+H)]
  train1 <- window(VNdata, start=c(1998,(j+1)), end=c(1998,j+fixed_length))
  
  obj <- apply(train1, 2, function(x){
    dfm <- tibble(x = x, time = rep(1:12, NROW(train1)/12)) 
    dfmm <- dfm %>%
      drop_na() %>%
      group_by(time) %>%
      summarise(mean = mean(x), .groups = "drop")
    dfv <- dfm %>%
      drop_na() %>%
      group_by(time) %>%
      summarise(var = var(x), .groups = "drop")
    n2 <- inner_join(dfm, dfmm, by = "time") %>% drop_na() %>%
      summarise(n2 = (x-mean)^2) %>% pull()
    fc_scale <- mean(abs(diff(x, 12)))
    list(v = sum(n2)/(length(n2)-12), m = pull(dfmm, var = "mean"),
         vs = pull(dfv, var = "var"), fc_scale = fc_scale)
  })
  
  mMatrix <- do.call("cbind", lapply(obj, function(x) x[["m"]]))
  var <- sapply(obj, function(x) x[["v"]])
  vMatrix <- do.call("cbind", lapply(obj, function(x) x[["vs"]]))
  fc_scale <- sapply(obj, function(x) x[["fc_scale"]])
  
  for(h in 1:sum(!is.na(month_test1))){
    mvdf <- mvdf %>% add_row("Date" = month_test1[h],
                             "Series" = colnames(mMatrix),
                             "Forecast Horizon" = h,
                             "Mean" = mMatrix[h,],
                             "SeasVar" = vMatrix[h,], 
                             "MeanVar" = var,
                             "fc_scale" = fc_scale,
                             "Training window_length" = fixed_length,
                             "Replication" = j+1,
                             "K" = 1)
  }
  
  
  means <- do.call(rbind, rep(list(mMatrix[,colnames(train1)]), NROW(train1)/NROW(mMatrix)))
  res_list <- list()
  res_list$k1 <- train1 - means
  
  residuals_mean[[j+1]] <- res_list
  
  time_end <- Sys.time()
  parJ <- end_traing-(j+1)
  par_time <- as.numeric(difftime(time_end, time_start_par, units = "secs"))
  tot_time <- as.numeric(difftime(time_end, time_start, units = "secs"))
  cat("Elapsed time: ", dhms(tot_time),"\nEstimated time remaining: ", 
      dhms(mean(par_time)*parJ),"\n", sep="")
}
save(mvdf, file = "./BaseForecasts/VN525_features.RData")
save(residuals_mean, file = "./BaseForecasts/VN525_features_residuals.RData")
