#' -----------------------------------------------------------------------------
#' VN525_features.R
#'
#' Creating an RData file with the residuals for the SA base forecasts.
#'
#' Input files: VN525.RData VN525_features.RData
#' Output files: VN525_features_residuals.RData
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
load("./BaseForecasts/VN525_features.RData")

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
residuals_mean <- list()

time_start <- Sys.time()
for(j in 0:end_traing){
  res_list <- list()
  
  month_test1 <- month_all[(j+fixed_length+1):(j+fixed_length+H)]
  train1 <- window(VNdata, start=c(1998,(j+1)), end=c(1998,j+fixed_length))
  
  means <- mvdf %>% filter(Replication==j+1, K == 1) %>% 
    select(Series, Mean, `Forecast Horizon`) %>% 
    pivot_wider(names_from = Series, values_from = Mean) %>%
    arrange(`Forecast Horizon`) %>% 
    select(-`Forecast Horizon`) %>% as.matrix()
  
  means <- do.call(rbind, rep(list(means[,colnames(train1)]), NROW(train1)/NROW(means)))
  
  res_list$k1 <- train1 - means
  
  residuals_mean[[j+1]] <- res_list
  cat(j+1, " ")
}

save(residuals_mean, file = "./BaseForecasts/VN525_features_residuals.RData")

