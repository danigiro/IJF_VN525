#' -----------------------------------------------------------------------------
#' VN555_base.R
#'
#' Creating an RData file with the VN555 ETS levels and SA base forecasts.
#'
#' Input files: DFbase.RData residuals.RData VN525_features.RData
#'              Ctools555.RData VN525_features_residuals.RData
#' Output files: VN555_base.RData
#'
#' This code is written by Daniele Girolimetto
#' Department of Statistics, University of Padua (Italy)
#' -----------------------------------------------------------------------------
library(tidyverse)
load("./BaseForecasts/ETSlev/DFbase.RData")
load("./BaseForecasts/ETSlev/residuals.RData")
load("./BaseForecasts/VN525_features.RData")
load("./Ctools555.RData")
load("./BaseForecasts/VN525_features_residuals.RData")
df_not_unique <- data.frame(new = not_unique_series, 
                            old = names(not_unique_series))
list_df_new <- apply(df_not_unique, 1,
                     function(x) DFbase %>% filter(Series == x[2]) %>% mutate(Series = x[1]))

DFnew <- do.call(rbind, list_df_new)
DFbase <- rbind(DFbase, DFnew)
names(DFbase$Series) <- NULL

DFbase <- DFbase %>% 
  mutate(Series = factor(Series, c(rownames(C555), colnames(C555)), 
                         ordered = TRUE))

list_df_new <- apply(df_not_unique, 1,
                     function(x) mvdf %>% filter(Series == x[2]) %>% mutate(Series = x[1]))

DFnew <- do.call(rbind, list_df_new)
mvdf <- rbind(mvdf, DFnew)
names(mvdf$Series) <- NULL

mvdf <- mvdf %>% 
  mutate(Series = factor(Series, c(rownames(C555), colnames(C555)), 
                         ordered = TRUE))

name555 <- c(rownames(C555), colnames(C555))
name525 <- c(rownames(C525), colnames(C525))

ser555 <- names(not_unique_series)
names(ser555) <- not_unique_series

num_series <- 1:525
names(num_series) <- name525
id_res <- recode(recode(name555, !!!ser555), !!!num_series)

for(i in 1:length(Residuals_ETS)){
  for(j in 1:length(Residuals_ETS[[i]])){
    Residuals_ETS[[i]][[j]] <- Residuals_ETS[[i]][[j]][,id_res]
    residuals_mean[[i]][[j]] <- residuals_mean[[i]][[j]][,id_res]
  }
}
save(mvdf, not_unique_series, DFbase, Residuals_ETS, C555, name555, name525, 
     residuals_mean, file = "./BaseForecasts/ETSlev/555/VN555_base.RData")

