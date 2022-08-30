all <- list()
obj <- c(obj, "all", "i")
i=0

### Mean ----
print("------- Mean -------")
try(expr = {
  load("./BaseForecasts/VN525_features.RData")
  DFmean <- mvdf
  DFmean <- DFmean %>% rename(Forecasts = Mean) %>%
    add_column(`F-method` = "Mean", `R-method` = "Mean", `R-comb` = "Mean") %>%
    mutate(Series = factor(Series, series, ordered = TRUE), 
           K = factor(K, c(12,6,4,3,2,1), ordered = TRUE))
  i = i + 1
  all[[i]] <- mcb_data(DFmean)
  rm(list=setdiff(ls(), obj))
}, silent = TRUE)

### ETSlev ----
## base ----
print("--------------------------")
try(expr = {
  load("./BaseForecasts/ETSlev/DFbase.RData")
  DFbase <- DFbase %>% mutate(Series = factor(Series, series, ordered = TRUE), 
                              K = factor(K, c(12,6,4,3,2,1), ordered = TRUE))
  i <- i + 1
  all[[i]] <- mcb_data(DFbase)
  rm(list=setdiff(ls(), obj))
}, silent = TRUE)

## base+ ----
print("--------------------------")
try(expr = {
  load("./BaseForecasts/ETSlev/DFbase.RData")
  DFplus <- DFbase %>% mutate(Series = factor(Series, series, ordered = TRUE), 
                              K = factor(K, c(12,6,4,3,2,1), ordered = TRUE))
  DFplus$`R-method` <- paste(DFplus$`R-method`, "+", sep="")
  DFplus$Forecasts[DFplus$Forecasts<0] <- 0
  i <- i + 1
  all[[i]] <- mcb_data(DFplus)
  rm(list=setdiff(ls(), obj))
}, silent = TRUE)

## hts ----
print("--------------------------")
try(expr = {
  load("./Reconciliation/ETSlev/VN525_htsrec.RData")
  DFhts <- DFhts %>% mutate(Series = factor(Series, series, ordered = TRUE), 
                            K = factor(K, c(12,6,4,3,2,1), ordered = TRUE))
  i <- i + 1
  all[[i]] <- mcb_data(DFhts)
  rm(list=setdiff(ls(), obj))
}, silent = TRUE)

## hts_mean ----
print("--------------------------")
try(expr = {
  load("./Reconciliation/ETSlev/VN525_htsrec_mean.RData")
  DFhts <- DFhts %>% mutate(Series = factor(Series, series, ordered = TRUE), 
                            K = factor(K, c(12,6,4,3,2,1), ordered = TRUE))
  i <- i + 1
  all[[i]] <- mcb_data(DFhts)
  rm(list=setdiff(ls(), obj))
}, silent = TRUE)

## csmix ----
print("--------------------------")
try(expr = {
  load("./Reconciliation/ETSlev/VN525_csmix.RData")
  DFmix <- DFmix %>%
    mutate(Series = factor(Series, series, ordered = TRUE), 
           K = factor(K, c(12,6,4,3,2,1), ordered = TRUE))
  i <- i + 1
  all[[i]] <- mcb_data(DFmix)
  rm(list=setdiff(ls(), obj))
}, silent = TRUE)

## cslccd ----
print("--------------------------")
try(expr = {
  load("./Reconciliation/ETSlev/VN555_cslccd.RData")
  DFcslev <- DFcslev %>% filter(Series %in% series) %>%
    mutate(Series = factor(Series, series, ordered = TRUE), 
           K = factor(K, c(12,6,4,3,2,1), ordered = TRUE))
  i <- i + 1
  all[[i]] <- mcb_data(DFcslev)
  rm(list=setdiff(ls(), obj))
}, silent = TRUE)

## mixCCC ----
print("--------------------------")
try(expr = {
  load("./Reconciliation/ETSlev/VN555_CCCmix.RData")
  DFmix <- DFmix %>% filter(Series %in% series) %>%
    mutate(Series = factor(Series, series, ordered = TRUE), 
           K = factor(K, c(12,6,4,3,2,1), ordered = TRUE))
  i <- i + 1
  all[[i]] <- mcb_data(DFmix)
  rm(list=setdiff(ls(), obj))
}, silent = TRUE)

## mixCCC endo----
print("--------------------------")
try(expr = {
  load("./Reconciliation/ETSlev/VN555_CCCmix_endo.RData")
  DFmix <- DFmix %>% filter(Series %in% series) %>%
    mutate(Series = factor(Series, series, ordered = TRUE), 
           K = factor(K, c(12,6,4,3,2,1), ordered = TRUE))
  i <- i + 1
  all[[i]] <- mcb_data(DFmix)
  rm(list=setdiff(ls(), obj))
}, silent = TRUE)

## mixCCCred ----
print("--------------------------")
try(expr = {
  load("./Reconciliation/ETSlev/VN555_CCCmix_red.RData")
  DFmix <- DFmix %>% filter(Series %in% series) %>%
    mutate(Series = factor(Series, series, ordered = TRUE), 
           K = factor(K, c(12,6,4,3,2,1), ordered = TRUE))
  i <- i + 1
  all[[i]] <- mcb_data(DFmix)
  rm(list=setdiff(ls(), obj))
}, silent = TRUE)

## mixCCCred_endo ----
print("--------------------------")
try(expr = {
  load("./Reconciliation/ETSlev/VN555_CCCmix_red_endo.RData")
  DFmix <- DFmix %>% filter(Series %in% series) %>%
    mutate(Series = factor(Series, series, ordered = TRUE), 
           K = factor(K, c(12,6,4,3,2,1), ordered = TRUE))
  i <- i + 1
  all[[i]] <- mcb_data(DFmix)
  rm(list=setdiff(ls(), obj))
}, silent = TRUE)

## cslccd_mean ----
print("--------------------------")
try(expr = {
  load("./Reconciliation/ETSlev/VN555_cslccd_mean.RData")
  DFcslcc <- DFcslcc %>% filter(Series %in% series) %>%
    mutate(Series = factor(Series, series, ordered = TRUE), 
           K = factor(K, c(12,6,4,3,2,1), ordered = TRUE))
  i <- i + 1
  all[[i]] <- mcb_data(DFcslcc)
  rm(list=setdiff(ls(), obj))
}, silent = TRUE)

## cslccd_mean_red ----
print("--------------------------")
try(expr = {
  load("./Reconciliation/ETSlev/VN555_cslccd_mean_red.RData")
  DFcslcc <- DFcslcc %>% filter(Series %in% series) %>%
    mutate(Series = factor(Series, series, ordered = TRUE), 
           K = factor(K, c(12,6,4,3,2,1), ordered = TRUE))
  i <- i + 1
  all[[i]] <- mcb_data(DFcslcc)
  rm(list=setdiff(ls(), obj))
}, silent = TRUE)

## cslccd_mean_red endo ----
print("--------------------------")
try(expr = {
  load("./Reconciliation/ETSlev/VN555_cslccd_mean_red_endo.RData")
  DFcslcc <- DFcslcc %>% filter(Series %in% series) %>%
    mutate(Series = factor(Series, series, ordered = TRUE), 
           K = factor(K, c(12,6,4,3,2,1), ordered = TRUE))
  i <- i + 1
  all[[i]] <- mcb_data(DFcslcc)
  rm(list=setdiff(ls(), obj))
}, silent = TRUE)

## cslccd_bCCCred ----
print("--------------------------")
try(expr = {
  load("./Reconciliation/ETSlev/VN555_cslccd_bCCCred.RData")
  DFcslcc <- DFcslcc %>% filter(Series %in% series) %>%
    mutate(Series = factor(Series, series, ordered = TRUE), 
           K = factor(K, c(12,6,4,3,2,1), ordered = TRUE))
  i <- i + 1
  all[[i]] <- mcb_data(DFcslcc)
  rm(list=setdiff(ls(), obj))
}, silent = TRUE)

## cslccd_bCCCred_endo ----
print("--------------------------")
try(expr = {
  load("./Reconciliation/ETSlev/VN555_cslccd_bCCCred_endo.RData")
  DFcslcc <- DFcslcc %>% filter(Series %in% series) %>%
    mutate(Series = factor(Series, series, ordered = TRUE), 
           K = factor(K, c(12,6,4,3,2,1), ordered = TRUE))
  i <- i + 1
  all[[i]] <- mcb_data(DFcslcc)
  rm(list=setdiff(ls(), obj))
}, silent = TRUE)

## cslccd_mLCC ----
print("--------------------------")
try(expr = {
  load("./Reconciliation/ETSlev/VN555_cslccd_mLCC.RData")
  DFcslev <- DFcslev %>% filter(Series %in% series) %>%
    mutate(Series = factor(Series, series, ordered = TRUE), 
           K = factor(K, c(12,6,4,3,2,1), ordered = TRUE))
  i <- i + 1
  all[[i]] <- mcb_data(DFcslev)
  rm(list=setdiff(ls(), obj))
}, silent = TRUE)

## cslccd_mLCCendo ----
print("--------------------------")
try(expr = {
  load("./Reconciliation/ETSlev/VN555_cslccd_mLCCendo.RData")
  DFcslev <- DFcslev %>% filter(Series %in% series) %>%
    mutate(Series = factor(Series, series, ordered = TRUE), 
           K = factor(K, c(12,6,4,3,2,1), ordered = TRUE))
  i <- i + 1
  all[[i]] <- mcb_data(DFcslev)
  rm(list=setdiff(ls(), obj))
}, silent = TRUE)

## cslccd_bLCC ----
print("--------------------------")
try(expr = {
  load("./Reconciliation/ETSlev/VN555_cslccd_bLCC.RData")
  DFcslev <- DFcslev %>% filter(Series %in% series) %>%
    mutate(Series = factor(Series, series, ordered = TRUE), 
           K = factor(K, c(12,6,4,3,2,1), ordered = TRUE))
  i <- i + 1
  all[[i]] <- mcb_data(DFcslev)
  rm(list=setdiff(ls(), obj))
}, silent = TRUE)

## cslccd_bLCC endo----
print("--------------------------")
try(expr = {
  load("./Reconciliation/ETSlev/VN555_cslccd_bLCCendo.RData")
  DFcslev <- DFcslev %>% filter(Series %in% series) %>%
    mutate(Series = factor(Series, series, ordered = TRUE), 
           K = factor(K, c(12,6,4,3,2,1), ordered = TRUE))
  i <- i + 1
  all[[i]] <- mcb_data(DFcslev)
  rm(list=setdiff(ls(), obj))
}, silent = TRUE)

### DATASET ----
if(length(all)>0){
  mae <- lapply(all, function(x) x[["mae"]])
  mae <- do.call(rbind, mae)
  
  mse <- lapply(all, function(x) x[["mse"]])
  mse <- do.call(rbind, mse)
  
  gmae <- lapply(all, function(x) x[["gmae"]])
  gmae <- do.call(rbind, gmae)
  
  gmse <- lapply(all, function(x) x[["gmse"]])
  gmse <- do.call(rbind, gmse)
  
  imae <- lapply(all, function(x) x[["imae"]])
  imae <- do.call(rbind, imae)
  
  imse <- lapply(all, function(x) x[["imse"]])
  imse <- do.call(rbind, imse)
}
