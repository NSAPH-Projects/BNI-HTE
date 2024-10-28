##### Load Libraries & Data #####
library(gdata)
library(stringr)
library(dplyr)
library(data.table)
library(xgboost)
library(randomForest)
library(ggplot2)
library(svMisc)
library(progress)
library(caret)
library(fst)
library(sf)
library(maptools)
library(devtools)

options(scipen=999)


##### Load data #####
load("../Data/processed_data.RData")

expit <- function(x) {
  return(exp(x)/(1+exp(x)))
}

## set power plant and zipcode covariates to be included
covs <- c("logPop", "PctUrban", "PctHighSchool", "PctPoor", "PctOccupied", "PctMovedIn5", #census
          "smokerate", #smoking
          "temp", "rh", #weather
          "totNumNOxControls", "logHeatInput", "logOpTime", "pctCapacity", "pctS_n_CR", "Phase2", #power plant characteristics
          "mean_age", "Female_rate", "PctNonwhite", #Medicare characteristics
          "HyADS") #HyADS
zip.covs <- c("logPop", "PctUrban", "PctHighSchool", "PctPoor", "PctOccupied", "PctMovedIn5", #census
              "smokerate", #smoking
              "temp", "rh", #weather
              "mean_age", "Female_rate", "PctNonwhite") #Medicare characteristics
pp.covs = c("totNumNOxControls", "logHeatInput", "logOpTime", "pctCapacity", "pctS_n_CR", "Phase2")  # - these are power plant covariates that are very predictive in PS model, but likely not confounders

downweight_ratio <- 0.6786576 ## mean of ratio of 2nd key-associated and key-associated HyADS for each zip code

################################################################
## Simulation with Key-Associated Plants
################################################################

process_data <- function(denom_counts, facility_var_2005, zip_var, hyads_var, zip.covs, pp.covs) {
  zips <- data.frame("ZIP" = zip_var$ZIP) ## list of zipcodes
  plants <- data.frame("FacID" = facility_var_2005$FacID) ## list of power plant IDs
  
  ## transform covariates
  zip_var$logPop <- log(zip_var$TotPop)
  zip_var$PctNonwhite <- 1 - zip_var$White_rate
  zip_var$ZIP <- as.factor(zip_var$ZIP)
  facility_var_2005$logHeatInput <- log(facility_var_2005$totHeatInput)
  facility_var_2005$logOpTime <- log(facility_var_2005$totOpTime)
  
  zip_var <- left_join(zip_var, denom_counts, by = "ZIP")
  zip_var$IHD_admissions_rate <- zip_var$IHD_admissions/zip_var$medicare_denom
  zip_var <- zip_var[is.finite(zip_var$IHD_admissions_rate),]
  
  zip_covariates <- select(zip_var, c("ZIP", all_of(zip.covs), "num_ppl", "medicare_denom", "IHD_admissions_rate", "latitude", "longitude")) ## select zipcode covariates to be used
  plant_covariates <- select(facility_var_2005, c("FacID", all_of(pp.covs), "Fac.Latitude", "Fac.Longitude")) ## select power plant covariates to be used
  
  
  
  
  
  ##### Without filtering for complete covariate data #####
  ## connect key-associated and 2nd key-associated plant to each outcome unit
  dat <- as.data.frame(t(sapply(hyads_var, function(x) head(row.names(hyads_var)[order(x, decreasing = TRUE)], 2))))
  setDT(dat, keep.rownames = "ZIP")
  colnames(dat) <- c("ZIP", "Key", "SecondKey")
  dat <- as.data.frame(dat)
  dat$Key <- as.character(dat$Key)
  dat$SecondKey <- as.character(dat$SecondKey)
  dat$Key_HyADS <- hyads_var[as.matrix(cbind(dat$Key, dat$ZIP))] ## hyads for key-associated plants
  dat$SecondKey_HyADS <- hyads_var[as.matrix(cbind(dat$SecondKey, dat$ZIP))]
  
  
  
  
  ## aggregate mean of covariates for each power plant for its key-associated zipcodes
  zip_agg <- dat
  zip_agg <- zip_agg %>% ## merge key-associated with zip covariates
    full_join(zip_covariates, by = "ZIP")
  zip_agg[zip_agg == -Inf] <- NA ## convert -Inf to NA
  
  ## create summary of zipcodes key-associated and second key-associated with plant
  key_covariates <- aggregate(zip_agg[, !names(zip_agg) %in% c("ZIP", "Key", "SecondKey", "Key_HyADS", "SecondKey_HyADS", "IHD_admissions_rate")], list(zip_agg$Key), mean, na.rm=T)
  colnames(key_covariates) <- c("Key", paste("Key", colnames(key_covariates)[-1], sep = "_"))
  secondkey_covariates <- aggregate(zip_agg[, !names(zip_agg) %in% c("ZIP", "Key", "SecondKey", "Key_HyADS", "SecondKey_HyADS", "IHD_admissions_rate")], list(zip_agg$SecondKey), mean, na.rm=T)
  colnames(secondkey_covariates) <- c("SecondKey", paste("SecondKey",colnames(secondkey_covariates)[-1], sep = "_"))
  secondkey_covariates[,-1] <- secondkey_covariates[,-1] * downweight_ratio ## downweight second key-associated
  
  plant_dat <- key_covariates %>% ## plant, key-associated covariates, second key-associated covariates, plant covariates
    full_join(secondkey_covariates, by = c("Key"="SecondKey")) %>%
    left_join(plant_covariates, by = c("Key"="FacID"))
  
  
  
  
  
  ##### Filter by complete covariate data #####
  complete_all <- plant_dat[complete.cases(plant_dat),]
  names(complete_all)[names(complete_all) == "Key"] <- "FacID"
  plants_complete <- complete_all$FacID
  
  
  
  ## get treatments
  treatments <- left_join(complete_all, facility_var_2005[,c(1,7)], by = "FacID")
  names(treatments)[names(treatments) == "Key"] <- "FacID"
  treatments$T <- treatments$ScrubbedFacility
  
  options(scipen = 10)
  treat.covs <- c("Key_logPop", "Key_PctUrban", "Key_PctHighSchool", "Key_PctPoor", "Key_PctOccupied", "Key_PctMovedIn5", "Key_smokerate", "Key_temp", "Key_rh", "Key_mean_age", "Key_Female_rate", "Key_PctNonwhite",
                  "SecondKey_logPop", "SecondKey_PctUrban", "SecondKey_PctHighSchool", "SecondKey_PctPoor", "SecondKey_PctOccupied", "SecondKey_PctMovedIn5", "SecondKey_smokerate", "SecondKey_temp", "SecondKey_rh", "SecondKey_mean_age", "SecondKey_Female_rate", "SecondKey_PctNonwhite",
                  "totNumNOxControls", "logHeatInput", "logOpTime", "pctCapacity", "pctS_n_CR", "Phase2")
  
  
  
  
  
  #############################################
  ##### propensity score estimation #####
  #############################################
  
  treat.formula <- paste(treat.covs, collapse=" + ")
  formulaT <- as.formula(paste('T~', treat.formula))
  boost.covs <- c(treat.covs, "T")
  datT <- select(treatments, all_of(boost.covs))
  
  treatments$T <- as.numeric(treatments$T)
  
  set.seed(2023)
  rf <- randomForest(formulaT, data = treatments)
  set.seed(NULL)
  treatments$prop_est <- rf$predicted
  
  
  
  
  
  
  
  
  
  
  
  
  
  ## merge zip-level data with key-associated plant covariates and covariate summaries
  dat_merge <- dat_complete %>%
    left_join(plant_covariates, by = c("Key"="FacID")) %>%
    left_join(key_covariates, by = "Key") %>%
    left_join(secondkey_covariates, by = "SecondKey")
  if(sum(is.na(dat_merge)) > 0) {print("warning: NAs in covariates")}
  
  ## Get key and upwind treatments using plant treatment mapping
  treat_mapping <- setNames(treatments$T, treatments$FacID)
  dat_merge$Z <- as.integer(treat_mapping[dat_merge$Key]) ## based on key-associated plant: Z=1 if treated, Z=0 if not
  dat_merge$G <- as.integer(treat_mapping[dat_merge$SecondKey])
  if(sum(is.na(dat_merge$Z)) + sum(is.na(dat_merge$G)) > 0) {print("warning: NAs produced in Z and/or G")}
  
  ## Map true treatment probabilities to zipcode level
  prop_est_mapping <- setNames(treatments$prop_est, treatments$FacID)
  dat_merge$Z_prop_est <- as.numeric(prop_est_mapping[dat_merge$Key])
  dat_merge$G_prop_est <- as.numeric(prop_est_mapping[dat_merge$SecondKey])
  
  dat_merge <- left_join(dat_merge, zip_covariates, by = "ZIP")
  dat_merge <- dat_merge[complete.cases(dat_merge),]
  
  hyads_covs <- dat_merge %>% select("Key_HyADS", "SecondKey_HyADS")
  dat_merge <- dat_merge %>% select(-contains("Key_"))
  dat_merge <- cbind(dat_merge, hyads_covs)
  
  ## Get facility-zip distances
  fac_long_mapping <- setNames(treatments$Fac.Longitude, treatments$FacID)
  fac_lat_mapping <- setNames(treatments$Fac.Latitude, treatments$FacID)
  dat_merge$Key_Fac.Longitude <- fac_long_mapping[dat_merge$Key]
  dat_merge$Key_Fac.Latitude <- fac_lat_mapping[dat_merge$Key]
  dat_merge$SecondKey_Fac.Longitude <- fac_long_mapping[dat_merge$SecondKey]
  dat_merge$SecondKey_Fac.Latitude <- fac_lat_mapping[dat_merge$SecondKey]
  
  library(geosphere)
  dat_merge <- dat_merge %>%
    rowwise() %>% 
    mutate(Key_distance = as.numeric(distm(c(longitude, latitude), c(Key_Fac.Longitude, Key_Fac.Latitude), fun = distHaversine)),
           SecondKey_distance = as.numeric(distm(c(longitude, latitude), c(SecondKey_Fac.Longitude, SecondKey_Fac.Latitude), fun = distHaversine)))
  
  
  dat_merge <- filter(dat_merge, ZIP < 80000)
  
  
  # truncate propensity scores
  quantileZ <- quantile(dat_merge$Z_prop_est, c(0.05,1))
  quantileG <- quantile(dat_merge$G_prop_est, c(0.05,1))
  dat_truncated <- dat_merge %>% mutate(Z_prop_trunc = case_when(
    Z_prop_est < quantileZ[1] ~ quantileZ[1],
    Z_prop_est > quantileZ[2] ~ quantileZ[2],
    TRUE ~ Z_prop_est),
    G_prop_trunc = case_when(
      G_prop_est < quantileG[1] ~ quantileG[1],
      G_prop_est > quantileG[2] ~ quantileG[2],
      TRUE ~ G_prop_est
    )
  )
  
  dat_truncated <- subset(dat_truncated, IHD_admissions_rate < quantile(dat_truncated$IHD_admissions_rate,0.98))
  
  ## Joint propensity scores
  dat_truncated$Z11 <- dat_truncated$Z_prop_trunc * dat_truncated$G_prop_trunc ## P(Z=1,G=1 | X_out)
  dat_truncated$Z10 <- dat_truncated$Z_prop_trunc * (1-dat_truncated$G_prop_trunc) ## P(Z=1,G=0 | X_out)
  dat_truncated$Z01 <- (1-dat_truncated$Z_prop_trunc) * dat_truncated$G_prop_trunc ## P(Z=0,G=1 | X_out)
  dat_truncated$Z00 <- (1-dat_truncated$Z_prop_trunc) * (1-dat_truncated$G_prop_trunc) ## P(Z=0,G=0 | X_out)
  
  quantileZ11 <- quantile(dat_truncated$Z11, c(0.25,1))
  quantileZ10 <- quantile(dat_truncated$Z10, c(0.1,1))
  quantileZ01 <- quantile(dat_truncated$Z01, c(0.1,1))
  
  dat_truncated <- dat_truncated %>% mutate(Z11 = case_when(
    Z11 < quantileZ11[1] ~ quantileZ11[1],
    TRUE ~ Z11),
    Z10 = case_when(
      Z10 < quantileZ10[1] ~ quantileZ10[1],
      TRUE ~ Z10),
    Z01 = case_when(
      Z01 < quantileZ01[1] ~ quantileZ01[1],
      TRUE ~ Z01)
  )
  
  
  dat_truncated_orig <- dat_truncated
  return(dat_truncated)
}
                                           
####### END DATA CLEANING & TREATMENT ASSIGNMENT #######



















##############################################################################################
##### ESTIMATOR FUNCTIONS
##############################################################################################


##### AIPW estimator function #####
aipw_po <- function(z, g, dat, stabilized = FALSE,
                    pred_Z1G1, pred_Z1G0, pred_Z0G1, pred_Z0G0) {
  
  
  if(stabilized == TRUE) {
    mean_weight_z11 <- mean((dat$Z == 1 & dat$G == 1)/dat$Z11)
    mean_weight_z10 <- mean((dat$Z == 1 & dat$G == 0)/dat$Z10)
    mean_weight_z01 <- mean((dat$Z == 0 & dat$G == 1)/dat$Z01)
    mean_weight_z00 <- mean((dat$Z == 0 & dat$G == 0)/dat$Z00)
  } else {
    mean_weight_z11 <- mean_weight_z10 <- mean_weight_z01 <- mean_weight_z00 <- 1
  }
  
  po = case_when(
    z == 1 & g == 1 ~ ((dat$Z == 1 & dat$G == 1)/dat$Z11/mean_weight_z11)*dat$IHD_admissions_rate +
      (1 - (dat$Z == 1 & dat$G == 1)/dat$Z11/mean_weight_z11)*pred_Z1G1,
    z == 1 & g == 0 ~ ((dat$Z == 1 & dat$G == 0)/dat$Z10/mean_weight_z10)*dat$IHD_admissions_rate +
      (1 - (dat$Z == 1 & dat$G == 0)/dat$Z10/mean_weight_z10)*pred_Z1G0,
    z == 0 & g == 1 ~ ((dat$Z == 0 & dat$G == 1)/dat$Z01/mean_weight_z01)*dat$IHD_admissions_rate +
      (1 - (dat$Z == 0 & dat$G == 1)/dat$Z01/mean_weight_z01)*pred_Z0G1,
    z == 0 & g == 0 ~ ((dat$Z == 0 & dat$G == 0)/dat$Z00/mean_weight_z00)*dat$IHD_admissions_rate +
      (1 - (dat$Z == 0 & dat$G == 0)/dat$Z00/mean_weight_z00)*pred_Z0G0)
  
  return(po)
  
}



gcomp <- function(z, g, dat,
                  pred_Z1G1, pred_Z1G0, pred_Z0G1, pred_Z0G0) {
  
  po = case_when(
    z == 1 & g == 1 ~ pred_Z1G1,
    z == 1 & g == 0 ~ pred_Z1G0,
    z == 0 & g == 1 ~ pred_Z0G1,
    z == 0 & g == 0 ~ pred_Z0G0)
  
  return(po)
  
}









#########################################################################################################################
##### Application
#########################################################################################################################

run_application <- function(dat_truncated, zip.covs) {
  Xzip.formula <- paste(zip.covs, collapse="+")
  .env = environment()
  
  formulaY <- as.formula(paste('IHD_admissions_rate ~ ', Xzip.formula), env = .env)
  
  subset_Z1G1 <- subset(dat_truncated, Z == 1 & G == 1)
  subset_Z1G0 <- subset(dat_truncated, Z == 1 & G == 0)
  subset_Z0G1 <- subset(dat_truncated, Z == 0 & G == 1)
  subset_Z0G0 <- subset(dat_truncated, Z == 0 & G == 0)
  
  out_mod_Z1G1_glm <- glm(formulaY, data = subset_Z1G1)
  out_mod_Z1G0_glm <- glm(formulaY, data = subset_Z1G0)
  out_mod_Z0G1_glm <- glm(formulaY, data = subset_Z0G1)
  out_mod_Z0G0_glm <- glm(formulaY, data = subset_Z0G0)
  pred_Z1G1 <- predict(out_mod_Z1G1_glm, dat_truncated, type = "response")
  pred_Z1G0 <- predict(out_mod_Z1G0_glm, dat_truncated, type = "response")
  pred_Z0G1 <- predict(out_mod_Z0G1_glm, dat_truncated, type = "response")
  pred_Z0G0 <- predict(out_mod_Z0G0_glm, dat_truncated, type = "response")
  
  
  ## G-comp
  dat_truncated$G_DE_G1 <- G_DE_G1 <- (gcomp(1,1,dat_truncated, pred_Z1G1, pred_Z1G0, pred_Z0G1, pred_Z0G0) - gcomp(0,1,dat_truncated, pred_Z1G1, pred_Z1G0, pred_Z0G1, pred_Z0G0))*10000
  dat_truncated$G_DE_G0 <- G_DE_G0 <- (gcomp(1,0,dat_truncated, pred_Z1G1, pred_Z1G0, pred_Z0G1, pred_Z0G0) - gcomp(0,0,dat_truncated, pred_Z1G1, pred_Z1G0, pred_Z0G1, pred_Z0G0))*10000
  dat_truncated$G_SE_Z1 <- G_SE_Z1 <- (gcomp(1,1,dat_truncated, pred_Z1G1, pred_Z1G0, pred_Z0G1, pred_Z0G0) - gcomp(1,0,dat_truncated, pred_Z1G1, pred_Z1G0, pred_Z0G1, pred_Z0G0))*10000
  dat_truncated$G_SE_Z0 <- G_SE_Z0 <- (gcomp(0,1,dat_truncated, pred_Z1G1, pred_Z1G0, pred_Z0G1, pred_Z0G0) - gcomp(0,0,dat_truncated, pred_Z1G1, pred_Z1G0, pred_Z0G1, pred_Z0G0))*10000
  
  ## AIPW
  dat_truncated$A_DE_G1 <- A_DE_G1 <- (aipw_po(1,1,dat_truncated,stabilized=F, pred_Z1G1, pred_Z1G0, pred_Z0G1, pred_Z0G0) - aipw_po(0,1,dat_truncated,stabilized=F, pred_Z1G1, pred_Z1G0, pred_Z0G1, pred_Z0G0))*10000
  dat_truncated$A_DE_G0 <- A_DE_G0 <- (aipw_po(1,0,dat_truncated,stabilized=F, pred_Z1G1, pred_Z1G0, pred_Z0G1, pred_Z0G0) - aipw_po(0,0,dat_truncated,stabilized=F, pred_Z1G1, pred_Z1G0, pred_Z0G1, pred_Z0G0))*10000
  dat_truncated$A_SE_Z1 <- A_SE_Z1 <- (aipw_po(1,1,dat_truncated,stabilized=F, pred_Z1G1, pred_Z1G0, pred_Z0G1, pred_Z0G0) - aipw_po(1,0,dat_truncated,stabilized=F, pred_Z1G1, pred_Z1G0, pred_Z0G1, pred_Z0G0))*10000
  dat_truncated$A_SE_Z0 <- A_SE_Z0 <- (aipw_po(0,1,dat_truncated,stabilized=F, pred_Z1G1, pred_Z1G0, pred_Z0G1, pred_Z0G0) - aipw_po(0,0,dat_truncated,stabilized=F, pred_Z1G1, pred_Z1G0, pred_Z0G1, pred_Z0G0))*10000
  
  ## SAIPW
  dat_truncated$S_DE_G1 <- S_DE_G1 <- (aipw_po(1,1,dat_truncated,stabilized=T, pred_Z1G1, pred_Z1G0, pred_Z0G1, pred_Z0G0) - aipw_po(0,1,dat_truncated,stabilized=T, pred_Z1G1, pred_Z1G0, pred_Z0G1, pred_Z0G0))*10000
  dat_truncated$S_DE_G0 <- S_DE_G0 <- (aipw_po(1,0,dat_truncated,stabilized=T, pred_Z1G1, pred_Z1G0, pred_Z0G1, pred_Z0G0) - aipw_po(0,0,dat_truncated,stabilized=T, pred_Z1G1, pred_Z1G0, pred_Z0G1, pred_Z0G0))*10000
  dat_truncated$S_SE_Z1 <- S_SE_Z1 <- (aipw_po(1,1,dat_truncated,stabilized=T, pred_Z1G1, pred_Z1G0, pred_Z0G1, pred_Z0G0) - aipw_po(1,0,dat_truncated,stabilized=T, pred_Z1G1, pred_Z1G0, pred_Z0G1, pred_Z0G0))*10000
  dat_truncated$S_SE_Z0 <- S_SE_Z0 <- (aipw_po(0,1,dat_truncated,stabilized=T, pred_Z1G1, pred_Z1G0, pred_Z0G1, pred_Z0G0) - aipw_po(0,0,dat_truncated,stabilized=T, pred_Z1G1, pred_Z1G0, pred_Z0G1, pred_Z0G0))*10000
  
  
  return(dat_truncated)
}
dat_truncated <- process_data(denom_counts, facility_var_2005, zip_var, hyads_var, zip.covs, pp.covs, perturb_hyads = F)
results <- run_application(dat_truncated, zip.covs)
summary(results)








#########################################################################################################################
##### Bootstrap
#########################################################################################################################

##### Residual bootstrap #####
boot_resid <- function(data, plants, m, reps = 10000) {
  all_resamples <- list()
  boot_results <- data.frame(matrix(ncol=12, nrow=reps))
  colnames(boot_results) <- c("G_DE_G1", "G_DE_G0", "G_SE_Z1", "G_SE_Z0", 
                              "A_DE_G1", "A_DE_G0", "A_SE_Z1", "A_SE_Z0",
                              "S_DE_G1", "S_DE_G0", "S_SE_Z1", "S_SE_Z0")
  
  pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull || Time remaining :eta]",
                         total = reps,
                         complete = "=",
                         incomplete = "-",
                         current = ">",
                         clear = FALSE,
                         width = 140)
  
  zip.covs <- c("logPop", "PctUrban", "PctHighSchool", "PctPoor", "PctOccupied", "PctMovedIn5", #census
                "smokerate", #smoking
                "temp", "rh", #weather
                "mean_age", "Female_rate", "PctNonwhite") #Medicare characteristics
  Xzip.formula <- paste(zip.covs, collapse="+")
  .env = environment()
  treat.covs <- c("Key_logPop", "Key_PctUrban", "Key_PctHighSchool", "Key_PctPoor", "Key_PctOccupied", "Key_PctMovedIn5", "Key_smokerate", "Key_temp", "Key_rh", "Key_mean_age", "Key_Female_rate", "Key_PctNonwhite",
                  "SecondKey_logPop", "SecondKey_PctUrban", "SecondKey_PctHighSchool", "SecondKey_PctPoor", "SecondKey_PctOccupied", "SecondKey_PctMovedIn5", "SecondKey_smokerate", "SecondKey_temp", "SecondKey_rh", "SecondKey_mean_age", "SecondKey_Female_rate", "SecondKey_PctNonwhite",
                  "totNumNOxControls", "logHeatInput", "logOpTime", "pctCapacity", "pctS_n_CR", "Phase2")
  treat.formula <- paste(treat.covs, collapse=" + ")
  formulaT <- as.formula(paste('T~', treat.formula))
  
  
  
  formulaY <- as.formula(paste('IHD_admissions_rate ~ ', Xzip.formula), env = .env)
  formulaY_new <- as.formula(paste('IHD_admissions_rate_new ~ ', Xzip.formula), env = .env)
  
  subset_Z1G1 <- subset(data, Z == 1 & G == 1)
  subset_Z1G0 <- subset(data, Z == 1 & G == 0)
  subset_Z0G1 <- subset(data, Z == 0 & G == 1)
  subset_Z0G0 <- subset(data, Z == 0 & G == 0)
  
  mod_Z1G1_orig <- glm(formulaY, data = subset_Z1G1)
  mod_Z1G0_orig <- glm(formulaY, data = subset_Z1G0)
  mod_Z0G1_orig <- glm(formulaY, data = subset_Z0G1)
  mod_Z0G0_orig <- glm(formulaY, data = subset_Z0G0)
  
  residuals_orig <- as.vector(c(mod_Z1G1_orig$residuals, mod_Z1G0_orig$residuals, mod_Z0G1_orig$residuals, mod_Z0G0_orig$residuals))
  
  
  for(r in 1:reps) {
    pb$tick()
    
    ## OUTCOME RESIDUAL BOOTSTRAP
    subset_Z1G1_new <- subset_Z1G1 %>% ungroup() %>% mutate(
      resid = sample(residuals_orig, nrow(subset_Z1G1), replace = T),
      IHD_admissions_rate_new = IHD_admissions_rate + resid
    )
    mod_Z1G1_new <- glm(formulaY_new, data = subset_Z1G1_new)
    pred_Z1G1_new <- predict(mod_Z1G1_new, data, type = "response")
    
    
    subset_Z1G0_new <- subset_Z1G0 %>% ungroup() %>% mutate(
      resid = sample(residuals_orig, nrow(subset_Z1G0), replace = T),
      IHD_admissions_rate_new = IHD_admissions_rate + resid
    )
    mod_Z1G0_new <- glm(formulaY_new, data = subset_Z1G0_new)
    pred_Z1G0_new <- predict(mod_Z1G0_new, data, type = "response")
    
    
    subset_Z0G1_new <- subset_Z0G1 %>% ungroup() %>% mutate(
      resid = sample(residuals_orig, nrow(subset_Z0G1), replace = T),
      IHD_admissions_rate_new = IHD_admissions_rate + resid
    )
    mod_Z0G1_new <- glm(formulaY_new, data = subset_Z0G1_new)
    pred_Z0G1_new <- predict(mod_Z0G1_new, data, type = "response")
    
    
    subset_Z0G0_new <- subset_Z0G0 %>% ungroup() %>% mutate(
      resid = sample(residuals_orig, nrow(subset_Z0G0), replace = T),
      IHD_admissions_rate_new = IHD_admissions_rate + resid
    )
    mod_Z0G0_new <- glm(formulaY_new, data = subset_Z0G0_new)
    pred_Z0G0_new <- predict(mod_Z0G0_new, data, type = "response")
    
    
    G_DE_G1 <- (gcomp(1,1,data_new, pred_Z1G1_new, pred_Z1G0_new, pred_Z0G1_new, pred_Z0G0_new) - gcomp(0,1,data_new, pred_Z1G1_new, pred_Z1G0_new, pred_Z0G1_new, pred_Z0G0_new))*10000
    G_DE_G0 <- (gcomp(1,0,data_new, pred_Z1G1_new, pred_Z1G0_new, pred_Z0G1_new, pred_Z0G0_new) - gcomp(0,0,data_new, pred_Z1G1_new, pred_Z1G0_new, pred_Z0G1_new, pred_Z0G0_new))*10000
    G_SE_Z1 <- (gcomp(1,1,data_new, pred_Z1G1_new, pred_Z1G0_new, pred_Z0G1_new, pred_Z0G0_new) - gcomp(1,0,data_new, pred_Z1G1_new, pred_Z1G0_new, pred_Z0G1_new, pred_Z0G0_new))*10000
    G_SE_Z0 <- (gcomp(0,1,data_new, pred_Z1G1_new, pred_Z1G0_new, pred_Z0G1_new, pred_Z0G0_new) - gcomp(0,0,data_new, pred_Z1G1_new, pred_Z1G0_new, pred_Z0G1_new, pred_Z0G0_new))*10000
    
    ## AIPW
    A_DE_G1 <- (aipw_po(1,1,data_new,stabilized=F, pred_Z1G1_new, pred_Z1G0_new, pred_Z0G1_new, pred_Z0G0_new) - aipw_po(0,1,data_new,stabilized=F, pred_Z1G1_new, pred_Z1G0_new, pred_Z0G1_new, pred_Z0G0_new))*10000
    A_DE_G0 <- (aipw_po(1,0,data_new,stabilized=F, pred_Z1G1_new, pred_Z1G0_new, pred_Z0G1_new, pred_Z0G0_new) - aipw_po(0,0,data_new,stabilized=F, pred_Z1G1_new, pred_Z1G0_new, pred_Z0G1_new, pred_Z0G0_new))*10000
    A_SE_Z1 <- (aipw_po(1,1,data_new,stabilized=F, pred_Z1G1_new, pred_Z1G0_new, pred_Z0G1_new, pred_Z0G0_new) - aipw_po(1,0,data_new,stabilized=F, pred_Z1G1_new, pred_Z1G0_new, pred_Z0G1_new, pred_Z0G0_new))*10000
    A_SE_Z0 <- (aipw_po(0,1,data_new,stabilized=F, pred_Z1G1_new, pred_Z1G0_new, pred_Z0G1_new, pred_Z0G0_new) - aipw_po(0,0,data_new,stabilized=F, pred_Z1G1_new, pred_Z1G0_new, pred_Z0G1_new, pred_Z0G0_new))*10000
    
    ## SAIPW
    S_DE_G1 <- (aipw_po(1,1,data_new,stabilized=T, pred_Z1G1_new, pred_Z1G0_new, pred_Z0G1_new, pred_Z0G0_new) - aipw_po(0,1,data_new,stabilized=T, pred_Z1G1_new, pred_Z1G0_new, pred_Z0G1_new, pred_Z0G0_new))*10000
    S_DE_G0 <- (aipw_po(1,0,data_new,stabilized=T, pred_Z1G1_new, pred_Z1G0_new, pred_Z0G1_new, pred_Z0G0_new) - aipw_po(0,0,data_new,stabilized=T, pred_Z1G1_new, pred_Z1G0_new, pred_Z0G1_new, pred_Z0G0_new))*10000
    S_SE_Z1 <- (aipw_po(1,1,data_new,stabilized=T, pred_Z1G1_new, pred_Z1G0_new, pred_Z0G1_new, pred_Z0G0_new) - aipw_po(1,0,data_new,stabilized=T, pred_Z1G1_new, pred_Z1G0_new, pred_Z0G1_new, pred_Z0G0_new))*10000
    S_SE_Z0 <- (aipw_po(0,1,data_new,stabilized=T, pred_Z1G1_new, pred_Z1G0_new, pred_Z0G1_new, pred_Z0G0_new) - aipw_po(0,0,data_new,stabilized=T, pred_Z1G1_new, pred_Z1G0_new, pred_Z0G1_new, pred_Z0G0_new))*10000
    
    
    boot_results$G_DE_G1[r] <- mean(G_DE_G1)
    boot_results$G_DE_G0[r] <- mean(G_DE_G0)
    boot_results$G_SE_Z1[r] <- mean(G_SE_Z1)
    boot_results$G_SE_Z0[r] <- mean(G_SE_Z0)
    boot_results$A_DE_G1[r] <- mean(A_DE_G1)
    boot_results$A_DE_G0[r] <- mean(A_DE_G0)
    boot_results$A_SE_Z1[r] <- mean(A_SE_Z1)
    boot_results$A_SE_Z0[r] <- mean(A_SE_Z0)
    boot_results$S_DE_G1[r] <- mean(S_DE_G1)
    boot_results$S_DE_G0[r] <- mean(S_DE_G0)
    boot_results$S_SE_Z1[r] <- mean(S_SE_Z1)
    boot_results$S_SE_Z0[r] <- mean(S_SE_Z0)
    
    all_resamples[[r]] <- plants_new
  }
  
  return(list(boot_results, all_resamples))
  
}



boot_resid_results <- boot_resid(dat_truncated, treatments, m = nrow(dat_truncated), reps = 10000)
sapply(boot_resid_results[[1]], function(y) {quantile(y, c(0.025, 0.975))})










##########################################################################################
## Misc analyses/items
##########################################################################################




##### plots #####
load("processed_data.RData")
ggplot(treatments, aes(prop_est, fill = ScrubbedFacility)) + 
  geom_histogram(alpha = 0.35, aes(y = ..density..), position = 'identity', bins = 15) +
  geom_density(alpha = 0.1, aes(color = ScrubbedFacility)) +
  xlab("Estimated propensity score") + ylab("Density") + labs(color = "Scrubber Installed", fill = "Scrubber Installed") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank())
ggsave("prop_overlap.png", width = 6, height = 4, dpi = 700)



quantileT <- quantile(treatments$prop_est, c(0.05,1))
treatments <- treatments %>% mutate(prop_trunc = case_when(
  prop_est < quantileT[1] ~ quantileT[1],
  prop_est > quantileT[2] ~ quantileT[2],
  TRUE ~ prop_est))
ggplot(treatments, aes(prop_trunc, fill = ScrubbedFacility)) + 
  geom_histogram(alpha = 0.35, aes(y = ..density..), position = 'identity', bins = 15) +
  geom_density(alpha = 0.1, aes(color = ScrubbedFacility)) +
  xlab("Estimated propensity score") + ylab("Density") + labs(color = "Scrubber Installed", fill = "Scrubber Installed") +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank())



