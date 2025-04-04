## simulations for AIPW/SAIPW estimators under various scenarios

## read command line arguments ##
args<-commandArgs(TRUE)
for (i in 1:length(args)) { eval (parse (text = args[[i]] )) }

## command line inputs are
## wd: working directory
## simnum: slurm job array number, used to set seed for data generation
## variance: outcome error variance (default 1)
## effsize: effect size (default -1)
## reps: number of reps per job
## boot_reps: number of boot reps per sim

## load packages ##
library(gdata)
library(stringr)
library(dplyr)
library(data.table)
library(ggplot2)
library(reshape2)
library(Rcpp)


setwd(wd)
options(scipen=5)

load("dat_corr_prop.rda")
dat_corr_prop <- dat_truncated_full[order(dat_truncated_full$ZIP),]  
load("dat_miss_prop.rda") 
dat_miss_prop <- dat_truncated_full[order(dat_truncated_full$ZIP),]  
data_orig <- left_join(dat_corr_prop, dat_miss_prop[,c("ZIP", "Z","G", "Z_prop_est","G_prop_est","Z11","Z10","Z01","Z00")], 
                       by = c("ZIP"), 
                       suffix = c("_corr", "_miss"))


dimnames = list(NULL, 
                c("DE_G1_2", "DE_G1_1", "DE_G1_0", "DE_G0_2", "DE_G0_1", "DE_G0_0", 
                  "SE_Z1_2", "SE_Z1_1", "SE_Z1_0", "SE_Z0_2", "SE_Z0_1", "SE_Z0_0"))
dimnames_len <- length(dimnames[[2]])

dimnames2 = list(NULL,
                 c("G_DE_G1", "G_DE_G0", "G_SE_Z1", "G_SE_Z0", 
                   "A_DE_G1", "A_DE_G0", "A_SE_Z1", "A_SE_Z0",
                   "S_DE_G1", "S_DE_G0", "S_SE_Z1", "S_SE_Z0"))
dimnames2_len <- length(dimnames2[[2]])
dimnames2_TE = list(NULL,
                    c("DE", "SE", 
                      "G_DE_G1", "G_DE_G0", "G_SE_Z1", "G_SE_Z0", 
                      "A_DE_G1", "A_DE_G0", "A_SE_Z1", "A_SE_Z0",
                      "S_DE_G1", "S_DE_G0", "S_SE_Z1", "S_SE_Z0"))
dimnames2_TE_len <- length(dimnames2_TE[[2]])


expit <- function(x) {
  return(exp(x)/(1+exp(x)))
}

generate_po <- function(data, var, t) {
  
  #####################################################
  ## Normal outcome model, heterogeneity in both DE & SE, same effect modifiers
  
  data <- data %>% mutate(DE = t*(ifelse(PctNonwhite_bin == 0, 0,
                                         ifelse(PctPoor_bin == 0, 1, 2))), ## White = 0, Nonwhite & Not Poor = 1, Nonwhite & Poor = 2
                          SE = t*(ifelse(PctNonwhite_bin == 0, 0,
                                         ifelse(PctPoor_bin == 0, 1, 2))))
  
  data <- data %>% mutate(mu = 2*logPop + 5*smokerate + 5*PctPoor + 10*PctNonwhite
                          + 10*PctNonwhite*smokerate)
  
  data <- data %>% mutate(Y00 = rnorm(mu, mu, var),
                          Y11 = Y00 + DE + SE,
                          Y10 = Y00 + DE,
                          Y01 = Y00 + SE)
  data <- data %>% mutate(Y = Y00*(1-Z)*(1-G) + Y10*Z*(1-G) + Y01*(1-Z)*G + Y11*Z*G)
  
  return(data)
}






##############################################################################################
##### Estimator functions
##############################################################################################


##### AIPW estimator function #####
aipw_po <- function(z, g, dat, formulaY, prop_spec, stabilized = FALSE) {
  
  if(prop_spec == "corr") {
    dat <- dat %>%
      mutate(Z11 = Z11_corr,
             Z10 = Z10_corr,
             Z01 = Z01_corr,
             Z00 = Z00_corr)
  } else if(prop_spec == "miss") {
    dat <- dat %>%
      mutate(Z11 = Z11_miss,
             Z10 = Z10_miss,
             Z01 = Z01_miss,
             Z00 = Z00_miss)
  } else {print("Improper propensity score specification")}
  
  
  subset_Z1G1 <- subset(dat, Z == 1 & G == 1)
  subset_Z1G0 <- subset(dat, Z == 1 & G == 0)
  subset_Z0G1 <- subset(dat, Z == 0 & G == 1)
  subset_Z0G0 <- subset(dat, Z == 0 & G == 0)
  
  out_mod_Z1G1 <- glm(formulaY, data = subset_Z1G1, family = gaussian())
  out_mod_Z1G0 <- glm(formulaY, data = subset_Z1G0, family = gaussian())
  out_mod_Z0G1 <- glm(formulaY, data = subset_Z0G1, family = gaussian())
  out_mod_Z0G0 <- glm(formulaY, data = subset_Z0G0, family = gaussian())
  
  if(stabilized == TRUE) {
    mean_weight_z11 <- mean((dat$Z == 1 & dat$G == 1)/dat$Z11)
    mean_weight_z10 <- mean((dat$Z == 1 & dat$G == 0)/dat$Z10)
    mean_weight_z01 <- mean((dat$Z == 0 & dat$G == 1)/dat$Z01)
    mean_weight_z00 <- mean((dat$Z == 0 & dat$G == 0)/dat$Z00)
  } else {
    mean_weight_z11 <- mean_weight_z10 <- mean_weight_z01 <- mean_weight_z00 <- 1
  }
  
  
  est_df <- dat %>%
    mutate(po = case_when(
      z == 1 & g == 1 ~ ((Z == 1 & G == 1)/Z11/mean_weight_z11)*Y + 
        (1 - (Z == 1 & G == 1)/Z11/mean_weight_z11)*predict(out_mod_Z1G1, newdata = dat, type = "response"),
      z == 1 & g == 0 ~ ((Z == 1 & G == 0)/Z10/mean_weight_z10)*Y + 
        (1 - (Z == 1 & G == 0)/Z10/mean_weight_z10)*predict(out_mod_Z1G0, newdata = dat, type = "response"),
      z == 0 & g == 1 ~ ((Z == 0 & G == 1)/Z01/mean_weight_z01)*Y + 
        (1 - (Z == 0 & G == 1)/Z01/mean_weight_z01)*predict(out_mod_Z0G1, newdata = dat, type = "response"),
      z == 0 & g == 0 ~ ((Z == 0 & G == 0)/Z00/mean_weight_z00)*Y + 
        (1 - (Z == 0 & G == 0)/Z00/mean_weight_z00)*predict(out_mod_Z0G0, newdata = dat, type = "response")
    )
    )
  
  return(est_df$po)
}


gcomp <- function(z, g, dat, formulaY) {
  
  subset_Z1G1 <- subset(dat, Z == 1 & G == 1)
  subset_Z1G0 <- subset(dat, Z == 1 & G == 0)
  subset_Z0G1 <- subset(dat, Z == 0 & G == 1)
  subset_Z0G0 <- subset(dat, Z == 0 & G == 0)
  
  out_mod_Z1G1 <- glm(formulaY, data = subset_Z1G1, family = gaussian())
  out_mod_Z1G0 <- glm(formulaY, data = subset_Z1G0, family = gaussian())
  out_mod_Z0G1 <- glm(formulaY, data = subset_Z0G1, family = gaussian())
  out_mod_Z0G0 <- glm(formulaY, data = subset_Z0G0, family = gaussian())
  
  est_df <- dat %>%
    mutate(po = case_when(
      z == 1 & g == 1 ~ predict(out_mod_Z1G1, newdata = dat, type = "response"),
      z == 1 & g == 0 ~ predict(out_mod_Z1G0, newdata = dat, type = "response"),
      z == 0 & g == 1 ~ predict(out_mod_Z0G1, newdata = dat, type = "response"),
      z == 0 & g == 0 ~ predict(out_mod_Z0G0, newdata = dat, type = "response")
      
    )
    )
  
  return(est_df$po)
}


## BOOT ESTIMATOR FUNCTIONS
aipw_po_boot <- function(z, g, dat, prop_spec, stabilized = FALSE) {
  
  if(prop_spec == "corr") {
    dat <- dat %>%
      mutate(Z11 = Z11_corr,
             Z10 = Z10_corr,
             Z01 = Z01_corr,
             Z00 = Z00_corr)
  } else if(prop_spec == "miss") {
    dat <- dat %>%
      mutate(Z11 = Z11_miss,
             Z10 = Z10_miss,
             Z01 = Z01_miss,
             Z00 = Z00_miss)
  } else {print("Improper propensity score specification")}
  
  if(stabilized == TRUE) {
    mean_weight_z11 <- mean((dat$Z == 1 & dat$G == 1)/dat$Z11)
    mean_weight_z10 <- mean((dat$Z == 1 & dat$G == 0)/dat$Z10)
    mean_weight_z01 <- mean((dat$Z == 0 & dat$G == 1)/dat$Z01)
    mean_weight_z00 <- mean((dat$Z == 0 & dat$G == 0)/dat$Z00)
  } else {
    mean_weight_z11 <- mean_weight_z10 <- mean_weight_z01 <- mean_weight_z00 <- 1
  }
  
  est_df <- dat %>%
    mutate(po = case_when(
      z == 1 & g == 1 ~ ((Z == 1 & G == 1)/Z11/mean_weight_z11)*Y_new +
        (1 - (Z == 1 & G == 1)/Z11/mean_weight_z11)*pred_Z1G1_new,
      z == 1 & g == 0 ~ ((Z == 1 & G == 0)/Z10/mean_weight_z10)*Y_new +
        (1 - (Z == 1 & G == 0)/Z10/mean_weight_z10)*pred_Z1G0_new,
      z == 0 & g == 1 ~ ((Z == 0 & G == 1)/Z01/mean_weight_z01)*Y_new +
        (1 - (Z == 0 & G == 1)/Z01/mean_weight_z01)*pred_Z0G1_new,
      z == 0 & g == 0 ~ ((Z == 0 & G == 0)/Z00/mean_weight_z00)*Y_new +
        (1 - (Z == 0 & G == 0)/Z00/mean_weight_z00)*pred_Z0G0_new
    )
    )
  
  return(est_df$po)
  
}



gcomp_boot <- function(z, g, dat) {
  
  est_df <- dat %>%
    mutate(po = case_when(
      z == 1 & g == 1 ~ pred_Z1G1_new,
      z == 1 & g == 0 ~ pred_Z1G0_new,
      z == 0 & g == 1 ~ pred_Z0G1_new,
      z == 0 & g == 0 ~ pred_Z0G0_new)
    )
  
  return(est_df$po)
  
}




##### Residual bootstrap #####

boot_resid <- function(data, formulaY, prop_spec, reps = 500) {
  
  boot_ATE <- data.frame(matrix(ncol=dimnames2_TE_len, nrow=reps,
                                dimnames=dimnames2_TE))
  boot_CATE_2 <- data.frame(matrix(ncol=dimnames2_TE_len, nrow=reps,
                                   dimnames=dimnames2_TE))
  boot_CATE_1 <- data.frame(matrix(ncol=dimnames2_TE_len, nrow=reps,
                                   dimnames=dimnames2_TE))
  boot_CATE_0 <- data.frame(matrix(ncol=dimnames2_TE_len, nrow=reps,
                                   dimnames=dimnames2_TE))
  
  
  subset_Z1G1 <- subset(data, Z == 1 & G == 1)
  subset_Z1G0 <- subset(data, Z == 1 & G == 0)
  subset_Z0G1 <- subset(data, Z == 0 & G == 1)
  subset_Z0G0 <- subset(data, Z == 0 & G == 0)
  
  mod_Z1G1_orig <- glm(formulaY, data = subset_Z1G1)
  mod_Z1G0_orig <- glm(formulaY, data = subset_Z1G0)
  mod_Z0G1_orig <- glm(formulaY, data = subset_Z0G1)
  mod_Z0G0_orig <- glm(formulaY, data = subset_Z0G0)
  
  residuals_Z1G1 <- mod_Z1G1_orig$residuals
  residuals_Z1G0 <- mod_Z1G0_orig$residuals
  residuals_Z0G1 <- mod_Z0G1_orig$residuals
  residuals_Z0G0 <- mod_Z0G0_orig$residuals
  residuals_all <- c(residuals_Z1G1, residuals_Z1G0, residuals_Z0G1, residuals_Z0G0)
  
  for(r in 1:reps) {
    set.seed(simnum*10000+r)
    data_new <- data
    
    formulaY_new <- update.formula(formulaY, Y_new ~ .)
    
    ## OUTCOME RESIDUAL BOOTSTRAP
    subset_Z1G1_new <- subset_Z1G1 %>% ungroup() %>% mutate(
      resid = sample(residuals_Z1G1, nrow(subset_Z1G1), replace = T),
      Y_new = Y + resid
    )
    mod_Z1G1_new <- glm(formulaY_new, data = subset_Z1G1_new)
    data_new$pred_Z1G1_new <- predict(mod_Z1G1_new, data, type = "response")
    
    
    subset_Z1G0_new <- subset_Z1G0 %>% ungroup() %>% mutate(
      resid = sample(residuals_Z1G0, nrow(subset_Z1G0), replace = T),
      Y_new = Y + resid
    )
    mod_Z1G0_new <- glm(formulaY_new, data = subset_Z1G0_new)
    data_new$pred_Z1G0_new <- predict(mod_Z1G0_new, data, type = "response")
    
    
    subset_Z0G1_new <- subset_Z0G1 %>% ungroup() %>% mutate(
      resid = sample(residuals_Z0G1, nrow(subset_Z0G1), replace = T),
      Y_new = Y + resid
    )
    mod_Z0G1_new <- glm(formulaY_new, data = subset_Z0G1_new)
    data_new$pred_Z0G1_new <- predict(mod_Z0G1_new, data, type = "response")
    
    
    subset_Z0G0_new <- subset_Z0G0 %>% ungroup() %>% mutate(
      resid = sample(residuals_Z0G0, nrow(subset_Z0G0), replace = T),
      Y_new = Y + resid
    )
    mod_Z0G0_new <- glm(formulaY_new, data = subset_Z0G0_new)
    data_new$pred_Z0G0_new <- predict(mod_Z0G0_new, data, type = "response")
    
    new_Y <- select(rbind(subset_Z1G1_new, subset_Z1G0_new, subset_Z0G1_new, subset_Z0G0_new), ZIP, resid, Y_new) 
    data_new <- left_join(data_new, new_Y, by = "ZIP")
    
    ## GCOMP
    G_DE_G1 <- (gcomp_boot(1,1,data_new) - gcomp_boot(0,1,data_new))
    G_DE_G0 <- (gcomp_boot(1,0,data_new) - gcomp_boot(0,0,data_new))
    G_SE_Z1 <- (gcomp_boot(1,1,data_new) - gcomp_boot(1,0,data_new))
    G_SE_Z0 <- (gcomp_boot(0,1,data_new) - gcomp_boot(0,0,data_new))
    
    ## AIPW
    A_DE_G1 <- (aipw_po_boot(1,1,data_new,prop_spec,stabilized=F) - aipw_po_boot(0,1,data_new,prop_spec,stabilized=F))
    A_DE_G0 <- (aipw_po_boot(1,0,data_new,prop_spec,stabilized=F) - aipw_po_boot(0,0,data_new,prop_spec,stabilized=F))
    A_SE_Z1 <- (aipw_po_boot(1,1,data_new,prop_spec,stabilized=F) - aipw_po_boot(1,0,data_new,prop_spec,stabilized=F))
    A_SE_Z0 <- (aipw_po_boot(0,1,data_new,prop_spec,stabilized=F) - aipw_po_boot(0,0,data_new,prop_spec,stabilized=F))
    
    ## SAIPW
    S_DE_G1 <- (aipw_po_boot(1,1,data_new,prop_spec,stabilized=T) - aipw_po_boot(0,1,data_new,prop_spec,stabilized=T))
    S_DE_G0 <- (aipw_po_boot(1,0,data_new,prop_spec,stabilized=T) - aipw_po_boot(0,0,data_new,prop_spec,stabilized=T))
    S_SE_Z1 <- (aipw_po_boot(1,1,data_new,prop_spec,stabilized=T) - aipw_po_boot(1,0,data_new,prop_spec,stabilized=T))
    S_SE_Z0 <- (aipw_po_boot(0,1,data_new,prop_spec,stabilized=T) - aipw_po_boot(0,0,data_new,prop_spec,stabilized=T))
    
    ## estimates
    boot_est <- data.frame(DE = data_new$DE, SE = data_new$SE,
                           DE_G1_G = G_DE_G1, DE_G0_G = G_DE_G0, SE_Z1_G = G_SE_Z1, SE_Z0_G = G_SE_Z0,
                           DE_G1_A = A_DE_G1, DE_G0_A = A_DE_G0, SE_Z1_A = A_SE_Z1, SE_Z0_A = A_SE_Z0,
                           DE_G1_S = S_DE_G1, DE_G0_S = S_DE_G0, SE_Z1_S = S_SE_Z1, SE_Z0_S = S_SE_Z0)
    
    ## overall means
    boot_ATE[r, ] <- colMeans(boot_est)
    
    ## subgroup means
    boot_CATE <- boot_est %>%
      group_by(DE) %>%
      summarize(across(everything(), mean, na.rm = T))
    boot_CATE_2[r,] <- boot_CATE[1,]
    boot_CATE_1[r,] <- boot_CATE[2,]
    boot_CATE_0[r,] <- boot_CATE[3,]
    
  }
  
  
  return(list(boot_ATE, boot_CATE_2, boot_CATE_1, boot_CATE_0))
  
}





#########################################################################################################################
##### SIMULATIONS
#########################################################################################################################


DE_G1_G <- DE_G0_G <- SE_Z1_G <- SE_Z0_G <- list()
DE_G1_A <- DE_G0_A <- SE_Z1_A <- SE_Z0_A <- list()
DE_G1_S <- DE_G0_S <- SE_Z1_S <- SE_Z0_S <- list()
ATE_scenA <- ATE_scenB <- ATE_scenC <- ATE_scenD <- data.frame(matrix(ncol=dimnames2_TE_len, nrow=reps, 
                                                                      dimnames=dimnames2_TE))
CATE_2_scenA <- CATE_2_scenB <- CATE_2_scenC <- CATE_2_scenD <- data.frame(matrix(ncol=dimnames2_TE_len, nrow=reps, 
                                                                                  dimnames=dimnames2_TE))
CATE_1_scenA <- CATE_1_scenB <- CATE_1_scenC <- CATE_1_scenD <- data.frame(matrix(ncol=dimnames2_TE_len, nrow=reps, 
                                                                                  dimnames=dimnames2_TE))
CATE_0_scenA <- CATE_0_scenB <- CATE_0_scenC <- CATE_0_scenD <- data.frame(matrix(ncol=dimnames2_TE_len, nrow=reps, 
                                                                                  dimnames=dimnames2_TE))

boot_ATE_scenA <- boot_ATE_scenB <- boot_ATE_scenC <- boot_ATE_scenD <- list()
boot_ATE_CI_scenA <- boot_ATE_CI_scenB <- boot_ATE_CI_scenC <- boot_ATE_CI_scenD <- list()
boot_ATE_SE_scenA <- boot_ATE_SE_scenB <- boot_ATE_SE_scenC <- boot_ATE_SE_scenD <- data.frame(matrix(ncol=dimnames2_len, nrow=reps,
                                                                                                      dimnames=dimnames2))
coverage_ATE_scenA <- coverage_ATE_scenB <- coverage_ATE_scenC <- coverage_ATE_scenD <- data.frame(matrix(ncol=dimnames2_len, nrow=reps,
                                                                                                          dimnames=dimnames2))
boot_CATE_2_scenA <- boot_CATE_2_scenB <- boot_CATE_2_scenC <- boot_CATE_2_scenD <- list()
boot_CATE_2_CI_scenA <- boot_CATE_2_CI_scenB <- boot_CATE_2_CI_scenC <- boot_CATE_2_CI_scenD <- list()
boot_CATE_2_SE_scenA <- boot_CATE_2_SE_scenB <- boot_CATE_2_SE_scenC <- boot_CATE_2_SE_scenD <- data.frame(matrix(ncol=dimnames2_len, nrow=reps,
                                                                                                                  dimnames=dimnames2))
coverage_CATE_2_scenA <- coverage_CATE_2_scenB <- coverage_CATE_2_scenC <- coverage_CATE_2_scenD <- data.frame(matrix(ncol=dimnames2_len, nrow=reps,
                                                                                                                      dimnames=dimnames2))
boot_CATE_1_scenA <- boot_CATE_1_scenB <- boot_CATE_1_scenC <- boot_CATE_1_scenD <- list()
boot_CATE_1_CI_scenA <- boot_CATE_1_CI_scenB <- boot_CATE_1_CI_scenC <- boot_CATE_1_CI_scenD <- list()
boot_CATE_1_SE_scenA <- boot_CATE_1_SE_scenB <- boot_CATE_1_SE_scenC <- boot_CATE_1_SE_scenD <- data.frame(matrix(ncol=dimnames2_len, nrow=reps,
                                                                                                                  dimnames=dimnames2))
coverage_CATE_1_scenA <- coverage_CATE_1_scenB <- coverage_CATE_1_scenC <- coverage_CATE_1_scenD <- data.frame(matrix(ncol=dimnames2_len, nrow=reps,
                                                                                                                      dimnames=dimnames2))
boot_CATE_0_scenA <- boot_CATE_0_scenB <- boot_CATE_0_scenC <- boot_CATE_0_scenD <- list()
boot_CATE_0_CI_scenA <- boot_CATE_0_CI_scenB <- boot_CATE_0_CI_scenC <- boot_CATE_0_CI_scenD <- list()
boot_CATE_0_SE_scenA <- boot_CATE_0_SE_scenB <- boot_CATE_0_SE_scenC <- boot_CATE_0_SE_scenD <- data.frame(matrix(ncol=dimnames2_len, nrow=reps,
                                                                                                                  dimnames=dimnames2))
coverage_CATE_0_scenA <- coverage_CATE_0_scenB <- coverage_CATE_0_scenC <- coverage_CATE_0_scenD <- data.frame(matrix(ncol=dimnames2_len, nrow=reps,
                                                                                                                      dimnames=dimnames2))

for(i in 1:reps) { 
  set.seed(simnum*1000+i)
  data <- generate_po(data_orig, var = variance, t = effsize) 
  
  zip.covs <- c("logPop", "smokerate", "PctPoor", "PctNonwhite")
  Xzip.formula_corr <- paste(zip.covs, collapse=" + ")
  formulaY_corr <- as.formula(paste('Y ~ ', Xzip.formula_corr, '+ PctNonwhite:smokerate'))

  Xzip.formula_miss <- paste(zip.covs, collapse=" + ")
  formulaY_miss <- as.formula(paste('Y ~ ', Xzip.formula_miss))
  
  DE_all <- data$DE
  SE_all <- data$SE
  
  ## Gcomp Estimator
  DE_G1_G[[i]] <- data.frame(DE = DE_all, 
                             est_scenA = gcomp(1,1,data,formulaY_corr) - gcomp(0,1,data,formulaY_corr),
                             est_scenB = gcomp(1,1,data,formulaY_corr) - gcomp(0,1,data,formulaY_corr),
                             est_scenC = gcomp(1,1,data,formulaY_miss) - gcomp(0,1,data,formulaY_miss),
                             est_scenD = gcomp(1,1,data,formulaY_miss) - gcomp(0,1,data,formulaY_miss)
  )
  DE_G0_G[[i]] <- data.frame(DE = DE_all, 
                             est_scenA = gcomp(1,0,data,formulaY_corr) - gcomp(0,0,data,formulaY_corr),
                             est_scenB = gcomp(1,0,data,formulaY_corr) - gcomp(0,0,data,formulaY_corr),
                             est_scenC = gcomp(1,0,data,formulaY_miss) - gcomp(0,0,data,formulaY_miss),
                             est_scenD = gcomp(1,0,data,formulaY_miss) - gcomp(0,0,data,formulaY_miss)
  )
  SE_Z1_G[[i]] <- data.frame(SE = SE_all, 
                             est_scenA = gcomp(1,1,data,formulaY_corr) - gcomp(1,0,data,formulaY_corr),
                             est_scenB = gcomp(1,1,data,formulaY_corr) - gcomp(1,0,data,formulaY_corr),
                             est_scenC = gcomp(1,1,data,formulaY_miss) - gcomp(1,0,data,formulaY_miss),
                             est_scenD = gcomp(1,1,data,formulaY_miss) - gcomp(1,0,data,formulaY_miss)
  )
  SE_Z0_G[[i]] <- data.frame(SE = SE_all, 
                             est_scenA = gcomp(0,1,data,formulaY_corr) - gcomp(0,0,data,formulaY_corr),
                             est_scenB = gcomp(0,1,data,formulaY_corr) - gcomp(0,0,data,formulaY_corr),
                             est_scenC = gcomp(0,1,data,formulaY_miss) - gcomp(0,0,data,formulaY_miss),
                             est_scenD = gcomp(0,1,data,formulaY_miss) - gcomp(0,0,data,formulaY_miss)
  )
  
  ## AIPW Estimator
  DE_G1_A[[i]] <- data.frame(DE = DE_all, 
                             est_scenA = aipw_po(1,1,data,formulaY_corr,prop_spec = "corr",stabilized = F) - aipw_po(0,1,data,formulaY_corr,prop_spec = "corr",stabilized = F),
                             est_scenB = aipw_po(1,1,data,formulaY_corr,prop_spec = "miss",stabilized = F) - aipw_po(0,1,data,formulaY_corr,prop_spec = "miss",stabilized = F),
                             est_scenC = aipw_po(1,1,data,formulaY_miss,prop_spec = "corr",stabilized = F) - aipw_po(0,1,data,formulaY_miss,prop_spec = "corr",stabilized = F),
                             est_scenD = aipw_po(1,1,data,formulaY_miss,prop_spec = "miss",stabilized = F) - aipw_po(0,1,data,formulaY_miss,prop_spec = "miss",stabilized = F)
  )
  DE_G0_A[[i]] <- data.frame(DE = DE_all, 
                             est_scenA = aipw_po(1,0,data,formulaY_corr,prop_spec = "corr",stabilized = F) - aipw_po(0,0,data,formulaY_corr,prop_spec = "corr",stabilized = F),
                             est_scenB = aipw_po(1,0,data,formulaY_corr,prop_spec = "miss",stabilized = F) - aipw_po(0,0,data,formulaY_corr,prop_spec = "miss",stabilized = F),
                             est_scenC = aipw_po(1,0,data,formulaY_miss,prop_spec = "corr",stabilized = F) - aipw_po(0,0,data,formulaY_miss,prop_spec = "corr",stabilized = F),
                             est_scenD = aipw_po(1,0,data,formulaY_miss,prop_spec = "miss",stabilized = F) - aipw_po(0,0,data,formulaY_miss,prop_spec = "miss",stabilized = F)
  )
  SE_Z1_A[[i]] <- data.frame(SE = SE_all, 
                             est_scenA = aipw_po(1,1,data,formulaY_corr,prop_spec = "corr",stabilized = F) - aipw_po(1,0,data,formulaY_corr,prop_spec = "corr",stabilized = F),
                             est_scenB = aipw_po(1,1,data,formulaY_corr,prop_spec = "miss",stabilized = F) - aipw_po(1,0,data,formulaY_corr,prop_spec = "miss",stabilized = F),
                             est_scenC = aipw_po(1,1,data,formulaY_miss,prop_spec = "corr",stabilized = F) - aipw_po(1,0,data,formulaY_miss,prop_spec = "corr",stabilized = F),
                             est_scenD = aipw_po(1,1,data,formulaY_miss,prop_spec = "miss",stabilized = F) - aipw_po(1,0,data,formulaY_miss,prop_spec = "miss",stabilized = F)
  )
  SE_Z0_A[[i]] <- data.frame(SE = SE_all, 
                             est_scenA = aipw_po(0,1,data,formulaY_corr,prop_spec = "corr",stabilized = F) - aipw_po(0,0,data,formulaY_corr,prop_spec = "corr",stabilized = F),
                             est_scenB = aipw_po(0,1,data,formulaY_corr,prop_spec = "miss",stabilized = F) - aipw_po(0,0,data,formulaY_corr,prop_spec = "miss",stabilized = F),
                             est_scenC = aipw_po(0,1,data,formulaY_miss,prop_spec = "corr",stabilized = F) - aipw_po(0,0,data,formulaY_miss,prop_spec = "corr",stabilized = F),
                             est_scenD = aipw_po(0,1,data,formulaY_miss,prop_spec = "miss",stabilized = F) - aipw_po(0,0,data,formulaY_miss,prop_spec = "miss",stabilized = F)
  )
  
  ## SAIPW Estimator
  DE_G1_S[[i]] <- data.frame(DE = DE_all, 
                             est_scenA = aipw_po(1,1,data,formulaY_corr,prop_spec = "corr",stabilized = T) - aipw_po(0,1,data,formulaY_corr,prop_spec = "corr",stabilized = T),
                             est_scenB = aipw_po(1,1,data,formulaY_corr,prop_spec = "miss",stabilized = T) - aipw_po(0,1,data,formulaY_corr,prop_spec = "miss",stabilized = T),
                             est_scenC = aipw_po(1,1,data,formulaY_miss,prop_spec = "corr",stabilized = T) - aipw_po(0,1,data,formulaY_miss,prop_spec = "corr",stabilized = T),
                             est_scenD = aipw_po(1,1,data,formulaY_miss,prop_spec = "miss",stabilized = T) - aipw_po(0,1,data,formulaY_miss,prop_spec = "miss",stabilized = T)
  )
  DE_G0_S[[i]] <- data.frame(DE = DE_all, 
                             est_scenA = aipw_po(1,0,data,formulaY_corr,prop_spec = "corr",stabilized = T) - aipw_po(0,0,data,formulaY_corr,prop_spec = "corr",stabilized = T),
                             est_scenB = aipw_po(1,0,data,formulaY_corr,prop_spec = "miss",stabilized = T) - aipw_po(0,0,data,formulaY_corr,prop_spec = "miss",stabilized = T),
                             est_scenC = aipw_po(1,0,data,formulaY_miss,prop_spec = "corr",stabilized = T) - aipw_po(0,0,data,formulaY_miss,prop_spec = "corr",stabilized = T),
                             est_scenD = aipw_po(1,0,data,formulaY_miss,prop_spec = "miss",stabilized = T) - aipw_po(0,0,data,formulaY_miss,prop_spec = "miss",stabilized = T)
  )
  SE_Z1_S[[i]] <- data.frame(SE = SE_all, 
                             est_scenA = aipw_po(1,1,data,formulaY_corr,prop_spec = "corr",stabilized = T) - aipw_po(1,0,data,formulaY_corr,prop_spec = "corr",stabilized = T),
                             est_scenB = aipw_po(1,1,data,formulaY_corr,prop_spec = "miss",stabilized = T) - aipw_po(1,0,data,formulaY_corr,prop_spec = "miss",stabilized = T),
                             est_scenC = aipw_po(1,1,data,formulaY_miss,prop_spec = "corr",stabilized = T) - aipw_po(1,0,data,formulaY_miss,prop_spec = "corr",stabilized = T),
                             est_scenD = aipw_po(1,1,data,formulaY_miss,prop_spec = "miss",stabilized = T) - aipw_po(1,0,data,formulaY_miss,prop_spec = "miss",stabilized = T)
  )
  SE_Z0_S[[i]] <- data.frame(SE = SE_all, 
                             est_scenA = aipw_po(0,1,data,formulaY_corr,prop_spec = "corr",stabilized = T) - aipw_po(0,0,data,formulaY_corr,prop_spec = "corr",stabilized = T),
                             est_scenB = aipw_po(0,1,data,formulaY_corr,prop_spec = "miss",stabilized = T) - aipw_po(0,0,data,formulaY_corr,prop_spec = "miss",stabilized = T),
                             est_scenC = aipw_po(0,1,data,formulaY_miss,prop_spec = "corr",stabilized = T) - aipw_po(0,0,data,formulaY_miss,prop_spec = "corr",stabilized = T),
                             est_scenD = aipw_po(0,1,data,formulaY_miss,prop_spec = "miss",stabilized = T) - aipw_po(0,0,data,formulaY_miss,prop_spec = "miss",stabilized = T)
  )
  
  
  
  
  ## summarize all estimates
  all_est_scenA <- data.frame(DE = DE_all, SE = SE_all,
                              DE_G1_G = DE_G1_G[[i]]$est_scenA, DE_G0_G = DE_G0_G[[i]]$est_scenA, SE_Z1_G = SE_Z1_G[[i]]$est_scenA, SE_Z0_G = SE_Z0_G[[i]]$est_scenA,
                              DE_G1_A = DE_G1_A[[i]]$est_scenA, DE_G0_A = DE_G0_A[[i]]$est_scenA, SE_Z1_A = SE_Z1_A[[i]]$est_scenA, SE_Z0_A = SE_Z0_A[[i]]$est_scenA,
                              DE_G1_S = DE_G1_S[[i]]$est_scenA, DE_G0_S = DE_G0_S[[i]]$est_scenA, SE_Z1_S = SE_Z1_S[[i]]$est_scenA, SE_Z0_S = SE_Z0_S[[i]]$est_scenA)
  all_est_scenB <- data.frame(DE = DE_all, SE = SE_all,
                              DE_G1_G = DE_G1_G[[i]]$est_scenB, DE_G0_G = DE_G0_G[[i]]$est_scenB, SE_Z1_G = SE_Z1_G[[i]]$est_scenB, SE_Z0_G = SE_Z0_G[[i]]$est_scenB,
                              DE_G1_A = DE_G1_A[[i]]$est_scenB, DE_G0_A = DE_G0_A[[i]]$est_scenB, SE_Z1_A = SE_Z1_A[[i]]$est_scenB, SE_Z0_A = SE_Z0_A[[i]]$est_scenB,
                              DE_G1_S = DE_G1_S[[i]]$est_scenB, DE_G0_S = DE_G0_S[[i]]$est_scenB, SE_Z1_S = SE_Z1_S[[i]]$est_scenB, SE_Z0_S = SE_Z0_S[[i]]$est_scenB)
  all_est_scenC <- data.frame(DE = DE_all, SE = SE_all,
                              DE_G1_G = DE_G1_G[[i]]$est_scenC, DE_G0_G = DE_G0_G[[i]]$est_scenC, SE_Z1_G = SE_Z1_G[[i]]$est_scenC, SE_Z0_G = SE_Z0_G[[i]]$est_scenC,
                              DE_G1_A = DE_G1_A[[i]]$est_scenC, DE_G0_A = DE_G0_A[[i]]$est_scenC, SE_Z1_A = SE_Z1_A[[i]]$est_scenC, SE_Z0_A = SE_Z0_A[[i]]$est_scenC,
                              DE_G1_S = DE_G1_S[[i]]$est_scenC, DE_G0_S = DE_G0_S[[i]]$est_scenC, SE_Z1_S = SE_Z1_S[[i]]$est_scenC, SE_Z0_S = SE_Z0_S[[i]]$est_scenC)
  all_est_scenD <- data.frame(DE = DE_all, SE = SE_all,
                              DE_G1_G = DE_G1_G[[i]]$est_scenD, DE_G0_G = DE_G0_G[[i]]$est_scenD, SE_Z1_G = SE_Z1_G[[i]]$est_scenD, SE_Z0_G = SE_Z0_G[[i]]$est_scenD,
                              DE_G1_A = DE_G1_A[[i]]$est_scenD, DE_G0_A = DE_G0_A[[i]]$est_scenD, SE_Z1_A = SE_Z1_A[[i]]$est_scenD, SE_Z0_A = SE_Z0_A[[i]]$est_scenD,
                              DE_G1_S = DE_G1_S[[i]]$est_scenD, DE_G0_S = DE_G0_S[[i]]$est_scenD, SE_Z1_S = SE_Z1_S[[i]]$est_scenD, SE_Z0_S = SE_Z0_S[[i]]$est_scenD)
  
  
  ## overall means
  ATE_scenA[i,] <- colMeans(all_est_scenA)
  ATE_scenB[i,] <- colMeans(all_est_scenB)
  ATE_scenC[i,] <- colMeans(all_est_scenC)
  ATE_scenD[i,] <- colMeans(all_est_scenD)
  
  
  ## subgroup means
  CATE_scenA <- all_est_scenA %>%
    group_by(DE) %>%
    summarize(across(everything(), mean, na.rm = T))
  CATE_scenB <- all_est_scenB %>%
    group_by(DE) %>%
    summarize(across(everything(), mean, na.rm = T))
  CATE_scenC <- all_est_scenC %>%
    group_by(DE) %>%
    summarize(across(everything(), mean, na.rm = T))
  CATE_scenD <- all_est_scenD %>%
    group_by(DE) %>%
    summarize(across(everything(), mean, na.rm = T))
  
  CATE_2_scenA[i,] <- CATE_scenA[1,]
  CATE_2_scenB[i,] <- CATE_scenB[1,]
  CATE_2_scenC[i,] <- CATE_scenC[1,]
  CATE_2_scenD[i,] <- CATE_scenD[1,]
  CATE_1_scenA[i,] <- CATE_scenA[2,]
  CATE_1_scenB[i,] <- CATE_scenB[2,]
  CATE_1_scenC[i,] <- CATE_scenC[2,]
  CATE_1_scenD[i,] <- CATE_scenD[2,]
  CATE_0_scenA[i,] <- CATE_scenA[3,]
  CATE_0_scenB[i,] <- CATE_scenB[3,]
  CATE_0_scenC[i,] <- CATE_scenC[3,]
  CATE_0_scenD[i,] <- CATE_scenD[3,]
  
  
  ## run bootstrap
  boot_scenA <- boot_resid(data, formulaY_corr, prop_spec = "corr", reps = boot_reps)
  boot_ATE_scenA[[i]] <- boot_scenA[[1]]
  boot_ATE_SE_scenA[i,] <- sapply(boot_ATE_scenA[[i]], FUN = sd)[-c(1:2)]
  boot_ATE_CI_scenA[[i]] <- sapply(boot_ATE_scenA[[i]], function(y) {quantile(y, c(0.025, 0.975))})[,-c(1:2)]
  coverage_ATE_scenA[i,] <- as.numeric(ATE_scenA$DE[i] >= boot_ATE_CI_scenA[[i]][1,] & ATE_scenA$DE[i] <= boot_ATE_CI_scenA[[i]][2,])
  
  boot_CATE_2_scenA[[i]] <- boot_scenA[[2]]
  boot_CATE_2_SE_scenA[i,] <- sapply(boot_CATE_2_scenA[[i]], FUN = sd)[-c(1:2)]
  boot_CATE_2_CI_scenA[[i]] <- sapply(boot_CATE_2_scenA[[i]], function(y) {quantile(y, c(0.025, 0.975))})[,-c(1:2)]
  coverage_CATE_2_scenA[i,] <- as.numeric(2*effsize >= boot_CATE_2_CI_scenA[[i]][1,] & 2*effsize <= boot_CATE_2_CI_scenA[[i]][2,])
  
  boot_CATE_1_scenA[[i]] <- boot_scenA[[3]]
  boot_CATE_1_SE_scenA[i,] <- sapply(boot_CATE_1_scenA[[i]], FUN = sd)[-c(1:2)]
  boot_CATE_1_CI_scenA[[i]] <- sapply(boot_CATE_1_scenA[[i]], function(y) {quantile(y, c(0.025, 0.975))})[,-c(1:2)]
  coverage_CATE_1_scenA[i,] <- as.numeric(1*effsize >= boot_CATE_1_CI_scenA[[i]][1,] & 1*effsize <= boot_CATE_1_CI_scenA[[i]][2,])
  
  boot_CATE_0_scenA[[i]] <- boot_scenA[[4]]
  boot_CATE_0_SE_scenA[i,] <- sapply(boot_CATE_0_scenA[[i]], FUN = sd)[-c(1:2)]
  boot_CATE_0_CI_scenA[[i]] <- sapply(boot_CATE_0_scenA[[i]], function(y) {quantile(y, c(0.025, 0.975))})[,-c(1:2)]
  coverage_CATE_0_scenA[i,] <- as.numeric(0*effsize >= boot_CATE_0_CI_scenA[[i]][1,] & 0*effsize <= boot_CATE_0_CI_scenA[[i]][2,])
  
  
  
  boot_scenB <- boot_resid(data, formulaY_corr, prop_spec = "miss", reps = boot_reps)
  boot_ATE_scenB[[i]] <- boot_scenB[[1]]
  boot_ATE_SE_scenB[i,] <- sapply(boot_ATE_scenB[[i]], FUN = sd)[-c(1:2)]
  boot_ATE_CI_scenB[[i]] <- sapply(boot_ATE_scenB[[i]], function(y) {quantile(y, c(0.025, 0.975))})[,-c(1:2)]
  coverage_ATE_scenB[i,] <- as.numeric(ATE_scenB$DE[i] >= boot_ATE_CI_scenB[[i]][1,] & ATE_scenB$DE[i] <= boot_ATE_CI_scenB[[i]][2,])
  
  boot_CATE_2_scenB[[i]] <- boot_scenB[[2]]
  boot_CATE_2_SE_scenB[i,] <- sapply(boot_CATE_2_scenB[[i]], FUN = sd)[-c(1:2)]
  boot_CATE_2_CI_scenB[[i]] <- sapply(boot_CATE_2_scenB[[i]], function(y) {quantile(y, c(0.025, 0.975))})[,-c(1:2)]
  coverage_CATE_2_scenB[i,] <- as.numeric(2*effsize >= boot_CATE_2_CI_scenB[[i]][1,] & 2*effsize <= boot_CATE_2_CI_scenB[[i]][2,])
  
  boot_CATE_1_scenB[[i]] <- boot_scenB[[3]]
  boot_CATE_1_SE_scenB[i,] <- sapply(boot_CATE_1_scenB[[i]], FUN = sd)[-c(1:2)]
  boot_CATE_1_CI_scenB[[i]] <- sapply(boot_CATE_1_scenB[[i]], function(y) {quantile(y, c(0.025, 0.975))})[,-c(1:2)]
  coverage_CATE_1_scenB[i,] <- as.numeric(1*effsize >= boot_CATE_1_CI_scenB[[i]][1,] & 1*effsize <= boot_CATE_1_CI_scenB[[i]][2,])
  
  boot_CATE_0_scenB[[i]] <- boot_scenB[[4]]
  boot_CATE_0_SE_scenB[i,] <- sapply(boot_CATE_0_scenB[[i]], FUN = sd)[-c(1:2)]
  boot_CATE_0_CI_scenB[[i]] <- sapply(boot_CATE_0_scenB[[i]], function(y) {quantile(y, c(0.025, 0.975))})[,-c(1:2)]
  coverage_CATE_0_scenB[i,] <- as.numeric(0*effsize >= boot_CATE_0_CI_scenB[[i]][1,] & 0*effsize <= boot_CATE_0_CI_scenB[[i]][2,])
  
  
  
  boot_scenC <- boot_resid(data, formulaY_miss, prop_spec = "corr", reps = boot_reps)
  boot_ATE_scenC[[i]] <- boot_scenC[[1]]
  boot_ATE_SE_scenC[i,] <- sapply(boot_ATE_scenC[[i]], FUN = sd)[-c(1:2)]
  boot_ATE_CI_scenC[[i]] <- sapply(boot_ATE_scenC[[i]], function(y) {quantile(y, c(0.025, 0.975))})[,-c(1:2)]
  coverage_ATE_scenC[i,] <- as.numeric(ATE_scenC$DE[i] >= boot_ATE_CI_scenC[[i]][1,] & ATE_scenC$DE[i] <= boot_ATE_CI_scenC[[i]][2,])
  
  boot_CATE_2_scenC[[i]] <- boot_scenC[[2]]
  boot_CATE_2_SE_scenC[i,] <- sapply(boot_CATE_2_scenC[[i]], FUN = sd)[-c(1:2)]
  boot_CATE_2_CI_scenC[[i]] <- sapply(boot_CATE_2_scenC[[i]], function(y) {quantile(y, c(0.025, 0.975))})[,-c(1:2)]
  coverage_CATE_2_scenC[i,] <- as.numeric(2*effsize >= boot_CATE_2_CI_scenC[[i]][1,] & 2*effsize <= boot_CATE_2_CI_scenC[[i]][2,])
  
  boot_CATE_1_scenC[[i]] <- boot_scenC[[3]]
  boot_CATE_1_SE_scenC[i,] <- sapply(boot_CATE_1_scenC[[i]], FUN = sd)[-c(1:2)]
  boot_CATE_1_CI_scenC[[i]] <- sapply(boot_CATE_1_scenC[[i]], function(y) {quantile(y, c(0.025, 0.975))})[,-c(1:2)]
  coverage_CATE_1_scenC[i,] <- as.numeric(1*effsize >= boot_CATE_1_CI_scenC[[i]][1,] & 1*effsize <= boot_CATE_1_CI_scenC[[i]][2,])
  
  boot_CATE_0_scenC[[i]] <- boot_scenC[[4]]
  boot_CATE_0_SE_scenC[i,] <- sapply(boot_CATE_0_scenC[[i]], FUN = sd)[-c(1:2)]
  boot_CATE_0_CI_scenC[[i]] <- sapply(boot_CATE_0_scenC[[i]], function(y) {quantile(y, c(0.025, 0.975))})[,-c(1:2)]
  coverage_CATE_0_scenC[i,] <- as.numeric(0*effsize >= boot_CATE_0_CI_scenC[[i]][1,] & 0*effsize <= boot_CATE_0_CI_scenC[[i]][2,])
  
  
  
  boot_scenD <- boot_resid(data, formulaY_miss, prop_spec = "miss", reps = boot_reps)
  boot_ATE_scenD[[i]] <- boot_scenD[[1]]
  boot_ATE_SE_scenD[i,] <- sapply(boot_ATE_scenD[[i]], FUN = sd)[-c(1:2)]
  boot_ATE_CI_scenD[[i]] <- sapply(boot_ATE_scenD[[i]], function(y) {quantile(y, c(0.025, 0.975))})[,-c(1:2)]
  coverage_ATE_scenD[i,] <- as.numeric(ATE_scenD$DE[i] >= boot_ATE_CI_scenD[[i]][1,] & ATE_scenD$DE[i] <= boot_ATE_CI_scenD[[i]][2,])
  
  boot_CATE_2_scenD[[i]] <- boot_scenD[[2]]
  boot_CATE_2_SE_scenD[i,] <- sapply(boot_CATE_2_scenD[[i]], FUN = sd)[-c(1:2)]
  boot_CATE_2_CI_scenD[[i]] <- sapply(boot_CATE_2_scenD[[i]], function(y) {quantile(y, c(0.025, 0.975))})[,-c(1:2)]
  coverage_CATE_2_scenD[i,] <- as.numeric(2*effsize >= boot_CATE_2_CI_scenD[[i]][1,] & 2*effsize <= boot_CATE_2_CI_scenD[[i]][2,])
  
  boot_CATE_1_scenD[[i]] <- boot_scenD[[3]]
  boot_CATE_1_SE_scenD[i,] <- sapply(boot_CATE_1_scenD[[i]], FUN = sd)[-c(1:2)]
  boot_CATE_1_CI_scenD[[i]] <- sapply(boot_CATE_1_scenD[[i]], function(y) {quantile(y, c(0.025, 0.975))})[,-c(1:2)]
  coverage_CATE_1_scenD[i,] <- as.numeric(1*effsize >= boot_CATE_1_CI_scenD[[i]][1,] & 1*effsize <= boot_CATE_1_CI_scenD[[i]][2,])
  
  boot_CATE_0_scenD[[i]] <- boot_scenD[[4]]
  boot_CATE_0_SE_scenD[i,] <- sapply(boot_CATE_0_scenD[[i]], FUN = sd)[-c(1:2)]
  boot_CATE_0_CI_scenD[[i]] <- sapply(boot_CATE_0_scenD[[i]], function(y) {quantile(y, c(0.025, 0.975))})[,-c(1:2)]
  coverage_CATE_0_scenD[i,] <- as.numeric(0*effsize >= boot_CATE_0_CI_scenD[[i]][1,] & 0*effsize <= boot_CATE_0_CI_scenD[[i]][2,])
  
}








save(all_est_scenA, all_est_scenB, all_est_scenC, all_est_scenD,
     ATE_scenA, ATE_scenB, ATE_scenC, ATE_scenD,
     CATE_2_scenA, CATE_2_scenB, CATE_2_scenC, CATE_2_scenD,
     CATE_1_scenA, CATE_1_scenB, CATE_1_scenC, CATE_1_scenD,
     CATE_0_scenA, CATE_0_scenB, CATE_0_scenC, CATE_0_scenD,
     boot_ATE_SE_scenA, boot_ATE_SE_scenB, boot_ATE_SE_scenC, boot_ATE_SE_scenD,
     boot_CATE_2_SE_scenA, boot_CATE_2_SE_scenB, boot_CATE_2_SE_scenC, boot_CATE_2_SE_scenD,
     boot_CATE_1_SE_scenA, boot_CATE_1_SE_scenB, boot_CATE_1_SE_scenC, boot_CATE_1_SE_scenD,
     boot_CATE_0_SE_scenA, boot_CATE_0_SE_scenB, boot_CATE_0_SE_scenC, boot_CATE_0_SE_scenD,
     boot_ATE_CI_scenA, boot_ATE_CI_scenB, boot_ATE_CI_scenC, boot_ATE_CI_scenD, 
     boot_CATE_2_CI_scenA, boot_CATE_2_CI_scenB, boot_CATE_2_CI_scenC, boot_CATE_2_CI_scenD,
     boot_CATE_1_CI_scenA, boot_CATE_1_CI_scenB, boot_CATE_1_CI_scenC, boot_CATE_1_CI_scenD,
     boot_CATE_0_CI_scenA, boot_CATE_0_CI_scenB, boot_CATE_0_CI_scenC, boot_CATE_0_CI_scenD,
     coverage_ATE_scenA, coverage_ATE_scenB, coverage_ATE_scenC, coverage_ATE_scenD,
     coverage_CATE_2_scenA, coverage_CATE_2_scenB, coverage_CATE_2_scenC, coverage_CATE_2_scenD, 
     coverage_CATE_1_scenA, coverage_CATE_1_scenB, coverage_CATE_1_scenC, coverage_CATE_1_scenD, 
     coverage_CATE_0_scenA, coverage_CATE_0_scenB, coverage_CATE_0_scenC, coverage_CATE_0_scenD, 
     file=paste0('results/sim','_reps',reps,'_',simnum,'.RData'))
