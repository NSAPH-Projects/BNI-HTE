library(gdata)
library(stringr)
library(dplyr)
library(data.table)
library(ggplot2)
library(svMisc)
library(progress)
library(caret)
library(reshape2)
library(Rcpp)

#install.packages("https://cran.r-project.org/src/contrib/Archive/randomForest/randomForest_4.6-12.tar.gz", repos=NULL, type="source") 
library(randomForest)




expit <- function(x) {
  return(exp(x)/(1+exp(x)))
}

generate_po <- function(dat, scenario, var, t) {
  
  n <- nrow(dat)
  
  if(scenario == 1) {
    #####################################################
    ## Normal outcome model, homogeneous treatment effect (no effect modification)
    
    dat$DE <- rep(-2, nrow(dat))
    dat$SE <- rep(-1, nrow(dat))
    
    dat <- dat %>% mutate(mu = 2*logPop + 5*smokerate + 5*PctPoor + 10*PctNonwhite
                          + 5*PctNonwhite*smokerate)
    
    dat <- dat %>% mutate(Y00 = rnorm(mu, mu, var))
  }
  
  if(scenario == 2) {
    #####################################################
    ## Normal outcome model, heterogeneity in DE only
    
    dat <- dat %>% mutate(DE = t*(ifelse(PctNonwhite_bin == 0, 0,
                                         ifelse(PctPoor_bin == 0, 1, 2))), ## White = 0, Nonwhite & Not Poor = 1, Nonwhite & Poor = 2
                          SE = t)
    
    dat <- dat %>% mutate(mu = 2*logPop + 5*smokerate + 5*PctPoor + 10*PctNonwhite
                          + 5*PctNonwhite*smokerate)
    
    dat <- dat %>% mutate(Y00 = rnorm(mu, mu, var))
  }
  
  if(scenario == 3) {
    ####################################################
    ## Normal outcome model, heterogeneity in SE only
    
    dat <- dat %>% mutate(DE = t,
                          SE = t*(ifelse(PctNonwhite_bin == 0, 0,
                                         ifelse(PctPoor_bin == 0, 1, 2)))) ## White = 0, Nonwhite & Not Poor = 1, Nonwhite & Poor = 2
    
    dat <- dat %>% mutate(mu = 2*logPop + 5*smokerate + 5*PctPoor + 10*PctNonwhite
                          + 5*PctNonwhite*smokerate)
    
    dat <- dat %>% mutate(Y00 = rnorm(mu, mu, var))
  }
  
  if(scenario == 4) {
    #####################################################
    ## Normal outcome model, heterogeneity in both DE & SE, same effect modifiers
    
    dat <- dat %>% mutate(DE = t*(ifelse(PctNonwhite_bin == 0, 0,
                                         ifelse(PctPoor_bin == 0, 1, 2))), ## White = 0, Nonwhite & Not Poor = 1, Nonwhite & Poor = 2
                          SE = t*(ifelse(PctNonwhite_bin == 0, 0,
                                         ifelse(PctPoor_bin == 0, 1, 2))))
    
    dat <- dat %>% mutate(mu = 2*logPop + 5*smokerate + 5*PctPoor + 10*PctNonwhite
                          + 10*PctNonwhite*smokerate)
    
    dat <- dat %>% mutate(Y00 = rnorm(mu, mu, var))
  }
  
  if(scenario == 5) {
    #####################################################
    ## Normal outcome model, heterogeneity in both DE & SE, different effect modifiers
    
    dat <- dat %>% mutate(DE = t*(ifelse(PctNonwhite_bin == 0, 0,
                                         ifelse(PctPoor_bin == 0, 1, 2))), ## White = 0, Nonwhite & Not Poor = 1, Nonwhite & Poor = 2
                          SE = t*(ifelse(PctNonwhite_bin == 0, 0, 1))) ## White = 0, Nonwhite = 1
    
    dat <- dat %>% mutate(mu = 2*logPop + 5*smokerate + 5*PctPoor + 10*PctNonwhite
                          + 5*PctNonwhite*smokerate)
    
    dat <- dat %>% mutate(Y00 = rnorm(mu, mu, var))
  }
  
  dat <- dat %>% mutate(Y11 = Y00 + DE + SE,
                        Y10 = Y00 + DE,
                        Y01 = Y00 + SE)
  dat <- dat %>% mutate(Y = Y00*(1-Z)*(1-G) + Y10*Z*(1-G) + Y01*(1-Z)*G + Y11*Z*G)
  
  
  return(dat)
}










##############################################################################################
##### AIPW ESTIMATION
##############################################################################################


##### AIPW estimator function #####
aipw_po <- function(z, g, dat, formulaY, stabilized = FALSE) {
  
  subset_Z1G1 <- subset(dat, Z == 1 & G == 1)
  subset_Z1G0 <- subset(dat, Z == 1 & G == 0)
  subset_Z0G1 <- subset(dat, Z == 0 & G == 1)
  subset_Z0G0 <- subset(dat, Z == 0 & G == 0)
  
  out_mod_Z1G1 <- glm(formulaY, data = subset_Z1G1, family = gaussian())
  out_mod_Z1G0 <- glm(formulaY, data = subset_Z1G0, family = gaussian())
  out_mod_Z0G1 <- glm(formulaY, data = subset_Z0G1, family = gaussian())
  out_mod_Z0G0 <- glm(formulaY, data = subset_Z0G0, family = gaussian())
  
  # out_mod_Z1G1 <- randomForest(formulaY, data = subset_Z1G1)
  # out_mod_Z1G0 <- randomForest(formulaY, data = subset_Z1G0)
  # out_mod_Z0G1 <- randomForest(formulaY, data = subset_Z0G1)
  # out_mod_Z0G0 <- randomForest(formulaY, data = subset_Z0G0)
  
  if(stabilized == TRUE) {
    mean_weight_z11 <- mean((dat$Z == 1 & dat$G == 1)/dat$Z11)
    mean_weight_z10 <- mean((dat$Z == 1 & dat$G == 0)/dat$Z10)
    mean_weight_z01 <- mean((dat$Z == 0 & dat$G == 1)/dat$Z01)
    mean_weight_z00 <- mean((dat$Z == 0 & dat$G == 0)/dat$Z00)
  } else {
    mean_weight_z11 <- mean_weight_z10 <- mean_weight_z01 <- mean_weight_z00 <- 1
  }
  
  
  aipw_df <- dat %>%
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
  
  return(aipw_df$po)
  
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
  
  # out_mod_Z1G1 <- randomForest(formulaY, data = subset_Z1G1)
  # out_mod_Z1G0 <- randomForest(formulaY, data = subset_Z1G0)
  # out_mod_Z0G1 <- randomForest(formulaY, data = subset_Z0G1)
  # out_mod_Z0G0 <- randomForest(formulaY, data = subset_Z0G0)
  
  aipw_df <- dat %>%
    mutate(po = case_when(
      z == 1 & g == 1 ~ predict(out_mod_Z1G1, newdata = dat, type = "response"),
      z == 1 & g == 0 ~ predict(out_mod_Z1G0, newdata = dat, type = "response"),
      z == 0 & g == 1 ~ predict(out_mod_Z0G1, newdata = dat, type = "response"),
      z == 0 & g == 0 ~ predict(out_mod_Z0G0, newdata = dat, type = "response")
      
      # z == 1 & g == 1 ~ predict(out_mod_Z1G1, dat),
      # z == 1 & g == 0 ~ predict(out_mod_Z1G0, dat),
      # z == 0 & g == 1 ~ predict(out_mod_Z0G1, dat),
      # z == 0 & g == 0 ~ predict(out_mod_Z0G0, dat)
    )
    )
  
  return(aipw_df$po)
  
}




#data <- generate_po(dat_truncated_full, scenario = 2, var = 1, t = -1)

#########################################################################################################################
##### SIMULATIONS
#########################################################################################################################

run_sims <- function(num_iter = 1, reps = 10000, sample_prop = 1, scenario = 4, variance = 1, effsize = -1,
                    prop_type, outcome_type, output_filename) {
  
  
  if(prop_type == "corr") { load("dat_corr_prop.rda") } 
  else if(prop_type == "miss") { load("dat_miss_prop.rda") }
  else { print("invalid propensity score specification")}

  zip.covs <- c("logPop", "smokerate", "PctPoor", "PctNonwhite")
  Xzip.formula <- paste(zip.covs, collapse=" + ")
  
  
  if(outcome_type == "corr") {
    formulaY <- as.formula(paste('Y ~ ', Xzip.formula, '+ PctNonwhite:smokerate'))
  }
  else if(outcome_type == "miss") {
    formulaY <- as.formula(paste('Y ~ ', Xzip.formula))
  }
  
  if(scenario == 2) { dimnames = list(NULL, 
                                      c("DE_G1_2", "DE_G1_1", "DE_G1_0", "DE_G0_2", "DE_G0_1", "DE_G0_0", 
                                        "SE_Z1_1", "SE_Z0_1")) }
  
  if(scenario == 4) { dimnames = list(NULL, 
                                      c("DE_G1_2", "DE_G1_1", "DE_G1_0", "DE_G0_2", "DE_G0_1", "DE_G0_0", 
                                        "SE_Z1_2", "SE_Z1_1", "SE_Z1_0", "SE_Z0_2", "SE_Z0_1", "SE_Z0_0")) }
  
  if(scenario == 5) { dimnames = list(NULL, 
                                      c("DE_G1_2", "DE_G1_1", "DE_G1_0", "DE_G0_2", "DE_G0_1", "DE_G0_0", 
                                        "SE_Z1_1", "SE_Z1_0", "SE_Z0_1", "SE_Z0_0")) }
  
  dimnames_len <- length(dimnames[[2]])
  
  ATE_G <- data.frame(matrix(ncol=dimnames_len, nrow=num_iter, 
                             dimnames=dimnames))
  
  ATE_A <- data.frame(matrix(ncol=dimnames_len, nrow=num_iter, 
                             dimnames=dimnames))
  
  ATE_S <- data.frame(matrix(ncol=dimnames_len, nrow=num_iter, 
                             dimnames=dimnames))
  
  
  
  pb <- progress_bar$new(format = "(:spin) [:bar] :percent [Elapsed time: :elapsedfull Est completion: :eta]",
                         total = reps,
                         complete = "=",   # Completion bar character
                         incomplete = "-", # Incomplete bar character
                         current = ">",    # Current bar character
                         clear = FALSE,    # If TRUE, clears the bar when finish
                         width = 140)      # Width of the progress bar
  
  for(j in 1:num_iter) {
    if(sample_prop != 1) { pb$tick() }
    
    dat_truncated <- slice_sample(dat_truncated_full, n = floor(sample_prop*nrow(dat_truncated_full)), replace = FALSE)
    
    
    DE_G1_G <- DE_G0_G <- SE_Z1_G <- SE_Z0_G <- list()
    DE_G1_A <- DE_G0_A <- SE_Z1_A <- SE_Z0_A <- list()
    DE_G1_S <- DE_G0_S <- SE_Z1_S <- SE_Z0_S <- list()
    
    for(i in 1:reps) {
      if(sample_prop == 1) { pb$tick() }
      
      data <- generate_po(dat_truncated, scenario = scenario, var = variance, t = effsize)
      
      
      ## Gcomp Estimator
      DE_G1_G[[i]] <- data.frame(est = gcomp(1,1,data,formulaY) - gcomp(0,1,data,formulaY))
      DE_G0_G[[i]] <- data.frame(est = gcomp(1,0,data,formulaY) - gcomp(0,0,data,formulaY))
      SE_Z1_G[[i]] <- data.frame(est = gcomp(1,1,data,formulaY) - gcomp(1,0,data,formulaY))
      SE_Z0_G[[i]] <- data.frame(est = gcomp(0,1,data,formulaY) - gcomp(0,0,data,formulaY))
      
      # AIPW Estimator
      DE_G1_A[[i]] <- data.frame(est = aipw_po(1,1,data,formulaY,stabilized = F) - aipw_po(0,1,data,formulaY,stabilized = F))
      DE_G0_A[[i]] <- data.frame(est = aipw_po(1,0,data,formulaY,stabilized = F) - aipw_po(0,0,data,formulaY,stabilized = F))
      SE_Z1_A[[i]] <- data.frame(est = aipw_po(1,1,data,formulaY,stabilized = F) - aipw_po(1,0,data,formulaY,stabilized = F))
      SE_Z0_A[[i]] <- data.frame(est = aipw_po(0,1,data,formulaY,stabilized = F) - aipw_po(0,0,data,formulaY,stabilized = F))
      
      # SAIPW Estimator
      DE_G1_S[[i]] <- data.frame(est = aipw_po(1,1,data,formulaY,stabilized = T) - aipw_po(0,1,data,formulaY,stabilized = T))
      DE_G0_S[[i]] <- data.frame(est = aipw_po(1,0,data,formulaY,stabilized = T) - aipw_po(0,0,data,formulaY,stabilized = T))
      SE_Z1_S[[i]] <- data.frame(est = aipw_po(1,1,data,formulaY,stabilized = T) - aipw_po(1,0,data,formulaY,stabilized = T))
      SE_Z0_S[[i]] <- data.frame(est = aipw_po(0,1,data,formulaY,stabilized = T) - aipw_po(0,0,data,formulaY,stabilized = T))
      
    }
    
    options(scipen=5)
    
    if(sample_prop != 1) {
      
      # G-comp
      mean_DE_G1_G <- as.vector(rowMeans(do.call(cbind, DE_G1_G)))
      mean_DE_G0_G <- as.vector(rowMeans(do.call(cbind, DE_G0_G)))
      mean_SE_Z1_G <- as.vector(rowMeans(do.call(cbind, SE_Z1_G)))
      mean_SE_Z0_G <- as.vector(rowMeans(do.call(cbind, SE_Z0_G)))
      
      df_G <- cbind(data, mean_DE_G1_G, mean_DE_G0_G, mean_SE_Z1_G, mean_SE_Z0_G)
      
      ate_DE_G1_G <- df_G %>%
        group_by(DE) %>%
        summarize(mean = mean(mean_DE_G1_G))
      ate_DE_G0_G <- df_G %>%
        group_by(DE) %>%
        summarize(mean = mean(mean_DE_G0_G))
      ate_SE_Z1_G <- df_G %>%
        group_by(SE) %>%
        summarize(mean = mean(mean_SE_Z1_G))
      ate_SE_Z0_G <- df_G %>%
        group_by(SE) %>%
        summarize(mean = mean(mean_SE_Z0_G))
      
      ATE_G[j,]$DE_G1_2 <- ate_DE_G1_G[1,]$mean
      ATE_G[j,]$DE_G1_1 <- ate_DE_G1_G[2,]$mean
      ATE_G[j,]$DE_G1_0 <- ate_DE_G1_G[3,]$mean
      ATE_G[j,]$DE_G0_2 <- ate_DE_G0_G[1,]$mean
      ATE_G[j,]$DE_G0_1 <- ate_DE_G0_G[2,]$mean
      ATE_G[j,]$DE_G0_0 <- ate_DE_G0_G[3,]$mean
      if(scenario == 2) {
        ATE_G[j,]$SE_Z1_1 <- ate_SE_Z1_G[1,]$mean
        ATE_G[j,]$SE_Z0_1 <- ate_SE_Z0_G[1,]$mean
      }
      if(scenario == 4) {
        ATE_G[j,]$SE_Z1_2 <- ate_SE_Z1_G[1,]$mean
        ATE_G[j,]$SE_Z1_1 <- ate_SE_Z1_G[2,]$mean
        ATE_G[j,]$SE_Z1_0 <- ate_SE_Z1_G[3,]$mean
        ATE_G[j,]$SE_Z0_2 <- ate_SE_Z0_G[1,]$mean
        ATE_G[j,]$SE_Z0_1 <- ate_SE_Z0_G[2,]$mean
        ATE_G[j,]$SE_Z0_0 <- ate_SE_Z0_G[3,]$mean
      }
      
      # AIPW
      mean_DE_G1_A <- as.vector(rowMeans(do.call(cbind, DE_G1_A)))
      mean_DE_G0_A <- as.vector(rowMeans(do.call(cbind, DE_G0_A)))
      mean_SE_Z1_A <- as.vector(rowMeans(do.call(cbind, SE_Z1_A)))
      mean_SE_Z0_A <- as.vector(rowMeans(do.call(cbind, SE_Z0_A)))
      
      df_A <- cbind(data, mean_DE_G1_A, mean_DE_G0_A, mean_SE_Z1_A, mean_SE_Z0_A)
      
      ate_DE_G1_A <- df_A %>%
        group_by(DE) %>%
        summarize(mean = mean(mean_DE_G1_A))
      ate_DE_G0_A <- df_A %>%
        group_by(DE) %>%
        summarize(mean = mean(mean_DE_G0_A))
      ate_SE_Z1_A <- df_A %>%
        group_by(SE) %>%
        summarize(mean = mean(mean_SE_Z1_A))
      ate_SE_Z0_A <- df_A %>%
        group_by(SE) %>%
        summarize(mean = mean(mean_SE_Z0_A))
      
      ATE_A[j,]$DE_G1_2 <- ate_DE_G1_A[1,]$mean
      ATE_A[j,]$DE_G1_1 <- ate_DE_G1_A[2,]$mean
      ATE_A[j,]$DE_G1_0 <- ate_DE_G1_A[3,]$mean
      ATE_A[j,]$DE_G0_2 <- ate_DE_G0_A[1,]$mean
      ATE_A[j,]$DE_G0_1 <- ate_DE_G0_A[2,]$mean
      ATE_A[j,]$DE_G0_0 <- ate_DE_G0_A[3,]$mean
      if(scenario == 2) {
        ATE_A[j,]$SE_Z1_1 <- ate_SE_Z1_A[1,]$mean
        ATE_A[j,]$SE_Z0_1 <- ate_SE_Z0_A[1,]$mean
      }
      if(scenario == 4) {
        ATE_A[j,]$SE_Z1_2 <- ate_SE_Z1_A[1,]$mean
        ATE_A[j,]$SE_Z1_1 <- ate_SE_Z1_A[2,]$mean
        ATE_A[j,]$SE_Z1_0 <- ate_SE_Z1_A[3,]$mean
        ATE_A[j,]$SE_Z0_2 <- ate_SE_Z0_A[1,]$mean
        ATE_A[j,]$SE_Z0_1 <- ate_SE_Z0_A[2,]$mean
        ATE_A[j,]$SE_Z0_0 <- ate_SE_Z0_A[3,]$mean
      }
      
      
      # SAIPW
      mean_DE_G1_S <- as.vector(rowMeans(do.call(cbind, DE_G1_S)))
      mean_DE_G0_S <- as.vector(rowMeans(do.call(cbind, DE_G0_S)))
      mean_SE_Z1_S <- as.vector(rowMeans(do.call(cbind, SE_Z1_S)))
      mean_SE_Z0_S <- as.vector(rowMeans(do.call(cbind, SE_Z0_S)))
      
      df_S <- cbind(data, mean_DE_G1_S, mean_DE_G0_S, mean_SE_Z1_S, mean_SE_Z0_S)
      
      ate_DE_G1_S <- df_S %>%
        group_by(DE) %>%
        summarize(mean = mean(mean_DE_G1_S))
      ate_DE_G0_S <- df_S %>%
        group_by(DE) %>%
        summarize(mean = mean(mean_DE_G0_S))
      ate_SE_Z1_S <- df_S %>%
        group_by(SE) %>%
        summarize(mean = mean(mean_SE_Z1_S))
      ate_SE_Z0_S <- df_S %>%
        group_by(SE) %>%
        summarize(mean = mean(mean_SE_Z0_S))
      
      ATE_S[j,]$DE_G1_2 <- ate_DE_G1_S[1,]$mean
      ATE_S[j,]$DE_G1_1 <- ate_DE_G1_S[2,]$mean
      ATE_S[j,]$DE_G1_0 <- ate_DE_G1_S[3,]$mean
      ATE_S[j,]$DE_G0_2 <- ate_DE_G0_S[1,]$mean
      ATE_S[j,]$DE_G0_1 <- ate_DE_G0_S[2,]$mean
      ATE_S[j,]$DE_G0_0 <- ate_DE_G0_S[3,]$mean
      if(scenario == 2) {
        ATE_S[j,]$SE_Z1_1 <- ate_SE_Z1_S[1,]$mean
        ATE_S[j,]$SE_Z0_1 <- ate_SE_Z0_S[1,]$mean
      }
      if(scenario == 4) {
        ATE_S[j,]$SE_Z1_2 <- ate_SE_Z1_S[1,]$mean
        ATE_S[j,]$SE_Z1_1 <- ate_SE_Z1_S[2,]$mean
        ATE_S[j,]$SE_Z1_0 <- ate_SE_Z1_S[3,]$mean
        ATE_S[j,]$SE_Z0_2 <- ate_SE_Z0_S[1,]$mean
        ATE_S[j,]$SE_Z0_1 <- ate_SE_Z0_S[2,]$mean
        ATE_S[j,]$SE_Z0_0 <- ate_SE_Z0_S[3,]$mean
      }
      
    }
    cat("Iteration", j)
    .POSIXct(Sys.time(), "EST")
  }
  
  
  
  
  
  if(sample_prop == 1) {
  
    ###################################################
    ## ATE boxplots (full pop)
    ###################################################
    #data <- generate_po(dat_truncated_full, scenario = 2, var = 1, t = -1)
    
    ## G - ONLY RUN THIS BLOCK AFTER LOADING G-COMP RESULTS
    ATE_G_full <- data.frame(matrix(ncol=dimnames_len, nrow=length(DE_G1_G),
                                    dimnames=dimnames))
    
    for(i in 1:length(DE_G1_G)) {
      comb_G <- data.frame(DE = DE_G1_G[[i]]$DE, SE = SE_Z1_G[[i]]$SE,
                         DE_G1 = DE_G1_G[[i]]$est, DE_G0 = DE_G0_G[[i]]$est, SE_Z1 = SE_Z1_G[[i]]$est, SE_Z0 = SE_Z0_G[[i]]$est)
      
      table_DE_G1_G <- comb_G %>%
        group_by(DE) %>%
        summarize(mean = mean(DE_G1), .groups = "drop")
      table_DE_G0_G <- comb_G %>%
        group_by(DE) %>%
        summarize(mean = mean(DE_G0), .groups = "drop")
      table_SE_Z1_G <- comb_G %>%
        group_by(SE) %>%
        summarize(mean = mean(SE_Z1), .groups = "drop")
      table_SE_Z0_G <- comb_G %>%
        group_by(SE) %>%
        summarize(mean = mean(SE_Z0), .groups = "drop")
      ATE_G_full[i,]$DE_G1_2 <- table_DE_G1_G[1,]$mean
      ATE_G_full[i,]$DE_G1_1 <- table_DE_G1_G[2,]$mean
      ATE_G_full[i,]$DE_G1_0 <- table_DE_G1_G[3,]$mean
      ATE_G_full[i,]$DE_G0_2 <- table_DE_G0_G[1,]$mean
      ATE_G_full[i,]$DE_G0_1 <- table_DE_G0_G[2,]$mean
      ATE_G_full[i,]$DE_G0_0 <- table_DE_G0_G[3,]$mean
      if(scenario == 2) {
        ATE_G_full[i,]$SE_Z1_1 <- table_SE_Z1_G[1,]$mean
        ATE_G_full[i,]$SE_Z0_1 <- table_SE_Z0_G[1,]$mean
      }
      if(scenario == 4) {
        ATE_G_full[i,]$SE_Z1_2 <- table_SE_Z1_G[1,]$mean
        ATE_G_full[i,]$SE_Z1_1 <- table_SE_Z1_G[2,]$mean
        ATE_G_full[i,]$SE_Z1_0 <- table_SE_Z1_G[3,]$mean
        ATE_G_full[i,]$SE_Z0_2 <- table_SE_Z0_G[1,]$mean
        ATE_G_full[i,]$SE_Z0_1 <- table_SE_Z0_G[2,]$mean
        ATE_G_full[i,]$SE_Z0_0 <- table_SE_Z0_G[3,]$mean
      }
      if(scenario == 5) { 
        ATE_G_full[i,]$SE_Z1_1 <- table_SE_Z1_G[1,]$mean
        ATE_G_full[i,]$SE_Z1_0 <- table_SE_Z1_G[2,]$mean
        ATE_G_full[i,]$SE_Z0_1 <- table_SE_Z0_G[1,]$mean
        ATE_G_full[i,]$SE_Z0_0 <- table_SE_Z0_G[2,]$mean
      }
    }
    
    
    ## A - ONLY RUN THIS BLOCK AFTER LOADING AIPW RESULTS
    ATE_A_full <- data.frame(matrix(ncol=dimnames_len, nrow=length(DE_G1_A),
                                    dimnames=dimnames))
    
    for(i in 1:length(DE_G1_A)) {
      comb_A <- data.frame(DE = DE_G1_A[[i]]$DE, SE = SE_Z1_A[[i]]$SE,
                           DE_G1 = DE_G1_A[[i]]$est, DE_G0 = DE_G0_A[[i]]$est, SE_Z1 = SE_Z1_A[[i]]$est, SE_Z0 = SE_Z0_A[[i]]$est)
      table_DE_G1_A <- comb_A %>%
        group_by(DE) %>%
        #filter(DE_G1 > quantile(DE_G1, 0.005) & DE_G1 < quantile(DE_G1, 0.995)) %>%
        summarize(mean = mean(DE_G1), .groups = "drop")
      table_DE_G0_A <- comb_A %>%
        group_by(DE) %>%
        summarize(mean = mean(DE_G0), .groups = "drop")
      table_SE_Z1_A <- comb_A %>%
        group_by(SE) %>%
        #filter(SE_Z1 > quantile(SE_Z1, 0.005) & SE_Z1 < quantile(SE_Z1, 0.995)) %>%
        summarize(mean = mean(SE_Z1), .groups = "drop")
      table_SE_Z0_A <- comb_A %>%
        group_by(SE) %>%
        summarize(mean = mean(SE_Z0), .groups = "drop")
      ATE_A_full[i,]$DE_G1_2 <- table_DE_G1_A[1,]$mean
      ATE_A_full[i,]$DE_G1_1 <- table_DE_G1_A[2,]$mean
      ATE_A_full[i,]$DE_G1_0 <- table_DE_G1_A[3,]$mean
      ATE_A_full[i,]$DE_G0_2 <- table_DE_G0_A[1,]$mean
      ATE_A_full[i,]$DE_G0_1 <- table_DE_G0_A[2,]$mean
      ATE_A_full[i,]$DE_G0_0 <- table_DE_G0_A[3,]$mean
      if(scenario == 2) {
        ATE_A_full[i,]$SE_Z1_1 <- table_SE_Z1_A[1,]$mean
        ATE_A_full[i,]$SE_Z0_1 <- table_SE_Z0_A[1,]$mean
      }
      if(scenario == 4) {
        ATE_A_full[i,]$SE_Z1_2 <- table_SE_Z1_A[1,]$mean
        ATE_A_full[i,]$SE_Z1_1 <- table_SE_Z1_A[2,]$mean
        ATE_A_full[i,]$SE_Z1_0 <- table_SE_Z1_A[3,]$mean
        ATE_A_full[i,]$SE_Z0_2 <- table_SE_Z0_A[1,]$mean
        ATE_A_full[i,]$SE_Z0_1 <- table_SE_Z0_A[2,]$mean
        ATE_A_full[i,]$SE_Z0_0 <- table_SE_Z0_A[3,]$mean
      }
      if(scenario == 5) { 
        ATE_A_full[i,]$SE_Z1_1 <- table_SE_Z1_A[1,]$mean
        ATE_A_full[i,]$SE_Z1_0 <- table_SE_Z1_A[2,]$mean
        ATE_A_full[i,]$SE_Z0_1 <- table_SE_Z0_A[1,]$mean
        ATE_A_full[i,]$SE_Z0_0 <- table_SE_Z0_A[2,]$mean
      }
    }
    
    
    ## S - ONLY RUN THIS BLOCK AFTER LOADING SAIPW RESULTS
    ATE_S_full <- data.frame(matrix(ncol=dimnames_len, nrow=length(DE_G1_S),
                                    dimnames=dimnames))
    
    for(i in 1:length(DE_G1_S)) {
      comb_S <- data.frame(DE = DE_G1_S[[i]]$DE, SE = SE_Z1_S[[i]]$SE,
                           DE_G1 = DE_G1_S[[i]]$est, DE_G0 = DE_G0_S[[i]]$est, SE_Z1 = SE_Z1_S[[i]]$est, SE_Z0 = SE_Z0_S[[i]]$est)
      table_DE_G1_S <- comb_S %>%
        group_by(DE) %>%
        #filter(DE_G1 > quantile(DE_G1, 0.005) & DE_G1 < quantile(DE_G1, 0.995)) %>%
        summarize(mean = mean(DE_G1), .groups = "drop")
      table_DE_G0_S <- comb_S %>%
        group_by(DE) %>%
        summarize(mean = mean(DE_G0), .groups = "drop")
      table_SE_Z1_S <- comb_S %>%
        group_by(SE) %>%
        #filter(SE_Z1 > quantile(SE_Z1, 0.005) & SE_Z1 < quantile(SE_Z1, 0.995)) %>%
        summarize(mean = mean(SE_Z1), .groups = "drop")
      table_SE_Z0_S <- comb_S %>%
        group_by(SE) %>%
        summarize(mean = mean(SE_Z0), .groups = "drop")
      ATE_S_full[i,]$DE_G1_2 <- table_DE_G1_S[1,]$mean
      ATE_S_full[i,]$DE_G1_1 <- table_DE_G1_S[2,]$mean
      ATE_S_full[i,]$DE_G1_0 <- table_DE_G1_S[3,]$mean
      ATE_S_full[i,]$DE_G0_2 <- table_DE_G0_S[1,]$mean
      ATE_S_full[i,]$DE_G0_1 <- table_DE_G0_S[2,]$mean
      ATE_S_full[i,]$DE_G0_0 <- table_DE_G0_S[3,]$mean
      if(scenario == 2) {
        ATE_S_full[i,]$SE_Z1_1 <- table_SE_Z1_S[1,]$mean
        ATE_S_full[i,]$SE_Z0_1 <- table_SE_Z0_S[1,]$mean
      }
      if(scenario == 4) {
        ATE_S_full[i,]$SE_Z1_2 <- table_SE_Z1_S[1,]$mean
        ATE_S_full[i,]$SE_Z1_1 <- table_SE_Z1_S[2,]$mean
        ATE_S_full[i,]$SE_Z1_0 <- table_SE_Z1_S[3,]$mean
        ATE_S_full[i,]$SE_Z0_2 <- table_SE_Z0_S[1,]$mean
        ATE_S_full[i,]$SE_Z0_1 <- table_SE_Z0_S[2,]$mean
        ATE_S_full[i,]$SE_Z0_0 <- table_SE_Z0_S[3,]$mean
      }
      if(scenario == 5) { 
        ATE_S_full[i,]$SE_Z1_1 <- table_SE_Z1_S[1,]$mean
        ATE_S_full[i,]$SE_Z1_0 <- table_SE_Z1_S[2,]$mean
        ATE_S_full[i,]$SE_Z0_1 <- table_SE_Z0_S[1,]$mean
        ATE_S_full[i,]$SE_Z0_0 <- table_SE_Z0_S[2,]$mean
      }
    }
    
    
    ## putting it all together
    
    ATE_G_full$DE_G1_2 <- abs(ATE_G_full$DE_G1_2 - 2*effsize)/abs(effsize)
    ATE_G_full$DE_G1_1 <- abs(ATE_G_full$DE_G1_1 - 1*effsize)/abs(effsize)
    ATE_G_full$DE_G1_0 <- abs(ATE_G_full$DE_G1_0 - 0*effsize)/abs(effsize)
    ATE_G_full$DE_G0_2 <- abs(ATE_G_full$DE_G0_2 - 2*effsize)/abs(effsize)
    ATE_G_full$DE_G0_1 <- abs(ATE_G_full$DE_G0_1 - 1*effsize)/abs(effsize)
    ATE_G_full$DE_G0_0 <- abs(ATE_G_full$DE_G0_0 - 0*effsize)/abs(effsize)
    ATE_G_full$SE_Z1_1 <- abs(ATE_G_full$SE_Z1_1 - 1*effsize)/abs(effsize)
    ATE_G_full$SE_Z0_1 <- abs(ATE_G_full$SE_Z0_1 - 1*effsize)/abs(effsize)
    
    if(scenario == 4) {
      ATE_G_full$SE_Z1_2 <- abs(ATE_G_full$SE_Z1_2 - 2*effsize)/abs(effsize)
      ATE_G_full$SE_Z1_0 <- abs(ATE_G_full$SE_Z1_0 - 0*effsize)/abs(effsize)
      ATE_G_full$SE_Z0_2 <- abs(ATE_G_full$SE_Z0_2 - 2*effsize)/abs(effsize)
      ATE_G_full$SE_Z0_0 <- abs(ATE_G_full$SE_Z0_0 - 0*effsize)/abs(effsize)
    }
    
    if(scenario == 5) {
      ATE_G_full$SE_Z1_0 <- abs(ATE_G_full$SE_Z1_0 - 0*effsize)/abs(effsize)
      ATE_G_full$SE_Z0_0 <- abs(ATE_G_full$SE_Z0_0 - 0*effsize)/abs(effsize)
    }
    
    ATE_A_full$DE_G1_2 <- abs(ATE_A_full$DE_G1_2 - 2*effsize)/abs(effsize)
    ATE_A_full$DE_G1_1 <- abs(ATE_A_full$DE_G1_1 - 1*effsize)/abs(effsize)
    ATE_A_full$DE_G1_0 <- abs(ATE_A_full$DE_G1_0 - 0*effsize)/abs(effsize)
    ATE_A_full$DE_G0_2 <- abs(ATE_A_full$DE_G0_2 - 2*effsize)/abs(effsize)
    ATE_A_full$DE_G0_1 <- abs(ATE_A_full$DE_G0_1 - 1*effsize)/abs(effsize)
    ATE_A_full$DE_G0_0 <- abs(ATE_A_full$DE_G0_0 - 0*effsize)/abs(effsize)
    ATE_A_full$SE_Z1_1 <- abs(ATE_A_full$SE_Z1_1 - 1*effsize)/abs(effsize)
    ATE_A_full$SE_Z0_1 <- abs(ATE_A_full$SE_Z0_1 - 1*effsize)/abs(effsize)
    
    if(scenario == 4) {
      ATE_A_full$SE_Z1_2 <- abs(ATE_A_full$SE_Z1_2 - 2*effsize)/abs(effsize)
      ATE_A_full$SE_Z1_0 <- abs(ATE_A_full$SE_Z1_0 - 0*effsize)/abs(effsize)
      ATE_A_full$SE_Z0_2 <- abs(ATE_A_full$SE_Z0_2 - 2*effsize)/abs(effsize)
      ATE_A_full$SE_Z0_0 <- abs(ATE_A_full$SE_Z0_0 - 0*effsize)/abs(effsize)
    }
    
    if(scenario == 5) {
      ATE_A_full$SE_Z1_0 <- abs(ATE_A_full$SE_Z1_0 - 0*effsize)/abs(effsize)
      ATE_A_full$SE_Z0_0 <- abs(ATE_A_full$SE_Z0_0 - 0*effsize)/abs(effsize)
    }
    
    ATE_S_full$DE_G1_2 <- abs(ATE_S_full$DE_G1_2 - 2*effsize)/abs(effsize)
    ATE_S_full$DE_G1_1 <- abs(ATE_S_full$DE_G1_1 - 1*effsize)/abs(effsize)
    ATE_S_full$DE_G1_0 <- abs(ATE_S_full$DE_G1_0 - 0*effsize)/abs(effsize)
    ATE_S_full$DE_G0_2 <- abs(ATE_S_full$DE_G0_2 - 2*effsize)/abs(effsize)
    ATE_S_full$DE_G0_1 <- abs(ATE_S_full$DE_G0_1 - 1*effsize)/abs(effsize)
    ATE_S_full$DE_G0_0 <- abs(ATE_S_full$DE_G0_0 - 0*effsize)/abs(effsize)
    ATE_S_full$SE_Z1_1 <- abs(ATE_S_full$SE_Z1_1 - 1*effsize)/abs(effsize)
    ATE_S_full$SE_Z0_1 <- abs(ATE_S_full$SE_Z0_1 - 1*effsize)/abs(effsize)
    
    if(scenario == 4) {
      ATE_S_full$SE_Z1_2 <- abs(ATE_S_full$SE_Z1_2 - 2*effsize)/abs(effsize)
      ATE_S_full$SE_Z1_0 <- abs(ATE_S_full$SE_Z1_0 - 0*effsize)/abs(effsize)
      ATE_S_full$SE_Z0_2 <- abs(ATE_S_full$SE_Z0_2 - 2*effsize)/abs(effsize)
      ATE_S_full$SE_Z0_0 <- abs(ATE_S_full$SE_Z0_0 - 0*effsize)/abs(effsize)
    }
    
    if(scenario == 5) {
      ATE_S_full$SE_Z1_0 <- abs(ATE_S_full$SE_Z1_0 - 0*effsize)/abs(effsize)
      ATE_S_full$SE_Z0_0 <- abs(ATE_S_full$SE_Z0_0 - 0*effsize)/abs(effsize)
    }
    
    
    
    if(scenario == 2) {
      spilloverlist_G <- c(ATE_G_full$SE_Z1_1,
                           ATE_G_full$SE_Z0_1)
      spilloverlist_A <- c(ATE_A_full$SE_Z1_1,
                           ATE_A_full$SE_Z0_1)
      spilloverlist_S <- c(ATE_S_full$SE_Z1_1,
                           ATE_S_full$SE_Z0_1)
    }
    if(scenario == 4) {
      spilloverlist_G <- c(ATE_G_full$SE_Z1_2, ATE_G_full$SE_Z1_1,
                           ATE_G_full$SE_Z1_0,
                           ATE_G_full$SE_Z0_2, ATE_G_full$SE_Z0_1, ATE_G_full$SE_Z0_0)
      spilloverlist_A <- c(ATE_A_full$SE_Z1_2, ATE_A_full$SE_Z1_1, 
                           ATE_A_full$SE_Z1_0,
                           ATE_A_full$SE_Z0_2, ATE_A_full$SE_Z0_1, ATE_A_full$SE_Z0_0)
      spilloverlist_S <- c(ATE_S_full$SE_Z1_2, ATE_S_full$SE_Z1_1, 
                           ATE_S_full$SE_Z1_0,
                           ATE_S_full$SE_Z0_2, ATE_S_full$SE_Z0_1, ATE_S_full$SE_Z0_0)
    }
    if(scenario == 5) {
      spilloverlist_G <- c(ATE_G_full$SE_Z1_1, ATE_G_full$SE_Z1_0,
                           ATE_G_full$SE_Z0_1, ATE_G_full$SE_Z0_0)
      spilloverlist_A <- c(ATE_A_full$SE_Z1_1, ATE_A_full$SE_Z1_0,
                           ATE_A_full$SE_Z0_1, ATE_A_full$SE_Z0_0)
      spilloverlist_S <- c(ATE_S_full$SE_Z1_1, ATE_S_full$SE_Z1_0,
                           ATE_S_full$SE_Z0_1, ATE_S_full$SE_Z0_0)
    }
    
    
    ATE_G_full_df_DE <- data.frame("Direct" = c(ATE_G_full$DE_G1_2, ATE_G_full$DE_G1_1, 
                                                ATE_G_full$DE_G1_0,
                                             ATE_G_full$DE_G0_2, ATE_G_full$DE_G0_1, ATE_G_full$DE_G0_0))
    ATE_G_full_df_DE$group <- rep('G-comp', nrow(ATE_G_full_df_DE))
    ATE_G_full_df_DE <- melt(ATE_G_full_df_DE, id.vars = 'group')
    ATE_G_full_df_SE <- data.frame("Spillover" = spilloverlist_G)
    ATE_G_full_df_SE$group <- rep('G-comp', nrow(ATE_G_full_df_SE))
    ATE_G_full_df_SE <- melt(ATE_G_full_df_SE, id.vars = 'group')
    ATE_G_full_df <- rbind(ATE_G_full_df_DE, ATE_G_full_df_SE)
    
    
    ATE_A_full_df_DE <- data.frame("Direct" = c(ATE_A_full$DE_G1_2, ATE_A_full$DE_G1_1, 
                                                ATE_A_full$DE_G1_0,
                                                ATE_A_full$DE_G0_2, ATE_A_full$DE_G0_1, ATE_A_full$DE_G0_0))
    ATE_A_full_df_DE$group <- rep('AIPW', nrow(ATE_A_full_df_DE))
    ATE_A_full_df_DE <- melt(ATE_A_full_df_DE, id.vars = 'group')
    ATE_A_full_df_SE <- data.frame("Spillover" = spilloverlist_A)
    ATE_A_full_df_SE$group <- rep('AIPW', nrow(ATE_A_full_df_SE))
    ATE_A_full_df_SE <- melt(ATE_A_full_df_SE, id.vars = 'group')
    ATE_A_full_df <- rbind(ATE_A_full_df_DE, ATE_A_full_df_SE)
    
    
    ATE_S_full_df_DE <- data.frame("Direct" = c(ATE_S_full$DE_G1_2, ATE_S_full$DE_G1_1, 
                                                ATE_S_full$DE_G1_0,
                                                ATE_S_full$DE_G0_2, ATE_S_full$DE_G0_1, ATE_S_full$DE_G0_0))
    ATE_S_full_df_DE$group <- rep('SAIPW', nrow(ATE_S_full_df_DE))
    ATE_S_full_df_DE <- melt(ATE_S_full_df_DE, id.vars = 'group')
    ATE_S_full_df_SE <- data.frame("Spillover" = spilloverlist_S)
    ATE_S_full_df_SE$group <- rep('SAIPW', nrow(ATE_S_full_df_SE))
    ATE_S_full_df_SE <- melt(ATE_S_full_df_SE, id.vars = 'group')
    ATE_S_full_df <- rbind(ATE_S_full_df_DE, ATE_S_full_df_SE)
    
    allData_full <- rbind(ATE_G_full_df, ATE_A_full_df, ATE_S_full_df)
    allData_full$group <- factor(allData_full$group, levels = c("G-comp", "AIPW", "SAIPW"))
    colnames(allData_full) <- c("group", "variable", "bias")
    
    #save(data, comb, ATE_A_full, ATE_A_full_df_DE, file = "test.rda")
    save(data, allData_full, ATE_G_full_df, ATE_A_full_df, ATE_S_full_df, file = output_filename)
  } else {
    save(ATE_G, ATE_A, ATE_S, file = output_filename)
  }
  print("Finished:")
  .POSIXct(Sys.time(), "EST")
}








###############################################################################################################
## Simulation parameters
## When compiling results from different scenarios... ONLY RUN ONE BLOCK AT A TIME
###############################################################################################################

num_iter <- 1
reps <- 10000
sample_prop <- 1





#############################################################################################################################################
## Misspecification Scenarios
#############################################################################################################################################

run_sims(num_iter, reps, sample_prop, scenario = 4, variance = 1, effsize = -1,
         prop_type = "corr", outcome_type = "corr", output_filename = "out_scenA_10k.Rda")
run_sims(num_iter, reps, sample_prop, scenario = 4, variance = 1, effsize = -1,
         prop_type = "miss", outcome_type = "corr", output_filename = "out_scenB_10k.Rda")
run_sims(num_iter, reps, sample_prop, scenario = 4, variance = 1, effsize = -1,
         prop_type = "corr", outcome_type = "miss", output_filename = "out_scenC_10k.Rda")
run_sims(num_iter, reps, sample_prop, scenario = 4, variance = 1, effsize = -1,
         prop_type = "miss", outcome_type = "miss", output_filename = "out_scenD_10k.Rda")


##### Different misspecification scenarios, var 1, effect size 0-2
## load scen A results
allData_scenA_DE <- allData_full[allData_full$variable == "Direct",]
allData_scenA_DE$variable <- rep("Scenario A", nrow(allData_scenA_DE))
allData_scenA_DE$variable <- as.factor(allData_scenA_DE$variable)
allData_scenA_SE <- allData_full[allData_full$variable == "Spillover",]
allData_scenA_SE$variable <- rep("Scenario A", nrow(allData_scenA_SE))
allData_scenA_SE$variable <- as.factor(allData_scenA_SE$variable)

## load scen B results
allData_scenB_DE <- allData_full[allData_full$variable == "Direct",]
allData_scenB_DE$variable <- rep("Scenario B", nrow(allData_scenB_DE))
allData_scenB_DE$variable <- as.factor(allData_scenB_DE$variable)
allData_scenB_SE <- allData_full[allData_full$variable == "Spillover",]
allData_scenB_SE$variable <- rep("Scenario B", nrow(allData_scenB_SE))
allData_scenB_SE$variable <- as.factor(allData_scenB_SE$variable)

## load scen C results
allData_scenC_DE <- allData_full[allData_full$variable == "Direct",]
allData_scenC_DE$variable <- rep("Scenario C", nrow(allData_scenC_DE))
allData_scenC_DE$variable <- as.factor(allData_scenC_DE$variable)
allData_scenC_SE <- allData_full[allData_full$variable == "Spillover",]
allData_scenC_SE$variable <- rep("Scenario C", nrow(allData_scenC_SE))
allData_scenC_SE$variable <- as.factor(allData_scenC_SE$variable)

## load scen D results
allData_scenD_DE <- allData_full[allData_full$variable == "Direct",]
allData_scenD_DE$variable <- rep("Scenario D", nrow(allData_scenD_DE))
allData_scenD_DE$variable <- as.factor(allData_scenD_DE$variable)
allData_scenD_SE <- allData_full[allData_full$variable == "Spillover",]
allData_scenD_SE$variable <- rep("Scenario D", nrow(allData_scenD_SE))
allData_scenD_SE$variable <- as.factor(allData_scenD_SE$variable)


all_combined_scen_DE <- rbind(allData_scenA_DE, allData_scenB_DE, allData_scenC_DE, allData_scenD_DE)
all_combined_scen_DE$bias <- 100*all_combined_scen_DE$bias
all_combined_scen_SE <- rbind(allData_scenA_SE, allData_scenB_SE, allData_scenC_SE, allData_scenD_SE)
all_combined_scen_SE$bias <- 100*all_combined_scen_SE$bias


ggplot(all_combined_scen_DE, aes(x = variable, y = bias, fill = group)) +
  geom_boxplot(outlier.size = 0.05, outlier.alpha = 0.05) +
  coord_cartesian(ylim=c(0,475)) +
  geom_hline(yintercept = 0 ,alpha=0.4) +
  facet_wrap( ~ variable, scales = 'free_x', ncol = 7) +
  theme_bw() + 
  theme(
    axis.title.x = element_text(vjust = 0, size = 12),
    axis.title.y = element_text(vjust = 2, size = 12),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  ) +
  xlab("Misspecification Scenario") + 
  ylab("Percent Absolute Bias (%)") +
  labs(fill = "Estimator") +
  ggtitle("Percent Absolute Bias of ATE Estimates (Direct Effects)", subtitle = "Under Various Misspecification Scenarios")
ggsave("bias_misspec_DE.png", width = 11, height = 7, dpi = 700)

ggplot(all_combined_scen_SE, aes(x = variable, y = bias, fill = group)) +
  geom_boxplot(outlier.size = 0.05, outlier.alpha = 0.05) +
  coord_cartesian(ylim=c(0,325)) +
  geom_hline(yintercept = 0 ,alpha=0.4) +
  facet_wrap( ~ variable, scales = 'free_x', ncol = 7) +
  theme_bw() + 
  theme(
    axis.title.x = element_text(vjust = 0, size = 12),
    axis.title.y = element_text(vjust = 2, size = 12),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  ) +
  xlab("Misspecification Scenario") + 
  ylab("Percent Absolute Bias (%)") +
  labs(fill = "Estimator") +
  ggtitle("Percent Absolute Bias of ATE Estimates (Spillover Effects)", subtitle = "Under Various Misspecification Scenarios")
ggsave("bias_misspec_SE.png", width = 11, height = 7, dpi = 700)





#############################################################################################################################################
## Outcome error variance
#############################################################################################################################################

run_sims(num_iter, reps, sample_prop, scenario = 4, variance = 0.2, effsize = -1,
         prop_type = "corr", outcome_type = "corr", output_filename = "out_var02_10k.Rda")
run_sims(num_iter, reps, sample_prop, scenario = 4, variance = 5, effsize = -1,
         prop_type = "corr", outcome_type = "corr", output_filename = "out_var5_10k.Rda")
run_sims(num_iter, reps, sample_prop, scenario = 4, variance = 10, effsize = -1,
         prop_type = "corr", outcome_type = "corr", output_filename = "out_var10_10k.Rda")


##### Different outcome variance, effect size 0-2
## var = 0.2
allData_scenA_var0_2_DE <- allData_full[allData_full$variable == "Direct",]
allData_scenA_var0_2_DE$variable <- rep("var = 0.2", nrow(allData_scenA_var0_2_DE))
allData_scenA_var0_2_DE$variable <- as.factor(allData_scenA_var0_2_DE$variable)
allData_scenA_var0_2_SE <- allData_full[allData_full$variable == "Spillover",]
allData_scenA_var0_2_SE$variable <- rep("var = 0.2", nrow(allData_scenA_var0_2_SE))
allData_scenA_var0_2_SE$variable <- as.factor(allData_scenA_var0_2_SE$variable)

## var = 1
allData_scenA_var1_DE <- allData_full[allData_full$variable == "Direct",]
allData_scenA_var1_DE$variable <- rep("var = 1", nrow(allData_scenA_var1_DE))
allData_scenA_var1_DE$variable <- as.factor(allData_scenA_var1_DE$variable)
allData_scenA_var1_SE <- allData_full[allData_full$variable == "Spillover",]
allData_scenA_var1_SE$variable <- rep("var = 1", nrow(allData_scenA_var1_SE))
allData_scenA_var1_SE$variable <- as.factor(allData_scenA_var1_SE$variable)

## var = 5
allData_scenA_var5_DE <- allData_full[allData_full$variable == "Direct",]
allData_scenA_var5_DE$variable <- rep("var = 5", nrow(allData_scenA_var5_DE))
allData_scenA_var5_DE$variable <- as.factor(allData_scenA_var5_DE$variable)
allData_scenA_var5_SE <- allData_full[allData_full$variable == "Spillover",]
allData_scenA_var5_SE$variable <- rep("var = 5", nrow(allData_scenA_var5_SE))
allData_scenA_var5_SE$variable <- as.factor(allData_scenA_var5_SE$variable)

## var = 10
allData_scenA_var10_DE <- allData_full[allData_full$variable == "Direct",]
allData_scenA_var10_DE$variable <- rep("var = 10", nrow(allData_scenA_var10_DE))
allData_scenA_var10_DE$variable <- as.factor(allData_scenA_var10_DE$variable)
allData_scenA_var10_SE <- allData_full[allData_full$variable == "Spillover",]
allData_scenA_var10_SE$variable <- rep("var = 10", nrow(allData_scenA_var10_SE))
allData_scenA_var10_SE$variable <- as.factor(allData_scenA_var10_SE$variable)


all_combined_scen_DE <- rbind(allData_scenA_var0_2_DE, allData_scenA_var1_DE, allData_scenA_var5_DE, allData_scenA_var10_DE)
all_combined_scen_DE$bias <- 100*all_combined_scen_DE$bias
all_combined_scen_SE <- rbind(allData_scenA_var0_2_SE, allData_scenA_var1_SE, allData_scenA_var5_SE, allData_scenA_var10_SE)
all_combined_scen_SE$bias <- 100*all_combined_scen_SE$bias


ggplot(all_combined_scen_DE, aes(x = variable, y = bias, fill = group)) +
  geom_boxplot(outlier.size = 0.05, outlier.alpha = 0.05) +
  coord_cartesian(ylim=c(0,475)) +
  geom_hline(yintercept = 0 ,alpha=0.4) +
  facet_wrap( ~ variable, scales = 'free_x', ncol = 7) +
  theme_bw() + 
  theme(
    axis.title.x = element_text(vjust = 0, size = 12),
    axis.title.y = element_text(vjust = 2, size = 12),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  ) +
  xlab("Outcome Model Error Variance") + 
  ylab("Percent Absolute Bias (%)") +
  labs(fill = "Estimator") +
  ggtitle("Percent Absolute Bias of ATE Estimates (Direct Effects)", subtitle = "Varying Outcome Model Error Variance")
ggsave("bias_var_DE.png", width = 11, height = 7, dpi = 700)

ggplot(all_combined_scen_SE, aes(x = variable, y = bias, fill = group)) +
  geom_boxplot(outlier.size = 0.05, outlier.alpha = 0.05) +
  coord_cartesian(ylim=c(0,475)) +
  geom_hline(yintercept = 0 ,alpha=0.4) +
  facet_wrap( ~ variable, scales = 'free_x', ncol = 7) +
  theme_bw() + 
  theme(
    axis.title.x = element_text(vjust = 0, size = 12),
    axis.title.y = element_text(vjust = 2, size = 12),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  ) +
  xlab("Outcome Model Error Variance") + 
  ylab("Percent Absolute Bias (%)") +
  labs(fill = "Estimator") +
  ggtitle("Percent Absolute Bias of ATE Estimates (Spillover Effects)", subtitle = "Varying Outcome Model Error Variance")
ggsave("bias_var_SE.png", width = 11, height = 7, dpi = 700)





#############################################################################################################################################
## Effect sizes
#############################################################################################################################################

run_sims(num_iter, reps, sample_prop, scenario = 4, variance = 1, effsize = -5,
         prop_type = "corr", outcome_type = "corr", output_filename = "out_eff5_10k.Rda")
run_sims(num_iter, reps, sample_prop, scenario = 4, variance = 1, effsize = -10,
         prop_type = "corr", outcome_type = "corr", output_filename = "out_eff10_10k.Rda")


##### Different effect sizes, var 1
## effect size = -1
allData_scenA_eff1_DE <- allData_full[allData_full$variable == "Direct",]
allData_scenA_eff1_DE$variable <- rep("effsize = -1", nrow(allData_scenA_eff1_DE))
allData_scenA_eff1_DE$variable <- as.factor(allData_scenA_eff1_DE$variable)
allData_scenA_eff1_SE <- allData_full[allData_full$variable == "Spillover",]
allData_scenA_eff1_SE$variable <- rep("effsize = -1", nrow(allData_scenA_eff1_SE))
allData_scenA_eff1_SE$variable <- as.factor(allData_scenA_eff1_SE$variable)

## effect size = -5
allData_scenA_eff5_DE <- allData_full[allData_full$variable == "Direct",]
allData_scenA_eff5_DE$variable <- rep("effsize = -5", nrow(allData_scenA_eff5_DE))
allData_scenA_eff5_DE$variable <- as.factor(allData_scenA_eff5_DE$variable)
allData_scenA_eff5_SE <- allData_full[allData_full$variable == "Spillover",]
allData_scenA_eff5_SE$variable <- rep("effsize = -5", nrow(allData_scenA_eff5_SE))
allData_scenA_eff5_SE$variable <- as.factor(allData_scenA_eff5_SE$variable)

## effect size = -10
allData_scenA_eff10_DE <- allData_full[allData_full$variable == "Direct",]
allData_scenA_eff10_DE$variable <- rep("effsize = -10", nrow(allData_scenA_eff10_DE))
allData_scenA_eff10_DE$variable <- as.factor(allData_scenA_eff10_DE$variable)
allData_scenA_eff10_SE <- allData_full[allData_full$variable == "Spillover",]
allData_scenA_eff10_SE$variable <- rep("effsize = -10", nrow(allData_scenA_eff10_SE))
allData_scenA_eff10_SE$variable <- as.factor(allData_scenA_eff10_SE$variable)


all_combined_scen_DE <- rbind(allData_scenA_eff1_DE, allData_scenA_eff5_DE, allData_scenA_eff10_DE)
all_combined_scen_DE$bias <- 100*all_combined_scen_DE$bias
all_combined_scen_SE <- rbind(allData_scenA_eff1_SE, allData_scenA_eff5_SE, allData_scenA_eff10_SE)
all_combined_scen_SE$bias <- 100*all_combined_scen_SE$bias


ggplot(all_combined_scen_DE, aes(x = variable, y = bias, fill = group)) +
  geom_boxplot(outlier.size = 0.05, outlier.alpha = 0.05) +
  coord_cartesian(ylim=c(0,125)) +
  geom_hline(yintercept = 0 ,alpha=0.4) +
  facet_wrap( ~ variable, scales = 'free_x', ncol = 7) +
  theme_bw() + 
  theme(
    axis.title.x = element_text(vjust = 0, size = 12),
    axis.title.y = element_text(vjust = 2, size = 12),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  ) +
  xlab("Population Average Treatment Effect") + 
  ylab("Percent Absolute Bias (%)") +
  labs(fill = "Estimator") +
  ggtitle("Percent Absolute Bias of ATE Estimates (Direct Effects)", subtitle = "Varying Population Average Treatment Effect Size")
ggsave("bias_effsize_DE.png", width = 11, height = 7, dpi = 700)

ggplot(all_combined_scen_SE, aes(x = variable, y = bias, fill = group)) +
  geom_boxplot(outlier.size = 0.05, outlier.alpha = 0.05) +
  coord_cartesian(ylim=c(0,125)) +
  geom_hline(yintercept = 0 ,alpha=0.4) +
  facet_wrap( ~ variable, scales = 'free_x', ncol = 7) +
  theme_bw() + 
  theme(
    axis.title.x = element_text(vjust = 0, size = 12),
    axis.title.y = element_text(vjust = 2, size = 12),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  ) +
  xlab("Population Average Treatment Effect") + 
  ylab("Percent Absolute Bias (%)") +
  labs(fill = "Estimator") +
  ggtitle("Percent Absolute Bias of ATE Estimates (Spillover Effects)", subtitle = "Varying Population Average Treatment Effect Size")
ggsave("bias_effsize_SE.png", width = 11, height = 7, dpi = 700)






ggplot(all_combined_scen_DE, aes(x = variable, y = bias, fill = group)) +
  geom_boxplot(outlier.size = 0.05, outlier.alpha = 0.05) +
  coord_cartesian(ylim=c(0,125)) +
  geom_hline(yintercept = 0 ,alpha=0.4) +
  facet_wrap( ~ variable, scales = 'free_x', ncol = 7) +
  theme_bw() + 
  theme(
    axis.title.x = element_text(vjust = 0, size = 12),
    axis.title.y = element_text(vjust = 2, size = 12),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  ) +
  xlab("Population Average Treatment Effect") + 
  ylab("Percent Absolute Bias (%)") +
  labs(fill = "Estimator") +
  ggtitle("Percent Absolute Bias of ATE Estimates (Direct Effects)", subtitle = "Varying Population Average Treatment Effect Size")































##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################















########################################################################################
## Sample size boxplots
########################################################################################

num_iter <- 100
reps <- 10000


run_sims(num_iter, reps, 0.005, scenario = 4, variance = 1, effsize = -1,
         prop_type = "corr", outcome_type = "corr", output_filename = "out_ss'005_1k_test.Rda")
run_sims(num_iter, reps, 0.01, scenario = 4, variance = 1, effsize = -1,
         prop_type = "corr", outcome_type = "corr", output_filename = "out_ss0'1_1k_test.Rda")
run_sims(num_iter, reps, 0.03, scenario = 4, variance = 1, effsize = -1,
         prop_type = "corr", outcome_type = "corr", output_filename = "out_ss0'3_1k_test.Rda")
run_sims(num_iter, reps, 0.05, scenario = 4, variance = 1, effsize = -1,
         prop_type = "corr", outcome_type = "corr", output_filename = "out_ss0'5_1k_test.Rda")
run_sims(num_iter, reps, 0.1, scenario = 4, variance = 1, effsize = -1,
         prop_type = "corr", outcome_type = "corr", output_filename = "out_ss1_1k_test.Rda")
run_sims(num_iter, reps, 0.2, scenario = 4, variance = 1, effsize = -1,
         prop_type = "corr", outcome_type = "corr", output_filename = "out_ss2_1k_test.Rda")
run_sims(num_iter, reps, 0.5, scenario = 4, variance = 1, effsize = -1,
         prop_type = "corr", outcome_type = "corr", output_filename = "out_ss5_1k_test.Rda")




load("out_ss'005_500_test.Rda")
ATE_G_05 <- ATE_G
ATE_A_05 <- ATE_A
ATE_S_05 <- ATE_S

load("out_ss0'1_500_test.Rda")
ATE_G_1 <- ATE_G
ATE_A_1 <- ATE_A
ATE_S_1 <- ATE_S

load("out_ss0'3_500_test.Rda")
ATE_G_3 <- ATE_G
ATE_A_3 <- ATE_A
ATE_S_3 <- ATE_S

load("out_ss0'5_500_test.Rda")
ATE_G_5 <- ATE_G
ATE_A_5 <- ATE_A
ATE_S_5 <- ATE_S

load("out_ss1_500_test.Rda")
ATE_G_10 <- ATE_G
ATE_A_10 <- ATE_A
ATE_S_10 <- ATE_S

load("out_ss2_500_test.Rda")
ATE_G_20 <- ATE_G
ATE_A_20 <- ATE_A
ATE_S_20 <- ATE_S

load("out_ss5_500_test.Rda")
ATE_G_50 <- ATE_G
ATE_A_50 <- ATE_A
ATE_S_50 <- ATE_S


## DE_G1_2

DEG1_2_G <- data.frame("0.5" = ATE_G_50$DE_G1_2, "0.2" = ATE_G_20$DE_G1_2, "0.1" = ATE_G_10$DE_G1_2,
                       "0.05" = ATE_G_5$DE_G1_2, "0.03" = ATE_G_3$DE_G1_2, "0.01" = ATE_G_1$DE_G1_2, "0.005" = ATE_G_05$DE_G1_2)
names(DEG1_2_G) <- sub("^X", "", names(DEG1_2_G))
DEG1_2_G$group <- rep('G-comp', nrow(DEG1_2_G))
DEG1_2_G <- melt(DEG1_2_G, id.vars = 'group')


DEG1_2_A <- data.frame("0.5" = ATE_A_50$DE_G1_2, "0.2" = ATE_A_20$DE_G1_2, "0.1" = ATE_A_10$DE_G1_2,
                       "0.05" = ATE_A_5$DE_G1_2, "0.03" = ATE_A_3$DE_G1_2, "0.01" = ATE_A_1$DE_G1_2, "0.005" = ATE_A_05$DE_G1_2)
names(DEG1_2_A) <- sub("^X", "", names(DEG1_2_A))
DEG1_2_A$group <- rep('AIPW', nrow(DEG1_2_A))
DEG1_2_A <- melt(DEG1_2_A, id.vars = 'group')

DEG1_2_S <- data.frame("0.5" = ATE_S_50$DE_G1_2, "0.2" = ATE_S_20$DE_G1_2, "0.1" = ATE_S_10$DE_G1_2,
                       "0.05" = ATE_S_5$DE_G1_2, "0.03" = ATE_S_3$DE_G1_2, "0.01" = ATE_S_1$DE_G1_2, "0.005" = ATE_S_05$DE_G1_2)
names(DEG1_2_S) <- sub("^X", "", names(DEG1_2_S))
DEG1_2_S$group <- rep('SAIPW', nrow(DEG1_2_S))
DEG1_2_S <- melt(DEG1_2_S, id.vars = 'group')


DEG1_2_allData <- rbind(DEG1_2_G, DEG1_2_A, DEG1_2_S)
DEG1_2_allData$group <- factor(DEG1_2_allData$group, levels = c("G-comp", "AIPW", "SAIPW"))

ggplot(DEG1_2_allData, aes(x = variable, y = value, fill = group)) +
  geom_boxplot() +
  coord_cartesian(ylim=c(-3.3,0.7)) +
  geom_hline(yintercept = -2 ,alpha=0.4) +
  facet_wrap( ~ variable, scales = 'free_x', ncol = 6) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  xlab("Sample Size") +
  ylab("ATE Estimate") +
  labs(fill = "Estimator")



## DE_G1_1

DEG1_1_G <- data.frame("0.5" = ATE_G_50$DE_G1_1, "0.2" = ATE_G_20$DE_G1_1, "0.1" = ATE_G_10$DE_G1_1,
                       "0.05" = ATE_G_5$DE_G1_1, "0.03" = ATE_G_3$DE_G1_1, "0.01" = ATE_G_1$DE_G1_1, "0.005" = ATE_G_05$DE_G1_1)
names(DEG1_1_G) <- sub("^X", "", names(DEG1_1_G))
DEG1_1_G$group <- rep('G-comp', nrow(DEG1_1_G))
DEG1_1_G <- melt(DEG1_1_G, id.vars = 'group')


DEG1_1_A <- data.frame("0.5" = ATE_A_50$DE_G1_1, "0.2" = ATE_A_20$DE_G1_1, "0.1" = ATE_A_10$DE_G1_1,
                       "0.05" = ATE_A_5$DE_G1_1, "0.03" = ATE_A_3$DE_G1_1, "0.01" = ATE_A_1$DE_G1_1, "0.005" = ATE_A_05$DE_G1_1)
names(DEG1_1_A) <- sub("^X", "", names(DEG1_1_A))
DEG1_1_A$group <- rep('AIPW', nrow(DEG1_1_A))
DEG1_1_A <- melt(DEG1_1_A, id.vars = 'group')

DEG1_1_S <- data.frame("0.5" = ATE_S_50$DE_G1_1, "0.2" = ATE_S_20$DE_G1_1, "0.1" = ATE_S_10$DE_G1_1,
                       "0.05" = ATE_S_5$DE_G1_1, "0.03" = ATE_S_3$DE_G1_1, "0.01" = ATE_S_1$DE_G1_1, "0.005" = ATE_S_05$DE_G1_1)
names(DEG1_1_S) <- sub("^X", "", names(DEG1_1_S))
DEG1_1_S$group <- rep('SAIPW', nrow(DEG1_1_S))
DEG1_1_S <- melt(DEG1_1_S, id.vars = 'group')


DEG1_1_allData <- rbind(DEG1_1_G, DEG1_1_A, DEG1_1_S)
DEG1_1_allData$group <- factor(DEG1_1_allData$group, levels = c("G-comp", "AIPW", "SAIPW"))

ggplot(DEG1_1_allData, aes(x = variable, y = value, fill = group)) +
  geom_boxplot() +
  coord_cartesian(ylim=c(-2.9,1.1)) +
  geom_hline(yintercept = -1 ,alpha=0.4) +
  facet_wrap( ~ variable, scales = 'free_x', ncol = 6) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  xlab("Sample Size") +
  ylab("ATE Estimate") +
  labs(fill = "Estimator")



## DE_G1_0

DEG1_0_G <- data.frame("0.5" = ATE_G_50$DE_G1_0, "0.2" = ATE_G_20$DE_G1_0, "0.1" = ATE_G_10$DE_G1_0,
                       "0.05" = ATE_G_5$DE_G1_0, "0.03" = ATE_G_3$DE_G1_0, "0.01" = ATE_G_1$DE_G1_0, "0.005" = ATE_G_05$DE_G1_0)
names(DEG1_0_G) <- sub("^X", "", names(DEG1_0_G))
DEG1_0_G$group <- rep('G-comp', nrow(DEG1_0_G))
DEG1_0_G <- melt(DEG1_0_G, id.vars = 'group')


DEG1_0_A <- data.frame("0.5" = ATE_A_50$DE_G1_0, "0.2" = ATE_A_20$DE_G1_0, "0.1" = ATE_A_10$DE_G1_0,
                       "0.05" = ATE_A_5$DE_G1_0, "0.03" = ATE_A_3$DE_G1_0, "0.01" = ATE_A_1$DE_G1_0, "0.005" = ATE_A_05$DE_G1_0)
names(DEG1_0_A) <- sub("^X", "", names(DEG1_0_A))
DEG1_0_A$group <- rep('AIPW', nrow(DEG1_0_A))
DEG1_0_A <- melt(DEG1_0_A, id.vars = 'group')

DEG1_0_S <- data.frame("0.5" = ATE_S_50$DE_G1_0, "0.2" = ATE_S_20$DE_G1_0, "0.1" = ATE_S_10$DE_G1_0,
                       "0.05" = ATE_S_5$DE_G1_0, "0.03" = ATE_S_3$DE_G1_0, "0.01" = ATE_S_1$DE_G1_0, "0.005" = ATE_S_05$DE_G1_0)
names(DEG1_0_S) <- sub("^X", "", names(DEG1_0_S))
DEG1_0_S$group <- rep('SAIPW', nrow(DEG1_0_S))
DEG1_0_S <- melt(DEG1_0_S, id.vars = 'group')


DEG1_0_allData <- rbind(DEG1_0_G, DEG1_0_A, DEG1_0_S)
DEG1_0_allData$group <- factor(DEG1_0_allData$group, levels = c("G-comp", "AIPW", "SAIPW"))

ggplot(DEG1_0_allData, aes(x = variable, y = value, fill = group)) +
  geom_boxplot() +
  coord_cartesian(ylim=c(-2.7,1.3)) +
  geom_hline(yintercept = 0 ,alpha=0.4) +
  facet_wrap( ~ variable, scales = 'free_x', ncol = 6) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  xlab("Sample Size") +
  ylab("ATE Estimate") +
  labs(fill = "Estimator")



## DE_G0_2

DEG0_2_G <- data.frame("0.5" = ATE_G_50$DE_G0_2, "0.2" = ATE_G_20$DE_G0_2, "0.1" = ATE_G_10$DE_G0_2,
                       "0.05" = ATE_G_5$DE_G0_2, "0.03" = ATE_G_3$DE_G0_2, "0.01" = ATE_G_1$DE_G0_2, "0.005" = ATE_G_05$DE_G0_2)
names(DEG0_2_G) <- sub("^X", "", names(DEG0_2_G))
DEG0_2_G$group <- rep('G-comp', nrow(DEG0_2_G))
DEG0_2_G <- melt(DEG0_2_G, id.vars = 'group')


DEG0_2_A <- data.frame("0.5" = ATE_A_50$DE_G0_2, "0.2" = ATE_A_20$DE_G0_2, "0.1" = ATE_A_10$DE_G0_2,
                       "0.05" = ATE_A_5$DE_G0_2, "0.03" = ATE_A_3$DE_G0_2, "0.01" = ATE_A_1$DE_G0_2, "0.005" = ATE_A_05$DE_G0_2)
names(DEG0_2_A) <- sub("^X", "", names(DEG0_2_A))
DEG0_2_A$group <- rep('AIPW', nrow(DEG0_2_A))
DEG0_2_A <- melt(DEG0_2_A, id.vars = 'group')

DEG0_2_S <- data.frame("0.5" = ATE_S_50$DE_G0_2, "0.2" = ATE_S_20$DE_G0_2, "0.1" = ATE_S_10$DE_G0_2,
                       "0.05" = ATE_S_5$DE_G0_2, "0.03" = ATE_S_3$DE_G0_2, "0.01" = ATE_S_1$DE_G0_2, "0.005" = ATE_S_05$DE_G0_2)
names(DEG0_2_S) <- sub("^X", "", names(DEG0_2_S))
DEG0_2_S$group <- rep('SAIPW', nrow(DEG0_2_S))
DEG0_2_S <- melt(DEG0_2_S, id.vars = 'group')


DEG0_2_allData <- rbind(DEG0_2_G, DEG0_2_A, DEG0_2_S)
DEG0_2_allData$group <- factor(DEG0_2_allData$group, levels = c("G-comp", "AIPW", "SAIPW"))

ggplot(DEG0_2_allData, aes(x = variable, y = value, fill = group)) +
  geom_boxplot() +
  coord_cartesian(ylim=c(-4,0)) +
  geom_hline(yintercept = -2 ,alpha=0.4) +
  facet_wrap( ~ variable, scales = 'free_x', ncol = 6) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  xlab("Sample Size") +
  ylab("ATE Estimate") +
  labs(fill = "Estimator")



## DE_G0_1

DEG0_1_G <- data.frame("0.5" = ATE_G_50$DE_G0_1, "0.2" = ATE_G_20$DE_G0_1, "0.1" = ATE_G_10$DE_G0_1,
                       "0.05" = ATE_G_5$DE_G0_1, "0.03" = ATE_G_3$DE_G0_1, "0.01" = ATE_G_1$DE_G0_1, "0.005" = ATE_G_05$DE_G0_1)
names(DEG0_1_G) <- sub("^X", "", names(DEG0_1_G))
DEG0_1_G$group <- rep('G-comp', nrow(DEG0_1_G))
DEG0_1_G <- melt(DEG0_1_G, id.vars = 'group')


DEG0_1_A <- data.frame("0.5" = ATE_A_50$DE_G0_1, "0.2" = ATE_A_20$DE_G0_1, "0.1" = ATE_A_10$DE_G0_1,
                       "0.05" = ATE_A_5$DE_G0_1, "0.03" = ATE_A_3$DE_G0_1, "0.01" = ATE_A_1$DE_G0_1, "0.005" = ATE_A_05$DE_G0_1)
names(DEG0_1_A) <- sub("^X", "", names(DEG0_1_A))
DEG0_1_A$group <- rep('AIPW', nrow(DEG0_1_A))
DEG0_1_A <- melt(DEG0_1_A, id.vars = 'group')

DEG0_1_S <- data.frame("0.5" = ATE_S_50$DE_G0_1, "0.2" = ATE_S_20$DE_G0_1, "0.1" = ATE_S_10$DE_G0_1,
                       "0.05" = ATE_S_5$DE_G0_1, "0.03" = ATE_S_3$DE_G0_1, "0.01" = ATE_S_1$DE_G0_1, "0.005" = ATE_S_05$DE_G0_1)
names(DEG0_1_S) <- sub("^X", "", names(DEG0_1_S))
DEG0_1_S$group <- rep('SAIPW', nrow(DEG0_1_S))
DEG0_1_S <- melt(DEG0_1_S, id.vars = 'group')


DEG0_1_allData <- rbind(DEG0_1_G, DEG0_1_A, DEG0_1_S)
DEG0_1_allData$group <- factor(DEG0_1_allData$group, levels = c("G-comp", "AIPW", "SAIPW"))

ggplot(DEG0_1_allData, aes(x = variable, y = value, fill = group)) +
  geom_boxplot() +
  coord_cartesian(ylim=c(-2.9,1.1)) +
  geom_hline(yintercept = -1 ,alpha=0.4) +
  facet_wrap( ~ variable, scales = 'free_x', ncol = 6) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  xlab("Sample Size") +
  ylab("ATE Estimate") +
  labs(fill = "Estimator")



## DE_G0_0

DEG0_0_G <- data.frame("0.5" = ATE_G_50$DE_G0_0, "0.2" = ATE_G_20$DE_G0_0, "0.1" = ATE_G_10$DE_G0_0,
                       "0.05" = ATE_G_5$DE_G0_0, "0.03" = ATE_G_3$DE_G0_0, "0.01" = ATE_G_1$DE_G0_0, "0.005" = ATE_G_05$DE_G0_0)
names(DEG0_0_G) <- sub("^X", "", names(DEG0_0_G))
DEG0_0_G$group <- rep('G-comp', nrow(DEG0_0_G))
DEG0_0_G <- melt(DEG0_0_G, id.vars = 'group')


DEG0_0_A <- data.frame("0.5" = ATE_A_50$DE_G0_0, "0.2" = ATE_A_20$DE_G0_0, "0.1" = ATE_A_10$DE_G0_0,
                       "0.05" = ATE_A_5$DE_G0_0, "0.03" = ATE_A_3$DE_G0_0, "0.01" = ATE_A_1$DE_G0_0, "0.005" = ATE_A_05$DE_G0_0)
names(DEG0_0_A) <- sub("^X", "", names(DEG0_0_A))
DEG0_0_A$group <- rep('AIPW', nrow(DEG0_0_A))
DEG0_0_A <- melt(DEG0_0_A, id.vars = 'group')

DEG0_0_S <- data.frame("0.5" = ATE_S_50$DE_G0_0, "0.2" = ATE_S_20$DE_G0_0, "0.1" = ATE_S_10$DE_G0_0,
                       "0.05" = ATE_S_5$DE_G0_0, "0.03" = ATE_S_3$DE_G0_0, "0.01" = ATE_S_1$DE_G0_0, "0.005" = ATE_S_05$DE_G0_0)
names(DEG0_0_S) <- sub("^X", "", names(DEG0_0_S))
DEG0_0_S$group <- rep('SAIPW', nrow(DEG0_0_S))
DEG0_0_S <- melt(DEG0_0_S, id.vars = 'group')


DEG0_0_allData <- rbind(DEG0_0_G, DEG0_0_A, DEG0_0_S)
DEG0_0_allData$group <- factor(DEG0_0_allData$group, levels = c("G-comp", "AIPW", "SAIPW"))

ggplot(DEG0_0_allData, aes(x = variable, y = value, fill = group)) +
  geom_boxplot() +
  coord_cartesian(ylim=c(-2.5,1.5)) +
  geom_hline(yintercept = 0 ,alpha=0.4) +
  facet_wrap( ~ variable, scales = 'free_x', ncol = 6) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  xlab("Sample Size") +
  ylab("ATE Estimate") +
  labs(fill = "Estimator")



## SE_Z1_2

SEZ1_2_G <- data.frame("0.5" = ATE_G_50$SE_Z1_2, "0.2" = ATE_G_20$SE_Z1_2, "0.1" = ATE_G_10$SE_Z1_2,
                       "0.05" = ATE_G_5$SE_Z1_2, "0.03" = ATE_G_3$SE_Z1_2, "0.01" = ATE_G_1$SE_Z1_2, "0.005" = ATE_G_05$SE_Z1_2)
names(SEZ1_2_G) <- sub("^X", "", names(SEZ1_2_G))
SEZ1_2_G$group <- rep('G-comp', nrow(SEZ1_2_G))
SEZ1_2_G <- melt(SEZ1_2_G, id.vars = 'group')


SEZ1_2_A <- data.frame("0.5" = ATE_A_50$SE_Z1_2, "0.2" = ATE_A_20$SE_Z1_2, "0.1" = ATE_A_10$SE_Z1_2,
                       "0.05" = ATE_A_5$SE_Z1_2, "0.03" = ATE_A_3$SE_Z1_2, "0.01" = ATE_A_1$SE_Z1_2, "0.005" = ATE_A_05$SE_Z1_2)
names(SEZ1_2_A) <- sub("^X", "", names(SEZ1_2_A))
SEZ1_2_A$group <- rep('AIPW', nrow(SEZ1_2_A))
SEZ1_2_A <- melt(SEZ1_2_A, id.vars = 'group')

SEZ1_2_S <- data.frame("0.5" = ATE_S_50$SE_Z1_2, "0.2" = ATE_S_20$SE_Z1_2, "0.1" = ATE_S_10$SE_Z1_2,
                       "0.05" = ATE_S_5$SE_Z1_2, "0.03" = ATE_S_3$SE_Z1_2, "0.01" = ATE_S_1$SE_Z1_2, "0.005" = ATE_S_05$SE_Z1_2)
names(SEZ1_2_S) <- sub("^X", "", names(SEZ1_2_S))
SEZ1_2_S$group <- rep('SAIPW', nrow(SEZ1_2_S))
SEZ1_2_S <- melt(SEZ1_2_S, id.vars = 'group')


SEZ1_2_allData <- rbind(SEZ1_2_G, SEZ1_2_A, SEZ1_2_S)
SEZ1_2_allData$group <- factor(SEZ1_2_allData$group, levels = c("G-comp", "AIPW", "SAIPW"))

ggplot(SEZ1_2_allData, aes(x = variable, y = value, fill = group)) +
  geom_boxplot() +
  coord_cartesian(ylim=c(-2.7,1.3)) +
  geom_hline(yintercept = -1 ,alpha=0.4) +
  facet_wrap( ~ variable, scales = 'free_x', ncol = 6) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  xlab("Sample Size") +
  ylab("ATE Estimate") +
  labs(fill = "Estimator")


## SE_Z1_1

SEZ1_1_G <- data.frame("0.5" = ATE_G_50$SE_Z1_1, "0.2" = ATE_G_20$SE_Z1_1, "0.1" = ATE_G_10$SE_Z1_1,
                       "0.05" = ATE_G_5$SE_Z1_1, "0.03" = ATE_G_3$SE_Z1_1, "0.01" = ATE_G_1$SE_Z1_1, "0.005" = ATE_G_05$SE_Z1_1)
names(SEZ1_1_G) <- sub("^X", "", names(SEZ1_1_G))
SEZ1_1_G$group <- rep('G-comp', nrow(SEZ1_1_G))
SEZ1_1_G <- melt(SEZ1_1_G, id.vars = 'group')


SEZ1_1_A <- data.frame("0.5" = ATE_A_50$SE_Z1_1, "0.2" = ATE_A_20$SE_Z1_1, "0.1" = ATE_A_10$SE_Z1_1,
                       "0.05" = ATE_A_5$SE_Z1_1, "0.03" = ATE_A_3$SE_Z1_1, "0.01" = ATE_A_1$SE_Z1_1, "0.005" = ATE_A_05$SE_Z1_1)
names(SEZ1_1_A) <- sub("^X", "", names(SEZ1_1_A))
SEZ1_1_A$group <- rep('AIPW', nrow(SEZ1_1_A))
SEZ1_1_A <- melt(SEZ1_1_A, id.vars = 'group')

SEZ1_1_S <- data.frame("0.5" = ATE_S_50$SE_Z1_1, "0.2" = ATE_S_20$SE_Z1_1, "0.1" = ATE_S_10$SE_Z1_1,
                       "0.05" = ATE_S_5$SE_Z1_1, "0.03" = ATE_S_3$SE_Z1_1, "0.01" = ATE_S_1$SE_Z1_1, "0.005" = ATE_S_05$SE_Z1_1)
names(SEZ1_1_S) <- sub("^X", "", names(SEZ1_1_S))
SEZ1_1_S$group <- rep('SAIPW', nrow(SEZ1_1_S))
SEZ1_1_S <- melt(SEZ1_1_S, id.vars = 'group')


SEZ1_1_allData <- rbind(SEZ1_1_G, SEZ1_1_A, SEZ1_1_S)
SEZ1_1_allData$group <- factor(SEZ1_1_allData$group, levels = c("G-comp", "AIPW", "SAIPW"))

ggplot(SEZ1_1_allData, aes(x = variable, y = value, fill = group)) +
  geom_boxplot() +
  coord_cartesian(ylim=c(-2.7,1.3)) +
  geom_hline(yintercept = -1 ,alpha=0.4) +
  facet_wrap( ~ variable, scales = 'free_x', ncol = 6) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  xlab("Sample Size") +
  ylab("ATE Estimate") +
  labs(fill = "Estimator")


## SE_Z1_0

SEZ1_0_G <- data.frame("0.5" = ATE_G_50$SE_Z1_0, "0.2" = ATE_G_20$SE_Z1_0, "0.1" = ATE_G_10$SE_Z1_0,
                       "0.05" = ATE_G_5$SE_Z1_0, "0.03" = ATE_G_3$SE_Z1_0, "0.01" = ATE_G_1$SE_Z1_0, "0.005" = ATE_G_05$SE_Z1_0)
names(SEZ1_0_G) <- sub("^X", "", names(SEZ1_0_G))
SEZ1_0_G$group <- rep('G-comp', nrow(SEZ1_0_G))
SEZ1_0_G <- melt(SEZ1_0_G, id.vars = 'group')


SEZ1_0_A <- data.frame("0.5" = ATE_A_50$SE_Z1_0, "0.2" = ATE_A_20$SE_Z1_0, "0.1" = ATE_A_10$SE_Z1_0,
                       "0.05" = ATE_A_5$SE_Z1_0, "0.03" = ATE_A_3$SE_Z1_0, "0.01" = ATE_A_1$SE_Z1_0, "0.005" = ATE_A_05$SE_Z1_0)
names(SEZ1_0_A) <- sub("^X", "", names(SEZ1_0_A))
SEZ1_0_A$group <- rep('AIPW', nrow(SEZ1_0_A))
SEZ1_0_A <- melt(SEZ1_0_A, id.vars = 'group')

SEZ1_0_S <- data.frame("0.5" = ATE_S_50$SE_Z1_0, "0.2" = ATE_S_20$SE_Z1_0, "0.1" = ATE_S_10$SE_Z1_0,
                       "0.05" = ATE_S_5$SE_Z1_0, "0.03" = ATE_S_3$SE_Z1_0, "0.01" = ATE_S_1$SE_Z1_0, "0.005" = ATE_S_05$SE_Z1_0)
names(SEZ1_0_S) <- sub("^X", "", names(SEZ1_0_S))
SEZ1_0_S$group <- rep('SAIPW', nrow(SEZ1_0_S))
SEZ1_0_S <- melt(SEZ1_0_S, id.vars = 'group')


SEZ1_0_allData <- rbind(SEZ1_0_G, SEZ1_0_A, SEZ1_0_S)
SEZ1_0_allData$group <- factor(SEZ1_0_allData$group, levels = c("G-comp", "AIPW", "SAIPW"))

ggplot(SEZ1_0_allData, aes(x = variable, y = value, fill = group)) +
  geom_boxplot() +
  coord_cartesian(ylim=c(-2.7,1.3)) +
  geom_hline(yintercept = -1 ,alpha=0.4) +
  facet_wrap( ~ variable, scales = 'free_x', ncol = 6) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  xlab("Sample Size") +
  ylab("ATE Estimate") +
  labs(fill = "Estimator")


## SE_Z0_2

SEZ0_2_G <- data.frame("0.5" = ATE_G_50$SE_Z0_2, "0.2" = ATE_G_20$SE_Z0_2, "0.1" = ATE_G_10$SE_Z0_2,
                       "0.05" = ATE_G_5$SE_Z0_2, "0.03" = ATE_G_3$SE_Z0_2, "0.01" = ATE_G_1$SE_Z0_2, "0.005" = ATE_G_05$SE_Z0_2)
names(SEZ0_2_G) <- sub("^X", "", names(SEZ0_2_G))
SEZ0_2_G$group <- rep('G-comp', nrow(SEZ0_2_G))
SEZ0_2_G <- melt(SEZ0_2_G, id.vars = 'group')


SEZ0_2_A <- data.frame("0.5" = ATE_A_50$SE_Z0_2, "0.2" = ATE_A_20$SE_Z0_2, "0.1" = ATE_A_10$SE_Z0_2,
                       "0.05" = ATE_A_5$SE_Z0_2, "0.03" = ATE_A_3$SE_Z0_2, "0.01" = ATE_A_1$SE_Z0_2, "0.005" = ATE_A_05$SE_Z0_2)
names(SEZ0_2_A) <- sub("^X", "", names(SEZ0_2_A))
SEZ0_2_A$group <- rep('AIPW', nrow(SEZ0_2_A))
SEZ0_2_A <- melt(SEZ0_2_A, id.vars = 'group')

SEZ0_2_S <- data.frame("0.5" = ATE_S_50$SE_Z0_2, "0.2" = ATE_S_20$SE_Z0_2, "0.1" = ATE_S_10$SE_Z0_2,
                       "0.05" = ATE_S_5$SE_Z0_2, "0.03" = ATE_S_3$SE_Z0_2, "0.01" = ATE_S_1$SE_Z0_2, "0.005" = ATE_S_05$SE_Z0_2)
names(SEZ0_2_S) <- sub("^X", "", names(SEZ0_2_S))
SEZ0_2_S$group <- rep('SAIPW', nrow(SEZ0_2_S))
SEZ0_2_S <- melt(SEZ0_2_S, id.vars = 'group')


SEZ0_2_allData <- rbind(SEZ0_2_G, SEZ0_2_A, SEZ0_2_S)
SEZ0_2_allData$group <- factor(SEZ0_2_allData$group, levels = c("G-comp", "AIPW", "SAIPW"))

ggplot(SEZ0_2_allData, aes(x = variable, y = value, fill = group)) +
  geom_boxplot() +
  coord_cartesian(ylim=c(-2.7,1.3)) +
  geom_hline(yintercept = -1 ,alpha=0.4) +
  facet_wrap( ~ variable, scales = 'free_x', ncol = 6) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  xlab("Sample Size") +
  ylab("ATE Estimate") +
  labs(fill = "Estimator")



## SE_Z0_1

SEZ0_1_G <- data.frame("0.5" = ATE_G_50$SE_Z0_1, "0.2" = ATE_G_20$SE_Z0_1, "0.1" = ATE_G_10$SE_Z0_1,
                       "0.05" = ATE_G_5$SE_Z0_1, "0.03" = ATE_G_3$SE_Z0_1, "0.01" = ATE_G_1$SE_Z0_1, "0.005" = ATE_G_05$SE_Z0_1)
names(SEZ0_1_G) <- sub("^X", "", names(SEZ0_1_G))
SEZ0_1_G$group <- rep('G-comp', nrow(SEZ0_1_G))
SEZ0_1_G <- melt(SEZ0_1_G, id.vars = 'group')


SEZ0_1_A <- data.frame("0.5" = ATE_A_50$SE_Z0_1, "0.2" = ATE_A_20$SE_Z0_1, "0.1" = ATE_A_10$SE_Z0_1,
                       "0.05" = ATE_A_5$SE_Z0_1, "0.03" = ATE_A_3$SE_Z0_1, "0.01" = ATE_A_1$SE_Z0_1, "0.005" = ATE_A_05$SE_Z0_1)
names(SEZ0_1_A) <- sub("^X", "", names(SEZ0_1_A))
SEZ0_1_A$group <- rep('AIPW', nrow(SEZ0_1_A))
SEZ0_1_A <- melt(SEZ0_1_A, id.vars = 'group')

SEZ0_1_S <- data.frame("0.5" = ATE_S_50$SE_Z0_1, "0.2" = ATE_S_20$SE_Z0_1, "0.1" = ATE_S_10$SE_Z0_1,
                       "0.05" = ATE_S_5$SE_Z0_1, "0.03" = ATE_S_3$SE_Z0_1, "0.01" = ATE_S_1$SE_Z0_1, "0.005" = ATE_S_05$SE_Z0_1)
names(SEZ0_1_S) <- sub("^X", "", names(SEZ0_1_S))
SEZ0_1_S$group <- rep('SAIPW', nrow(SEZ0_1_S))
SEZ0_1_S <- melt(SEZ0_1_S, id.vars = 'group')


SEZ0_1_allData <- rbind(SEZ0_1_G, SEZ0_1_A, SEZ0_1_S)
SEZ0_1_allData$group <- factor(SEZ0_1_allData$group, levels = c("G-comp", "AIPW", "SAIPW"))

ggplot(SEZ0_1_allData, aes(x = variable, y = value, fill = group)) +
  geom_boxplot() +
  coord_cartesian(ylim=c(-2.7,1.3)) +
  geom_hline(yintercept = -1 ,alpha=0.4) +
  facet_wrap( ~ variable, scales = 'free_x', ncol = 6) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  xlab("Sample Size") +
  ylab("ATE Estimate") +
  labs(fill = "Estimator")


## SE_Z0_0

SEZ0_0_G <- data.frame("0.5" = ATE_G_50$SE_Z0_0, "0.2" = ATE_G_20$SE_Z0_0, "0.1" = ATE_G_10$SE_Z0_0,
                       "0.05" = ATE_G_5$SE_Z0_0, "0.03" = ATE_G_3$SE_Z0_0, "0.01" = ATE_G_1$SE_Z0_0, "0.005" = ATE_G_05$SE_Z0_0)
names(SEZ0_0_G) <- sub("^X", "", names(SEZ0_0_G))
SEZ0_0_G$group <- rep('G-comp', nrow(SEZ0_0_G))
SEZ0_0_G <- melt(SEZ0_0_G, id.vars = 'group')


SEZ0_0_A <- data.frame("0.5" = ATE_A_50$SE_Z0_0, "0.2" = ATE_A_20$SE_Z0_0, "0.1" = ATE_A_10$SE_Z0_0,
                       "0.05" = ATE_A_5$SE_Z0_0, "0.03" = ATE_A_3$SE_Z0_0, "0.01" = ATE_A_1$SE_Z0_0, "0.005" = ATE_A_05$SE_Z0_0)
names(SEZ0_0_A) <- sub("^X", "", names(SEZ0_0_A))
SEZ0_0_A$group <- rep('AIPW', nrow(SEZ0_0_A))
SEZ0_0_A <- melt(SEZ0_0_A, id.vars = 'group')

SEZ0_0_S <- data.frame("0.5" = ATE_S_50$SE_Z0_0, "0.2" = ATE_S_20$SE_Z0_0, "0.1" = ATE_S_10$SE_Z0_0,
                       "0.05" = ATE_S_5$SE_Z0_0, "0.03" = ATE_S_3$SE_Z0_0, "0.01" = ATE_S_1$SE_Z0_0, "0.005" = ATE_S_05$SE_Z0_0)
names(SEZ0_0_S) <- sub("^X", "", names(SEZ0_0_S))
SEZ0_0_S$group <- rep('SAIPW', nrow(SEZ0_0_S))
SEZ0_0_S <- melt(SEZ0_0_S, id.vars = 'group')


SEZ0_0_allData <- rbind(SEZ0_0_G, SEZ0_0_A, SEZ0_0_S)
SEZ0_0_allData$group <- factor(SEZ0_0_allData$group, levels = c("G-comp", "AIPW", "SAIPW"))

ggplot(SEZ0_0_allData, aes(x = variable, y = value, fill = group)) +
  geom_boxplot() +
  coord_cartesian(ylim=c(-2.7,1.3)) +
  geom_hline(yintercept = -1 ,alpha=0.4) +
  facet_wrap( ~ variable, scales = 'free_x', ncol = 6) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  xlab("Sample Size") +
  ylab("ATE Estimate") +
  labs(fill = "Estimator")


## combined DE/SE plots

DEG1_2_G$bias <- abs(DEG1_2_G$value + 2)
DEG1_1_G$bias <- abs(DEG1_1_G$value + 1)
DEG1_0_G$bias <- abs(DEG1_0_G$value + 0)
DEG0_2_G$bias <- abs(DEG0_2_G$value + 2)
DEG0_1_G$bias <- abs(DEG0_1_G$value + 1)
DEG0_0_G$bias <- abs(DEG0_0_G$value + 0)
SEZ1_2_G$bias <- abs(SEZ1_2_G$value + 2)
SEZ1_1_G$bias <- abs(SEZ1_1_G$value + 1)
SEZ1_0_G$bias <- abs(SEZ1_0_G$value + 0)
SEZ0_2_G$bias <- abs(SEZ0_2_G$value + 2)
SEZ0_1_G$bias <- abs(SEZ0_1_G$value + 1)
SEZ0_0_G$bias <- abs(SEZ0_0_G$value + 0)

DEG1_2_A$bias <- abs(DEG1_2_A$value + 2)
DEG1_1_A$bias <- abs(DEG1_1_A$value + 1)
DEG1_0_A$bias <- abs(DEG1_0_A$value + 0)
DEG0_2_A$bias <- abs(DEG0_2_A$value + 2)
DEG0_1_A$bias <- abs(DEG0_1_A$value + 1)
DEG0_0_A$bias <- abs(DEG0_0_A$value + 0)
SEZ1_2_A$bias <- abs(SEZ1_2_A$value + 2)
SEZ1_1_A$bias <- abs(SEZ1_1_A$value + 1)
SEZ1_0_A$bias <- abs(SEZ1_0_A$value + 0)
SEZ0_2_A$bias <- abs(SEZ0_2_A$value + 2)
SEZ0_1_A$bias <- abs(SEZ0_1_A$value + 1)
SEZ0_0_A$bias <- abs(SEZ0_0_A$value + 0)

DEG1_2_S$bias <- abs(DEG1_2_S$value + 2)
DEG1_1_S$bias <- abs(DEG1_1_S$value + 1)
DEG1_0_S$bias <- abs(DEG1_0_S$value + 0)
DEG0_2_S$bias <- abs(DEG0_2_S$value + 2)
DEG0_1_S$bias <- abs(DEG0_1_S$value + 1)
DEG0_0_S$bias <- abs(DEG0_0_S$value + 0)
SEZ1_2_S$bias <- abs(SEZ1_2_S$value + 2)
SEZ1_1_S$bias <- abs(SEZ1_1_S$value + 1)
SEZ1_0_S$bias <- abs(SEZ1_0_S$value + 0)
SEZ0_2_S$bias <- abs(SEZ0_2_S$value + 2)
SEZ0_1_S$bias <- abs(SEZ0_1_S$value + 1)
SEZ0_0_S$bias <- abs(SEZ0_0_S$value + 0)

DE_G <- rbind(DEG1_2_G, DEG1_1_G, DEG1_0_G, DEG0_2_G, DEG0_1_G, DEG0_0_G)
DE_A <- rbind(DEG1_2_A, DEG1_1_A, DEG1_0_A, DEG0_2_A, DEG0_1_A, DEG0_0_A)
DE_S <- rbind(DEG1_2_S, DEG1_1_S, DEG1_0_S, DEG0_2_S, DEG0_1_S, DEG0_0_S)
SE_G <- rbind(SEZ1_2_G, SEZ1_1_G, SEZ1_0_G, SEZ0_2_G, SEZ0_1_G, SEZ0_0_G)
SE_A <- rbind(SEZ1_2_A, SEZ1_1_A, SEZ1_0_A, SEZ0_2_A, SEZ0_1_A, SEZ0_0_A)
SE_S <- rbind(SEZ1_2_S, SEZ1_1_S, SEZ1_0_S, SEZ0_2_S, SEZ0_1_S, SEZ0_0_S)


DE_allData <- rbind(DE_G, DE_A, DE_S)
DE_allData$group <- factor(DE_allData$group, levels = c("G-comp", "AIPW", "SAIPW"))
DE_allData$value <- NULL
DE_allData$bias <- 100*DE_allData$bias

SE_allData <- rbind(SE_G, SE_A, SE_S)
SE_allData$group <- factor(SE_allData$group, levels = c("G-comp", "AIPW", "SAIPW"))
SE_allData$value <- NULL
SE_allData$bias <- 100*SE_allData$bias

## by sample size
ggplot(DE_allData, aes(x = variable, y = bias, fill = group)) +
  geom_boxplot(outlier.size = 0.05, outlier.alpha = 0.05) +
  coord_cartesian(ylim=c(0,300)) +
  geom_hline(yintercept = 0 ,alpha=0.4) +
  facet_wrap( ~ variable, scales = 'free_x', ncol = 7) +
  theme_bw() + 
  theme(
    axis.title.x = element_text(vjust = 0, size = 12),
    axis.title.y = element_text(vjust = 2, size = 12),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  ) +
  xlab("Sample Proportion Size") + 
  ylab("Percent Absolute Bias (%)") +
  labs(fill = "Estimator") +
  ggtitle("Percent Absolute Bias of ATE Estimates (Direct Effects)", subtitle = "Varying Sample Proportion Size (of Population)")
ggsave("bias_samplesize_DE.png", width = 14, height = 7, dpi = 700)

ggplot(SE_allData, aes(x = variable, y = bias, fill = group)) +
  geom_boxplot(outlier.size = 0.05, outlier.alpha = 0.05) +
  coord_cartesian(ylim=c(0,300)) +
  geom_hline(yintercept = 0 ,alpha=0.4) +
  facet_wrap( ~ variable, scales = 'free_x', ncol = 7) +
  theme_bw() + 
  theme(
    axis.title.x = element_text(vjust = 0, size = 12),
    axis.title.y = element_text(vjust = 2, size = 12),
    strip.background = element_blank(),
    strip.text.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
  ) +
  xlab("Sample Proportion Size") + 
  ylab("Percent Absolute Bias (%)") +
  labs(fill = "Estimator") +
  ggtitle("Percent Absolute Bias of ATE Estimates (Spillover Effects)", subtitle = "Varying Sample Proportion Size (of Population)")
ggsave("bias_samplesize_SE.png", width = 14, height = 7, dpi = 700)






## RMSE
rmse_DE_G <- aggregate(DE_G$bias, by = list(DE_G$variable), FUN = function(x) sqrt(mean(x^2)))
rmse_DE_A <- aggregate(DE_A$bias, by = list(DE_A$variable), FUN = function(x) sqrt(mean(x^2)))
rmse_DE_S <- aggregate(DE_S$bias, by = list(DE_S$variable), FUN = function(x) sqrt(mean(x^2)))

rmse_SE_G <- aggregate(SE_G$bias, by = list(SE_G$variable), FUN = function(x) sqrt(mean(x^2)))
rmse_SE_A <- aggregate(SE_A$bias, by = list(SE_A$variable), FUN = function(x) sqrt(mean(x^2)))
rmse_SE_S <- aggregate(SE_S$bias, by = list(SE_S$variable), FUN = function(x) sqrt(mean(x^2)))

rmse_full_G <- aggregate(ATE_G_full_df$value, by = list(ATE_G_full_df$variable), FUN = function(x) sqrt(mean(x^2)))
rmse_full_A <- aggregate(ATE_A_full_df$value, by = list(ATE_A_full_df$variable), FUN = function(x) sqrt(mean(x^2)))
rmse_full_S <- aggregate(ATE_S_full_df$value, by = list(ATE_S_full_df$variable), FUN = function(x) sqrt(mean(x^2)))

