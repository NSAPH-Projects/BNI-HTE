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
library(sf)
library(maptools)
library(tmap)
library(choroplethrZip)

load(".../HyADSmat.Rda")
load(".../zip_covs.rda")
load(".../facilities_for_analysis.RDa")

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
pp.covs = c("totNumNOxControls", "logHeatInput", "logOpTime", "pctCapacity", "pctS_n_CR", "Phase2", "ScrubbedFacility")  # - these are power plant covariates that are very predictive in PS model, but likely not confounders

downweight_ratio <- 0.6786576 ## mean of ratio of 2nd key-associated and key-associated HyADS for each zip code

################################################################
## Simulation with Key-Associated Plants
################################################################

##### Data Preparation and Cleaning #####
hyads_var <- as.data.frame(sumMAP.2005.fac) ## hyads matrix
zip_var <- data ## zipcode covariates
facility_var_2005 <- dat_facility_2005 ## power plant covariates
rm("sumMAP.2005.fac")
rm("data")
rm("dat_facility_2005")


zips <- data.frame("ZIP" = zip_var$ZIP) ## list of zipcodes
plants <- data.frame("FacID" = facility_var_2005$FacID) ## list of power plant IDs

## transform covariates
zip_var$logPop <- log(zip_var$TotPop)
zip_var$PctNonwhite <- 1 - zip_var$White_rate
facility_var_2005$logHeatInput <- log(facility_var_2005$totHeatInput)
facility_var_2005$logOpTime <- log(facility_var_2005$totOpTime)

zip_covariates <- select(zip_var, c("ZIP", all_of(zip.covs))) ## select zipcode covariates to be used
plant_covariates <- select(facility_var_2005, c("FacID", all_of(pp.covs))) ## select power plant covariates to be used


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
key_covariates <- aggregate(zip_agg[,-c(1:5)], list(zip_agg$Key), mean, na.rm=T)
colnames(key_covariates) <- c("Key", paste("Key", colnames(key_covariates)[-1], sep = "_"))
secondkey_covariates <- aggregate(zip_agg[,-c(1:5)], list(zip_agg$SecondKey), mean, na.rm=T)
colnames(secondkey_covariates) <- c("SecondKey", paste("SecondKey",colnames(secondkey_covariates)[-1], sep = "_"))
secondkey_covariates[,-1] <- secondkey_covariates[,-1] * downweight_ratio ## downweight second key-associated
plant_dat <- key_covariates %>% ## plant, key-associated covariates, second key-associated covariates, plant covariates
  full_join(secondkey_covariates, by = c("Key"="SecondKey")) %>%
  left_join(plant_covariates, by = c("Key"="FacID"))






##### Filter by complete covariate data #####
  complete_all <- plant_dat[complete.cases(plant_dat),]
names(complete_all)[names(complete_all) == "Key"] <- "FacID"
plants_complete <- complete_all$FacID

hyads_mat_complete <- hyads_var[plants_complete,] ## hyads matrix based on plants with complete data
## connect key-associated and 2nd key-associated plant to each outcome unit
dat_complete <- as.data.frame(t(sapply(hyads_mat_complete, function(x) head(row.names(hyads_mat_complete)[order(x, decreasing = TRUE)], 2))))
setDT(dat_complete, keep.rownames = "ZIP")
colnames(dat_complete) <- c("ZIP", "Key", "SecondKey")
dat_complete <- as.data.frame(dat_complete)
dat_complete$Key <- as.character(dat_complete$Key)
dat_complete$SecondKey <- as.character(dat_complete$SecondKey)
dat_complete$Key_HyADS <- hyads_mat_complete[as.matrix(cbind(dat_complete$Key, dat_complete$ZIP))] ## hyads for key-associated plants
dat_complete$SecondKey_HyADS <- hyads_mat_complete[as.matrix(cbind(dat_complete$SecondKey, dat_complete$ZIP))]

changes <- anti_join(dat, dat_complete, by=c("ZIP", "Key", "SecondKey")) ## look at changes

dat_complete_original <- dat_complete

zips_df <- data.frame(region = dat_complete$ZIP,
                      value = rep(1, length(dat_complete$ZIP)))

zip_choropleth(zips_df) + scale_fill_discrete(na.value = "gray")


plot(dat_complete$Key_HyADS,dat_complete$SecondKey_HyADS)




## generate treatments
treatments <- left_join(complete_all, facility_var_2005[,c(1,7)], by = "FacID")
names(treatments)[names(treatments) == "Key"] <- "FacID"
treatments <- treatments %>% mutate(x = 0.1*Key_logPop
                                    - 1.5*Key_logPop*Key_PctUrban + 0.05*logOpTime^2)
treatments$prop <- expit(treatments$x)
treatments$T <- rbinom(nrow(treatments),1,treatments$prop)
options(scipen = 10)




#############################################
##### propensity score estimation #####
#############################################

treat.formula_corr <- paste(treat.covs_corr, collapse=" + ")
treat.formula_miss <- paste(treat.covs_miss, collapse=" + ")
formulaT_corr <- as.formula(paste('T~ Key_logPop + Key_logPop:Key_PctUrban + I(logOpTime^2)'))
formulaT_miss <- as.formula(paste('T~ Key_logPop + Key_PctUrban + logOpTime'))


## run this for correct specification
#############################################
treat.mod_corr <- glm(formulaT_corr, treatments, family = binomial(link = "logit"))
treatments$prop_est <- treat.mod_corr$fitted.values


## run this for misspecification
#############################################
treat.mod_miss <- glm(formulaT_miss, treatments, family = binomial(link = "logit"))
treatments$prop_est <- treat.mod_miss$fitted.values




## merge zip-level data with key-associated plant covariates and covariate summaries
dat_merge <- dat_complete %>%
left_join(plant_covariates, by = c("Key"="FacID")) %>%
left_join(key_covariates, by = "Key") %>%
left_join(secondkey_covariates, by = "SecondKey")
if(sum(is.na(dat_merge)) > 0) {print("warning: NAs in covariates")}
dat_merge <- filter(dat_merge, Key_HyADS > summary(dat_merge$Key_HyADS)[2]) ## Drop bottom 25th percentile of zip key HyADS scores
dat_merge <- filter(dat_merge, ZIP < 80000)


## Get key and upwind treatments using plant treatment mapping
treat_mapping <- setNames(treatments$T, treatments$FacID)
dat_merge$Z <- as.integer(treat_mapping[dat_merge$Key]) ## based on key-associated plant: Z=1 if treated, Z=0 if not
dat_merge$G <- as.integer(treat_mapping[dat_merge$SecondKey])
if(sum(is.na(dat_merge$Z)) + sum(is.na(dat_merge$G)) > 0) {print("warning: NAs produced in Z and/or G")}

## Map true treatment probabilities to zipcode level
prop_mapping <- setNames(treatments$prop, treatments$FacID)
prop_est_mapping <- setNames(treatments$prop_est, treatments$FacID)
dat_merge$Z_prop <- as.numeric(prop_mapping[dat_merge$Key])
dat_merge$G_prop <- as.numeric(prop_mapping[dat_merge$SecondKey])
dat_merge$Z_prop_est <- as.numeric(prop_est_mapping[dat_merge$Key])
dat_merge$G_prop_est <- as.numeric(prop_est_mapping[dat_merge$SecondKey])

dat_merge <- left_join(dat_merge, zip_covariates, by = "ZIP")
dat_merge <- dat_merge[complete.cases(dat_merge),]

hyads_scores <- dat_merge %>% select("Key_HyADS","SecondKey_HyADS")
dat_merge <- dat_merge %>% select(-contains("Key_"))
dat_merge <- cbind(dat_merge, hyads_scores)

effect_mod <- paste(c("PctNonwhite", "PctPoor", "PctUrban"), collapse=" + ")

dat_merge <- dat_merge %>% mutate(PctNonwhite_bin = ifelse(PctNonwhite > quantile(PctNonwhite, 0.33), 1, 0))
dat_merge <- dat_merge %>% mutate(PctPoor_bin = ifelse(PctPoor > median(PctPoor), 1, 0))
dat_merge <- dat_merge %>% mutate(PctUrban_bin = ifelse(PctUrban > median(PctUrban), 1, 0))


Xzip.formula <- paste(zip.covs, collapse=" + ")
.env <- environment() ## identify the environment of cv.step
# formulaZ <- as.formula(paste('Z~', Xzip.formula), env = .env)
# formulaG <- as.formula(paste('G~Z+', Xzip.formula), env = .env)






## truncate est propensity score
quantileZ <- quantile(dat_merge$Z_prop_est, c(0.1,1))
quantileG <- quantile(dat_merge$G_prop_est, c(0.1,1))
dat_merge$Z_prop_est <- ifelse(dat_merge$Z_prop_est < quantileZ[1], quantileZ[1], dat_merge$Z_prop_est)
dat_merge$G_prop_est <- ifelse(dat_merge$G_prop_est < quantileG[1], quantileG[1], dat_merge$G_prop_est)

## Joint propensity scores
dat_merge$Z11 <- dat_merge$Z_prop_est * dat_merge$G_prop_est ## P(Z=1,G=1 | X_out)
dat_merge$Z10 <- dat_merge$Z_prop_est * (1-dat_merge$G_prop_est) ## P(Z=1,G=0 | X_out)
dat_merge$Z01 <- (1-dat_merge$Z_prop_est) * dat_merge$G_prop_est ## P(Z=0,G=1 | X_out)
dat_merge$Z00 <- (1-dat_merge$Z_prop_est) * (1-dat_merge$G_prop_est) ## P(Z=0,G=0 | X_out)








