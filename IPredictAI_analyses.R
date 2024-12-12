## PACKAGES INSTALLATION
install.packages("lme4") 
install.packages("plyr") 
install.packages("tidyr") 
install.packages("openxlsx", dependencies = TRUE)
install.packages("survival")
install.packages("xlsx")
install.packages("readxl")

## LOADING PACKAGES
library("readxl")
library("lme4")
library(MASS)
library("plyr")
library(tidyr)
library(openxlsx)
library(survival)

"
DEFINITIONS
best case scenario = individual mean exposure + max iref use (100) + no injury (0)
worst case scenario = individual mean exposure + max iref use (100) + injury (1)
"

## IMPORT DATA
path <- "./IPREDICTAI/CODE/OUTPUT/"
date_database <- "20240623"

data <- read_excel(sprintf("%sIPredictAI_results_main %s.xlsx", path, date_database)) # The main dataset
data_best <- read_excel(sprintf("%sIPredictAI_results_main_best %s.xlsx", path, date_database)) # The best case scenario dataset
data_worst <- read_excel(sprintf("%sIPredictAI_results_main_worst %s.xlsx", path, date_database)) # The worst case scenario dataset
data_tte <- read_excel(sprintf("%sIPredictAI_results_tte_total_exposure %s.xlsx", path, date_database)) # The time to event dataset

## VARIABLES INITIALISATION
# NUMBER OF ATHLETES BY RESPONSE RATE
athletes_response_rate <- c()
rr <- c()

for (i in 0:100){
  rr_subset <- subset(data, data$response_rate_day_iref >= i)
  if (nrow(rr_subset) > 0){
    athletes_response_rate[[i+1]] <- rr_subset$u_id
    rr[i+1] <- i
  } else {
    athletes_response_rate[[i+1]] <- NA
    rr[i+1] <- i
  }
}

max_length <- max(sapply(athletes_response_rate, length))
 
padded_lists <- lapply(athletes_response_rate, function(x) {
  if(length(x) < max_length) {
    c(x, rep(NA, max_length - length(x)))} 
  else {x}})

athletes_rr <- as.data.frame(do.call(rbind, padded_lists))
athletes_rr <- as.data.frame(t(athletes_rr))
colnames(athletes_rr) <- c(lapply(rr, function(x) paste0(x, "%")))
write.xlsx(athletes_rr, sprintf("%sathletes_response_rate.xlsx", path))

# OUTPUT DATAFRAME
df_dict <- list(
  goal = c(), # The goal of the analysis
  analysis_name = c(), # The name of the analysis
  model_name = c(), # The name of the model used in the analysis
  model_formula = c(), # The formula of the model used in the analysis
  iref = c(), # The level of I-REF use in the analysis when perfoming sensitivity analysis
  injury = c(), # The injury status in the analysis when perfoming sensitivity analysis
  exposure = c(), # The exposure status in the analysis when perfoming sensitivity analysis
  coeff = c(), # The coefficient of the level of I-REF use from the model
  CI_lower = c(), # The lower bound of the confidence interval of the coefficient
  CI_upper = c(), # The upper bound of the confidence interval of the coefficient
  pvalue = c() # The p-value of the association between factors and the outcome
)

# TIME TO EVENT OBJECT TO USE FOR THE COX REGRESSION
survObject <- Surv(data_tte$exposure_athletics, data_tte$event_main_category_bin)

## FUNCTIONS 
# FUNCTION TO OUTPUT THE RESULTS OF THE ANALYSES
func_output <- function(name, model, df_dict, formula_str){
  # goal
  if (grepl("1_", name)) {
    df_dict$goal <- c(df_dict$goal, "primary")
  } else if (grepl("2_", name)) {
    df_dict$goal <- c(df_dict$goal, "secondary")
  } else if (grepl("3_", name)) {
    df_dict$goal <- c(df_dict$goal, "tertiary")
  }
  # analysis
  df_dict$analysis_name <- c(df_dict$analysis_name, name)
  # model_name
  df_dict$model_name <- c(df_dict$model_name, class(model)[1])
  # model_formula
  df_dict$model_formula <- c(df_dict$model_formula, formula_str)
  # iref / injury / exposure
  if (grepl("best", name)) {
    df_dict$iref <- c(df_dict$iref, 100)
    df_dict$injury <- c(df_dict$injury, 0)
    df_dict$exposure <- c(df_dict$exposure, "mean(indiv)")
  } else if (grepl("worst", name)) {
    df_dict$iref <- c(df_dict$iref, 100)
    df_dict$injury <- c(df_dict$injury, 1)
    df_dict$exposure <- c(df_dict$exposure, "mean(indiv)")
  } else {
    df_dict$iref <- c(df_dict$iref, NA)
    df_dict$injury <- c(df_dict$injury, NA)
    df_dict$exposure <- c(df_dict$exposure, NA)
  }
  # coefficient / CI_lower / CI_upper
  var_coeff <- summary(model)$coefficients[2,1]
  var_ci_lower <- confint.default(model)[2,1]
  var_ci_upper <- confint.default(model)[2,2]

  if (class(model)[1]=="negbin" || class(model)[1]=="glm"){
    df_dict$coeff <- c(df_dict$coeff, round(exp(var_coeff), 4))
    df_dict$CI_lower <- c(df_dict$CI_lower, round(exp(var_ci_lower), 4))
    df_dict$CI_upper <- c(df_dict$CI_upper, round(exp(var_ci_upper), 4))
    df_dict$pvalue <- c(df_dict$pvalue, round(summary(model)$coefficients[2,4], 4))
  } else if (class(model)[1]=="coxph") {
    df_dict$coeff <- c(df_dict$coeff, round(exp(summary(model)$coefficients[1,1]), 4))
    df_dict$CI_lower <- c(df_dict$CI_lower, round(exp(confint.default(model)[1,1]), 4))
    df_dict$CI_upper <- c(df_dict$CI_upper, round(exp(confint.default(model)[1,2]), 4))
    df_dict$pvalue <- c(df_dict$pvalue, round(summary(model)$coefficients[1,5], 4))
  } else {
    df_dict$coeff <- c(df_dict$coeff, round(var_coeff, 4))
    df_dict$CI_lower <- c(df_dict$CI_lower, round(var_ci_lower, 4))
    df_dict$CI_upper <- c(df_dict$CI_upper, round(var_ci_upper, 4))
    df_dict$pvalue <- c(df_dict$pvalue, round(summary(model)$coefficients[2,4], 4))
  }
  return(df_dict)
}

# FUNCTION TO PERFORM THE LOGISTIC REGRESSION
func_model_glm <- function(model_name, df, outcome, var, df_dict){
    formula_str  <- paste(outcome, "~",  var, "+ ages + gender + injury_history_20122")
    model <- glm(formula=formula_str, data=df, family="binomial")
    df_dict <- func_output(name=model_name, model=model, df_dict=df_dict, formula_str=formula_str)
    return(df_dict)
}

# FUNCTION TO PERFORM THE LOGISTIC REGRESSION WITH RESPONSE RATE
func_model_response_rate_glm <- function(df, outcome, var, athletes_response_rate, file_name, path){
  rr <- c()
  association <- c()
  pvalue <- c()
  part_included <- c()
  confint_low <- c()
  confint_high <- c()
  error <- c()
  
  for (i in 0:100){
    rr_subset <- subset(df, df$u_id %in% athletes_response_rate[[i+1]])
    if (nrow(rr_subset) > 0){
      formula_str  <- paste(outcome, "~",  var, "+ ages + gender + injury_history_20122")
      tryCatch({
        model <- glm(formula=formula_str, data=rr_subset, family="binomial")
        association[i+1] <- coefficients(model)[2]
        pvalue[i+1] <- summary(model)$coefficients[2,4]
        confint_low[i+1] <- confint.default(model)[2,1]
        confint_high[i+1] <- confint.default(model)[2,2]
        part_included[i+1] <- nrow(rr_subset)
        rr[i+1] <- i
        }, 
        error = function(e) {
          association[i+1] <<- NA
          pvalue[i+1] <<- NA
          confint_low[i+1] <<- NA
          confint_high[i+1] <<- NA
          part_included[i+1] <<- nrow(rr_subset)
          rr[i+1] <<- i
        }    
      )
    } else {
      association[i+1] <- NA
      pvalue[i+1] <- NA
      confint_low[i+1] <- NA
      confint_high[i+1] <- NA
      part_included[i+1] <- 0
      rr[i+1] <- i
    }
  }

  association <- exp(association)
  confint_low <- exp(confint_low)
  confint_high <- exp(confint_high)

  x <- data.frame(rr, association, pvalue, confint_low, confint_high, part_included)
  x <- as.data.frame(sapply(x, function(e) replace(e, is.nan(e), NA)))
  xlsx_name <- sprintf("%s%s.xlsx", path, file_name)
  write.xlsx(x, xlsx_name)
}

# FUNCTION TO PERFORM THE NEGATIVE BINOMIAL REGRESSION
func_model_nb <- function(model_name, df, outcome, var, var_offset, df_dict){
    formula_str  <- paste(outcome, "~",  var, "+ ages + gender + injury_history_20122 + offset(log(", var_offset, "))")
    model <- glm.nb(formula=formula_str, data=df)
    df_dict <- func_output(name=model_name, model=model, df_dict=df_dict, formula_str=formula_str)
    return(df_dict)
}

# FUNCTION TO PERFORM THE NEGATIVE BINOMIAL REGRESSION WITH RESPONSE RATE
func_model_response_rate_nb <- function(df, outcome, var, var_offset, athletes_response_rate, file_name, path){
  rr <- c()
  association <- c()
  pvalue <- c()
  part_included <- c()
  confint_low <- c()
  confint_high <- c()
  error <- c()
  
  for (i in 0:100){
    rr_subset <- subset(df, df$u_id %in% athletes_response_rate[[i+1]])
    if (nrow(rr_subset) > 0){
      formula_str  <- paste(outcome, "~",  var, "+ ages + gender + injury_history_20122 + offset(log(", var_offset, "))")
      tryCatch({
        model <- glm.nb(formula=formula_str, data=rr_subset)
        association[i+1] <- coefficients(model)[2]
        pvalue[i+1] <- summary(model)$coefficients[2,4]
        confint_low[i+1] <- confint.default(model)[2,1]
        confint_high[i+1] <- confint.default(model)[2,2]
        part_included[i+1] <- nrow(rr_subset)
        rr[i+1] <- i
        }, 
        error = function(e) {
          association[i+1] <<- NA
          pvalue[i+1] <<- NA
          confint_low[i+1] <<- NA
          confint_high[i+1] <<- NA
          part_included[i+1] <<- nrow(rr_subset)
          rr[i+1] <<- i
          }    
      )
    } else {
      association[i+1] <- NA
      pvalue[i+1] <- NA
      confint_low[i+1] <- NA
      confint_high[i+1] <- NA
      part_included[i+1] <- 0
      rr[i+1] <- i
    }
  }

  association <- exp(association)
  confint_low <- exp(confint_low)
  confint_high <- exp(confint_high)

  x <- data.frame(rr, association, pvalue, confint_low, confint_high, part_included)
  x <- as.data.frame(sapply(x, function(e) replace(e, is.nan(e), NA)))
  xlsx_name <- sprintf("%s%s.xlsx", path, file_name)
  write.xlsx(x, xlsx_name)
}

print("PROCESSING...")
## ---- PRIMARY GOAL ----
# Relationship between the level of I-REF use and the ICPR burden during an athletics season.
df_dict <- func_model_nb(model_name='1_burden_main', df=data, outcome="number_day_injury_athletics", var="lvl_iref_use_generated_view", var_offset="exposure_athletics", df_dict=df_dict)
df_dict <- func_model_nb(model_name='1_burden_best', df=data_best, outcome="number_day_injury_athletics", var="lvl_iref_use_generated_view", var_offset="exposure_athletics", df_dict=df_dict)
df_dict <- func_model_nb(model_name='1_burden_worst', df=data_worst, outcome="number_day_injury_athletics", var="lvl_iref_use_generated_view", var_offset="exposure_athletics", df_dict=df_dict)

func_model_response_rate_nb(df=data, outcome='number_day_injury_athletics', var="lvl_iref_use_generated_view", var_offset="exposure_athletics", athletes_response_rate=athletes_response_rate, file_name='1_burden', path=path)
func_model_response_rate_nb(df=data_best, outcome='number_day_injury_athletics', var="lvl_iref_use_generated_view", var_offset="exposure_athletics", athletes_response_rate=athletes_response_rate, file_name='1_burden_best', path=path)
func_model_response_rate_nb(df=data_worst, outcome='number_day_injury_athletics', var="lvl_iref_use_generated_view", var_offset="exposure_athletics", athletes_response_rate=athletes_response_rate, file_name='1_burden_worst', path=path)

## ---- SECONDARY GOAL ----
# Association between the level of I-REF use and the percentage of athletes with at least one ICPR during an athletics season.
df_dict <- func_model_glm(model_name='2_prevalence_main', df=data, outcome="prevalence_individual", var="lvl_iref_use_generated_view", df_dict=df_dict)
df_dict <- func_model_glm(model_name='2_prevalence_best', df=data_best, outcome="prevalence_individual", var="lvl_iref_use_generated_view", df_dict=df_dict)
df_dict <- func_model_glm(model_name='2_prevalence_worst', df=data_worst, outcome="prevalence_individual", var="lvl_iref_use_generated_view", df_dict=df_dict)

func_model_response_rate_glm(df=data, outcome='prevalence_individual', var="lvl_iref_use_generated_view", athletes_response_rate=athletes_response_rate, file_name='2_prevalence', path=path)
func_model_response_rate_glm(df=data_best, outcome='prevalence_individual', var="lvl_iref_use_generated_view", athletes_response_rate=athletes_response_rate, file_name='2_prevalence_best', path=path)
func_model_response_rate_glm(df=data_worst, outcome='prevalence_individual', var="lvl_iref_use_generated_view", athletes_response_rate=athletes_response_rate, file_name='2_prevalence_worst', path=path)

# Association between the level of I-REF use and the time to the first ICPR during an athletics season.
res_cox_2 <- coxph(survObject ~ data_tte$lvl_iref_use_generated_view + ages + gender + injury_history_20122, data=data_tte)
df_dict <- func_output(name="2_tte", model=res_cox_2, df_dict=df_dict, formula_str="")

# Association between the level of I-REF use and the number of ICPR per 1000 hours of athletics activity during an athletics season
df_dict <- func_model_nb(model_name='2_incidence_main', df=data, outcome="number_injury_athletics", var="lvl_iref_use_generated_view", var_offset="exposure_athletics", df_dict=df_dict)
df_dict <- func_model_nb(model_name='2_incidence_best', df=data_best, outcome="number_injury_athletics", var="lvl_iref_use_generated_view", var_offset="exposure_athletics", df_dict=df_dict)
df_dict <- func_model_nb(model_name='2_incidence_worst', df=data_worst, outcome="number_injury_athletics", var="lvl_iref_use_generated_view", var_offset="exposure_athletics", df_dict=df_dict)

func_model_response_rate_nb(df=data, outcome='number_injury_athletics', var="lvl_iref_use_generated_view", var_offset="exposure_athletics", athletes_response_rate=athletes_response_rate, file_name='2_incidence', path=path)
func_model_response_rate_nb(df=data_best, outcome='number_injury_athletics', var="lvl_iref_use_generated_view", var_offset="exposure_athletics", athletes_response_rate=athletes_response_rate, file_name='2_incidence_best', path=path)
func_model_response_rate_nb(df=data_worst, outcome='number_injury_athletics', var="lvl_iref_use_generated_view", var_offset="exposure_athletics", athletes_response_rate=athletes_response_rate, file_name='2_incidence_worst', path=path)

## ---- TERTIARY GOAL ----
# Association between the frequency of I-REF view and the ICPR burden during an athletics season.
df_dict <- func_model_nb(model_name='3_burden_main', df=data, outcome="number_day_injury_athletics", var="iref_view_generated", var_offset="exposure_athletics", df_dict=df_dict)
df_dict <- func_model_nb(model_name='3_burden_best', df=data_best, outcome="number_day_injury_athletics", var="iref_view_generated", var_offset="exposure_athletics", df_dict=df_dict)
df_dict <- func_model_nb(model_name='3_burden_worst', df=data_worst, outcome="number_day_injury_athletics", var="iref_view_generated", var_offset="exposure_athletics", df_dict=df_dict)

func_model_response_rate_nb(df=data, outcome='number_day_injury_athletics', var="iref_view_generated", var_offset="exposure_athletics", athletes_response_rate=athletes_response_rate, file_name='3_burden', path=path)
func_model_response_rate_nb(df=data_best, outcome='number_day_injury_athletics', var="iref_view_generated", var_offset="exposure_athletics", athletes_response_rate=athletes_response_rate, file_name='3_burden_best', path=path)
func_model_response_rate_nb(df=data_worst, outcome='number_day_injury_athletics', var="iref_view_generated", var_offset="exposure_athletics", athletes_response_rate=athletes_response_rate, file_name='3_burden_worst', path=path)

# Association between the frequency of I-REF view and the percentage of athletes with at least one ICPR during an athletics season.
df_dict <- func_model_glm(model_name='3_prevalence_main', df=data, outcome="prevalence_individual", var="iref_view_generated", df_dict=df_dict)
df_dict <- func_model_glm(model_name='3_prevalence_best', df=data_best, outcome="prevalence_individual", var="iref_view_generated", df_dict=df_dict)
df_dict <- func_model_glm(model_name='3_prevalence_worst', df=data_worst, outcome="prevalence_individual", var="iref_view_generated", df_dict=df_dict)

func_model_response_rate_glm(df=data, outcome='prevalence_individual', var="iref_view_generated", athletes_response_rate=athletes_response_rate, file_name='3_prevalence', path=path)
func_model_response_rate_glm(df=data_best, outcome='prevalence_individual', var="iref_view_generated", athletes_response_rate=athletes_response_rate, file_name='3_prevalence_best', path=path)
func_model_response_rate_glm(df=data_worst, outcome='prevalence_individual', var="iref_view_generated", athletes_response_rate=athletes_response_rate, file_name='3_prevalence_worst', path=path)

# Association between the frequency of I-REF view and the time to the first ICPR during an athletics season.
res_cox_3 <- coxph(survObject ~ data_tte$iref_view_generated + ages + gender + injury_history_20122, data=data_tte)
df_dict <- func_output(name="3_tte", model=res_cox_3, df_dict=df_dict, formula_str="")

# Association between the frequency of I-REF view and the number of ICPR per 1000 hours of athletics activity during an athletics season
df_dict <- func_model_nb(model_name='3_incidence_main', df=data, outcome="number_injury_athletics", var="iref_view_generated", var_offset="exposure_athletics", df_dict=df_dict)
df_dict <- func_model_nb(model_name='3_incidence_best', df=data_best, outcome="number_injury_athletics", var="iref_view_generated", var_offset="exposure_athletics", df_dict=df_dict)
df_dict <- func_model_nb(model_name='3_incidence_worst', df=data_worst, outcome="number_injury_athletics", var="iref_view_generated", var_offset="exposure_athletics", df_dict=df_dict)

func_model_response_rate_nb(df=data, outcome='number_injury_athletics', var="iref_view_generated", var_offset="exposure_athletics", athletes_response_rate=athletes_response_rate, file_name='3_incidence', path=path)
func_model_response_rate_nb(df=data_best, outcome='number_injury_athletics', var="iref_view_generated", var_offset="exposure_athletics", athletes_response_rate=athletes_response_rate, file_name='3_incidence_best', path=path)
func_model_response_rate_nb(df=data_worst, outcome='number_injury_athletics', var="iref_view_generated", var_offset="exposure_athletics", athletes_response_rate=athletes_response_rate, file_name='3_incidence_worst', path=path)

## EXPORT
df <- as.data.frame(df_dict)
write.xlsx(df, sprintf("%sIPredictAI_results_statistics %s.xlsx", path, date_database))
print("DONE")