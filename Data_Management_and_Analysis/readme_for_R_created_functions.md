# R Package

This R package provides a suite of functions for epidemiological analyses. The package includes functions for descriptive statistics, Cox proportional hazards models, stratified analyses, competing risks models, and utilities for data preprocessing and visualization.

---
## Installation
To use this package, copy the code into an R script or R package environment. Install any necessary dependencies using:
```r
install.packages(c("here", "tidyverse", "tidymodels","survival", "rlang", "gtsummary", "purrr", "dplyr", "broom", "metafor", "gt", "survminer", "ggplot2"))
```
Below is a list of the functions in this package with descriptions of their purpose and usage.

## List of Functions

1. [Descriptive Table Generator](#1-descriptive-table-generator)
2. [Cox Proportional Hazard Assumption Check](#2-cox-proportional-hazard-assumption-check)
3. [Cox Models](#3-cox-models)
4. [Stratified Analyses](#4-stratified-analyses)
5. [Competing Risk Models](#5-competing-risk-models)
6. [E-value Calculation](#6-e-value-calculation)
7. [Forest Plot Preparation](#7-forest-plot-preparation)
8. [Interaction and Trend Tests](#8-interaction-and-trend-tests)
9. [Data Cleaning Functions](#9-data-cleaning-functions)


---

## 1. Descriptive Table Generator
### Function: `descript_table`
Generates a descriptive statistics table for specified variables.

**Usage:**

```r
descript_table(var_name, var_label, by = NULL, data)
```
**Arguments:**

`var_name`: Variable name(s) to summarize.

`var_label`: Label(s) for the variable(s).

`by`: Grouping variable (optional).

`data`: Data frame.

**Details:** Supports continuous and categorical variables with customizable statistics.

**Value:** A formatted tbl_summary object.

**Example:**

```r
var_label = list(
  "ENTRY_AGE" ~ "Age at baseline",
  "sex" ~ "Sex",
  "race_modif" ~ "Race"
  )

var_name = c("ENTRY_AGE", "sex", "race_modif")

descript_covar_by_night_shift = descript_table(var_name = var_name,
                                var_label = var_label,
                                by = "night_ev_all", #Night shift status (Yes/No)
                                data = night_shift_popu) 

```
## 2. Cox Proportional Hazard Assumption Check

### Function: `cox_models_assumption`
Checks the Cox proportional hazard assumption for **categorical** variables using time-dependent models.

**Usage**:
```r
cox_models_assumption(data, categories, response_formula)
```
**Arguments**:

`data`: Data frame containing survival data.

`categories`: Character vector of categorical variables to evaluate.

`response_formula`: Formula for the survival model (e.g., Surv(time, status)).

**Details:** 
https://bookdown.org/rwnahhas/RMPH/survival-phassumption.html

**Value:** 
Prints model summaries for variables including the time-dependent term.

**Example**:

```r
categories = c("night_ev_all", "night_perm_ev_all", "night_rotating_ev_all",
                 "night_perm_dur_all_cat", "night_perm_intensity_all_cat", 
                 "night_perm_cumulative_all_cat")
 
cox_models_assumption(data=night_shift_popu,
                        categories=categories,
                        response_formula = "Surv(ENTRY_AGE, EXIT_AGE, MOR_ALLCAUSE)")
```
## 3. Cox Models
### Function: `cox_model`
Performs Cox proportional hazards modeling.

**Usage:**

```r
cox_model(surv, var_name, var_label, covar, data, model_name = NULL)
```
**Arguments:**

`surv`: Survival formula.

`var_name`: Predictor variable(s).

`var_label`: Label(s) for predictor(s).

`covar`: Covariates for adjustment.

`data`: Data frame.

`model_name`: Optional model name for table header.

**Details:** Fits Cox models and outputs a formatted ready-to-publish regression table.

**Value:** A styled table with hazard ratios (HR) and confidence intervals (CI).

**Example:**

```r
var_label = list(
  "night_ev_all" ~ "Night shift work status"    ,
  "night_ev_all_detailed" ~ "Night shift work status (detailed)"
)

var_name = c("night_ev_all", "night_ev_all_detailed")

covariates = c(
  "sex" ~ "Sex",
  "race_modif" ~ "Race",
)

cox_model(surv = "Surv(ENTRY_AGE, EXIT_AGE, MOR_ALLCAUSE)",
                             var_name = var_name,
                             var_label = var_label,
                             covar = covariates,
                             data=night_shift_popu)
```

## 4. Stratified Analyses
### Functions: `stratification_wide_for_a_list`, `stratification_long_cat`

Perform stratified survival analyses across categorical or continuous variables.

**Usage:**

```r
stratification_wide_for_a_list(strat_list, surv, var_name, var_label, covar, add_p_interaction = FALSE, data)
stratification_long_for_a_list(strat_list, surv, var_name, var_label, covar, add_p_interaction = FALSE, data)
```
**Arguments:**

`strat_list`: Stratification variables.

`surv`: Survival formula.

`var_name`: Predictor variable(s).

`var_label`: Label(s) for predictor(s).

`covar`: Covariates for adjustment.

`data`: Data frame.

**Details:** These functions support both wide and long formats for stratified analyses.

**Value:** A merged or stacked table of stratified results.

**Example 1:**

```r
strat_list = list(
  "bmi_cat" ~ "BMI",
  "smk_status" ~ "Smoking status",
  "morning_chronotype" ~ "Morning chronotype"
)

covariates = c(
  "sex" ~ "Sex",
  "race_modif" ~ "Race",
)

stratification_wide_for_a_list(var_name="night_ev_all",
                               var_label=list("night_ev_all" ~ "Night shift work status"),
                               strat_list=strat_list,
                               add_p_interaction = TRUE,
                               surv="Surv(ENTRY_AGE, EXIT_AGE, MOR_ALLCAUSE)",
                               covar = covariates,
                               data=night_shift_popu)
```

**Example 2:**

```r
strat_list = list(
  "an_bmi_cat" ~ "Baseline BMI",
  "ses_edu_level" ~ "Personal attained education level",
  "ses_income" ~ "Annual household income")

strat_height_age10 = stratification_long_for_a_list(var_name = "an_height_age10",
                                                    var_label = list("an_height_age10" ~ "Height relative to peers at age 10"),
                                                    strat_list =  strat_list,
                                                    add_p_interaction = TRUE,
                                                    surv = "Surv(ident_ageexact_bl, ident_EOF, ident_DTC)",
                                                    covar = covariates_childhood,
                                                    data = earlylife_popu) 
```
## 5. Fine-Grey Competing Risk Models
### Function: `competing_risk_model`
Implements Fine-Gray competing risk models.

**Usage:**

```r
competing_risk_model(surv, var_name, var_label, covar, data, model_name = NULL)
```
**Arguments:** See Cox Models `cox_model`.

**Details:** Uses the Fine-Gray subdistribution hazard model for competing risks.

**Value:** A regression table similar to Cox models.

**Example:**

```r
 earlylife_popu$ident_compete_factor = factor(earlylife_popu$ident_compete)
   
 cox_cmprsk_childhood = competing_risk_model(
     surv = "Surv(ident_ageexact_bl, ident_EOF, ident_compete_factor)",
     var_name = var_name_childhood_sup,
     var_label = var_label_childhood_sup,
     covar = covariates_childhood,
     data = earlylife_popu
   )
   ```
## 6. E-value Calculation
### Functions: `evalues.HR`, `add_evalues`
**Usage:**

```r
evalues.HR(est, lo, hi, rare = NA, true = 1)
add_evalues(data)
```
**Arguments**

`est`: Point estimate of the hazard ratio.

`lo`: Lower confidence limit.

`hi`: Upper confidence limit.

`rare`: Boolean for rare outcomes.

`true`: True HR for comparison.

**Returns**
E-values or a data frame with added E-values.

**Details:** Calculates and appends E-values for hazard ratio estimates.

**Example:**

```r
multivariable_childhood = cox_model(surv = "Surv(ident_ageexact_bl, ident_EOF, ident_DTC)",
                             var_name = var_name_childhood,
                             var_label = var_label_childhood,
                             covar = covariates_childhood,
                             data=earlylife_popu)

Evalues = cbind(multivariable_childhood$`_data`$estimate,
                multivariable_childhood$`_data`$conf.low,
                multivariable_childhood$`_data`$conf.high)

colnames(Evalues) = c("estimate","conf.low","conf.high")

Evalues = round(add_evalues(Evalues),2)

Evalues = bind(multivariable_evalues$`_data`$label,Evalues)|>
            as.data.frame()
```
## 7. Forest Plot Preparation
### Function: `forest_plot_child`
Prepares forest plots for childhood stratified analyses.

**Usage:**

```r
forest_plot_child(dat, cat, label1, label2)
```

**Example:**

```r
strat_list = list(
  "an_bmi_cat" ~ "Baseline BMI",
  "ses_edu_level" ~ "Personal attained education level",
  "ses_income" ~ "Annual household income")

strat_height_age10 = stratification_long_for_a_list(var_name = "an_height_age10",
                                                    var_label = list("an_height_age10" ~ "Height relative to peers at age 10"),
                                                    strat_list =  strat_list,
                                                    add_p_interaction = TRUE,
                                                    surv = "Surv(ident_ageexact_bl, ident_EOF, ident_DTC)",
                                                    covar = covariates_childhood,
                                                    data = earlylife_popu) 
png(filename=here::here(childhood_height_age10_strat.png"))
strat_height_age10 |>forest_plot_child(cat = "Taller",label1 = "Height relative to peers at age 10: ",label2 = "Taller vs same height")
dev.off()
```
## 8. Interaction and Trend Tests
### Functions: `p_interaction`, `p_trend`
**Usage:**

```r
p_interaction(surv, var_name1, var_name2, covar, data)
p_trend(surv, var_name1, covar, data)
```
**Details:** Performs an F-test for nested models.
## 9. Miscellaneous Utilities 
### Functions:
`replace_values_with_na`: Replaces specific values with NA.

`set_unknown_level`: Assigns "Unknown" to missing values in factors.

`remove_numbering_in_levels`: Removes numeric prefixes from factor levels.

`factored_case_when`: A helper for creating ordered factors.

## Dependencies
`survival`: For survival analysis.

`gtsummary`: For descriptive and regression tables.

`purrr`, `dplyr`: For data manipulation.

`metafor`: For meta-analysis.