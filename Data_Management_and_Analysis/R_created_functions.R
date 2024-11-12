
###Create functions that I'll use in the codes in the Night Shift Project

#0.1 Check Cox proportional hazard assumption:
#https://bookdown.org/rwnahhas/RMPH/survival-phassumption.html
cox_models_assumption <- function(data, categories, response_formula) {
  
  # Helper function to create binary variables for each level in a category
  create_binary_vars <- function(df, category) {
    levels <- unique(df[[category]])
    for (level in levels) {
      var_name <- paste0(category, "_test_", gsub(" ", "_", tolower(level)))
      df[[var_name]] <- ifelse(df[[category]] == level, 1, 0)
    }
    return(df)
  }
  
  # Main loop through each category
  for (category in categories) {
    
    # Create binary variables for the current category
    data <- create_binary_vars(data, category)
    binary_vars <- grep(paste0("^", category, "_test_"), names(data), value = TRUE)
    binary_vars_all = paste0(binary_vars,collapse =  " + ")
    # Loop through each binary variable created for the category
    for (var in binary_vars) {
      
      # Time-dependent Cox model
      td_model <- coxph(as.formula(paste(response_formula, " ~ ", binary_vars_all," + tt(", var, ")")), 
                        data = data, tt = function(x, t, ...) x * t)
      cat("\nTime-dependent model for:", var, "\n")
      print(summary(td_model))
    }
  }
}
#For example instead of doing this: 
# night_shift_popu_test=night_shift_popu|>mutate(night_ev_all_test_no= case_when(night_ev_all=="No"~1,.default = 0),
# night_ev_all_test_yes= case_when(night_ev_all=="Yes"~1,.default = 0),
# night_ev_all_test_unknown= case_when(night_ev_all=="Unknown status"~1,.default = 0))
# 
# coxph(Surv(ENTRY_AGE, EXIT_AGE, MOR_ALLCAUSE)~night_ev_all_test_no + night_ev_all_test_yes + night_ev_all_test_unknown + tt(night_ev_all_test_yes),data=night_shift_popu_test,tt=function(x,t,...) x*t)
# coxph(Surv(ENTRY_AGE, EXIT_AGE, MOR_ALLCAUSE)~night_ev_all_test_no + night_ev_all_test_yes + night_ev_all_test_unknown + tt(night_ev_all_test_no),data=night_shift_popu_test,tt=function(x,t,...) x*t)
# coxph(Surv(ENTRY_AGE, EXIT_AGE, MOR_ALLCAUSE)~night_ev_all_test_no + night_ev_all_test_yes + night_ev_all_test_unknown + tt(night_ev_all_test_unknown),data=night_shift_popu_test,tt=function(x,t,...) x*t)
# just run this code: cox_models_assumption(data=night_shift_popu_test,categories="night_ev_all",response_formula = "Surv(ENTRY_AGE, EXIT_AGE, MOR_ALLCAUSE)")



#1. Function for descriptive info
descript_table = function(var_name, var_label,by=NULL,data)
{
  Descriptive_table = data[,-1]   |>
    tbl_summary(by=by,
      include = var_name,
      type = all_continuous() ~ "continuous",
      statistic = list(
        all_continuous()
        ~ c("{median} ({p25}, {p75})"),
        all_categorical() ~ "{n} ({p})"
      ),
      label = var_label,
      digits = list(
        all_continuous() ~ 1,
        all_categorical() ~ c(0,1)
      ),
      missing = "ifany"
    )  |>
    
    modify_footnote(all_stat_cols() ~ NA)
  return(Descriptive_table)
}


#2. Function for cox models
cox_model = function(surv, var_name, var_label, covar, data, model_name=NULL)
{
  results_cox_function = map2(
    .x = var_name,
    .y = var_label,
    ~ reformulate(c(.x, covar), response = surv)
       |>
      coxph(method = "breslow", data = data) |>
      tbl_regression(
        exponentiate = TRUE,
        conf.int = TRUE,
        conf.level = 0.95,
        add_estimate_to_reference_rows = TRUE,
        label = .y,
        include = .x
      ) |>
      modify_fmt_fun(estimate ~ function(x)
        ifelse(is.na(x), NA, "—"),
        rows = n_event == "0") |>
      modify_fmt_fun(ci ~ function(x)
        ifelse(is.na(x), NA, "—"),
        rows = n_event == "0")
    
  ) |>
    tbl_stack() |>
    modify_column_hide(p.value) |>
    modify_header(n_event ~ "**Number of events**",
                  exposure ~ "**Person-years**",
                  estimate ~ "**HR**") |>
    modify_table_body(~ .x |> dplyr::relocate(exposure, .before = estimate))
  
  if (!is.null(model_name)) {
    results_cox_function = results_cox_function |>   modify_spanning_header(c(n_event, exposure, estimate, ci) ~
                                                                              paste0("**", model_name, "**"))
  }
  return(results_cox_function)
  rm(results_cox_function)
}

#3 Stratified analyses (wide)

stratification_wide = function(strat_var,
                          surv,
                          var_name,
                          var_label,
                          covar,
                          data)
{
  data = data |> 
    filter({{strat_var}}!="Unknown status")|>droplevels()
  
  # Extract unique levels from the stratification variable
  unique_levels = data |> select({{strat_var}}) |>
    summarise_each(list( ~ levels(.))) |>
    pull() |> as.factor()
  
  
  header = unique_levels
  header = paste("**", header, "**",sep="")
  
  
  # Use purrr::map to create a list of stratified data frames and apply the function cox_model on them
  results = unique_levels |>
    map(
      ~ data |>
        filter({{strat_var}} == .x) |>
        cox_model(surv = surv,
                  var_name = var_name,
                  var_label = var_label,
                  covar = covar
        )
    ) |>
    tbl_merge(tab_spanner = header)
  return(results)
  rm(results)
}

stratification_wide_for_a_list = function(var_name,
                                          var_label,
                                          strat_list,
                                          add_p_interaction = FALSE,
                                          surv,
                                          covar,
                                          data) {
  results_list <- list()
  for (strat in strat_list) {
    strat=as.character(strat)
    strat_label=strat[3]
    strat_var=as.name(strat[2])
    # Call the stratification_wide function with the current strat_var
    result <- stratification_wide(
      strat_var = {
        {
          strat_var
        }
      },
      surv = surv,
      var_name = var_name,
      var_label = var_label,
      covar = covar,
      data = data
    ) |> modify_table_body(~ .x  |>
                             dplyr::add_row(label = ifelse(add_p_interaction,"p-interaction",""))|>
                             dplyr::mutate(p_interaction = p_interaction(
                               surv = surv,
                               var_name1 = as.character(strat_var),
                               var_name2 = {{var_name}},
                               covar = covar,
                               data = data
                             )$`Pr(>|Chi|)`[2]))
    # Store the result in the results_list
    results_list[[as.character(strat_var)]] <- result
  }
  
  final_result = tbl_merge(results_list,tab_spanner = FALSE)
  return(final_result)
}


#4.1 Stratified analyses (long for continuous variables)


 stratification_long = function(strat_var,
                                surv,
                                var_name,
                                var_label,
                                covar,
                                data)
 {
   data = data |> 
     filter({{strat_var}}!="Unknown status")|>droplevels()
   
   # Extract unique levels from the stratification variable
   unique_levels = data |> select({{strat_var}}) |>
     summarise_each(list( ~ levels(.))) |>
     pull() |> as.factor()
   
   
   header = unique_levels
   header = paste("**", header, "**",sep="")
   
   
   # Use purrr::map to create a list of stratified data frames and appy the function cox_model on them
   results = unique_levels |>
     map(
       ~ data |>
         filter({{strat_var}} == .x) |>
         cox_model(surv = surv,
                   var_name = var_name,
                   var_label = var_label,
                   covar = covar
         )) |>
     tbl_stack()|>modify_table_body(~ .x |> dplyr::mutate(label = unique_levels))
   return(results)
   rm(results)
 }


#4.2 Stratified analyses (long for categorical variables)

stratification_long_cat = function(strat_var,
                                   surv,
                                   var_name,
                                   var_label,
                                   cat=NULL,
                                   covar,
                                   data)
{
  data = data |> 
    filter({{strat_var}}!="Unknown status")|>droplevels()
  
  # Extract unique levels from the stratification variable
  unique_levels = data |> select({{strat_var}}) |>
    summarise_each(list( ~ levels(.))) |>
    pull() |> as.factor()
  
  
  header = unique_levels
  header = paste("**", header, "**",sep="")
  
  
  # Use purrr::map to create a list of stratified data frames and appy the function cox_model on them
  results = unique_levels |>
    map(~ {
      # Filter data based on the current level of the categorical variable
      label = as.character(.x)
      filtered_data <- data |>
        filter({{strat_var}} == .x)
      
      # Apply cox_model function to the filtered data
      cox_model_results <- cox_model(filtered_data,
                                     surv = surv,
                                     var_name = var_name,
                                     var_label = var_label,
                                     covar = covar)
      
      # Modify table body based on the cox_model results
      cox_model_results$table_body$label[1]= label
      
      return(cox_model_results)
    }) |>
    tbl_stack()|> 
    #modify_table_body(filter, is.na(n_event) | label %in% cat)|>
    modify_column_hide(exposure)
  return(results)
  rm(results)
}

stratification_long_for_a_list = function(var_name,
                                          var_label,
                                          strat_list,
                                          add_p_interaction = FALSE,
                                          surv,
                                          covar,
                                          data) {
  results_list <- list()
  for (strat in strat_list) {
    strat=as.character(strat)
    strat_label=strat[3]
    strat_var=as.name(strat[2])
    # Call the stratification_long_cat function with the current strat_var
    result <- stratification_long_cat(
      strat_var = {
        {
          strat_var
        }
      },
      surv = surv,
      var_name = var_name,
      var_label = var_label,
      covar = covar,
      data = data
    ) |> modify_table_body(~ .x  |>
                             dplyr::add_row(label = strat_label,
                                            .before = 0)|>
                             dplyr::add_row(label = ifelse(add_p_interaction,"p-interaction",""))|>
                             dplyr::mutate(p_interaction = p_interaction(
                               surv = surv,
                               var_name1 = as.character(strat_var),
                               var_name2 = {{var_name}},
                               covar = covar,
                               data = data
                             )$`Pr(>|Chi|)`[2]))
    # Store the result in the results_list
    results_list[[as.character(strat_var)]] <- result
  }
  
  final_result = tbl_stack(results_list)
  return(final_result)
}

#5. Function for competing risk models

competing_risk_model = function(surv, var_name, var_label, covar, data, model_name = NULL)
{
  competing = function(x, y)
  {
    formula = reformulate(c(x, covar), response = "Surv(fgstart, fgstop, fgstatus)")
    finegray_data = reformulate(c(x, covar), response = surv) |>
      finegray(id = PSID, data = data)
    fine_gray = coxph(formula, weight = fgwt, finegray_data) |>
      tbl_regression(
        exponentiate = TRUE,
        conf.int = TRUE,
        conf.level = 0.95,
        add_estimate_to_reference_rows = TRUE,
        label = y,
        include = x
      )|>
      modify_fmt_fun(estimate ~ function(x)
        ifelse(is.na(x), NA, "—"),
        rows = n_event == "0") |>
      modify_fmt_fun(ci ~ function(x)
        ifelse(is.na(x), NA, "—"),
        rows = n_event == "0") 
  }
  
  results_competing_risk_function <-
    map2(.x = var_name, .y = var_label, competing) |>
    tbl_stack() |>
    modify_column_hide(p.value) |>
    modify_header(n_event ~ "**Number of events**",
                  exposure ~ "**Person-years**",
                  estimate ~ "**HR**") |>
    modify_table_body(~ .x |> dplyr::relocate(exposure, .before = estimate))
  
  if (!is.null(model_name)) {
    results_competing_risk_function = results_competing_risk_function |>   modify_spanning_header(c(n_event, exposure, estimate, ci) ~
                                                                              paste0("**", model_name, "**"))
  }
  return(results_competing_risk_function)
  rm(results_competing_risk_function)
}

#6 E-value
evalues.HR <- function(est, lo = NA, hi = NA, rare = NA, true = 1, ...) {
  
  # Sanity checks
  if (est < 0) stop("HR cannot be negative")
  
  if (is.na(rare)) rare <- NULL # for compatibility w/ HR constructor
  
  if (!inherits(est, "HR")) est <- HR(est, rare = rare)
  if (!is.na(lo) && !inherits(lo, "HR")) lo <- HR(lo, rare = attr(est, "rare"))
  if (!is.na(hi) && !inherits(hi, "HR")) hi <- HR(hi, rare = attr(est, "rare"))
  if (!inherits(true, "HR")) true <- HR(true, rare = attr(est, "rare"))
  
  est <- toRR(est)
  if (!is.na(lo)) lo <- toRR(lo)
  if (!is.na(hi)) hi <- toRR(hi)
  true <- toRR(true)
  
  return(evalues.RR(est = est, lo = lo, hi = hi))
}

add_evalues <- function(data) {
  evalues_output <- apply(data, 1, function(row) {
    if (!is.na(row["conf.high"]) && row["conf.high"] > 0) {  
      evalues_result <- as.data.frame(evalues.HR(est = row["estimate"], lo = row["conf.low"], hi = row["conf.high"], rare = TRUE))
      evalues_result <- evalues_result[-1, ]
      c(point = evalues_result$point, lower = evalues_result$lower, upper = evalues_result$upper)
    } else {
      c(point = NA, lower = NA, upper = NA)
    }
  })
  
  evalues_output <- as.data.frame(t(evalues_output))
  evalues_output$Evalues_confidence = ifelse(is.na(evalues_output$lower),evalues_output$upper,evalues_output$lower)
  # Add new columns to the data frame
  data <- cbind(data, evalues_output$point,evalues_output$Evalues_confidence)
  colnames(data)[(ncol(data) - 1):(ncol(data))] <- c("Evalues_estimate", "Evalues_confidence")
  
  return(data)
}

#7 Prepare for forest plot (strat long) for childhood paper
change_numbers_next_value=function(numbers) {
  for (i in 1:(length(numbers))) {
    if (i == length(numbers)) {numbers[i] = NA}
    else {numbers[i] <- numbers[i + 1]}
  }
  return(numbers)
}

change_all_values_next_value <- function(df, columns) {
  for (column in columns) {
    df[[column]] <- change_numbers_next_value(df[[column]])
  }
  return(df)
}

columns_to_change <- c("n_obs", "n_event", "exposure", "estimate", "std.error", 
                       "statistic", "nevent", "conf.low", "conf.high", "ci", "p.value","cases","n_event_ref")

forest_plot_child = function(dat,cat,label1,label2){
  dat = dat$table_body|>group_by(N)|>
    mutate(n_event_ref = ifelse(or(is.na(n_event),row_number()<=2),NA, n_event[2]))|>
    ungroup()|> 
    mutate(cases=case_when(is.na(n_event) ~ "",
                           .default=paste0(n_event,"/",n_event_ref)))|>filter(is.na(n_event) |label == cat)|>
    change_all_values_next_value(columns_to_change)|>
    mutate(cases=case_when(label=="p-interaction" ~ as.character(round(p_interaction,2)),
                           .default=cases))|>
    dplyr::filter(or(label != cat,cases!=""))|>sort_childhood_not_gtsummary()|>
    mutate(ln_estimate=ifelse(n_event==0,NA,log(estimate)))|>
    mutate(ln_original_standard_error=ifelse(n_event==0,NA,(log(conf.high)-log(conf.low))/3.92))|>
    mutate(group_index = as.integer(factor(p_interaction, levels = unique(p_interaction))))|>
    group_by(group_index)|> 
    group_modify( ~ add_row(.x,label="",.before = 0))|>ungroup()|>
    filter(row_number() !=1)
  
  forest_meta_random=metagen(TE = ln_estimate,seTE = ln_original_standard_error,studlab=paste(label),
                             data=dat,comb.fixed = FALSE,
                             method.tau="DL",hakn=FALSE,prediction = FALSE,sm="HR",
                             random = FALSE, 
                             fixed = FALSE)
  plot=forest(forest_meta_random,print.tau2 = FALSE,col.diamond = "blue",
              label.right = "Risk higher",label.left = "Risk lower",
              colgap.left = "3mm",prediction = FALSE,
              leftcols=c("label","cases"),leftlabs = c("Potential modifying factors",paste0(label1,"\n",label2," (cases)")),
              header.line=TRUE,
              rightcols = c("estimate","ci" ),
              rightlabs = c("HR","95%CI"),colgap.right = "5mm",colgap.forest.right = "5mm",
              resid.hetstat = FALSE,col.by = "black",subgroup = FALSE,
              addrow.overall=TRUE,hetstat = TRUE,overall.hetstat = FALSE,
              overall = FALSE,xlim = c(0.3,3.5), weight.study="same", plotwidth="2inch",
              spacing = 1.3,
              addrow = TRUE,
              fontsize = 14,squaresize =0.7)
  return(plot)
}

#8.1 p-interaction

p_interaction = function(surv, var_name1, var_name2, covar, data, model_name=NULL)
{
  results_cox1 = reformulate(c(var_name1, var_name2, covar), response = surv)|>
    coxph(method = "breslow", data = data) 
  
  results_cox2 = reformulate(c(var_name1, var_name2,paste0(var_name1,"*",var_name2), covar), response = surv)|>
    coxph(method = "breslow", data = data) 
  
  p_interaction = anova(results_cox1,results_cox2)
  return(p_interaction)
}

#8.1 p-trend

p_trend = function(surv, var_name1, covar, data, model_name=NULL)
{
  data=data|>mutate(var_name1=as.numeric(var_name1))
    results_cox1 = reformulate(c(var_name1, covar), response = surv)|>
    coxph(method = "breslow", data = data) 
  
  results_cox2 = reformulate(c( covar), response = surv)|>
    coxph(method = "breslow", data = data) 
  
  p_trend = anova(results_cox1,results_cox2)
  return(p_trend)
}


#9. Replace letter A, B, C, D, K, N, Q, S and . to missing value
replace_values_with_na <- function(dat) {
  values_to_replace <- c("A", "B", "C", "D", "K", "N", "Q", "S", ".")
  dat[] <- lapply(dat, function(x) ifelse(x %in% values_to_replace, NA, x))
  return(dat)
}


#10 Assign "Unknown" to missing value
set_unknown_level <- function(my_factor, unknown_label = "Unknown") {
  # Add "Unknown" to the factor levels if it's not already there
  levels(my_factor) <- c(levels(my_factor), unknown_label)
  
  # Replace NA values with the unknown label
  my_factor[is.na(my_factor)] <- unknown_label
  
  # Return the modified factor
  return(my_factor)
}

#11 Polish tables for in night shift paper
polish_night_shift = function(data, group_name = FALSE)
{
  #polish the table
  if (group_name==TRUE)
  {
    data = data |>
      modify_table_body(
        ~ .x  |>
          dplyr::add_row(
            label = "p-trend",
            .before = grep("night_perm_cumulative_all_cat", data$table_body$variable)
          )  |>
          dplyr::add_row(
            label = "p-trend",
            .before = grep("night_perm_intensity_all_cat", data$table_body$variable) + 1
          )  |>
          dplyr::add_row(
            label = "p-trend",
            .before = grep("night_perm_age_start", data$table_body$variable) + 2
          )
      )
    
  }
   data = modify_footnote(data,
                          abbreviation = TRUE)|>
     as_gt() |> tab_footnote(footnote = "Model 1 was adjusted for sex, race, and the calendar year of employment. Model 2 was additionally adjusted for BMI, smoking status, and levels of physical activity")
  return(data)
}


#6. function to sort ui_birth_weight properly
#sort_birth_weight = function(table) {table|>
#  modify_table_body( ~.x |>
#                       dplyr::mutate(sort = 1:nrow(.x)) |>
#                        dplyr::mutate(sort = if_else(label == "Between 2500 and 3999 g",
#                                                     sort + 1,
#                                                     if_else(label =="< 2500 g", sort - 1, sort)
#                        )) |> arrange(sort))
# }
#7 function to sort properly rep_br_dev_age_cat, rep_age_menarche_cat, an_weight_age10, an_height_age10, an_weight_teen, ses_income_child
# sort_childhood = function(table) {table|>
#     modify_table_body( ~.x |>
#                          dplyr::mutate(sort = 1:nrow(.x)) |>
#                          dplyr::mutate(sort = if_else(label == "11-13 years of age",  sort + 1,
#                                                       if_else(label =="Less than 11 years of age", sort - 1, sort)))|>
#                          dplyr::mutate(sort = if_else(label == "12-14 years of age",  sort + 1,
#                                                       if_else(label =="Less than 12 years of age", sort - 1, sort)))|>
#                          dplyr::mutate(sort = if_else(label == "Same weight",  sort + 1,
#                                                       if_else(label =="Lighter", sort - 1, sort)))|>
#                          dplyr::mutate(sort = if_else(label == "Same height",  sort + 1,
#                                                       if_else(label =="Shorter", sort - 1, sort)))|>
#                          dplyr::mutate(sort = if_else(label == "Middle income",  sort + 1,
#                                                       if_else(label =="Well off", sort - 1, sort))) |> 
#                          arrange(sort))
# }
# 
# sort_childhood_not_gtsummary = function(table) {table|>
#                          dplyr::mutate(sort = 1:nrow(table)) |>
#                          dplyr::mutate(sort = if_else(label == "11-13 years of age",  sort + 1,
#                                                       if_else(label =="Less than 11 years of age", sort - 1, sort)))|>
#                          dplyr::mutate(sort = if_else(label == "12-14 years of age",  sort + 1,
#                                                       if_else(label =="Less than 12 years of age", sort - 1, sort)))|>
#                          dplyr::mutate(sort = if_else(label == "Same weight",  sort + 1,
#                                                       if_else(label =="Lighter", sort - 1, sort)))|>
#                          dplyr::mutate(sort = if_else(label == "Same height",  sort + 1,
#                                                       if_else(label =="Shorter", sort - 1, sort)))|>
#                          dplyr::mutate(sort = if_else(label == "Middle income",  sort + 1,
#                                                       if_else(label =="Well off", sort - 1, sort))) |> 
#                          arrange(sort)
# }
# 
# #9 Polish tables for childhood paper
# polish_childhood = function(data, group_name = FALSE) 
# {
#   #polish the table
#   if (group_name==TRUE)
#   {  
#     data = data |>
#     modify_table_body(
#       ~ .x  |>
#         dplyr::add_row(
#           label = "Growth and reproductive factors",
#           .before = grep("an_weight_age10", data$table_body$variable)
#         )  |>
#         dplyr::add_row(
#           label = "Lifestyle factors",
#           .before = grep("ph_MET_wk_before20_cat", data$table_body$variable) + 1
#         )  |>
#         dplyr::add_row(
#           label = "Socioeconomic factors",
#           .before = grep("ses_income_child", data$table_body$variable) + 2
#         ) |>
#         dplyr::filter(!label %in% c("Hormonal birth control duration under age 20","Started after 20","Unknown status",
#                                     "Total years of smoking before 20","Smoking before 20 (Pack-years)","Never smoked")))
#   data=data|>  modify_table_body(~ .x  |>    dplyr::mutate(sort = 1:nrow(.x)) |>
#         dplyr::mutate(sort = if_else(and(variable=="sm_age_start_smok_cat",label == "20 years of age or older"),  sort + 4,
#                                      if_else(label =="Unknown smoking status", sort + 3, sort)))|>arrange(sort))
#      
#   
#   data = data |>
#     modify_table_styling(
#       column = "label",
#       rows = data$table_body$label %in% c(
#         "Participant characteristics at baseline",
#         "Growth and reproductive factors",
#         "Lifestyle factors",
#         "Socioeconomic factors"
#       ),
#       text_format = "bold"
#     )
#   }
#   data = modify_footnote(data, abbreviation = TRUE) |>
#     as_gt() |> tab_footnote(
#       footnote = "Multivariable models were adjusted for attained age (timescale), and self-identified race/ethnicity"
#     )
#   
# }




#12. Remove the numbering
remove_numbering_in_levels <- function(dat) {
  for (i in colnames(dat)) {
    if (is.factor(dat[, i])) {
      levels(dat[, i]) <- gsub("^\\d+\\.?\\s*\\)\\s*", "", levels(dat[, i]), perl = TRUE)
    }
  }
  return(dat)
}
#13. factored_case_when
factored_case_when<- function(...) {
args <- list2(...)
rhs <- map(args, f_rhs)

cases <- case_when(
  !!!args
)

exec(fct_relevel, cases, !!!rhs)
}

