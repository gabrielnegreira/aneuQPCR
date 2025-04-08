#to do
#1) fix the function plot_model() and display on the plot the linear phase as well as the inflection point.
#2) make a function to plot the somies.
#3) simplify the scrip for other users to use it.

#libraries####
# List of required packages
required_packages <- c("readr", "xlsx", "dplyr", "tidyr", "tibble",
                       "drc", "scales", "ggrepel")

# Check which packages are not installed
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

# Install missing packages
if (length(missing_packages) > 0) {
  install.packages(missing_packages)
} else {
  message("All required packages are already installed.")
}

# Load all packages
invisible(lapply(required_packages, library, character.only = TRUE))

#create_qPCRobj
create_qPCRobj <- function(amp_curves, well_meta, ncycles, melting_curve = FALSE){
  #format well_meta
  rownames(well_meta) <- well_meta$well
  
  #format the amp_curves
  amp_curves <- as.matrix(amp_curves)
  #remove wells not present in the well meta
  nwells_before <- ncol(amp_curves)
  amp_curves <- amp_curves[,rownames(well_meta)]
  removed_wells <- nwells_before - ncol(amp_curves) 
  if(removed_wells > 0){
    message(paste("removed", removed_wells, "wells which were not present in well meta"))  
  }
  
  #ensure rownames represent cycles
  rownames(amp_curves) <- c(1:nrow(amp_curves))
  
  if(melting_curve){
    melting_curves <- amp_curves[(ncycles+1):nrow(amp_curves),] #assumes melting curves are the subsequent "cycles" after the PCR is done.
    amp_curves <- amp_curves[c(1:ncycles),]
    raw_data <- list(amp_curves = amp_curves, melting_curves = melting_curves)
  }else{
    amp_curves <- amp_curves[c(1:ncycles),]
    raw_data <- list(amp_curves = amp_curves)
  }
  
  qPCRobj <- list(data = raw_data, metadata = list(well_meta = well_meta))
  
  return(qPCRobj)
}

#fit_model_LL5####
#this function will fit five-parameter logistic (5PL) model to each PCR well.
#these models are used by subsequent functions to calculate Ct and amplification efficiency.
fit_model_LL5 <- function(qPCRobj){
  df <- qPCRobj$data$amp_curves
  well_meta <- qPCRobj$metadata$well_meta
  well_meta$LL5_model_Rsquared <- NA
  well_meta$linear_part_lm_Rsquared <- NA
  well_meta$efficiency <- NA
  well_meta$extrapolated_start <- NA
  well_meta$inflection_point <- NA
  well_meta$start_linear_part <- NA
  well_meta$end_linear_part <- NA
  models <- list()
  for(i in c(1:ncol(df))){
    well <- colnames(df[,i, drop = FALSE])
    to_model <- data.frame(cycle = as.numeric(rownames(df)), fluorescence = as.numeric(df[,i]))
    #log2 transformation makes it easier to find the linear phase
    to_model$fluorescence <- log2(to_model$fluorescence)
    model <- try(drm(fluorescence ~ cycle, data = to_model, fct = LL.5()), silent = TRUE)
    
    #if model fails, mark it (usually happens with negative samples)
    if(class(model) == "try-error"){
      model <- "not defined"
      efficiency <- NA
      inflection_point <- NA
    }else{
      #determine the R2 of the model
      residuals <- resid(model)
      sst <- sum((to_model$fluorescence - mean(to_model$fluorescence))^2)
      ssr <- sum(residuals^2)
      r_squared <- 1 - (ssr / sst)
      
      #determine the linear part based on derivatives
      start_cycle <- as.numeric(rownames(model$deriv1)[which(model$deriv[,1] == max(model$deriv[,1]))]) #cycle with the maximum derivative value
      end_cycle <- as.numeric(rownames(model$deriv1)[which(model$deriv[,1] == min(model$deriv[,1]))]) #cycle with the lowest derivative value
      
      model$linear_part <- to_model %>%
        filter(cycle >= start_cycle & cycle <= end_cycle)
      inflection_point <- coef(model)[4]
      
      #calculate the efficiency
      if(nrow(model$linear_part) < 2){
        efficiency <- NA
        extrapolated_start <- NA
      }else{
        if(inflection_point < model$linear_part$cycle[2]){
          efficiency <- NA
          extrapolated_start <- NA
        }else{
          lmodel <- lm(fluorescence ~ cycle, filter(model$linear_part, cycle < inflection_point))
          slope <- coef(lmodel)[2] #this is the efficiency in the linear amplification phase
          intercept <- coef(lmodel)[1] #this is the initial fluorescence of the target DNA
          efficiency <- 2^slope #since the LL5 model was build at log2 scale, we need to exponentiate it
          extrapolated_start <- 2^intercept #same here
          model$lmodel <- lmodel
          well_meta[well,]$linear_part_lm_Rsquared <- summary(lmodel)$r.squared
        }
        
      }
      
    }
    models[[well]] <- model
    well_meta[well,]$LL5_model_Rsquared <- r_squared 
    well_meta[well,]$efficiency <- efficiency
    well_meta[well,]$extrapolated_start <- extrapolated_start
    well_meta[well,]$inflection_point <- as.numeric(inflection_point)
    well_meta[well,]$start_linear_part <- as.numeric(start_cycle)
    well_meta[well,]$end_linear_part <- as.numeric(end_cycle)
  }
  qPCRobj$metadata$well_meta <- well_meta
  qPCRobj$models$LL5 <- models
  return(qPCRobj)
}

#calc_Cts####
#this function will calculate Ct values given a threshold
calc_Cts <- function(qPCRobj, threshold, model = "LL5"){
  #sub functions
  #this function simply return the predicted values minus the threshold. This is used in the optimize function to find
  #the cycle that minimizes the distance between the fitted values and the threshold.
  
  objective_function <- function(cycle, threshold, model) {
    return(abs(predict(model, newdata = data.frame(cycle = cycle)) - threshold))
  }
  
  if(!model %in% names(qPCRobj$model)){
    stop(paste("Could not find model named", model, "in the qPCR object provided"))
  }
  
  well_meta <- qPCRobj$metadata$well_meta
  well_meta$Ct_value <- NA
  
  for(well in names(qPCRobj$models[[model]])){
    model_to_use <- qPCRobj$models[[model]][[well]]
    #if for some reason there is no model for the well, return NA
    if(class(model_to_use) != "drc"){
      Ct_value <- NA
    }else{
      #if the maximum value in the model doens't reach the threshold, also return NA
      if(suppressWarnings(max(predict(model_to_use), na.rm = TRUE)) < log2(threshold)){
        Ct_value <- NA
      }else{
        cycle_range <- model_to_use$data$cycle
        Ct_value <- suppressWarnings(optimize(objective_function, 
                                              interval = cycle_range, 
                                              threshold = threshold, 
                                              model = model_to_use)$minimum)
      }
    }
    well_meta[well,]$Ct_value <- Ct_value
  }
  
  #warn if it is not possible to calculate Ct value for any well
  if(sum(!is.na(well_meta$Ct_value)) == 0){
    warning("Could not determine the Ct value for any well. Is the threhold too high?")
  }
  
  #warn if a Ct value is bigger than the inflection point in the qPCR
  test <- well_meta %>%
    filter(!is.na(Ct_value)) %>%
    filter(LL5_model_Rsquared > 0.99)
  
  if(sum(test$Ct_value > test$inflection_point) > 0){
    warning("Some Ct values are higher than their respective inflection points. Consider lowering the threshold.")
  }
  if(sum(test$Ct_value < test$start_linear_part) > 0){
    warning("Some Ct values are lower than the beggining of the linear amplification phase. Consider increasing the threshold.")
  }
  
  qPCRobj$metadata$well_meta <- well_meta
  return(qPCRobj)
}

#calc_copy_number####
calc_copy_number <- function(qPCRobj, ref_target = "chr36", ref_target_copy = 2, remove_outliers = TRUE){
  library(dplyr)
  
  if(!"LL5" %in% names(qPCRobj$models)){
    qPCRobj <- fit_model_LL5(qPCRobj)
  }
  
  required <- c("Ct_value", "inflection_point", "efficiency")
  
  if(sum(!required %in% colnames(qPCRobj$metadata$well_meta)) != 0){
    qPCRobj <- calc_Cts(qPCRobj)
  }
  
  well_meta <- qPCRobj$metadata$well_meta
  
  sample_meta <- well_meta %>%
    group_by(rep_of) %>%
    filter(!is.na(rep_of))%>%
    mutate(out = Ct_value - median(Ct_value, na.rm = TRUE),
           is_out = abs(out) > 0.2) %>%
    {if(remove_outliers) filter(., !is_out) else .} %>%
    mutate(mean_Ct = mean(Ct_value),
           sd_Ct = sd(Ct_value),
           mean_inflection_point = mean(inflection_point),
           sd_inflection_point = sd(inflection_point),
           mean_extrapolated_start = mean(extrapolated_start),
           sd_extrapolated_start = sd(extrapolated_start),
           mean_efficiency = mean(efficiency),
           sd_efficiency = sd(efficiency),
           rep_group = rep_of,
           n_rep = n(),
           n_out = sum(is_out)) %>%
    ungroup() %>%
    dplyr::select(c("sample", "target", "sample_type", "rep_group", "n_rep", "n_out", contains(c("mean", "sd_")))) %>%
    distinct() 
  
  reference <- sample_meta %>%
    filter(target == ref_target) %>%
    mutate(sample_to_compare = sample) %>%
    dplyr::select(-"sample")
  
  sample_meta$reference_extrapolated_start <- reference[match(sample_meta$sample, reference$sample_to_compare),]$mean_extrapolated_start
  sample_meta$reference_Ct <- reference[match(sample_meta$sample, reference$sample_to_compare),]$mean_Ct
  sample_meta$reference_efficiency <- reference[match(sample_meta$sample, reference$sample_to_compare),]$mean_efficiency
  
  sample_meta <- sample_meta %>%
    mutate(copy_number_by_Ct = ref_target_copy*((reference_efficiency^reference_Ct)/(mean_efficiency^mean_Ct))) %>%
    mutate(copy_number_by_extrapolation = ref_target_copy * (mean_extrapolated_start/reference_extrapolated_start))
  
  qPCRobj$metadata$sample_meta <- sample_meta
  return(qPCRobj)
}


#plot_amps####
#this function simply plots the amplification curves
plot_amps <- function(qPCRobj, threshold_line = NULL, group_by = NULL, split_by = NULL, scale = c("original", "log2")){
  
  scale <- match.arg(scale)
  
  wells_info <- qPCRobj$metadata$well_meta
  raw_data <- qPCRobj$data$amp_curves
  
  to_plot <- raw_data %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("well") %>%
    cbind(wells_info[colnames(raw_data),-which(tolower(colnames(wells_info)) == "well")]) %>%
    pivot_longer(cols = any_of(rownames(raw_data)), names_to = "cycle", values_to = "fluorescence") %>%
    mutate(cycle = as.integer(cycle), fluorescence = as.numeric(fluorescence))
  
  if(is.null(group_by)){
    group_by <- "well"
  }
  
  plot <- to_plot %>%
    ggplot(aes(x = cycle, y = fluorescence, color = .data[[group_by]], group = well))+
    geom_line()+
    {if(scale == "log2")scale_y_continuous(transform = "log2")}+
    {if(is.numeric(threshold_line))geom_hline(yintercept=threshold_line, linetype = "dashed")}+
    {if(group_by == "well")guides(color = "none")}+
    {if(!is.null(split_by))facet_wrap(vars(.data[[split_by]]))}
  
  return(plot)
}

plot_model <- function(model, scale = c("log2", "original")){
  scale <- match.arg(scale)
  
  to_plot <- model$data[,c(1,2)]
  model_fit <-fitted(model)
  coefficients <- coef(model)
  #this is here just to make my life easier
  names(coefficients) <- c("slope", "lower_asymptope", "upper_asymptope", "inflection_point", "assymetry")
  
  #get the linear part and the model for that part
  linear_part <- model$linear_part
  lmodel <- model$lmodel
  
  if(!is.null(lmodel)){
    extrapolated <- data.frame(cycle = c(0:coefficients["inflection_point"])) %>%
      mutate(fluorescence = predict(lmodel, data.frame(cycle))) %>%
      mutate(label = round(fluorescence, digits = 3))
      if(scale == "original"){
        extrapolated <- extrapolated %>%
          mutate(fluorescence = 2^fluorescence) %>%
          mutate(label = format(fluorescence, scientific = TRUE, digits = 3))
      } 
  }
  
  if(scale == "original"){
    to_plot$fluorescence <- 2^to_plot$fluorescence
    linear_part$fluorescence <- 2^linear_part$fluorescence
    model_fit <- 2^model_fit
    coefficients[c("lower_asymptope", "upper_asymptope")] <- 2^coefficients[c("lower_asymptope", "upper_asymptope")]
  }
  
  plot <- to_plot %>%
    ggplot(aes(x = cycle, y = fluorescence))+
    geom_line(aes(y = model_fit), linewidth = 0.5)+
    {if(!is.null(lmodel))geom_line(data = extrapolated, color = "blue", linewidth = 1, linetype = "dashed")}+
    geom_line(data = linear_part, color = "red", linewidth = 1)+
    geom_point()+
    {if(!is.null(lmodel))geom_label_repel(data = extrapolated[which(extrapolated$cycle == min(extrapolated$cycle)),], aes(label = label), color = "blue", nudge_y = diff(range(to_plot$fluorescence))*0.3, angle = 90)}+
    labs(subtitle = "red = linear amplification part; dashed blue = extrapolation of linear amplification, dashed green = inflection point")+
    geom_hline(yintercept = c(coefficients["lower_asymptope"], coefficients["upper_asymptope"]), color = "red", linetype = "dashed")+
    geom_vline(xintercept = c(coefficients["inflection_point"]), color = "green", linetype = "dashed")+
    {if(scale == "log2")labs(y = "Log2(Fluorescence)")}
  return(plot)
}


