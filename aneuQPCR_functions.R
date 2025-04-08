#to do
#1) fix the function plot_model() and display on the plot the linear phase as well as the inflection point.
#2) make a function to plot the somies.
#3) simplify the scrip for other users to use it.

#libraries####
library(readr)
library(readxl)
library(dplyr)
library(tidyr)
library(tibble)
library(drc)
library(scales)
library(ggrepel)

#create_qPCRobj
create_qPCRobj <- function(amp_curves, well_meta, ncycles, melting_curve = FALSE){
  #format well_meta
  rownames(well_meta) <- well_meta$well
  #remove wells not present in the well meta
  nwells_before <- ncol(amp_curves)
  amp_curves <- amp_curves[,rownames(well_meta)]
  message(paste("removed", nwells_before - ncol(amp_curves), "wells which were not present in well meta"))
  
  if(melting_curve){
    melting_curves <- amp_curves[ncycles+1:nrow(amp_curves),] #assumes melting curves are the subsequent "cycles" after the PCR is done.
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
  well_meta$efficiency <- NA
  well_meta$extrapolated_start <- NA
  well_meta$inflection_point <- NA
  models <- list()
  for(i in c(1:ncol(df))){
    well <- colnames(df[,i, drop = FALSE])
    to_model <- data.frame(cycle = as.numeric(rownames(df)), fluorescence = as.numeric(df[,i]))
    #log2 transformation makes it easier to find the linear phase
    to_model$fluorescence <- log2(to_model$fluorescence)
    model <- try(drm(fluorescence ~ cycle, data = to_model, fct = LL.5()))
    
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
        }
        
      }
      
    }
    models[[well]] <- model
    well_meta[well,]$LL5_model_Rsquared <- r_squared 
    well_meta[well,]$efficiency <- efficiency
    well_meta[well,]$extrapolated_start <- extrapolated_start
    well_meta[well,]$inflection_point <- as.numeric(inflection_point)
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
  
  well_info <- qPCRobj$metadata$well_meta
  well_info$Ct_value <- NA
  
  for(well in names(qPCRobj$models[[model]])){
    model_to_use <- qPCRobj$models[[model]][[well]]
    #if for some reason there is no model for the well, return NA
    if(class(model_to_use) != "drc"){
      Ct_value <- NA
    }else{
      #if the maximum value in the model doens't reach the threshold, also return NA
      if(max(predict(model_to_use), na.rm = TRUE) < threshold){
        Ct_value <- NA
      }else{
        cycle_range <- model_to_use$data$cycle
        Ct_value <- optimize(objective_function, 
                             interval = cycle_range, 
                             threshold = threshold, 
                             model = model_to_use)$minimum
      }
    }
      well_info[well,]$Ct_value <- Ct_value
  }
  
  #warn if a Ct value is bigger than the inflection point in the qPCR
  test <- well_info %>%
    filter(!is.na(Ct_value)) %>%
    mutate(test = Ct_value > inflection_point) %>%
    pull("test") %>%
    sum()
  
  if(test > 0){
    warning("Some Ct values are higher then their respective inflection points. Consider lowering the threshold.")
  }
  
  qPCRobj$metadata$well_meta <- well_info
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
  
  well_info <- qPCRobj$metadata$well_meta
  
  sample_meta <- well_info %>%
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
plot_amps <- function(qPCRobj, threshold_line = NULL, group_by = "target", split_by = "sample"){
  library(dplyr)
  library(ggplot2)
  
  wells_info <- qPCRobj$metadata$well_meta
  raw_data <- qPCRobj$data$amp_curves
  
  
 plot <- raw_data %>%
    t() %>%
    as.data.frame() %>%
    rownames_to_column("well") %>%
    cbind(wells_info[colnames(raw_data),-which(tolower(colnames(wells_info)) == "well")]) %>%
    pivot_longer(cols = any_of(rownames(raw_data)), names_to = "cycle", values_to = "fluorescence") %>%
    mutate(cycle = as.integer(cycle), fluorescence = as.numeric(fluorescence)) %>%
    ggplot(aes(x = cycle, y = fluorescence, color = .data[[group_by]], group = well))+
    geom_line()+
   {if(is.numeric(threshold_line))geom_hline(yintercept=threshold_line, linetype = "dashed")}+
   facet_wrap(vars(.data[[split_by]]))
 
 return(plot)
}

plot_model <- function(model, scale = c("log2", "original")){
  scale <- match.arg(scale)
  
  to_plot <- model$data[,c(1,2)]
  
  coefficients <- coef(model)
  #this is here just to make my life easier
  names(coefficients) <- c("slope", "lower_asymptope", "upper_asymptope", "inflection_point", "assymetry")
  
  #get the linear part and the model for that part
  linear_part <- model$linear_part
  lmodel <- model$lmodel
  
  #now to extrapolate the linear model, we need to find the cycles where the linear model crosses both asymptopes
  #this subfunction basically calculates the distance between the predicted fluorescence value in the linear `model` and the fluorecence in the asymptope (threshold)
  objective_function <- function(cycle, threshold, model) {
    return(abs(predict(model, newdata = data.frame(cycle = cycle)) - threshold)) 
  }
  #find the cycle where the linear model croses the lower asymptope
  min_cycle <- optimize(objective_function, 
                              interval = model$data$cycle, 
                              threshold = coefficients["lower_asymptope"], 
                              model = lmodel)$minimum
  #find the cycle where the linear model croses the upper asymptope
  max_cycle <- optimize(objective_function, 
                        interval = model$data$cycle, 
                        threshold = coefficients["upper_asymptope"], 
                        model = lmodel)$minimum
  
  #once we have the cycles, we can predict the values inside those cycles to display the extrapolated linear model 
  extrapolated <- data.frame(cycle = c(min_cycle, max_cycle)) %>%
    mutate(fluorescence = predict(lmodel, data.frame(cycle)))
 
  plot <- to_plot %>%
  ggplot(aes(x = cycle, y = fluorescence))+
    geom_line(aes(y = fitted(model)), linewidth = 0.5)+
    {if(scale == "log2")geom_line(data = extrapolated, color = "blue", linewidth = 1, linetype = "dashed")}+
    geom_line(data = linear_part, color = "red", linewidth = 1)+
    geom_point()+
    {if(scale == "original")geom_label_repel(data = extrapolated[which(extrapolated$cycle == min(extrapolated$cycle)),], aes(label = round(2^fluorescence, 3)), color = "blue", nudge_y = 2, angle = 90)}+
    {if(scale == "original")scale_y_continuous(transform = scales::trans_new("reverse_log2", transform = function(x)2^x, inverse = function(x)log2(x)),)}+
    labs(subtitle = "red = linear amplification part; dashed blue = extrapolation of linear amplification, dashed green = inflection point")+
    geom_hline(yintercept = c(coefficients["lower_asymptope"], coefficients["upper_asymptope"]), color = "red", linetype = "dashed")+
    geom_vline(xintercept = c(coefficients["inflection_point"]), color = "green", linetype = "dashed")
  return(plot)
}

