#This is a script that will use the aneuQPCR functions to calculate the efficiency and extrapolate initial fluoresceces for each well in a qPCR.
#This expects a tsv file exported from the LightCycler480 software. 

#set working directory to the place where the script is located (only works with rstudio)
if (!requireNamespace("rstudioapi", quietly = TRUE)) {
  install.packages("rstudioapi")
}

if (rstudioapi::isAvailable()) {
  script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
  setwd(script_dir)
}

#load main functions
source("aneuQPCR_functions.R")

#get amplification curves
message("Select the amplification curves file.")
amp_curves <- file.choose()
amp_curves <- as.data.frame(read_delim(amp_curves,  locale = locale(decimal_mark = ",")))

#set number of amplification cycles
ncycles <- 40
threshold <- 1

#format the table to convert it to a cycle X Well matrix
rownames(amp_curves) <- amp_curves$Index #set the cycle number to the row name
amp_curves <- amp_curves[,!grepl("^X\\.\\.\\..*|^Index$|^Text$", colnames(amp_curves))] #remove unecessary columns
colnames(amp_curves) <- gsub(":.*", "", colnames(amp_curves)) #keep colnames as just the well positions
#amp_curves <- apply(amp_curves, 2, function(x) as.numeric(gsub(",", ".", x))) #replace comma by dot and convert values to numeric

#create a dummy well_meta dataframe as it is required by the qPCR object 
well_meta <- data.frame(well = colnames(amp_curves))

#build the qPCR object
qPCRobj <- create_qPCRobj(amp_curves, well_meta, ncycles, melting_curve = TRUE)

#calculate efficiency and extrapolation of initial fluorescence with LL5 modeling. 
qPCRobj <- suppressWarnings(fit_model_LL5(qPCRobj))

#calculate Ct value with simple threshold line
qPCRobj <- calc_Cts(qPCRobj, threshold = threshold)

#check if threshold line is well set
plot_amps(qPCRobj, threshold_line = threshold, scale = "log2")

#save the excel file
xlsx::write.xlsx(qPCRobj$metadata$well_meta, file = "per_well_calculations.xlsx", row.names = FALSE)
message(paste("results stored at", paste0(getwd(), "/per_well_calculations.xlsx")))

#plot example
plot <-plot_model(qPCRobj$models$LL5$F1)+ggtitle("Log2 scale", subtitle = NULL)
ggsave("figures/LL5_model_log2_scale.png", plot = plot, width = 6, height = 5)
#plot example
plot <-plot_model(qPCRobj$models$LL5$F1, scale = "original")+ggtitle("Original scale", subtitle = NULL)
ggsave("figures/LL5_model_original_scale.png", plot = plot, width = 6, height = 5)