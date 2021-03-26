#--------------------------------------------------------------------------#
# Script for calculating Sobol indices using model M3 as well as
# a Latin hypercube design for the parameter values using a log scale
# Authors:      Kai Budde
# Created:      2021-02-23
# Last changed: 2021-03-26
# Used version of Sessl2R: 0.1.3
#--------------------------------------------------------------------------#

# Preface ##################################################################
remove(list = ls())

library(sensitivity)
library(lhs)
library(ggplot2)

# Install the R package for using Sessl output data
devtools::install_github("SFB-ELAINE/Sessl2R", ref = "v0.1.3")
require(Sessl2R)

directory_of_R_script <- rstudioapi::getSourceEditorContext()$path
directory_of_R_script <- gsub(pattern = "ExecuteSesslandPlotResults_M3.R",
                              replacement = "",
                              x = directory_of_R_script)
old_directory <- getwd()
setwd(directory_of_R_script)

# Parameter set ############################################################

# Start taking the time
time_begin <- Sys.time()

# Subdirectory with SESSL script
dir_sessl_script <- "experiments"

# Calculate design of experiment
n <- 350
# Number of design points for simulation: n*(number of parameters + 2)

parameter_names <- c("kLRAss", "ke_nonraft", "ke_raft")

set.seed(42)
LHS_design1 <- lhs::create_oalhs(n = n, k = length(parameter_names), bChooseLargerDesign = TRUE, bverbose = FALSE)
LHS_design1 <- data.frame(LHS_design1)
names(LHS_design1) <- parameter_names
# Rearrange samples for correct results
LHS_design1 <- LHS_design1[sample(1:length(LHS_design1[[1]]), replace = FALSE),]

set.seed(1234)
LHS_design2 <- lhs::create_oalhs(n = n, k = length(parameter_names), bChooseLargerDesign = TRUE, bverbose = FALSE)
LHS_design2 <- data.frame(LHS_design2)
names(LHS_design2) <- parameter_names
# Rearrange samples for correct results
LHS_design2 <- LHS_design2[sample(1:length(LHS_design2[[1]]), replace = FALSE),]

# LHS design with min and max values
LHS_design1$kLRAss <- qunif(LHS_design1$kLRAss, min = 1e08, max=1e10)
LHS_design2$kLRAss <- qunif(LHS_design2$kLRAss, min = 1e08, max=1e10)

LHS_design1$ke_nonraft <- qunif(LHS_design1$ke_nonraft, min = 0.01, max=0.5)
LHS_design2$ke_nonraft <- qunif(LHS_design2$ke_nonraft, min = 0.01, max=0.5)

LHS_design1$ke_raft <- qunif(LHS_design1$ke_raft, min = 0.01, max=0.5)
LHS_design2$ke_raft <- qunif(LHS_design2$ke_raft, min = 0.01, max=0.5)

# Call Sessl script with the parameter values of x$X
run_sessl_experiment <- function(design) {
  
  setwd(dir_sessl_script)
  
  number_of_design_points <- nrow(design)
  
  write.csv(x = design, file = "designOfExperiment.csv", row.names = FALSE)
  
  # Run Sessl script
  if(.Platform$OS.type == "unix") {
    system(command = "./run.sh", wait = TRUE)
  } else {
    system(command = "./run.bat", wait = TRUE)
  }
  
  # Find most recent results directory
  
  subdirectories <- list.dirs(recursive = FALSE)
  experiment_result_dirs <- grep(pattern = "results",
                                 x = subdirectories,
                                 value = TRUE)
  recent_result <- experiment_result_dirs[length(experiment_result_dirs)]
  
  # Load data from experiments
  # Save a dataframe with the experiment result
  print(paste0("Folder with most recent SESSL results: ", recent_result))
  
  setwd(recent_result)
  Sessl2R::getData(input_dir = getwd())
  
  # Load the data frame and copy loaded df to default name
  df_files <- list.files()
  df_files <- df_files[grepl(pattern = "rda", x = df_files)]
  
  # Go through every data frame saved in the directory and store the results
  # in df_complete
  for(j in 1:length(df_files)){
    
    load(df_files[j])
    loaded_dataframe <- gsub("\\.rda", "", df_files[j])
    
    df <- get(loaded_dataframe)
    
    if(j > 1){
      df_complete <- rbind(df_complete, df)
    }else{
      df_complete <- df
    }
    
    rm(loaded_dataframe)
  }
  
  setwd(directory_of_R_script)
  
  df_complete <- df_complete[order(df_complete$config),]
  df_complete <- df_complete[order(df_complete$run),]
  
  return(df_complete)
}


#Sensitivity analysis
soboljansen_result <- sensitivity::soboljansen(
  model = NULL, X1 = LHS_design1, X2 = LHS_design2,
  nboot = 100, conf = 0.95)

df_complete <- run_sessl_experiment(soboljansen_result$X)

simulation_result <- (df_complete$Lrp6 + df_complete$Lrp6Axin) / 4000
# (4000 is the number of total LRP6 in the membrane in the beginning)

sensitivity::tell(soboljansen_result, simulation_result)
sensitivity_plot <- ggplot2::ggplot(soboljansen_result)

# Save results
dir.create("results", showWarnings = FALSE)
setwd("results")
filename <- paste(gsub(pattern = " ", replacement = "_", x = Sys.time()),
                  "_SoboljansenResult_n", (dim(LHS_design1)[2]+2)*n, ".Rda", sep="")
save(soboljansen_result, file = filename)
ggplot2::ggsave(filename = gsub(pattern = "\\.Rda", replacement = "\\.pdf",
                                x = filename), plot = sensitivity_plot, width = 12,
                height = 6, device = "pdf")

# Stop taking the time
time_end <- Sys.time()
time_difference <- round(time_end-time_begin, digits = 1)
print(paste(time_difference, " ", units(time_difference), " for ",
            (dim(LHS_design1)[2]+2)*n, " experiments.", sep=""))

setwd(old_directory)
