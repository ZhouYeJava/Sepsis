#!/usr/bin/Rscript

###########################################
###Calculate Statistics For Measurements###
###########################################

###Zhou Ye###
###07/10/2014###

args <- commandArgs(T)
feature_file <- args[1]
patient_data <- args[2]
feature <- as.matrix(read.table(feature_file))
files <- list.files(path=patient_data, pattern="*.csv", full.names=T)
print(paste("Number of Patients = ", length(files), sep=""))

suppressMessages(library(foreach))
suppressMessages(library(doMC))
registerDoMC(50)

all_features <- foreach(f=files, .combine="rbind") %dopar% {
    data <- read.csv(f, header=F)
    colMeans(data[,5:ncol(data)], na.rm=T)
}

feature_mean <- apply(all_features, 2, function(x) {mean(x, na.rm=T)})
feature_sd <- apply(all_features, 2, function(x) {sd(x, na.rm=T)})
statistics <- data.frame(rbind(feature_mean, feature_sd))
colnames(statistics) <- feature
print(statistics)
