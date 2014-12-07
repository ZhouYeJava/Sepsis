#!/usr/bin/Rscript

args <- commandArgs(T)
feature_file <- args[1]
patient_data <- args[2]
feature <- as.matrix(read.table(feature_file))
files <- list.files(path=patient_data, pattern="*.csv", full.names=T)

for (f in files) {
    data <- read.csv(f)
    colnames(data) <- c("time", feature)
    write.csv(x=data, file=f, row.names=F)
}
