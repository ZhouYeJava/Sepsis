#!/usr/bin/Rscript

##################
###ROC Analysis###
##################

###Zhou Ye###
###06/16/2014###

args <- commandArgs(T)
result_path <- args[1]
suppressMessages(library(foreach))
suppressMessages(library(doMC))
suppressMessages(library(pROC))
suppressMessages(library(verification))
registerDoMC(50)

###Combine Files###
combine_files <- function(files) {
	data <- foreach (f=files, .combine="rbind") %dopar% {
		read.csv(file=f, header=F)
	}
	data <- data.frame(data[,2:ncol(data)])
	colnames(data) <- c("prediction", "label")
	return(data)
}

###Obtain AUC###
obtain_auc <- function(data) {
	return(roc.area(obs=data$label, pred=data$prediction)$A)
}

###Evaluation###
auc <- vector(length=5)
for (i in 1:5) {
	print(i)
	files <- list.files(path=result_path, pattern=paste("fold_", i, "_[0-9]+.csv", sep=""), full.names=T)
	auc[i] <- obtain_auc(combine_files(files))
}
print(mean(auc))
print(sd(auc))
