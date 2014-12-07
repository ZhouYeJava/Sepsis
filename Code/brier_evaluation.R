#!/usr/bin/Rscript

####################
###Brier Analysis###
####################

###Zhou Ye###
###06/16/2014###

args <- commandArgs(T)
result_path <- args[1]
#num_group <- as.integer(args[2])
suppressMessages(library(foreach))
suppressMessages(library(doMC))
#suppressMessages(library(ResourceSelection))
registerDoMC(50)

###Combine Files###
combine_files <- function(files) {
	data <- foreach (f = files, .combine="rbind") %dopar% {
		read.csv(file=f, header=F)
	}
	data <- data.frame(data[,2:ncol(data)])
	colnames(data) <- c("prediction", "label")
	return(data)
}

###Obtain Brier Score###
obtain_brier <- function(data) {
	#hl <- hoslem.test(x=data$label, y=data$prediction, g=num_group)
	#return(hl$statistic)
	return(mean((data$label-data$prediction)^2))
}

###Evaluation###
brier <- vector(length=5)
for (i in 1:5) {
	print(i)
	files <- list.files(path=result_path, pattern=paste("fold_", i, "_[0-9]+.csv", sep=""), full.names=T)
	brier[i] <- obtain_brier(combine_files(files))
}
print(mean(brier))
print(sd(brier))
