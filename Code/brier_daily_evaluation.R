#!/usr/bin/Rscript

##########################
###Brier Daily Analysis###
##########################

###Zhou Ye###
###07/10/2014###

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
	data <- data.frame(data)
	colnames(data) <- c("time", "prediction", "label")
	return(data)
}

###Obtain Brier Score###
obtain_brier <- function(data) {
	#hl <- hoslem.test(x=data$label, y=data$prediction, g=num_group)
	#return(hl$statistic)
	return(mean((data$label-data$prediction)^2))
}

###Evaluation###
brier <- matrix(NA, 5, 5) #five days and five folds
for (i in 1:5) {
        print(i)
        files <- list.files(path=result_path, pattern=paste("fold_", i, "_[0-9]+.csv", sep=""), full.names=T)
        data <- combine_files(files)
        for (j in 1:5) {
                if (j<5) {
                        d <- subset(data, data$time>=24*(j-1)+1 & data$time<=24*j)
                }
                else {
                        d <- subset(data, data$time>=24*(j-1)+1)
                }
                brier[i,j] <- obtain_brier(d)
        }
}
print(apply(brier, 2, mean))
print(apply(brier, 2, sd))
