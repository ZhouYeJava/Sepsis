#!/usr/bin/Rscript

###################
###Time Analysis###
###################

###Zhou Ye###
###06/16/2014###

args <- commandArgs(T)
result_path <- args[1]
time <- as.integer(args[2])
threshold <- as.numeric(args[3])
suppressMessages(library(foreach))
suppressMessages(library(doMC))
registerDoMC(50)

###Obtain Threshold###
obtain_threshold <- function(files) {
    data <- foreach (f = files, .combine="rbind") %dopar% {
        read.csv(file=f, header=F)
    }
    data <- data.frame(data[,2:ncol(data)])
    colnames(data) <- c("prediction", "label")
    threshold <- seq(min(from=data$prediction), to=max(data$prediction), length=100)
    acc <- foreach (i = threshold, .combine="c") %dopar% {
	pred <- as.numeric(data$prediction>i)
	sum(pred==data$label)/nrow(data)
    }
    return(threshold[which.max(acc)])
}

###Compute Sensitivity/Specificity###
compute_stat <- function(files, threshold) {
    data <- foreach (f = files, .combine="rbind") %dopar% {
        read.csv(file=f, header=F)
    }
    data <- data.frame(data[,2:ncol(data)])
    colnames(data) <- c("pred", "label")
    data$pred <- as.numeric(data$pred>threshold)
    sen <- length(which(data$pred==1 & data$label==1))/length(which(data$label==1))
    spe <- length(which(data$pred==0 & data$label==0))/length(which(data$label==0))
    return(c(sen, spe))
}

###Compute Time###
compute_time <- function(files, threshold) {
	time_to_event <- foreach (f = files, .combine="c") %dopar% {
		data <- read.csv(file=f, header=F)
		pred <- as.numeric(data[,2]>threshold)
		if (any(data[,3]==1)) {
		    if (nrow(data)<=time) {
                        detect <- 0
                        for (t in (nrow(data)-1):0) {
                            if (pred[nrow(data)-t]==1) {
                                detect <- t
                                break
                            }
                        }
                        detect
		    }
	       	    else {
            	        detect <- 0
           		for (t in time:0) {
                	    if (pred[nrow(data)-t]==1) {
                       	        detect <- t
                                break
                	    }
                 	}
            	        detect
		    }
		}
		else {
			NA
		}
	}
	return(mean(time_to_event, na.rm=T))
}

###Evaluation###
detect_time <- vector(length=5)
sensitivity <- vector(length=5)
specificity <- vector(length=5)
for (i in 1:5) {
	print(i)
	files <- list.files(path=result_path, pattern=paste("fold_", i, "_[0-9]+.csv", sep=""), full.names=T)
	temp <- compute_stat(files, threshold)
	sensitivity[i] <- temp[1]
	specificity[i] <- temp[2]
	detect_time[i] <- compute_time(files, threshold)
}
print(mean(sensitivity))
print(sd(sensitivity))
print(mean(specificity))
print(sd(specificity))
print(mean(detect_time))
print(sd(detect_time))
