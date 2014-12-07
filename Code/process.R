#!/usr/bin/Rscript

#########################
###Process Sepsis Data###
#########################

###Zhou Ye###
###06/13/2014###

rm(list=ls())
setwd("/home/zhouye/Data/")
data_path <- "/udata/PatientNew/"
data_path_shock <- "/udata/ClusteredShockPatientNew/"
data_path_normal <- "/udata/ClusteredNormalPatientNew/"
suppressMessages(library(foreach))
suppressMessages(library(doMC))
registerDoMC(50)
files <- list.files(path=data_path, pattern="*.csv")
feature <- as.matrix(read.table("feature.txt"))
feature_name <- c("time", "sirs", "sepsis", "shock", feature)
feature_new <- c("time", feature)
interval <- 60 #cluster every 5 minutes
threshold <- 0 #minimum sample points

build_new_data <- function(data1, data0) {
	first_time <- data1[1,]$time
	last_time <- data1[nrow(data1),]$time
	num_time <- ceiling((last_time-first_time)/interval)
	time_new <- seq(from=first_time, length=num_time, by=interval)
	data_new <- data.frame(matrix(NA, length(time_new), ncol(data0)-3))
	colnames(data_new) <- feature_new
	data_new$time <- time_new
	for (i in 1:nrow(data_new)) {
		d <- subset(data0, data0$time>=data_new[i,]$time-interval & data0$time<=data_new[i,]$time)
		if (nrow(d)>0) {
			data_new[i,2:ncol(data_new)] <- colMeans(d[,5:ncol(data0)], na.rm=T)
		}
		if (i>1) {
			index_missing <- which(is.na(data_new[i,]))
			data_new[i,index_missing] <- data_new[i-1,index_missing]
		}
	}
	indicator <- rowSums(data_new)
	index_missing <- which(is.na(indicator))
	data_new <- data_new[-index_missing,]
	if (nrow(data_new)==0) {
		return(NULL)
	}
	data_new$time <- 1:nrow(data_new)
	return(data_new)
}

result <- foreach (f = files) %dopar% {
	data <- read.csv(paste(data_path, f, sep=""), header=F)
	colnames(data) <- feature_name
	data_new <- NULL
	if (any(data$sirs==1)) {
		if (any(data$shock==1)) {
			cat(paste(f, " is Septic Shock", "\n", sep=""))
			index_sirs <- which(data$sirs==1)[1] #first time SIRS
			index_shock <- which(data$shock==1)[1] #first time Septic Shock
			data_shock <- data[index_sirs:index_shock,]
			data_new <- build_new_data(data_shock, data)
			if (!is.null(data_new)) {
				if (nrow(data_new)>threshold) {
					write.csv(x=data_new, file=paste(data_path_shock, f, sep=""), row.names=F)
				}
			}
		}
		else {
			cat(paste(f, " is Normal", "\n", sep=""))
			index_sirs <- which(data$sirs==1)[1] #first time SIRS
			data_normal <- data[index_sirs:nrow(data),]
			data_new <- build_new_data(data_normal, data)
			if (!is.null(data_new)) {
				if (nrow(data_new)>threshold) {
					write.csv(x=data_new, file=paste(data_path_normal, f, sep=""), row.names=F)
				}
			}
		}
	}
	data_new
}
