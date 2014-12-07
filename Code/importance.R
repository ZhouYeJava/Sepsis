#!/usr/bin/Rscript

########################
###Importance Ranking###
########################

###Zhou Ye###
###07/17/2014###

###Configuration###
args <- commandArgs(T)
shock_path <- args[1] #location of shock patients data
normal_path <- args[2] #location of normal patients data
output_path <- args[3] #location of output profiles
time <- as.integer(args[4]) #time to septic shock
suppressMessages(library(foreach))
suppressMessages(library(doMC))
suppressMessages(library(randomForest))
suppressMessages(library(stringr))
suppressMessages(library(verification))
registerDoMC(50)

###Read Data###
read_data <- function(data_path, flag) {
	files <- list.files(path=data_path, pattern="*.csv")
	data_all <- foreach (f = files) %dopar% {
		identifier <- str_extract(f, "\\d+") #obtain identifier
		list(identifier=identifier, flag=flag, data=read.csv(paste(data_path, f, sep="")))
	}
	return(data_all)
}
shock_patients <- read_data(shock_path, 0)
normal_patients <- read_data(normal_path, 1)
patients <- c(shock_patients, normal_patients) 
print("Finish Processing Data")

###Random Forest###
construct_train <- function(shock_train, normal_train) {
	shock_train_frame <- foreach (patient = shock_train, .combine="rbind") %dopar% {
		data <- patient$data
		label <- vector(length=nrow(data))
		patient_time <- data[,1]
		shock_time <- patient_time[length(patient_time)]
		label[patient_time<shock_time-time] <- 0
		label[patient_time>=shock_time-time] <- 1
		if (nrow(data)>1) {
			cbind(as.matrix(data[,2:ncol(data)]), label)
		}
		else {
			c(as.matrix(data[,2:ncol(data)]), label)
		}
	}
	normal_train_frame <- foreach (patient = normal_train, .combine="rbind") %dopar% {
		data <- patient$data
		if (nrow(data)>1) {
			cbind(as.matrix(data[2:ncol(data)]), rep(0, nrow(data)))
		}
		else {
			c(as.matrix(data[2:ncol(data)]), 0)
		}
	}
	training_frame <- rbind(shock_train_frame, normal_train_frame)
	training_0 <- subset(training_frame, training_frame[,ncol(training_frame)]==0)
	training_1 <- subset(training_frame, training_frame[,ncol(training_frame)]==1)
	training_set <- rbind(training_1, training_0[sample(1:nrow(training_0), nrow(training_1)),])
	return(training_set)
}
training <- construct_train(shock_patients, normal_patients)
model <- randomForest(x=training[,-ncol(training)], y=as.factor(training[,ncol(training)]), ntree=100, importance=T)
imp <- importance(x=model, type=1)
pdf(paste(output_path, time, ".pdf", sep=""))
varImpPlot(model)
dev.off()


