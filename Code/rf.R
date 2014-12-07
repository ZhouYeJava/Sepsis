#!/usr/bin/Rscript

###################
###Random Forest###
###################

###Zhou Ye###
###06/24/2014###

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

###Train/Test###
#construct training set
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
#test evaluation
test_evaluation <- function(shock_test, normal_test, model) {
	test <- foreach (patient = c(shock_test, normal_test), .combine="rbind") %dopar% {
		data <- patient$data
		prediction <- as.data.frame(predict(object=model, newdata=data[,-1], type="prob"))
		value <- unlist(prediction["1"])
		patient_time <- data[,1]
		label <- NULL
		if (patient$flag==1) {
			label <- rep(0, nrow(data))
		}
		else {
			shock_time <- patient_time[length(patient_time)]
			label <- vector(length=nrow(data))
			label[patient_time<shock_time-time] <- 0
			label[patient_time>=shock_time-time] <- 1
		}
		output <- NULL
		if (nrow(data)>1) {
			output <- cbind(value, label)
		}
		else {
			output <- t(c(value, label))
		}
		output
	}
	return(roc.area(obs=test[,2], pred=test[,1])$A)
}
#model selection
model_selection <- function(shock_patients, normal_patients, candicates) {
	result <- vector(length=length(candicates))
	num_shock <- length(shock_patients)
	num_normal <- length(normal_patients)
	shock_train <- shock_patients[1:ceiling(0.7*num_shock)]
	shock_test <- shock_patients[(ceiling(0.7*num_shock)+1):num_shock]
	normal_train <- normal_patients[1:ceiling(0.7*num_normal)]
	normal_test <- normal_patients[(ceiling(0.7*num_normal)+1):num_normal]
	training <- construct_train(shock_train, normal_train)
	for (j in 1:length(candicates)) {
		model <- randomForest(x=training[,-ncol(training)], y=as.factor(training[,ncol(training)]), ntree=candicates[j])
		result[j] <- test_evaluation(shock_test, normal_test, model)
	}
	return(candicates[which.max(result)])
}
#construct test output
construct_test <- function(patient, model, output_path, fold_num) {
	data <- patient$data
	prediction <- as.data.frame(predict(object=model, newdata=data[,-1], type="prob"))
	value <- unlist(prediction["1"])
	patient_time <- data[,1]
	label <- NULL
	if (patient$flag==1) {
		label <- rep(0, nrow(data))
	}
	else {
		shock_time <- patient_time[length(patient_time)]
		label <- vector(length=nrow(data))
		label[patient_time<shock_time-time] <- 0
		label[patient_time>=shock_time-time] <- 1
	}
	output <- NULL
	if (nrow(data)>1) {
		output <- cbind(patient_time, value, label)
	}
	else {
		output <- t(c(patient_time, value, label))
	}
	write.table(x=output, file=paste(output_path, "fold_", fold_num, "_", patient$identifier, ".csv", sep=""), sep=",", row.names=F, col.names=F)
}
#evaluation
num_fold <- 5
num_tree <- c(100, 200, 300, 400, 500)
num_shock <- length(shock_patients)
num_normal <- length(normal_patients)
system(paste("rm -f ", output_path, "*", sep=""))
for (i in 1:num_fold) {
	print(paste("fold", i, sep="_"))
	num <- ceiling(num_shock/num_fold)
	start <- num*(i-1)+1
	end <- ifelse(num*i<num_shock, num*i, num_shock)
	shock_test <- shock_patients[start:end]
	shock_train <- shock_patients[-(start:end)]
	num <- ceiling(num_normal/num_fold)
	start <- num*(i-1)+1
	end <- ifelse(num*i<num_normal, num*i, num_normal)
	normal_test <- normal_patients[start:end]
	normal_train <- normal_patients[-(start:end)]
	parameter <- model_selection(shock_patients, normal_patients, num_tree)
	print(paste("parameter = ", parameter, sep=""))
	training <- construct_train(shock_train, normal_train)
	model <- randomForest(x=training[,-ncol(training)], y=as.factor(training[,ncol(training)]), ntree=parameter)
	test <- foreach (patient = c(shock_test, normal_test)) %dopar% {
		construct_test(patient, model, output_path, i)
	}
}
print("Finish Evaluation")

