#!/usr/bin/Rscript

#####################
###Single Time Cox###
#####################

###Zhou Ye###
###05/14/2014###

suppressMessages(library(foreach))
suppressMessages(library(doMC))
suppressMessages(library(survival))
registerDoMC(50)
args <- commandArgs(TRUE)
shock_path <- args[1]
normal_path <- args[2]
t0 <- as.integer(args[3])
t1 <- as.integer(args[4])

###Read Data###
shock_files <- list.files(shock_path, pattern="*.csv")
normal_files <- list.files(normal_path, pattern="*.csv")
shock_data <- foreach(i = 1:length(shock_files), .combine = "rbind") %dopar% {
	data <- read.csv(paste(shock_path, shock_files[i], sep=""))
	if (t0<nrow(data)) {
		line <- vector(length=ncol(data)+1)
		line[1] <- nrow(data)-t0
		line[2] <- 1
		line[3:length(line)] <- as.numeric(data[t0, 2:ncol(data)])
		names(line) <- c("time", "status", paste("feature", 1:(ncol(data)-1), sep=""))
		line		
	}
	else {
		NULL
	}
}
normal_data <- foreach(i = 1:length(normal_files), .combine = "rbind") %dopar% {
	data <- read.csv(paste(normal_path, normal_files[i], sep=""))
	if (t0<nrow(data)) {
		line <- vector(length=ncol(data)+1)
		line[1] <- nrow(data)-t0
		line[2] <- 0
		line[3:length(line)] <- as.numeric(data[t0, 2:ncol(data)])
		names(line) <- c("time", "status", paste("feature", 1:(ncol(data)-1), sep=""))
		line		
	}
	else {
		NULL
	}
}
all_data <- as.data.frame(rbind(shock_data, normal_data))

###10 Fold Cross Validation###
num_fold = 5
all_data <- all_data[sample(1:nrow(all_data), nrow(all_data)),]
num <- ceiling(nrow(all_data)/num_fold)
brier <- vector(length=num_fold)
accuracy <- vector(length=num_fold)
for (i in 1:num_fold) {
	start <- num*(i-1)+1
	end <- ifelse(num*i>nrow(all_data), nrow(all_data), num*i)
	train <- all_data[-(start:end),]
	test <- all_data[start:end,]
	model <- coxph(formula=Surv(time, status)~., data=train)
	prob <- survfit(formula=model, newdata=test)$surv[t1,]
	pred <- as.numeric(prob<0.5)
	label <- as.numeric(test$time<(t0+t1))
	index <- which(test$status==1)
	brier[i] <- mean((label[index]-prob[index])^2)
	accuracy[i] <- sum(label[index]==pred[index])/length(label[index])
}
print(t1)
print(paste("Brier Score:", mean(brier), sd(brier), sep=" "))
print(paste("Accuracy:", mean(accuracy), sd(accuracy), sep=" "))
