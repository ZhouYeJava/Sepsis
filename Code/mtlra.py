#!/usr/bin/python

######################
###Single Time MTLR###
######################

###Zhou Ye###
###07/21/2014###

import os
import sys
import glob
import random
import math
import re

###Patient Class###
class Patient:
	def __init__(self, data_file, flag):
		self.identifier = int(re.findall(r'\d+', data_file)[-1])
		reader = open(data_file, "r")
		lines = reader.readlines() #first line is column names
		num_time = len(lines)-1 #total number of sample points
		self.shock_time = int(lines[num_time].strip().split(",")[0]) #time of septic shock
		self.flag = flag
		self.data = {}
		for i in range(1, num_time):
			line = lines[i].strip().split(",")
			curr_time = int(line[0])
			feature = ""
			for j in range(1, len(line)):
				if float(line[j])!=0:
					feature += str(j)+":"+line[j]+" "
			self.data[curr_time] = feature
	def __eq__(self, other):
		if not isinstance(other, self.__class__):
			return False
		else:
			if self.identifier==other.identifier:
				return True
			else:
				return False
	def single_select(self, time):
		if time not in self.data:
			return ""
		output = str(self.shock_time-time)+" "+str(self.flag)+" "+self.data[time]+"\n"
		return output

###Read Data###
def read_data(data_path, flag):
	files = glob.glob(data_path+"/*.csv")
	data = [Patient(data_file, flag) for data_file in files]
	return data

###Construct MTLR Files###
def construct_single_file(data, time):
	data_string = [d.single_select(time) for d in data]
	return "".join(data_string)

###Parse Result###
def parse_result(result_file):
	brier = None
	accuracy = None
	reader = open(result_file, "r")
	lines = reader.readlines()
	count = 0
	for line in lines:
		if line.startswith("#"):
			if count==0:
				try:
					accuracy = float(line.split(": ")[-1])
				except ValueError:
					pass
				count += 1
			elif count==1:
				try:
					brier = float(line.split(": ")[-1])
				except ValueError:
					pass
				count += 1
			else:
				break
	reader.close()
	return (brier, accuracy)

###Compute Statistics###
def mean(array):
	return sum(array)/len(array)
def sd(array):
	m = mean(array)
	diff = [(x-m)**2 for x in array]
	return math.sqrt(sum(diff)/len(diff))

if __name__=="__main__":
	###Configuration###
	data_shock = sys.argv[1] #location of shock patients data
	data_normal = sys.argv[2] #location of normal patients data
	mtlr = sys.argv[3] #location of multi-task logistic regression program
	t0 = int(sys.argv[4]) #current time points (features)
	t1 = int(sys.argv[5]) #next time points (labels)
	p0 = sys.argv[6] #penalty 0
	p1 = sys.argv[7] #penalty 1

	###Read Data###
	#0 for septic shock patients; 1 for normal patients
	shock_patients = read_data(data_shock, 0)
	normal_patients = read_data(data_normal, 1)
	data = shock_patients+normal_patients
	print "Finish Reading Data"

	###Cross Validation###
	num_fold = 5
	data = random.sample(data, len(data))
	num = int(math.ceil(len(data)/num_fold))
	os.chdir(mtlr)
	os.system("make")
	t_max = t1+1
	brier = []
	accuracy = []
	for i in range(1, num_fold+1):
		print i
		start = num*i
		end = num*(i+1) if num*(i+1)<len(data) else len(data)
		test = data[start:end]
		train = [d for d in data if d not in test]
		train_file = construct_single_file(train, t0)
		test_file = construct_single_file(test, t0)
		writer_train = open("fold_"+str(i)+".train", "w")
		writer_test = open("fold_"+str(i)+".test", "w")
		writer_train.write(train_file)
		writer_test.write(test_file)
		writer_train.close()
		writer_test.close()
		os.system("./mtlr_train -c "+p0+" -d "+p1+" -m "+str(t_max)+" -u 1 -i "+"fold_"+str(i)+".train"+" -o "+"fold_"+str(i)+".model")
		os.system("./mtlr_test -l class -t "+str(t1)+" -m "+str(t_max)+" -i "+"fold_"+str(i)+".test"+" -o "+"fold_"+str(i)+".model > fold_"+str(i)+".result")
		score = parse_result("fold_"+str(i)+".result")
		brier.append(score[0])
		accuracy.append(score[1])	
	os.system("make clean")
	print "Finish Training and Testing!"
	print t1
	print "Accuracy:"+str(mean(accuracy))+" "+str(sd(accuracy))
	print "Brier Score: "+str(mean(brier))+" "+str(sd(brier))
	os.system("rm -f fold*")

