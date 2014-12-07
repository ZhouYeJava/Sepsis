#!/usr/bin/python

############################
###Get Data From Database###
############################

###Zhou Ye###
###04/29/2014###

import os

begin = "psql -d MIMIC2 -t -A -F',' -c " #log in the database
setPath = "SET search_path TO workspace,public" #set path in the database

def getData(icuID, feature):
	query = "SELECT hosp_time, sirs_intp, severe_sepsis, septic_shock"
	for f in feature:
		query += ", "+f
	query += " FROM trend_vars WHERE trend_vars.icustay_id="+str(icuID)+" ORDER BY hosp_time"
	command = begin+"'"+setPath+";"+query+"'"+" > "+"/udata/PatientNew/ICU"+str(icuID)+".csv"
	os.system(command)

if __name__=="__main__":
	reader1 = open("id.txt", "r")
	reader2 = open("feature_new.txt", "r")
	allPatient = []
	allFeature = []

	for line in reader1.readlines():
		line = line.strip()
		allPatient.append(int(line))

	for line in reader2.readlines():
		line = line.strip()
		allFeature.append(line)

	reader1.close()
	reader2.close()

	for p in allPatient:
		print(p)
		getData(p, allFeature)
