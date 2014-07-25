#!/usr/bin/python
import os

#30 directories plus the zero one
np = 31
#truncation
top = 11

d = []

#make list of the working directories
for j in range(0,8):
	for i in range(0,10):
		dir = 'e' + str(j) + '-map0' + str(i)
		d.append(dir)

	for i in range(10,np):
		dir = 'e' + str(j) + '-map' + str(i)
		d.append(dir)

currentpath = os.getcwd()
for j in range(0,8):
	for i in range(0,top):
		zz = i + j*np 
		workingpath = currentpath + '/' + d[zz]
		os.chdir(workingpath)
		os.system('qsub submit.pbs')
