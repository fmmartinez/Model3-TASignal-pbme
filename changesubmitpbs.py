#!/usr/bin/python
import os

pname = 'TA152001'

#30 directories plus the zero one
np = 31

d = []
l = []

#make list of the working directories
for j in range(0,8):
	for i in range(0,10):
		dir = 'e' + str(j) + '-map0' + str(i)
		d.append(dir)

	for i in range(10,np):
		dir = 'e' + str(j) + '-map' + str(i)
		d.append(dir)

#original structure of the pbs file
l.append('#!/bin/bash -l\r\n')
l.append('#PBS -S /bin/bash\r\n')
l.append('#PBS -N ' + pname)
l.append('#PBS -l walltime=50:00:00\r\n')
l.append('cd $PBS_O_WORKDIR\r\n')
l.append('time ./a.out < map.in\r\n')

for k in range(0,8):
	if (k==0):
		p0name = pname + '_000-'
	if (k==1):
		p0name = pname + '_010-'
	if (k==2):
		p0name = pname + '_001-'
	if (k==3):
		p0name = pname + '_011-'
	if (k==4):
		p0name = pname + '_100-'
	if (k==5):
		p0name = pname + '_110-'
	if (k==6):
		p0name = pname + '_101-'
	if (k==7):
		p0name = pname + '_111-'
	
	for i in range(0,np):
		ii = i + k*np
		fs = open('./' + d[ii] + '/submit.pbs','w')
		l[2] = '#PBS -N ' + p0name + str(i) + '\r\n'
		for j in range(0,6):
			fs.write(l[j])
		fs.close()

   #os.system('qsub ./' + d[i] + '/submit.pbs')
