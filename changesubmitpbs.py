#!/usr/bin/python
import os

pname = 'TA152001'
wallt = '50:00:00'

#30 directories plus the zero one
np = 31

d = []

#make list of the working directories
for j in range(0,8):
	for i in range(0,10):
		dir = 'e' + str(j) + '-map0' + str(i)
		d.append(dir)

	for i in range(10,np):
		dir = 'e' + str(j) + '-map' + str(i)
		d.append(dir)

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
		l = '''#!/bin/bash -l
#PBS -S /bin/bash
#PBS -N %s%s
#PBS -l walltime=%s

cd $PBS_O_WORKDIR
time ./a.out < map.in''' %(p0name,str(i),wallt)

		fs.write(l)
		fs.close()

   #os.system('qsub ./' + d[i] + '/submit.pbs')
