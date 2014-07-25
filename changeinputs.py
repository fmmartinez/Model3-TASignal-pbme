#!/usr/bin/python
import os
import shutil

#run value
deltahere = 0.0

#new values
dt    = '5d-5'
nmds  = '954'
step  = '128'

nmcs = '10000'

#original values
#double check that dt, nmds and step are related with these!!!
#for example halving dt doubles nmds and step
nmds0 = '477'
dt0   = '1.d-4'
step0 = '64'

#number of runs per signal, n+1 because 0 is included
np = 31

d = []
l = []
m = []

#make list of the working directories
for j in range(0,8):
	for i in range(0,10):
		dir = 'e' + str(j) + '-map0' + str(i)
		d.append(dir)

	for i in range(10,np):
		dir = 'e' + str(j) + '-map' + str(i)
		d.append(dir)

#Generate Folders if not present
for i in range(0,np*8):
   if not(os.path.exists(d[i])):
      os.mkdir(d[i])

#original structure of map.in file in map00:
l.append('Np\tDELTA\tKONDO\tNOSC\tOME_MAX\t\r\n')
l.append('3\t' + str(deltahere) + '\t0.09\t20\t50\r\n')
l.append('NMCS\tNMDS\tseed\tDT\tLUMDA_D\r\n')
l.append(nmcs + '\t' + nmds  + '\t90\t' + dt  + '\t10\r\n')
l.append('Eg\tEb\tEd\tmu\tE0\tE1\tbeta\tvib_omega\r\n')
l.append('0\t240\t240\t1\t0.5\t0.0\t0.24\t0.9549d0\r\n')
l.append('TAU\tOMEGA\ttime3\tstep1\tstep2\r\n')
l.append('0.045\t260\t0.15\t' + step  + '\t0.04\r\n')
l.append('BATH(0:B EQ 1:T EQ)\tINIT\tNFILE\r\n')
l.append('0\t\t\t3\t0\r\n')
l.append('basispc\tg\tb\td\r\n')
l.append('0\t5\t10\t1\r\n')
l.append('pi\pj\pk\r\n')
l.append('0\t0\t0\r\n')

#original structure of intensity.in file in map00:
m.append('E1\tNMDS\ttau\tomega\tdelay\r\n')
m.append('0.0\t477\t0.045\t260\t0.15\r\n')

for k in range(0,8):
	e0 = 1.5
	e1 = 0.15
	
	if (k == 0):
		p1 = 0
		p2 = 0
		p3 = 0
	if (k == 1):
		p1 = 0
		p2 = 1
		p3 = 0
	if (k == 2):
		p1 = 0
		p2 = 0
		p3 = 1
	if (k == 3):
		p1 = 0
		p2 = 1
		p3 = 1
	if (k == 4):
		p1 = 1
		p2 = 0
		p3 = 0
	if (k == 5):
		p1 = 1
		p2 = 1
		p3 = 0
	if (k == 6):
		p1 = 1
		p2 = 0
		p3 = 1
	if (k == 7):
		p1 = 1
		p2 = 1
		p3 = 1
	
	for i in range(0,np):
		ii = i + k*np
		mapfile = open('./' + d[ii]  + '/map.in','w')
		intfile = open('./' + d[ii]  + '/intensity.in','w')
		
#		print 'cp map ../'+d[ii]+'/.'
		
		l[5] = '0\t240\t240\t1\t' + str(e0) + '\t' + str(e1) + '\t0.24\t37.7d0\r\n'
		l[9] = '0\t\t\t3\t' + str(i) + '\r\n'
		l[13]= str(p1) + '\t' + str(p2) + '\t' + str(p3) + '\r\n'
		
		for j in range(0,14):
			mapfile.write(l[j])

		nmdsstep = str(int(nmds) + i*int(step))
		delay = str(0.15 + i*0.04)
		m[1] = str(e1) + '\t' + nmdsstep  + '\t0.045\t260\t' + delay + '\r\n'
		for j in range(0,2):
			intfile.write(m[j])
		
		intfile.close()
		mapfile.close()

#copy executable
for i in range(0,8*np):
	shutil.copy2('a.out',d[i])
