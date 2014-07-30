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
#nmds0 = '477'
#dt0   = '1.d-4'
#step0 = '64'

#number of basis functions per center
bg = 20
bb = 25
bd = 25

#number of runs per signal, n+1 because 0 is included
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

#Generate Folders if not present
for i in range(0,np*8):
   if not(os.path.exists(d[i])):
      os.mkdir(d[i])


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

#write each map.in file
		l='''Np	DELTA	KONDO	NOSC	OME_MAX
3	%s	0.09	20	50
NMCS	NMDS	seed	DT	LUMDA_D
%s	%s	90	%s	10
Eg	Eb	Ed	mu	E0	E1	beta	vib_omega
0	240	240	1	%s	%s	0.24	37.7
TAU	OMEGA	time3	step1	step2
0.045	260	0.15	%s	0.04
BATH(0:B EQ 1:T EQ)	INIT	NFILE
0	3	%s
basisfg	b	d
%s	%s	%s
pi	pj	pk
%s	%s	%s''' % (str(deltahere),nmcs,nmds,dt,e0,e1,step,str(i),str(bg),
					str(bb),str(bd),str(p1),str(p2),str(p3))

		mapfile.write(l)

		nmdsstep = str(int(nmds) + i*int(step))
		delay = str(0.15 + i*0.04)

#write each intensity.in file
		m='''E1	NMDS	tau	omega	delay
%s	%s	0.045	260	%s''' % (str(e1),nmdsstep,delay)
		
		intfile.write(m)
		
		intfile.close()
		mapfile.close()

#copy executable
for i in range(0,8*np):
	shutil.copy2('a.out',d[i])
