#!/usr/bin/env python3
import os, sys, glob, shutil, multiprocessing, subprocess

#-------------------------------------------------------------------
offset = 0      # Number of MD steps for initial thermalization
increment = 100
first = offset
last = first + 1000
#-------------------------------------------------------------------

dir_list = sorted(glob.glob('0*'))

procs = []
logfiles = []

workdir = os.getcwd()	

# copy data and analyze
for datadir in dir_list:

	mdstep=int(datadir)

	first_file=first+mdstep
	last_file=last+mdstep

	print(f'\n\nworkdir = {workdir}\ndatadir = {datadir}')
	print(f'first_file,last_file,increment: {first_file} {last_file} {increment}')
	for n in range(first_file,last_file,increment):
		file_list=glob.glob(f'../../DAT/{n:09d}.xyz')
		for f in file_list:
			print(f'copying {f} to {datadir}')
			shutil.copy(f,datadir)

	print(workdir, datadir)
	os.chdir(datadir)

	p = subprocess.Popen(['../a.out','.'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	procs.append(p)

	os.chdir(workdir)

for p in procs:
    p.wait()

# make plots
procs = []

for datadir in dir_list:
	os.chdir(datadir)
	print(f'making plots in {datadir}')
	p = subprocess.Popen(['python','../plot.py'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
	procs.append(p)
	os.chdir(workdir)

for p in procs:
    p.wait()
