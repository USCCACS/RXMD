#!/usr/bin/env python3 

import glob, os
import pandas as pd

dirs = sorted(glob.glob("000*"))

keys = ['1-distance','H-H(gr)','O-H(gr)','H-O(gr)','O-O(gr)','H-H(nr)','H-O(nr)','O-H(nr)','O-O(nr)']

def sortdata(data):

	d = data[keys[0]].to_frame().join(data[keys[1]].to_frame())
	for k in keys[2:]:
		#print('key', k)
		d = d.join(data[k].to_frame())

	print(d.info())
	return d


# initialize DataFrame by the first dir. The rest of data will be appended to grdata.
filename = os.path.join(dirs[0],'gr.dat')
grdata = grdata0 = pd.read_csv(filename)
grdata.rename(columns=lambda x: x.strip(), inplace=True)
grdata = sortdata(grdata)

count=1
for d in dirs[1:]:
    filename = os.path.join(d,'gr.dat')

    print('\n',d,filename)

    # check if missing gr file
    if os.path.isfile(filename): 
        with open(filename, 'r') as fin:
            data = pd.read_csv(filename)
            data.rename(columns=lambda x: x.strip(), inplace=True)
            data = sortdata(data)

            grdata += data

            count += 1

grdata /= count
grdata.to_csv('./gr.dat', index=False)
print('# of frames', count, dirs)
