#!/usr/bin/env python3

import pandas as pd
import glob, os, sys
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter

plt.rcParams.update({'font.size': 24})
plt.rcParams["font.family"] = "Times New Roman"

def ba_plot(dirs, num_ba=2):
    
    for dir in dirs:
        for infile in glob.glob(dir+'/ba-*.dat'):

            dat = pd.read_csv(infile)
            dat = dat.iloc[:, :-1] # drop last column
            columns = dat.columns
            names = []
            for c in columns:
                names.append(c.lstrip())
            print(names)
            dat.columns = names
            dat = dat.reindex(sorted(dat.columns), axis=1)
            columns = dat.columns
            print(columns)
            
            num_rows = (len(columns)-1)/num_ba
            print('===== %s %d ===='%(infile, num_rows))
             
            fig = plt.figure(figsize=(num_rows*10, num_ba*8))

            for c in range(len(columns)-1):
                subplots = fig.add_subplot(num_ba, num_rows, c+1)
                plt.setp(subplots.get_yticklabels(),visible=False)
                
                #subplots.set_ylim([0,1e-6])
                subplots.set_xlim([0,180])
                subplots.set_xticks((0, 30, 60, 90, 120, 150, 180))
                subplots.plot(dat[columns[0]], dat[columns[c+1]], 
                              label=columns[c+1].lstrip(), linewidth=2)
                subplots.legend(loc='lower left')
                subplots.grid()
            outfile = infile+'.png'
            plt.savefig(outfile,bbox_inches='tight')
            print('------ %s -----'%(outfile))
    return

def sq_plot(dirs, num_sq=2):

    num_cols = num_sq*num_sq + 1
    num_pairs = num_sq*num_sq
    print('num_sq,num_cols,num_pairs: ', num_sq,num_cols,num_pairs)

    for dir in dirs:
        for infile in glob.glob(dir+'/sq.dat'):
            
            print('===== %s ===='%(infile))
            
            dat = pd.read_csv(infile)
            
            columns = dat.columns
            names = []
            for c in columns:
                names.append(c.lstrip())
            print('name: ', names)

            dat.columns = names
            dat = dat.iloc[:,0:num_cols+1]
            dat = dat.reindex(sorted(dat.columns), axis=1)
            print('dat.columns: ', dat.columns)
            
            columns = dat.columns
            print(columns)

            # Sn(q)
            fig = plt.figure(figsize=(10*num_sq, 8*num_sq))
            subplots = fig.add_subplot(1,1,1)
            subplots.plot(dat[columns[-1]], dat[columns[-2]], \
                label=columns[-2].lstrip(), linewidth=3)
            outfile = os.path.join(dir,'snq.png')
            subplots.legend()
            subplots.grid()
            plt.savefig(outfile,bbox_inches='tight')

            # S(q) partials
            fig = plt.figure(figsize=(10*num_sq, 8*num_sq))
            for c in range(num_pairs):
                subplots = fig.add_subplot(num_sq,num_sq,c+1)
                #subplots.set_ylim([0,6])
                subplots.plot(dat[columns[-1]], dat[columns[c]], 
                              label=columns[c].lstrip(), linewidth=3)
                subplots.legend()
                subplots.grid()
            outfile = os.path.join(dir,'sq.png')
            plt.savefig(outfile,bbox_inches='tight')
    
            print('------ %s -----'%(outfile))
            
    return 
                
def gr_plot(dirs, num_gr=2):

    num_cols = num_gr*num_gr*2 + 2
    num_pairs = num_gr*num_gr
    print('num_gr,num_cols,num_pairs: ', num_gr,num_cols,num_pairs)

    for dir in dirs:
        for infile in glob.glob(dir+'/gr.dat'):
            
            fig = plt.figure(figsize=(10*num_gr, 8*num_gr))
            
            print('===== %s ===='%(infile))
            
            dat = pd.read_csv(infile)
            
            columns = dat.columns
            names = []
            for c in columns:
                names.append(c.lstrip())
            print(names)
            dat.columns = names
            dat = dat.iloc[:,0:num_cols]
            dat = dat.reindex(sorted(dat.columns), axis=1)
            print(dat.columns)
            
            columns = dat.columns
            print(columns)

            for c in range(num_pairs):
                subplots = fig.add_subplot(num_gr,num_gr,c+1)
                subplots.set_ylim([0,6])
                #subplots.set_xlim([0,10])
                subplots.plot(dat[columns[0]], dat[columns[2*c+1]], 
                              label=columns[2*c+1].lstrip(), linewidth=3)
                subplots.plot(dat[columns[0]], dat[columns[2*c+2]],
                              label=columns[2*c+2].lstrip(), linewidth=3)
                subplots.legend()
                subplots.grid()
            outfile = os.path.join(dir,'gr.png')
            plt.savefig(outfile,bbox_inches='tight')
    
            print('------ %s -----'%(outfile))
            
    return 


dirs=['.']
gr_plot(dirs)
ba_plot(dirs,2)
sq_plot(dirs)
