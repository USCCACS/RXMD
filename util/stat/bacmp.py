#!/usr/bin/env python3
import sys, os, glob, multiprocessing
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt

plt.figure(figsize=(10,8))
plt.rcParams.update({'font.size': 18})

y_labels = ['H-O-H','Na-O-H']
x_label = '1-angle'

def gauss(x,p):
    return 1.0/(p[1]*np.sqrt(2*np.pi))*np.exp(-(x-p[0])**2/(2*p[1]**2))

def fit_gaussian(X,Y,p0):

    errfunc = lambda p, x, y: gauss(x,p) - y

    p1, success = opt.leastsq(errfunc, p0[:], args=(X,Y))
    fit_mu, fit_stdev = p1
    FWHM = 2*np.sqrt(2*np.log(2))*fit_stdev
    
    return fit_mu, fit_stdev, FWHM

def make_plot(qmdfile, nnmdfile, datadir):

    mdstep = os.path.basename(datadir)

    qmd = pd.read_csv(qmdfile)
    qmd.rename(columns=lambda x: x.strip(), inplace=True)
    qmd = qmd.loc[:,~qmd.columns.duplicated()]
    nnmd = pd.read_csv(nnmdfile)
    nnmd.rename(columns=lambda x: x.strip(), inplace=True)
    nnmd = nnmd.loc[:,~nnmd.columns.duplicated()]
    
    # normalize bond angle
    for y_label in y_labels:

        if not y_label in qmd: 
            print(f'WARNING: y_label ({y_label}) does not exist in {qmdfile}. no result saved.')
            continue

        qmd[y_label] /= qmd[y_label].sum()
        qmd_ = qmd[qmd[y_label]>0.0]

        if len(qmd_[x_label]) == 0:
            print(f'WARNING: y_label ({y_label}) is empty after normalization. no result saved.')
            continue

        minval = qmd_[x_label].iloc[0]
        maxval = qmd_[x_label].iloc[-1]
        
        nnmd[y_label] /= nnmd[y_label].sum()
        nnmd_ = nnmd[nnmd[y_label]>0.0]
        minval = nnmd_[x_label].iloc[0]
        maxval = nnmd_[x_label].iloc[-1]
    
        fig,ax = plt.subplots()
        ax2=ax.twinx()
        ax.set_xlim(0,180)
        ax.set_xticks([0,30,60,90,120,150,180])
        ax.set_title(f'{y_label} : {mdstep}')
        ax.plot(qmd[x_label],qmd[y_label],label=f'QMD',marker='.',color='red')
        ax2.plot(nnmd[x_label], nnmd[y_label], label=f'NNQMD',marker='.',color='blue')
        ax.legend()
        ax2.legend()
        #plt.show()
        print(f"saving {os.path.join(datadir,y_label+'.qmd.png')}")
        fig.savefig(os.path.join(datadir,'ba-'+y_label+'.qmd.png'),bbox_inches='tight',dpi=300)


qmdfile = sys.argv[1]
dataroot = sys.argv[2]

# make plots
if os.path.isfile(dataroot):
    nnmdfile = dataroot
    make_plot(qmdfile, nnmdfile, os.getcwd())
else:
   dir_list = glob.glob(dataroot+"/0*")
   procs = []
   for datadir in dir_list:
       print(f'making plots for {datadir}')
   
       nnmdfile = os.path.join(datadir,os.path.basename(qmdfile))
       p = multiprocessing.Process(target=make_plot, args=(qmdfile, nnmdfile, datadir))
       p.start()
       procs.append(p)
   
   for p in procs:
       p.join()
