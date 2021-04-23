#!/usr/bin/env python3
import sys, os, glob, multiprocessing
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt

labels = ['O-O(gr)','O-Na(gr)','O-H(gr)','Na-O(gr)','Na-Na(gr)','Na-H(gr)','H-O(gr)','H-Na(gr)','H-H(gr)']


def gauss(x,p):
    return 1.0/(p[1]*np.sqrt(2*np.pi))*np.exp(-(x-p[0])**2/(2*p[1]**2))

def fit_gaussian(X,Y,p0):

    errfunc = lambda p, x, y: gauss(x,p) - y

    p1, success = opt.leastsq(errfunc, p0[:], args=(X,Y))
    fit_mu, fit_stdev = p1
    FWHM = 2*np.sqrt(2*np.log(2))*fit_stdev
    
    return fit_mu, fit_stdev, FWHM

def make_plot(qmdfile, basepath, nnmdfile):

    #mdstep = os.path.basename(datadir)
    #nnmdfile = os.path.join(datadir,'gr.dat')

    qmd = pd.read_csv(qmdfile)
    qmd.rename(columns=lambda x: x.strip(), inplace=True)
    nnmd = pd.read_csv(nnmdfile)
    nnmd.rename(columns=lambda x: x.strip(), inplace=True)

    for label in labels:
        qmd_x = nnmd_x = '1-distance'
        qmd_y = nnmd_y = label
        qlabel=nlabel=label
    
        plt.figure(figsize=(10,8))
        plt.rcParams.update({'font.size': 14})
        plt.ylim(0,10)
        plt.xlim(0,6)
        plt.title(f'basepath: {basepath}')

        plt.plot(qmd[qmd_x],qmd[qmd_y],label=f'QMD : {qlabel}',marker='.')
        plt.plot(nnmd[nnmd_x], nnmd[nnmd_y], label=f'NNMD : {nlabel}',marker='.')
        plt.legend()
        #plt.show()
        plt.savefig(os.path.join(basepath,'gr-'+label+'gr.qmd.png'),bbox_inches='tight')
        plt.close('all')


qmdfile = sys.argv[1]
dataroot = sys.argv[2]

if os.path.isfile(dataroot):
    basepath = '.'
    nnmdfile = dataroot
    make_plot(qmdfile, basepath, nnmdfile)
else:
# make plots
    dir_list = glob.glob(os.path.join(dataroot,'0*'))
    procs = []
    for datadir in dir_list:

        print(f'making plots for {datadir}')
        basepath = os.path.basename(datadir)
        nnmdfile = os.path.join(datadir,'gr.dat')

        p = multiprocessing.Process(target=make_plot, args=(qmdfile, basepath, nnmdfile))
        p.start()

    for p in procs:
        p.join()

