#!/usr/bin/env python3

import pandas as pd
import glob, os, sys
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import MultipleLocator, AutoMinorLocator

plt.rcParams.update({'font.size': 24})
#plt.rcParams["font.family"] = "Times New Roman"
plt.rc('grid', linestyle=':', color='red', linewidth=1)

gr_xlabel='1-distance'
gr_ylabels=['H-H','H-N','N-H','N-N']

ba_xlabel='1-angle'
ba_ylabels=['N-H-N','N-N-N']

sq_xlabel='wave_number'
sq_ylabels=['H-H','H-N','N-H','N-N']


def ba_plot(dirs, num_columns=2):

    num_plots = len(ba_ylabels)
    num_rows = int(num_plots/num_columns)

    for dir in dirs:
        for infile in glob.glob(dir+'/ba.dat'):

            dat = pd.read_csv(infile)
            dat.columns = dat.columns.str.strip()

            fig = plt.figure(figsize=(num_rows*10, num_columns*10))

            for c in range(num_plots):
                ax = fig.add_subplot(num_rows,num_columns,c+1)

                ba_y = ba_ylabels[c]
                ax.plot(dat[ba_xlabel], dat[ba_y], label=ba_y, linewidth=2)

                #subplots.set_ylim([0,1e-6])
                ax.set_xlim([0,180])
                ax.set_xticks((0, 30, 60, 90, 120, 150, 180))
                ax.yaxis.set_major_locator(ticker.NullLocator())
                ax.legend(loc='lower left')

            outfile = 'ba.png'
            plt.savefig(outfile,bbox_inches='tight')

            print('=== BA ===================================')
            print(f'infile {infile}       outfile {outfile}')
            print(f'dat.cloumns {dat.columns}')
            print('==========================================')
    return

def sq_plot(dirs, num_columns=2):

    num_plots = len(sq_ylabels)
    num_rows = int(num_plots/num_columns)

    for dir in dirs:
        for infile in glob.glob(dir+'/sq.dat'):
            
            dat = pd.read_csv(infile)
            dat.columns = dat.columns.str.strip()

            # Sn(q) neutron
            fig = plt.figure(figsize=(10, 10))
            ax = fig.add_subplot()
            ax.plot(dat[sq_xlabel], dat['Snq'], label='Sn(q)', linewidth=3)
            plt.savefig(os.path.join(dir,'snq.png'),bbox_inches='tight')

            # S(q) partials
            fig = plt.figure(figsize=(10*num_rows, 10*num_columns))

            for c in range(num_plots):
                ax = fig.add_subplot(num_rows,num_columns,c+1)
                sq_y = sq_ylabels[c]
                ax.plot(dat[sq_xlabel], dat[sq_y], label=sq_y, linewidth=3)
                ax.grid(which='major', color='red', linewidth=2.0)
                ax.grid(which='minor', color='green', linestyle=':', linewidth=1.0)
                ax.legend()

            outfile = os.path.join(dir,'sq.png')
            plt.savefig(outfile,bbox_inches='tight')

            print('=== S(q) =================================')
            print(f'infile {infile}       outfile {outfile}')
            print(f'dat.cloumns {dat.columns}')
            print('==========================================')
            
    return 
                
def gr_plot(dirs, num_columns=2):

    num_plots = len(gr_ylabels)
    num_rows = int(num_plots/num_columns)

    for dir in dirs:
        for infile in glob.glob(dir+'/gr.dat'):
            
            fig = plt.figure(figsize=(10*num_rows, 10*num_columns))
            
            print('===== %s ===='%(infile))
            
            dat = pd.read_csv(infile)
            dat.columns = dat.columns.str.strip()

            for c in range(num_plots):
                ax= fig.add_subplot(num_rows,num_columns,c+1)
                ax.set_ylim([0,6])
                ax.set_xlim([0,5])
                #subplots.xaxis.set_minor_locator(AutoMinorLocator(10))
                ax.minorticks_on()

                gr_y = gr_ylabels[c]+'(gr)'
                nr_y = gr_ylabels[c]+'(nr)'

                ax.plot(dat[gr_xlabel], dat[gr_y], label=gr_y, linewidth=3)
                ax.plot(dat[gr_xlabel], dat[nr_y], label=nr_y, linewidth=3)
                ax.grid(which='major', color='red', linewidth=1.0)
                ax.grid(which='minor', color='green', linestyle=':', linewidth=0.8)
                ax.legend()

            outfile = os.path.join(dir,'gr.png')
            plt.savefig(outfile,bbox_inches='tight')
    
            print('=== g(r) =================================')
            print(f'gr: infile {infile}       outfile {outfile}')
            print(f'dat.cloumns {dat.columns}')
            print('==========================================')
            
    return 


if __name__ == "__main__":

    dirs=['.']

    gr_plot(dirs)
    sq_plot(dirs)
    ba_plot(dirs)
