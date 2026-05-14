import sys 
import numpy as np 
from glob import glob
import os 
import yaml 
import argparse

import matplotlib.pyplot as plt 
import matplotlib as mpl
mpl.rcParams['lines.linewidth'] = 0.5
mpl.rcParams['font.size'] = 10
mpl.rcParams['xtick.labelsize']=10
mpl.rcParams['ytick.labelsize']=10
mpl.rcParams['axes.labelsize']=10
mpl.rcParams['legend.fontsize'] = 8
mpl.rcParams['legend.fontsize'] = 8
mpl.rcParams['savefig.bbox'] = 'tight'

from fwat.const import PARAM_FILE

def main():
    parser = argparse.ArgumentParser(description='Plot histogram of noise traces for all events and stations.')
    parser.add_argument('--model', type=str, required=True,help='Model name (e.g., M03)')
    parser.add_argument("--model0", type=str, default="M00", help='Initial model name for comparison (default: M00)')
    parser.add_argument("--path", type=str, default="../../", help='Path to the working directory (default: ../../)')
    args = parser.parse_args()
    
    # set directory
    path = args.path
    seisdir= "noise_hist/"
    misfits = f"{path}/misfits"
    paramfile = f"{path}/{PARAM_FILE}"

    #### stop here

    # read model name
    M1 = args.model
    M0 = args.model0

    os.makedirs(seisdir,exist_ok=True)

    # load paramfile
    with open(paramfile,"r") as f:
        pdict = yaml.safe_load(f)['measure']['noise']
    nbands = len(pdict['FILTER_BANDS'])

    for ib in range(nbands):
        data = np.array([])
        data0 = np.array([])
        Tmin,Tmax = pdict['FILTER_BANDS'][ib]
        band="T%03g_T%03g" %(Tmin,Tmax)

        # load misfits data
        chi = 0.
        chi0 = 0.
        filenames = glob(f"{misfits}/{M1}/*_{band}_noise_window_chi")
        for i in range(len(filenames)):
            temp = np.loadtxt(filenames[i],usecols=[5,-1],ndmin=2)
            idx = np.where(abs(temp[:,0]) > 1.0e-5)[0]
            data =  np.append(data,temp[idx,0])
            chi += np.sum(temp[:,-1])
        filenames = glob(f"{misfits}/{M0}/*_{band}_noise_window_chi")
        for i in range(len(filenames)):
            temp = np.loadtxt(filenames[i],usecols=[5,-1],ndmin=2)
            idx = np.where(abs(temp[:,0]) > 1.0e-5)[0]
            data0 =  np.append(data0,temp[idx,0])
            chi0 += np.sum(temp[:,-1])
        print(f"{band}: {len(data)} traces, chi = {chi:.2f}, chi0 = {chi0:.2f}")

        # create figures
        fig,ax = plt.subplots(1,1,figsize=(12,5))
        ax.hist(data0,bins=30,range=(-5.5,5.5),label=f'{M0}',color='blue',alpha=0.5)
        ax.hist(data,bins=30,range=(-5.5,5.5),label=f'{M1}',color='gray',alpha=0.5)
        ax.legend()
        fig.savefig(f"{seisdir}/{M1}.{band}.jpg")
        fig.clear()


if __name__ == "__main__":
    main()
