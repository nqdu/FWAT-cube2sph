import sys 
import numpy as np 
from glob import glob
import os 
import yaml 

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

def main():
    if len(sys.argv) != 2:
        print("Usage: ./compare_tele.py model(M03)")
        exit(1)
    
    # set directory
    from fwat.const import PARAM_FILE
    path = "../.."
    seisdir= "noise_hist/"
    misfits = f"{path}/misfits"
    paramfile = f"{path}/{PARAM_FILE}"

    #### stop here

    # read model name
    M1 = sys.argv[1]

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
        filenames = glob(f"{misfits}/{M1}/*_{band}_noise_window_chi")
        for i in range(len(filenames)):
            temp = np.loadtxt(filenames[i],usecols=[14])
            data =  np.append(data,temp)
        filenames = glob(f"{misfits}/M00/*_{band}_noise_window_chi")
        for i in range(len(filenames)):
            temp = np.loadtxt(filenames[i],usecols=[14])
            data0 =  np.append(data0,temp)

        # create figures
        fig,ax = plt.subplots(1,1,figsize=(12,5))
        ax.hist(data0,bins=30,range=[-1,1],label='M00')
        ax.hist(data,bins=30,range=[-1,1],label=f'{M1}',color='gray',alpha=0.5)
        ax.legend()
        fig.savefig(f"{seisdir}/{M1}.{band}.jpg")
        fig.clear()
        print(np.max(abs(data)),np.max(abs(data0)))


if __name__ == "__main__":
    main()
