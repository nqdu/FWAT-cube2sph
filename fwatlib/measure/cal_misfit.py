import numpy as np 
from glob import glob 
import sys 
import yaml

def compute_misfits(evts:list,files:list):
    # loop every  events
    sumf = 0.
    sumn = 0
    for ievt in range(len(evts)):
        for i in range(len(files)):
            # load misfit file 
            d = np.loadtxt(files[i],usecols=[0,28],dtype=str,ndmin=2)
            # print(d.shape,files[i])

            # select matched events
            idx = d[:,0] == evts[ievt]
            sumn += d[idx,0].shape[0]
            sumf += np.sum(np.float64(d[idx,1]))
    
    return sumf,sumn

def main():
    # check input args
    if len(sys.argv) !=4 and len(sys.argv) !=3: 
        print("please run like: python ./cal_misfit.py M00 simu_type step_index")
        exit(1)

    # get 
    mod = sys.argv[1]
    simu_type = sys.argv[2]
    do_ls = False
    if len(sys.argv) == 4:
        do_ls = True
        step_index = int(sys.argv[3])
        assert(step_index in [0,1])
    
    # get period band used
    all_bands = []
    with open(f"fwat_params/FWAT.PAR.yaml","r") as f:
        pdict = yaml.safe_load(f)['measure'][f'{simu_type}']
    if simu_type == 'rf':
        # find F0
        bands = pdict['GAUSS_F0']
        for i in range(len(bands)):
            all_bands.append("%03g" % bands[i])
    else:
        bands_short = [x[0] for x in pdict['FILTER_BANDS']]
        bands_high = [x[1] for x in pdict['FILTER_BANDS']]

        for i in range(len(bands_high)):
            band = "T%03g_T%03g" %(float(bands_short[i]),float(bands_high[i]))
            all_bands.append(band)

    # now loop to compute misfits 
    sumf = 0.
    sumn = 0
    for i in range(len(all_bands)):
        band = all_bands[i]
        if do_ls:
            evts = np.loadtxt(f"src_rec/sources.dat.ls.{simu_type}",usecols=0,dtype=str)
            files = glob(f"./misfits/{mod}*.ls_{band}_{simu_type}_window_chi")
            if step_index == 0:
                tempfiles = glob(f"./misfits/{mod}.*_{band}_{simu_type}_window_chi")
                files = []
                for f in tempfiles:
                    if '.ls' not in f:
                        files.append(f)
            
            # compute misfit
            chi,n = compute_misfits(evts,files)
            sumf += chi 
            sumn += n 
        else:
            evts = np.loadtxt(f"src_rec/sources.dat.{simu_type}",usecols=0,dtype=str)
            files = []
            tempfiles = glob(f"./misfits/{mod}.*_{band}_{simu_type}_window_chi")
            for f in tempfiles:
                if '.ls' not in f:
                    files.append(f)
            chi,n = compute_misfits(evts,files)
            sumf += chi 
            sumn += n

    print(sumf,sumn) 

main()