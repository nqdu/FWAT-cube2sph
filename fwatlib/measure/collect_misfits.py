import numpy as np
from glob import glob 
import sys 
import yaml

def main():
    # check input args
    if len(sys.argv) !=4 and len(sys.argv) !=3: 
        print("please run like: python ./collect_misfits.py M00 simu_type [do_ls=0]")
        exit(1)

    mod = sys.argv[1]
    simu_type = sys.argv[2]
    do_ls = False
    if len(sys.argv) == 4:
        do_ls = (int(sys.argv[3]) == 1)
    print(do_ls)

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
    
    for i in range(len(all_bands)):
        # get all files
        band = all_bands[i]
        if do_ls:
            tempfiles = glob(f"./misfits/{mod}*.ls_{band}_{simu_type}_window_chi")
            newfile = f"{mod}.ls_{band}_{simu_type}_window_chi"
            files = []
            for f in tempfiles:
                if newfile not in f:
                    files.append(f) 
        else:
            tempfiles = glob(f"./misfits/{mod}.*_{band}_{simu_type}_window_chi")
            files = []
            for f in tempfiles:
                if '.ls' not in f:
                    files.append(f) 
            newfile = f"{mod}_{band}_{simu_type}_window_chi"
        
        # write new files
        if len(files) >= 1:
            fio = open(newfile,"w")
            for i in range(len(files)):
                d = np.loadtxt(files[i],usecols=[0,1,2,3,28],dtype=str,ndmin=2)
                for j in range(d.shape[0]):
                    fio.write(" ".join(d[j,:].tolist()) + "\n")
            fio.close()

main()    
