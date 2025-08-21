import numpy as np 
import os 
import sys 
import yaml

def run(argv):
    from fwat.const import PARAM_FILE,SRC_REC,MISFIT

    # check input args
    if len(argv) !=2:
        print("Usage: fwat misfit MODEL(M00) simu_type")
        exit(1)

    mod = argv[0]
    simu_type = argv[1]

    # get period band used
    all_bands = []
    with open(f"{PARAM_FILE}","r") as f:
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

    # now loop each event to compute misfits 
    sumf = 0.
    sumn = 0
    evts = np.loadtxt(f"{SRC_REC}/sources.dat.{simu_type}",usecols=0,dtype=str)
    for ievt in range(len(evts)):
        for ib in range(len(all_bands)):
            band = all_bands[ib]
            filename = f"./{MISFIT}/{mod}/{evts[ievt]}_{band}_{simu_type}_window_chi"
            if os.path.exists(filename):
                d = np.loadtxt(filename,usecols=28,ndmin=2)
                chi = np.sum(d[:,0])
                # print(band,evts[ievt],d.shape[0])
                if np.isnan(chi):
                    print(f"error in {evts[ievt]} {filename}")
                n = d.shape[0]

                # accumulate
                sumf += chi 
                sumn += n 

    print(sumf,sumn)

def main():
    # check input args
    if len(sys.argv) !=3:
        print("please run like: python ./cal_misfit.py M00 simu_type")
        exit(1)

    # run 
    run(sys.argv[1:])


if __name__ == "__main__":
    main()
