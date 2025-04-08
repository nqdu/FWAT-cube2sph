import numpy as np
from obspy.io.sac import SACTrace
import sys 
import os 

from utils import read_fwat_params
from scipy.signal import convolve

def main():
    if len(sys.argv) != 3:
        print("Usage: python syn_data.py iter evtid")
        exit(1)
    
    # get parma
    iter = int(sys.argv[1])
    evtid = sys.argv[2]

    # read station coordinates
    stationfile = 'solver/M%02d' %(iter) + f'.set{evtid}/DATA/STATIONS_FILTERED'
    statxt = np.loadtxt(stationfile,dtype=str,ndmin=2)

    # read fwat parameters
    pdict = read_fwat_params('solver/M%02d' %(iter) + f'.set{evtid}/DATA/FWAT.PAR.yaml')['measure']['tele']
    CCODE = "." + pdict['CH_CODE']
    
    # synthetic paramters
    syndir = 'solver/M%02d' %(iter) + f'.set{evtid}/OUTPUT_FILES'
    name = statxt[0,1] + "." + statxt[0,0] + CCODE + "Z.sac"
    syn_z_hd = SACTrace.read(syndir + name,headonly=True)
    npt_syn = syn_z_hd.npts
    dt_syn = syn_z_hd.delta

    # get frequency band
    Tmin_list = [x[0] for x in pdict['FILTER_BANDS']]
    Tmax_list = [x[1] for x in pdict['FILTER_BANDS']]

    # sum for stf for each period band if required
    ncomp = len(pdict['COMPS'])
    components = pdict['COMPS']
    stf = np.zeros((ncomp,npt_syn))
    nsta = statxt.shape[0]
    for ib in range(len(Tmax_list)):
        bandname = 'T%03g_T%03g' %(Tmin_list[ib],Tmax_list[ib])
        for ic in range(ncomp):
            ch = components[ic]
            tr = SACTrace(f'src_rec/stf_{ch}.sac.{bandname}_{id}').data 
            stf[ic,:] += tr.data 
    
    # synthetic by using summed stf
    outdir = syndir + '/SYN/'
    os.makedirs(outdir,exist_ok=True)
    for i in range(nsta):
        for ic in range(ncomp):
            ch = components[ic]
            name = statxt[i,1] + "." + statxt[i,0] + CCODE + ch + ".sac"
            syn_tr = SACTrace.read(syndir + "/" + name )

            # synthetic
            data = dt_syn * convolve(syn_tr.data,stf[ic,:],'same')

            # write data to outdir
            syn_tr.data = data 
            syn_tr.write(outdir + '/' + name)

    return 0

if __name__ == "__main__":
    main()