import numpy as np
from obspy.io.sac import SACTrace
import sys 
import os 

from utils import get_fwat_params

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
    
    # synthetic paramters
    syndir = 'solver/M%02d' %(iter) + f'.set{evtid}/OUTPUT_FILES'
    name = statxt[0,1] + "." + statxt[0,0] + ".BXZ.sac"
    syn_z_hd = SACTrace.read(syndir + name,headonly=True)
    npt_syn = syn_z_hd.npts
    dt_syn = syn_z_hd.delta

    # read fwat parameters
    pdict = get_fwat_params('solver/M%02d' %(iter) + f'.set{evtid}/DATA/FWAT.PAR')

    # get frequency band
    Tmin_list = list(map(lambda x:float(x),pdict['SHORT_P'].split()))
    Tmax_list = list(map(lambda x:float(x),pdict['LONG_P'].split()))
    assert(len(Tmin_list) == len(Tmax_list))

    # sum for stf for each period band if required
    ncomp = int(pdict['NRCOMP'])
    components = pdict['RCOMPS'].split()
    stf = np.zeros((ncomp,npt_syn))
    nsta = statxt.shape[0]
    for ib in range(len(Tmax_list)):
        bandname = 'T%03d_T%03d' %(Tmin_list[ib],Tmax_list[ib])
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
            name = statxt[i,1] + "." + statxt[i,0] + ".BX" + ch + ".sac"
            syn_tr = SACTrace.read(syndir + "/" + name )

            # synthetic
            data = dt_syn * np.convolve(syn_tr.data,stf[ic,:],'same')

            # write data to outdir
            syn_tr.data = data 
            syn_tr.write(outdir + '/' + name)

    return 0

if __name__ == "__main__":
    main()