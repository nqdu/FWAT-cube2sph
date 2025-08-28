import numpy as np 
from typing import Final

def measure_adj(t0_inp,dt_inp,npt_inp,
                t0_syn,dt_syn,npt_syn,
                tstart,tend,imeas:int,
                tlong,tshort,verbose:bool,
                obs_data,syn_data,
                compute_adj_source = True,
                run_bandpass=False,display_details=False,
                output_measure_files = False,
                tshift_min = -4.5,tshift_max=4.5,
                dlna_min = -1.5,dlna_max=1.5,
                cc_min = 0.8,err_type = 1,
                dt_sigma_min=1.,dlna_sigma_min=0.5,
                itaper = 1,wtr=0.02,npi=2.5,
                dt_fac=2.,err_fac=2.5,
                dt_max_scale = 3.5,
                ncyle_in_window=1.5,
                use_physical_disp=False):
    """
    MEASURE ADJ wrapper
    """
    from .lib import libmeas
    
    if imeas == 5:
        itaper = 2

    tr,amp, window,adj =  \
        libmeas.measure(
            t0_inp,dt_inp,npt_inp,
            t0_syn,dt_syn,npt_syn,
            tstart,tend,imeas,
            tlong,tshort,verbose,
            obs_data,syn_data,run_bandpass,
            display_details,output_measure_files,
            compute_adj_source,tshift_min,
            tshift_max,dlna_min,dlna_max,
            cc_min,err_type,dt_sigma_min,
            dlna_sigma_min,itaper,wtr,npi,
            dt_fac,err_fac,dt_max_scale,
            ncyle_in_window,use_physical_disp
        )
    
    return tr,amp,window,adj
