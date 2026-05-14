import numpy as np 
from typing import Final
from fwat.adjoint.MeasureStats import MeasureStats

def measure_adj(t0_inp: float,dt_inp: float,npt_inp: int,
                t0_syn: float,dt_syn: float,npt_syn: int,
                tstart: float,tend: float,imeas:int,
                tlong: float,tshort: float,verbose:bool,
                obs_data: np.ndarray,syn_data: np.ndarray,
                compute_adj_source: bool = True,
                run_bandpass: bool = False,display_details: bool = False,
                output_measure_files: bool = False,
                tshift_min: float = -4.5,tshift_max: float = 4.5,
                dlna_min: float = -1.5,dlna_max: float = 1.5,
                cc_min: float = 0.8,err_type: int = 1,
                dt_sigma_min: float = 1.,dlna_sigma_min: float = 0.5,
                itaper: int = 1,wtr: float = 0.02,npi: float = 2.5,
                dt_fac: float = 2.,err_fac: float = 2.5,
                dt_max_scale: float = 3.5,
                ncyle_in_window: float = 1.5,
                use_physical_disp: bool = False) -> tuple[MeasureStats,np.ndarray]:
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
    
    stats = MeasureStats(
        adj_type = str(imeas),
        misfit=tr,
        tstart=tstart,
        tend=tend,
        tr_chi=tr,
        am_chi=amp,
        tshift = window[6]
    )
    
    return stats,adj
