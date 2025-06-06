import sys 
import os 

from mpi4py import MPI
import numpy as np
from obspy.io.sac import SACTrace

from utils import preprocess,read_params,interpolate_syn
from utils import get_simu_info,dif1,cumtrapz1
from utils import rotate_RT_to_NE,get_source_loc
from tele.tele import compute_ak135_time,get_injection_time

def main():
    if len(sys.argv) != 4:
        print("Usage: python preprocess_sks.py iter evtid run_opt")
        exit(1)
    
    # get input parameter
    iter = int(sys.argv[1])
    evtid = sys.argv[2]
    run_opt = int(sys.argv[3])
    evtname = evtid
    mdir = 'M%02d' %iter
    if run_opt == 2:
        mdir = mdir + '.ls'

    # initialize mpi
    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()
    nprocs = comm.Get_size()

    # load paramfile as dictionary
    pdict = read_params(f'solver/{mdir}/{evtname}/DATA/FWAT.PAR.yaml')['measure']['sks']
    
    # print log
    verbose = pdict['VERBOSE_MODE']
    CCODE = "." + pdict['CH_CODE']

    # for SKS R/T components are required
    components = ['R','T']
    ncomp = 2

    # get frequency band
    Tmin_list = [x[0] for x in pdict['FILTER_BANDS']]
    Tmax_list = [x[1] for x in pdict['FILTER_BANDS']]

    # read injection starttime 
    t_inj = get_injection_time(evtid)

    # read station coordinates
    stationfile = f'src_rec/STATIONS_{evtid}_globe'
    statxt = np.loadtxt(stationfile,dtype=str,ndmin=2)
    nsta = statxt.shape[0]

    # read source loc
    sourcefile = "./src_rec/sources.dat.sks"
    evla,evlo,evdp = get_source_loc(evtid,sourcefile)

    # compute ak135 theoretical travel time for each station
    if '_' not in evtid:
        phase = 'SKS'
    else:
        phase = evtid.split('_')[0]
    tref = compute_ak135_time(evla,evlo,evdp,statxt,phase)
    
    # synthetic data parameters
    syndir = f'solver/{mdir}/{evtname}/OUTPUT_FILES/'
    name = statxt[0,1] + "." + statxt[0,0] + CCODE + f"T.sac"
    _,dt_syn,npt_syn = get_simu_info(syndir + name)
    t0_syn = t_inj * 1

    # get time window 
    win_tb,win_te = pdict['TIME_WINDOW']
    npt2 = int((win_te + win_tb) / dt_syn)
    if npt2 // 2 * 2 != npt2:
        npt2 += 1

    # if run_opt == 1, save synthetic data to SYN
    if run_opt == 1:
        # synthetic by using summed stf
        outdir =  f'fwat_data/{evtname}/'
        os.makedirs(outdir,exist_ok=True)
        if myrank == 0: print("Synthetic Observed Data")
        for i in range(myrank,nsta,nprocs):
            ncomp = len(components)
            for ic in range(ncomp):
                ch = components[ic]
                name = statxt[i,1] + "." + statxt[i,0]  + CCODE + ch + ".sac"
                syn_tr = SACTrace.read(syndir + "/" + name )

                # write data to outdir
                syn_tr.b = t_inj - tref[i]
                syn_tr.write(outdir + '/' + name)
        
        # finish here
        comm.Barrier()
        return 0

    # now loop every band to compute adjoint source
    nbands = len(Tmin_list)
    for ib in range(nbands):
        freqmin = 1. / Tmax_list[ib]
        freqmax = 1. / Tmin_list[ib]
        bandname = 'T%03g_T%03g' %(Tmin_list[ib],Tmax_list[ib])
        
        # mkdir to save all files
        os.makedirs(syndir + f"{bandname}",exist_ok=True)
        os.makedirs(syndir + f"{bandname}/OUTPUT_FILES",exist_ok=True)
        
        # loop each station to compute SI and adjoint source
        misfits = np.zeros((nsta,3))
        nsta = statxt.shape[0]
        for i in range(myrank,nsta,nprocs):
            # check time window
            name = statxt[i,1] + "." + statxt[i,0]  +  CCODE + f"T.sac"
            head = SACTrace.read(syndir + "/" + name,headonly=True)
            if tref[i] +  win_te > t0_syn + head.npts * head.delta:
                print("time window exceed data range")
                print("tref,win_te,t0_syn,npts,delta = ",tref[i],win_te,t0_syn,head.npts,head.delta)
                exit(1)
            
            # syn/obs
            syn_data = np.zeros((2,npt_syn))
            obs_data = np.zeros((2,npt_syn))

            # loop every component
            for ic in range(ncomp):
                ch = components[ic]
                name = statxt[i,1] + "." + statxt[i,0] + CCODE + ch + ".sac"
                obs_tr = SACTrace.read(f'fwat_data/{evtid}/' + name )
                syn_tr = SACTrace.read(syndir + "/" + name )
                
                # paramters for syn
                t0_syn = t_inj

                # parameters for obs
                t0_obs = obs_tr.b + tref[i] # theoretical time is 0 in data 
                dt_obs = obs_tr.delta
                npt_obs = obs_tr.npts

                # preprocessing
                obs_tr.data = preprocess(obs_tr.data,dt_obs,freqmin,freqmax)
                syn_tr.data = preprocess(syn_tr.data,dt_syn,freqmin,freqmax)
                
                # interpolate the obs/syn data to the same sampling of syn data
                t0_inp = tref[i] - win_tb
                u1 = interpolate_syn(obs_tr.data,t0_obs,dt_obs,npt_obs,t0_inp,dt_syn,npt2)
                w1 = interpolate_syn(syn_tr.data,t0_syn,dt_syn,npt_syn,t0_inp,dt_syn,npt2)
                u = interpolate_syn(u1,t0_inp,dt_syn,npt2,t0_syn,dt_syn,npt_syn)
                w = interpolate_syn(w1,t0_inp,dt_syn,npt2,t0_syn,dt_syn,npt_syn)

                # copy to global data
                obs_data[ic,:] = u * 1.
                syn_data[ic,:] = w * 1.
            
            # compute SI for synthetic and obs data
            Robs_dot = dif1(obs_data[0,:],dt_syn)
            norm = cumtrapz1(Robs_dot**2,dt_syn)[-1]
            if norm < 1.0e-15:
                norm = 0.
            else:
                norm = 1. / norm
            Tobs = obs_data[1,:]
            Tsyn = syn_data[1,:]
            SI_obs = -2. * cumtrapz1(Robs_dot*Tobs,dt_syn)[-1] * norm
            SI_syn = -2. * cumtrapz1(Robs_dot*Tsyn,dt_syn)[-1] * norm
            L = 0.5 * (SI_syn - SI_obs)**2
            misfits[i,:] = [SI_obs,SI_syn,L]

            # save obs and syn data 
            baz = 0
            t0 = tref[i] - t_inj
            for ic in range(ncomp):
                ch = components[ic]
                name = statxt[i,1] + "." + statxt[i,0] + CCODE + ch + ".sac"
                sac_head = SACTrace.read(syndir + "/" + name,headonly=True)
                name = syndir + f"{bandname}" + "/" +name 
                sac_head.data = obs_data[ic,:]
                sac_head.b = 0.
                sac_head.t0 = t0
                sac_head.write(name + ".obs")
                sac_head.data = syn_data[ic,:]
                sac_head.write(name + ".syn")

                if ic == 0:
                    baz = sac_head.baz
            
            # rotate adjoint source
            adjsrc_T = -2. * norm * (SI_syn - SI_obs) * Robs_dot
            adjsrc_R = adjsrc_T * 0.
            adjsrc_Z = adjsrc_T * 0.
            adjsrc_N,adjsrc_E = rotate_RT_to_NE(adjsrc_R,adjsrc_T,baz)

            # save adjoint source
            data = np.zeros((3,npt_syn,2))
            data[:,:,0] = np.arange(npt_syn) * dt_syn
            data[0,:,1] = adjsrc_N
            data[1,:,1] = adjsrc_E 
            data[2,:,1] = adjsrc_Z
            for ic,ch in enumerate(['N','E','Z']):
                filename = f"{syndir}/{bandname}/OUTPUT_FILES/" + statxt[i,1] + "." + statxt[i,0] + CCODE + f"{ch}.adj"
                np.savetxt(filename,data[ic,:,:],fmt="%g")
        
        # gather misfits
        tmp = comm.allreduce(misfits); misfits = tmp * 1.

        # write window_chi
        if myrank == 0:
            fio = open(f"{syndir}/{bandname}/window_chi","w")
            for i in range(nsta):
                net = statxt[i,1]
                sta = statxt[i,0]
                ch = pdict['CH_CODE'] + 'T'
                temp = np.zeros((20),'f4')
                fio.write(
                    f"{evtid} {sta} {net} {ch} {i} 0 %f %f {' '.join(map(str,temp))} %f %f %f 0. 0.\n" %(
                        -win_tb + tref[i] - t_inj,win_te + tref[i] - t_inj, 
                        misfits[i,-1],misfits[i,0],misfits[i,1]
                    )
                )
            fio.close()
    
    comm.Barrier()

if __name__ == "__main__":
    main()
