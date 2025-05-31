import numpy as np
from obspy import Trace
from obspy.io.sac import SACTrace
import sys 
import os 
from mpi4py import MPI

from utils import interpolate_syn,read_params,preprocess
from utils import get_simu_info,get_source_loc
from tele.tele import get_average_amplitude,compute_ak135_time
from tele.tele import get_injection_time
from scipy.signal import convolve,correlate

def shift_data(u,dt,t0):
    u_spec = np.fft.fft(u)
    freqs = np.fft.fftfreq(len(u),dt)

    # shift in freq domain
    u_spec *= np.exp(-1j * np.pi * 2 * freqs * t0)

    u_out = np.fft.ifft(u_spec).real 

    return u_out 

def seis_pca(stf_collect:np.ndarray):
    _,nt = stf_collect.shape 

    # remove average
    rec = stf_collect - np.mean(stf_collect,axis=1,keepdims=True)

    # compute eigenvalue/eigenvector
    _,v = np.linalg.eig(rec @ rec.T / nt)

    # first Principle component is what we need
    stf = (rec.T @ v[:,0]).real
    
    # make sure stf has same polarity as the average one
    stf_avg = np.mean(stf_collect,axis=0)
    sig1 = np.sign(stf_avg[np.argmax(abs(stf_avg))])
    sig2 = np.sign(stf[np.argmax(abs(stf))])

    # rescale
    stf = stf * sig1 * sig2 

    return stf

def compute_stf(glob_syn,glob_obs,dt_syn,freqmin,freqmax,components):
    # load stf function
    from tele.deconit import time_decon

    # get dimension
    nsta,ncomp,nt = glob_syn.shape

    # compute stf station by station 
    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()
    nprocs = comm.Get_size()
    stf_collect = np.zeros((nsta,nt))
    for i in range(myrank,nsta,nprocs):
        for ic in range(ncomp):
            ch = components[ic]
            if ch != 'Z': continue

            # get syn/obs data 
            u = glob_obs[i,ic,:] * 1. 
            w = glob_syn[i,ic,:] * 1.

            # time deconvolution
            stf1 = time_decon(u,w,dt_syn)
            #stf1 = compute_stf_freq(u,w,dt_syn)
            #stf1 = deconit(u,w,dt_syn,0.,1.5)
            tr = Trace(stf1); tr.stats.delta = dt_syn
            tr.filter("bandpass",freqmin=freqmin,freqmax=freqmax,zerophase = True,corners=4)
            stf1 = np.float32(tr.data) 

            # save stf
            stf_collect[i,:] = stf1 * 1.
    
    # all gather
    tmp = comm.allreduce(stf_collect); stf_collect = tmp * 1.

    # now get stf  by using PCA/ average    
    stf = np.zeros((2,glob_syn.shape[2]),'f4')
    temp = seis_pca(stf_collect)
    # temp = np.mean(stf_collect,axis=0)
    #temp = temp / np.max(np.abs(temp))
    stf[0,:] = temp
    stf[1,:] = temp * 1.
    
    # compute time shift between syn and obs
    time_shift = np.zeros((nsta))
    for i in range(myrank,nsta,nprocs):
        for ic in range(ncomp):
            if components[ic] != 'Z': continue
            out = dt_syn * convolve(glob_syn[i,ic,:],stf[ic,:],'same')
            corr = correlate(glob_obs[i,ic,:],out,"full")
            time_shift[i] = (np.argmax(abs(corr)) - nt + 1) * dt_syn 
    tmp = comm.allreduce(time_shift); time_shift = tmp * 1.
    t0 = np.mean(time_shift)

    # shift stf to remove tshift residulas
    for i in range(ncomp):
        stf[i,:] = shift_data(stf[i,:],dt_syn,t0)
    
    # recover data
    recp_syn = glob_syn * 0
    for i in range(myrank,nsta,nprocs):
        for ic in range(ncomp):
            recp_syn[i,ic,:] = dt_syn * convolve(glob_syn[i,ic,:],stf[ic,:],'same')
    tmp = comm.allreduce(recp_syn); recp_syn = tmp * 1.
    
    # calculate amplitude scale factor 
    avgamp = np.zeros((2))
    if myrank == 0: print('\ncalculate ampl scale factors ...')
    for ic in range(ncomp):
        avgarr = np.zeros((nsta))
        ch = components[ic]
        if myrank == 0: print(f'Initial amp factors: {ch}')
        sumamp = 0.
        j = 0
        for i in range(myrank,nsta,nprocs):
            if np.max(np.abs(recp_syn[i,ic,:])) > 0:
                avgarr[i] = np.dot(glob_obs[i,ic,:],recp_syn[i,ic,:]) / np.sum(recp_syn[i,ic,:]**2)
                sumamp += avgarr[i]
                j += 1 
            else :
                avgarr[i] = -9999.
        tmp = comm.allreduce(avgarr); avgarr = tmp * 1.
        j_all = comm.allreduce(j,MPI.SUM); j = j_all 
        sumamp_all = comm.allreduce(sumamp,MPI.SUM); sumamp = sumamp_all * 1.
        avgamp0 = sumamp / j 

        if j == 0:
            if myrank == 0: print(f"skip ic {ch}")
            continue 
        if myrank == 0:
            print(f"There are {j} non zero rec found")
            print("Initial average amp: %f" % avgamp0)
            print("selecting amp factors |A - A0| <=0.2")
        sumamp = 0. 
        j = 0
        for i in range(myrank,nsta,nprocs):
            if np.abs(avgarr[i] - avgamp0) <=0.2 and avgarr[i] != -9999.:
                sumamp += avgarr[i]
                j += 1
        j_all = comm.allreduce(j,MPI.SUM); j = j_all 
        sumamp_all = comm.allreduce(sumamp,MPI.SUM); sumamp = sumamp_all * 1.

        # check usage
        if j == 0 or j <= int(0.1 * nsta):
            print("Error: too little data satisfy |A-A0| <=0.2")
            exit(1)

        avgamp[ic] = sumamp / j 
        if myrank == 0: print("final average amp: %f" %(avgamp[ic]))

        # normalize stf
        stf[ic,:] *= avgamp[ic]

    return stf

def main():
    if len(sys.argv) != 4:
        print("Usage: python preprocess_tele.py iter evtid run_opt")
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
    pdict = read_params(f'solver/{mdir}/{evtname}/DATA/FWAT.PAR.yaml')['measure']['tele']
    
    # print log
    verbose = pdict['VERBOSE_MODE']
    CCODE = "." + pdict['CH_CODE']

    # get frequency band
    Tmin_list = [x[0] for x in pdict['FILTER_BANDS']]
    Tmax_list = [x[1] for x in pdict['FILTER_BANDS']]

    # get components
    ncomp = len(pdict['COMPS'])
    components = pdict['COMPS']

    # read injection starttime 
    t_inj = get_injection_time(evtid)
    
    # read station coordinates
    stationfile = f'src_rec/STATIONS_{evtid}_globe'
    statxt = np.loadtxt(stationfile,dtype=str,ndmin=2)
    nsta = statxt.shape[0]

    # read source loc
    sourcefile = "./src_rec/sources.dat.tele"
    evla,evlo,evdp = get_source_loc(evtid,sourcefile)

    # compute ak135 theoretical travel time for each station
    tref = compute_ak135_time(evla,evlo,evdp,statxt,'P')
    
    # synthetic data parameters
    syndir = f'solver/{mdir}/{evtname}/OUTPUT_FILES/'
    name = statxt[0,1] + "." + statxt[0,0] + CCODE + f"{components[0]}.sac"
    _,dt_syn,npt_syn = get_simu_info(syndir + name)
    t0_syn = t_inj # starttime is injection time

    # get time window 
    win_tb,win_te = pdict['TIME_WINDOW']
    npt2 = int((win_te + win_tb) / dt_syn)
    if npt2 // 2 * 2 != npt2:
        npt2 += 1

    # if run_opt == 1, save synthetic data to SYN
    if run_opt == 1:
        # sum stf from all bands
        stf = np.zeros((ncomp,npt_syn))
        nsta = statxt.shape[0]
        for ib in range(len(Tmax_list)):
            bandname = 'T%03g_T%03g' %(Tmin_list[ib],Tmax_list[ib])
            for ic in range(ncomp):
                ch = components[ic]
                tr = SACTrace.read(f'src_rec/stf_{ch}.sac.{bandname}_{evtid}')
                stf[ic,:] += tr.data 
        
        # synthetic by using summed stf
        outdir =  f'fwat_data/{evtname}/'
        os.makedirs(outdir,exist_ok=True)
        if myrank == 0: print("Synthetic Observed Data")
        for i in range(myrank,nsta,nprocs):
            for ic in range(ncomp):
                ch = components[ic]
                name = statxt[i,1] + "." + statxt[i,0]  + CCODE + ch + ".sac"
                syn_tr = SACTrace.read(syndir + "/" + name )

                # synthetic
                data = dt_syn * convolve(syn_tr.data,stf[ic,:],'same')

                # write data to outdir
                syn_tr.data = data 
                syn_tr.b = t_inj - tref[i]
                syn_tr.write(outdir + '/' + name)
        
        # finish here
        comm.Barrier()
        return 0


    # now loop every band to estimate the stf 
    nbands = len(Tmin_list)
    for ib in range(nbands):
        freqmin = 1. / Tmax_list[ib]
        freqmax = 1. / Tmin_list[ib]
        bandname = 'T%03g_T%03g' %(Tmin_list[ib],Tmax_list[ib])
        
        # mkdir to save all files
        os.makedirs(syndir + f"{bandname}",exist_ok=True)

        # deconvolution to get stf 
        nsta = statxt.shape[0]
        glob_obs = np.zeros((nsta,ncomp,npt_syn))
        glob_syn = np.zeros((nsta,ncomp,npt_syn))
        
        # allocate tasks
        for i in range(myrank,nsta,nprocs):
            # check time window
            name = statxt[i,1] + "." + statxt[i,0]  +  CCODE + f"{components[0]}.sac"
            head = SACTrace.read(syndir + "/" + name,headonly=True)
            if tref[i] +  win_te > t0_syn + head.npts * head.delta:
                print("time window exceed data range")
                print("tref,win_te,t0_syn,npts,delta = ",tref[i],win_te,t0_syn,head.npts,head.delta)
                exit(1)

            # loop every component
            for ic in range(ncomp):
                ch = components[ic]
                name = statxt[i,1] + "." + statxt[i,0] + CCODE + ch + ".sac"
                syn_tr = SACTrace.read(syndir + "/" + name )
                obs_tr = SACTrace.read(f'fwat_data/{evtid}/' + name )
                
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

                # save to  global array   
                glob_obs[i,ic,:] = u
                glob_syn[i,ic,:] = w 

        # all reduce
        tmp = comm.allreduce(glob_obs); glob_obs = tmp * 1.
        tmp = comm.allreduce(glob_syn); glob_syn = tmp * 1.
        
        # get source time function 
        # check if sac file exists
        stf_names = [f'src_rec/stf_{ch}.sac.' + f"{bandname}" + f"_{evtid}" for ch in components]
        has_stf_flag = True 
        for ic in range(ncomp):
            has_stf_flag = has_stf_flag and os.path.isfile(stf_names[ic])
        if has_stf_flag: # read stf from directory
            stf = np.zeros((ncomp,npt_syn),'f4')
            for ic in range(ncomp):
                stf[ic,:] = SACTrace.read(stf_names[ic]).data
        else: # estimate stf
            stf = compute_stf(glob_syn,glob_obs,dt_syn,freqmin,freqmax,components)

            # save to stfnames
            sac_head = SACTrace(b=0,npts=npt_syn,delta=dt_syn)
            for ic in range(ncomp):
                sac_head.data = stf[ic,:]
                if myrank ==0 : sac_head.write(stf_names[ic])

        # save stf to syndir
        if myrank == 0:
            path = syndir + f'{bandname}/'
            sac_head = SACTrace(b=0,npts=npt_syn,delta=dt_syn)
            for ic in range(ncomp):
                ch = components[ic]
                sac_head.data = stf[ic,:]
                sac_head.write(path + f'stf_{ch}.sac')

        # wait all jobs finished
        comm.Barrier()

        # save semdata 
        if verbose:
            for i in range(myrank,nsta,nprocs):
                for ic in range(ncomp):
                    ch = components[ic]
                    name = statxt[i,1] + "." + statxt[i,0] + CCODE + ch + ".sac"
                    sac_head = SACTrace.read(syndir + "/" + name,headonly=True)
                    name = syndir + f"{bandname}" + "/" +name 
                    sac_head.data = glob_syn[i,ic,:]
                    sac_head.write(name + ".org")
                    sac_head.b = 0
        
        # estimate averagemap
        avgamp = get_average_amplitude(glob_obs,0)
        if myrank == 0: 
            print("average amplitude for data gather: %g" %(avgamp))
            f = open(syndir + f"{bandname}" + "/average_amp.dat","w")
            f.write("%g\n%g\n" %(avgamp,dt_syn))
            f.close()

        # save obs and recovered data for measure adj
        for i in range(myrank,nsta,nprocs):
            t0 = tref[i] - t_inj
            for ic in range(ncomp):
                ch = components[ic]
                name = statxt[i,1] + "." + statxt[i,0] + CCODE + ch + ".sac"
                sac_head = SACTrace.read(syndir + "/" + name,headonly=True)
                name = syndir + f"{bandname}" + "/" +name 
                sac_head.data = glob_obs[i,ic,:]
                sac_head.b = 0.
                sac_head.t0 = t0
                sac_head.write(name + ".obs")
                data = dt_syn * convolve(glob_syn[i,ic,:],stf[ic,:],'same')
                sac_head.data = data
                sac_head.write(name + ".syn")

        # write measure_adj window
        if myrank == 0:
            f = open(syndir + f"{bandname}" + "/MEASUREMENT.WINDOWS","w")
            f.write("%d\n" %(nsta*ncomp))

            for i in range(nsta):
                # time window
                tstart = tref[i] - t_inj - win_tb 
                tend = tref[i] - t_inj + win_te

                for ic in range(ncomp):
                    ch = components[ic]
                    name = statxt[i,1] + "." + statxt[i,0] + CCODE + ch + ".sac"
                    f.write("%s\n" %(name + '.obs'))
                    f.write("%s\n" %(name + '.syn'))
                    f.write("1\n")
                    f.write("%f %f\n" %(tstart,tend))

            f.close()

            # write MEASUREMENTS.PAR
            from utils import measure_adj_file
            imeas = pdict['IMEAS']
            ccode = pdict['CH_CODE']
            txt = measure_adj_file(0,dt_syn,npt_syn,0.99*np.min(Tmin_list),1.01*np.max(Tmax_list),imeas,ccode)
            f = open(f"solver/{mdir}/{evtname}/DATA/MEASUREMENT.PAR","w")
            f.write(txt)
            f.close()
    
    comm.Barrier()

if __name__ == "__main__":
    main()
