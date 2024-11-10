import numpy as np
from obspy import Trace
from obspy.io.sac import SACTrace
from obspy.taup import TauPyModel
import sys 
import os 
from mpi4py import MPI

from utils import interpolate_syn,get_fwat_params,get_average_amplitude
from utils import preprocess,cumtrapz1
from tele.deconit import time_decon
from scipy.signal import fftconvolve

def compute_ak135_time(evtid,statxt,do_ls = False):
    nsta = statxt.shape[0]
    t_ref = np.zeros((nsta))
    sourcefile = 'src_rec/sources.dat.tele'
    if do_ls:
        sourcefile = 'src_rec/sources.dat.ls.tele'
    temp = np.loadtxt(sourcefile,dtype=str,ndmin=2)

    # create taup model
    model = TauPyModel("ak135")

    # loop every source to find travel time
    find_source = False
    evla,evlo,evdp = 0,0,0
    for i in range(temp.shape[0]):
        if temp[i,0] == evtid:
            evla,evlo,evdp = float(temp[i,1]),float(temp[i,2]),float(temp[i,3])
            find_source = True
    
    if not find_source:
        print(f"cannot find source, please check {evtid} and {sourcefile}")
        exit(1)
    
    for i in range(nsta):
        stla = float(statxt[i,2])
        stlo = float(statxt[i,3])
        t_ref[i] = model.get_travel_times_geo(evdp,evla,evlo,stla,stlo,['P'])[0].time
    
    return t_ref

def seis_pca(glob_syn,glob_obs,stf_collect,dt_syn,components):
    from tele.pca.libpca import PCA

    nsta = glob_syn.shape[0]
    stf = np.zeros((2,glob_syn.shape[2]),'f4')
    temp = PCA(np.float32(stf_collect))
    # temp = np.mean(stf_collect,axis=0)
    #temp = temp / np.max(np.abs(temp))

    # = xnew[:,0]
    stf[0,:] = temp
    stf[1,:] = temp * 1.
    avgamp = np.zeros((2),'f4')

    # recover data
    ncomp = len(components)
    recp_syn = glob_syn.copy()
    for i in range(nsta):
        for ic in range(ncomp):
            recp_syn[i,ic,:] = dt_syn * fftconvolve(glob_syn[i,ic,:],stf[ic,:],'same')
    
    # calculate amplitude scale factor 
    print('\ncalculate ampl scale factors ...')
    avgarr = np.zeros((nsta),'f4')
    for ic in range(ncomp):
        ch = components[ic]
        print(f'Initial amp factors: {ch}')
        sumamp = 0.
        j = 0
        for i in range(nsta):
            if np.max(np.abs(recp_syn[i,ic,:])) > 0:
                avgarr[i] = np.dot(glob_obs[i,ic,:],recp_syn[i,ic,:]) / np.sum(recp_syn[i,ic,:]**2)
                sumamp += avgarr[i]
                j += 1 
            else :
                avgarr[i] = -9999.
        
        print(f"There are {j} non zero rec found")
        if j == 0:
            print(f"skip ic {ch}")
            continue 
        avgamp0 = sumamp / j 
        print("Initial average amp: %f" % avgamp0)
        print("selecting amp factors |A - A0| <=0.2")
        sumamp = 0. 
        j = 0
        for i in range(nsta):
            if np.abs(avgarr[i] - avgamp0) <=0.2 and avgarr[i] != -9999.:
                sumamp += avgarr[i]
                j += 1
        if j == 0 or j <= int(0.1 * nsta):
            print("Error: too little data satisfy |A-A0| <=0.2")
            exit(1)
        avgamp[ic] = sumamp / j 
        print("final average amp: %f" %(avgamp[ic]))

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
    pdict = get_fwat_params(f'solver/{mdir}/{evtname}/DATA/FWAT.PAR')
    
    # print log
    verbose = False
    if pdict['VERBOSE_MODE'].split()[0] == '.true.':
        verbose = True

    # read injection starttime 
    temp = np.loadtxt('src_rec/injection_time',dtype=str,ndmin=2)
    find_src = False
    for i in range(temp.shape[0]):
        if temp[i,0] == evtid:
            t_inj = float(temp[i,1])
            find_src = True
            break
    if not find_src:
        print(f'please check {evtid} in src_rec/injection_time')
    
    # read station coordinates
    stationfile = f'src_rec/STATIONS_{evtid}_globe'
    statxt = np.loadtxt(stationfile,dtype=str,ndmin=2)
    nsta = statxt.shape[0]

    # compute ak135 theoretical travel time for each station
    do_ls = False
    if run_opt == 2: do_ls=True
    tref = compute_ak135_time(evtid,statxt,do_ls)
    
    # synthetic data parameters
    syndir = f'solver/{mdir}/{evtname}/OUTPUT_FILES/'
    name = statxt[0,1] + "." + statxt[0,0] + ".BXZ.sac"
    syn_z_hd = SACTrace.read(syndir + name,headonly=True)
    npt_syn = syn_z_hd.npts
    dt_syn = syn_z_hd.delta
    t0_syn =  t_inj # starttime is injection time
    sac_head = syn_z_hd.copy()

    # get time window 
    win_tb = np.float32(pdict['TW_BEFORE'])
    win_te = np.float32(pdict['TW_AFTER'])
    npt2 = int((win_te + win_tb) / dt_syn)
    if npt2 // 2 * 2 != npt2:
        npt2 += 1

    # get frequency band
    Tmin_list = list(map(lambda x:float(x),pdict['SHORT_P'].split()))
    Tmax_list = list(map(lambda x:float(x),pdict['LONG_P'].split()))
    assert(len(Tmin_list) == len(Tmax_list))

    # get components
    ncomp = int(pdict['NRCOMP'])
    components = pdict['RCOMPS'].split()

    # if run_opt == 1, save synthetic data to SYN
    if run_opt == 1:
        # sum stf from all bands
        stf = np.zeros((ncomp,npt_syn))
        nsta = statxt.shape[0]
        for ib in range(len(Tmax_list)):
            bandname = 'T%03d_T%03d' %(Tmin_list[ib],Tmax_list[ib])
            for ic in range(ncomp):
                ch = components[ic]
                tr = SACTrace.read(f'src_rec/stf_{ch}.sac.{bandname}_{evtid}')
                stf[ic,:] += tr.data 
        
        # synthetic by using summed stf
        outdir = syndir + '/SYN/'
        os.makedirs(outdir,exist_ok=True)
        for i in range(myrank,nsta,nprocs):
            for ic in range(ncomp):
                ch = components[ic]
                name = statxt[i,1] + "." + statxt[i,0] + ".BX" + ch + ".sac"
                syn_tr = SACTrace.read(syndir + "/" + name )

                # synthetic
                data = dt_syn * fftconvolve(syn_tr.data,stf[ic,:],'same')

                # write data to outdir
                syn_tr.data = data 
                syn_tr.data = t_inj - tref[i]
                syn_tr.write(outdir + '/' + name)
        
        #
        return 0

    # now loop every band to estimate the stf 
    nbands = len(Tmin_list)
    for ib in range(nbands):
        freqmin = 1. / Tmax_list[ib]
        freqmax = 1. / Tmin_list[ib]
        bandname = 'T%03d_T%03d' %(Tmin_list[ib],Tmax_list[ib])
        
        # mkdir to save all files
        os.makedirs(syndir + f"{bandname}",exist_ok=True)

        # deconvolution to get stf 
        nsta = statxt.shape[0]
        stf_collect = np.zeros((nsta,npt_syn),'f4')
        glob_obs = np.zeros((nsta,ncomp,npt_syn),'f4')
        glob_syn = np.zeros((nsta,ncomp,npt_syn),'f4')
        
        # allocate tasks
        for i in range(myrank,nsta,nprocs):
            # check time window
            name = statxt[i,1] + "." + statxt[i,0] + ".BXZ" + ".sac"
            head = SACTrace.read(syndir + "/" + name,headonly=True)
            if tref[i] +  win_te > t0_syn + head.npts * head.delta:
                print("time window exceed data range")
                print("tref,win_te,t0_syn,npts,delta = ",tref[i],win_te,t0_syn,head.npts,head.delta)
                exit(1)

            # loop every component
            for ic in range(ncomp):
                ch = components[ic]
                name = statxt[i,1] + "." + statxt[i,0] + ".BX" + ch + ".sac"
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

                if ch == 'Z':
                    # deconv
                    stf1 = time_decon(u,w,dt_syn)
                    #stf1 = compute_stf_freq(u,w,dt_syn)
                    #stf1 = deconit(u,w,dt_syn,0.,1.5)
                    tr = Trace(stf1); tr.stats.delta = dt_syn
                    tr.filter("bandpass",freqmin=freqmin,freqmax=freqmax,zerophase = True,corners=4)
                    stf1 = np.float32(tr.data) 

                    # save field
                    stf_collect[i,:] = stf1
                    
                glob_obs[i,ic,:] = u
                glob_syn[i,ic,:] = w 

        # all reduce
        tmp = comm.allreduce(glob_obs); glob_obs = tmp
        tmp = comm.allreduce(glob_syn); glob_syn = tmp
        tmp = comm.allreduce(stf_collect); stf_collect = tmp
        
        # compute PCA analysis to get primary STF
        # PCA to get optimal stf
        stf = None 
        if myrank == 0:
            # check if sac file exists
            stf_names = [f'src_rec/stf_{ch}.sac.' + f"{bandname}" + f"_{evtid}" for ch in components]
            flag = True 
            for ic in range(ncomp):
                flag = flag and os.path.isfile(stf_names[ic])
            if flag:
                stf = np.zeros((ncomp,npt_syn),'f4')
                for ic in range(ncomp):
                    stf[ic,:] = SACTrace.read(stf_names[ic]).data
            else:
                stf = seis_pca(glob_syn,glob_obs,stf_collect,dt_syn,components)
                for ic in range(ncomp):
                    sac_head.data = stf[ic,:]
                    sac_head.write(stf_names[ic])
            
            # write stf to syndir
            path = syndir + f'{bandname}/'
            for ic in range(ncomp):
                ch = components[ic]
                sac_head.data = stf[ic,:]
                sac_head.b = 0
                sac_head.delta = dt_syn
                sac_head.write(path + f'stf_{ch}.sac')
        stf = comm.bcast(stf)

        # save semdata 
        if verbose:
            for i in range(myrank,nsta,nprocs):
                for ic in range(ncomp):
                    ch = components[ic]
                    name = statxt[i,1] + "." + statxt[i,0] + ".BX" + ch + ".sac"
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
                name = statxt[i,1] + "." + statxt[i,0] + ".BX" + ch + ".sac"
                sac_head = SACTrace.read(syndir + "/" + name,headonly=True)
                name = syndir + f"{bandname}" + "/" +name 
                sac_head.data = glob_obs[i,ic,:]
                sac_head.b = 0.
                sac_head.t0 = t0
                sac_head.write(name + ".obs")
                data = dt_syn * fftconvolve(glob_syn[i,ic,:],stf[ic,:],'same')
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
                    name = statxt[i,1] + "." + statxt[i,0] + ".BX" + ch + ".sac"
                    f.write("%s\n" %(name + '.obs'))
                    f.write("%s\n" %(name + '.syn'))
                    f.write("1\n")
                    f.write("%f %f\n" %(tstart,tend))

            f.close()
    
    comm.Barrier()

if __name__ == "__main__":
    main()
