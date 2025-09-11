import numpy as np 

def get_injection_time(evtid:str):
    temp = np.loadtxt('src_rec/injection_time',dtype=str,ndmin=2)
    find_src = False
    t_inj = 0.
    for i in range(temp.shape[0]):
        if temp[i,0] == evtid:
            t_inj = float(temp[i,1])
            find_src = True
            break
    if not find_src:
        print(f'please check {evtid} in src_rec/injection_time')
        exit(1)
    
    return t_inj 

def compute_ak135_time(evla:float,evlo:float,evdp:float,
                       stla:np.ndarray,stlo:np.ndarray,phase='P'):
    from obspy.taup import TauPyModel
    from mpi4py import MPI

    # get mpi info 
    comm = MPI.COMM_WORLD
    myrank = comm.Get_rank()
    nprocs = comm.Get_size()
    
    nsta = len(stla)
    t_ref = np.zeros((nsta))

    # create taup model
    model = TauPyModel("ak135")
    
    for i in range(myrank,nsta,nprocs):
        t_ref[i] = model.get_travel_times_geo(evdp,evla,evlo,stla[i],stlo[i],[phase])[0].time
    
    # gather all
    tmp = comm.allreduce(t_ref); t_ref = tmp * 1.

    return t_ref

def get_average_amplitude(glob_obs,icomp):
    """
    get the average amplitude for icomp-th component 

    Parameters:
        glob_obs: np.ndarray
            observation data in global, shape(nsta,ncomp,nt)
        icomp: int
            icomp-th componnet
    """
    nsta = glob_obs.shape[0]
    avgamp0 = 0.
    for i in range(nsta):
        avgamp0 += np.max(np.abs(glob_obs[i,icomp,:]))
    avgamp0 /= nsta 
    igood =0 
    avgamp = 0. 
    for i in range(nsta):
        if np.max(np.abs(glob_obs[i,icomp,:])) - avgamp0 < 0.2 * avgamp0:
            igood += 1 
            avgamp += np.max(np.abs(glob_obs[i,icomp,:])) 
    avgamp /= igood 

    return avgamp


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
    w,v = np.linalg.eig(rec @ rec.T / nt)

    # sort eigenvalues
    w = np.real(w) 
    idx = np.argsort(w)[::-1]
    w = w[idx]
    v = v[:,idx]

    # first Principle component is what we need
    stf = (rec.T @ v[:,0]).real
    
    # make sure stf has same polarity as the average one
    stf_avg = np.mean(stf_collect,axis=0)
    sig1 = np.sign(stf_avg[np.argmax(abs(stf_avg))])
    sig2 = np.sign(stf[np.argmax(abs(stf))])

    # rescale
    stf = stf * sig1 * sig2 

    return stf,w

def compute_stf(glob_syn,glob_obs,dt_syn,freqmin,freqmax,components):
    # load stf function
    from fwat.measure.tele.deconit import time_decon
    from mpi4py import MPI
    from obspy import Trace
    from scipy.signal import convolve,correlate

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
            stf1 = np.asarray(tr.data,dtype='f4')

            # save stf
            stf_collect[i,:] = stf1 * 1.
    
    # all gather
    tmp = comm.allreduce(stf_collect); stf_collect = tmp * 1.

    # now get stf  by using PCA/ average    
    stf = np.zeros((2,glob_syn.shape[2]),'f4')
    temp,w = seis_pca(stf_collect)
    # temp = np.mean(stf_collect,axis=0)
    #temp = temp / np.max(np.abs(temp))
    stf[0,:] = temp
    stf[1,:] = temp * 1.

    if myrank == 0:
        print("\nPCA first component ratio = %f" % (w[0] / np.sum(w)))
    
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