import numpy as np 

def get_injection_time(evtid:str):
    temp = np.loadtxt('src_rec/injection_time',dtype=str,ndmin=2)
    find_src = False
    for i in range(temp.shape[0]):
        if temp[i,0] == evtid:
            t_inj = float(temp[i,1])
            find_src = True
            break
    if not find_src:
        print(f'please check {evtid} in src_rec/injection_time')
        exit(1)
    
    return t_inj 


def compute_ak135_time(evla:float,evlo:float,evdp:float,statxt,phase='P'):
    from obspy.taup import TauPyModel
    
    nsta = statxt.shape[0]
    t_ref = np.zeros((nsta))

    # create taup model
    model = TauPyModel("ak135")
    
    for i in range(nsta):
        stla = float(statxt[i,2])
        stlo = float(statxt[i,3])
        t_ref[i] = model.get_travel_times_geo(evdp,evla,evlo,stla,stlo,[phase])[0].time
    

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
        avgamp0 += np.max(np.abs(glob_obs[i,0,:]))
    avgamp0 /= nsta 
    igood =0 
    avgamp = 0. 
    for i in range(nsta):
        if np.max(np.abs(glob_obs[i,0,:])) - avgamp0 < 0.2 * avgamp0:
            igood += 1 
            avgamp += np.max(np.abs(glob_obs[i,0,:])) 
    avgamp /= igood 

    return avgamp