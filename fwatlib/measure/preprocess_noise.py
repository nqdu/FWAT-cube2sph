from obspy.io.sac import SACTrace
import sys 
import numpy as np
from mpi4py import MPI
import os 

from utils import read_fwat_params,interpolate_syn
from utils import preprocess,dif1,cumtrapz1


def main():
    if len(sys.argv) != 4:
        print("Usage: python preprocess_noise.py iter evtid run_opt")
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
    pdict = read_fwat_params(f'solver/{mdir}/{evtname}/DATA/FWAT.PAR.yaml')['measure']['noise']

    # print log
    verbose = pdict['VERBOSE_MODE']
    CCODE = "." + pdict['CH_CODE']

    # read station coordinates
    stationfile = f'solver/{mdir}/' + f'{evtname}/DATA/STATIONS_FILTERED'
    statxt = np.loadtxt(stationfile,dtype=str,ndmin=2)
    nsta = statxt.shape[0]
    
    # synthetic parameters
    syndir = f'solver/{mdir}/' + f'{evtname}/OUTPUT_FILES/'
    name = statxt[0,1] + "." + statxt[0,0] + CCODE + "Z.sac"
    syn_z_hd = SACTrace.read(syndir + name,headonly=True)
    npt_syn = syn_z_hd.npts
    dt_syn = syn_z_hd.delta
    t0_syn = syn_z_hd.b

    # get frequency band/ group velocity
    Tmin_list = [x[0] for x in pdict['FILTER_BANDS']]
    Tmax_list = [x[1] for x in pdict['FILTER_BANDS']]
    vmin_list = [x[0] for x in pdict['GROUPVEL_WIN']]
    vmax_list = [x[1] for x in pdict['GROUPVEL_WIN']]

    # get components
    ncomp = len(pdict['COMPS'])
    components = pdict['COMPS']

    # write synthetic data to SYN dir
    if run_opt == 1:
        outdir =  f'fwat_data/{evtname}/'
        if myrank == 0: print("Synthetic Observed Data")
        os.makedirs(outdir,exist_ok=True)
        for i in range(myrank,nsta,nprocs):
            for icomp in range(ncomp):
                ch = components[icomp]
            
                name = statxt[i,1] + "." + statxt[i,0] + CCODE + ch + ".sac"
                syn_tr = SACTrace.read(syndir + "/" + name )

                # EGF or CCF
                if pdict['USE_EGF'] == False:
                    if myrank == 0: print("EGF => CCF ...")
                    data = syn_tr.data 
                    data1 = -cumtrapz1(data,syn_tr.delta)
                    syn_tr.data = data1 
                
                # write  to fwat_data/
                syn_tr.write(outdir + "/" + name )
        
        comm.Barrier()
        return 0

    # loop around each band
    nbands = len(Tmin_list)
    for ib in range(nbands):
        freqmin = 1. / Tmax_list[ib]
        freqmax = 1. / Tmin_list[ib]
        bandname = 'T%03g_T%03g' %(Tmin_list[ib],Tmax_list[ib])
        if myrank == 0 and run_opt != 1:
            print(f"preprocessing for band {bandname} ...")

        # new dt/nt for interpolate
        dt_inp = 0.01
        t0_inp = -10.
        
        # mkdir to save all files
        os.makedirs(syndir + f"{bandname}",exist_ok=True)

        # loop every station to do preprocessing
        for i in range(myrank,nsta,nprocs):
            for icomp in range(ncomp):
                ch = components[icomp]
                name = statxt[i,1] + "." + statxt[i,0] + CCODE + ch + ".sac"
                obs_tr = SACTrace.read(f'fwat_data/{evtid}/' + name )
                syn_tr = SACTrace.read(syndir + "/" + name )

                # parameters for obs
                t0_obs = obs_tr.b
                dt_obs = obs_tr.delta
                npt_obs = obs_tr.npts

                # new npts/t0 for interpolate
                t1_inp = t0_syn + (npt_syn - 1) * dt_syn
                npt1_inp = int((npt_obs - 1) * dt_obs / dt_inp)
                npt_cut = int((t1_inp - t0_inp) / dt_inp) + 1
                
                # interp obs/syn
                dat_inp1 = interpolate_syn(obs_tr.data,t0_obs,dt_obs,npt_obs,
                                        t0_obs + dt_inp,dt_inp,npt1_inp)
                syn_inp = interpolate_syn(syn_tr.data,t0_syn,dt_syn,npt_syn,
                                         t0_inp,dt_inp,npt_cut)

                # compute time derivative 
                if pdict['USE_EGF'] == False:
                    dat_inp1 = -dif1(dat_inp1,dt_inp)
                    if myrank == 0: print("CCFs => EGFs ...")
                
                # cut 
                dat_inp = interpolate_syn(dat_inp1,t0_obs + dt_inp,dt_inp,npt1_inp,
                                         t0_inp,dt_inp,npt_cut)

                # preprocess
                dat_inp = preprocess(dat_inp,dt_inp,freqmin,freqmax)
                syn_inp = preprocess(syn_inp,dt_inp,freqmin,freqmax)

                # normalize my amplitude of the window
                dist = obs_tr.dist
                win_b = np.floor((dist / vmax_list[ib] - Tmax_list[ib] / 2. + t0_inp) / dt_inp)
                win_e = np.floor((dist / vmin_list[ib] + Tmax_list[ib] / 2. - t0_inp) / dt_inp)
                win_b = int(max(win_b,0))
                win_e = int(min(win_e,len(dat_inp)-1))
                dat_inp *= np.max(np.abs(syn_inp[win_b:win_e])) / np.max(np.abs(dat_inp[win_b:win_e]))

                # write sac for measure adj
                filename = syndir + f"{bandname}/" + name + '.obs'
                obs_tr.delta = dt_inp; obs_tr.b = t0_inp 
                obs_tr.data = dat_inp 
                obs_tr.write(filename)
                filename = syndir + f"{bandname}/" + name + '.syn'
                syn_tr.delta = dt_inp; syn_tr.b = t0_inp 
                syn_tr.data = syn_inp 
                syn_tr.write(filename)

        comm.Barrier()

        # write window for measure adj
        if myrank == 0:
            f = open(syndir + f"{bandname}" + "/MEASUREMENT.WINDOWS","w")
            f.write("%d\n" %(nsta*ncomp))

            for i in range(nsta):
                for ic in range(ncomp):
                    ch = components[ic]
                    name = statxt[i,1] + "." + statxt[i,0] + CCODE + ch + ".sac"
                    head = SACTrace.read(syndir + f"{bandname}/" + name + '.obs',headonly=True)
                    dist = head.dist
                    tstart = dist / vmax_list[ib] - Tmax_list[ib] * 0.5 
                    tend = dist / vmin_list[ib] + Tmax_list[ib] * 0.5 
                    tstart = max(tstart,t0_inp)
                    tend = min(tend,head.e)
                    f.write("%s\n" %(name + '.obs'))
                    f.write("%s\n" %(name + '.syn'))
                    f.write("1\n")
                    f.write("%f %f\n" %(tstart,tend))
            f.close()
    comm.Barrier()            


if __name__ == "__main__":
    main()