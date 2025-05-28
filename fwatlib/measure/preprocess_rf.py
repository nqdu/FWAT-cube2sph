import numpy as np
from obspy.io.sac import SACTrace
import sys 
import os 
from mpi4py import MPI

from utils import interpolate_syn,read_params
from tele.tele import get_average_amplitude
from tele.deconit import deconit,gauss_filter,apply_gaussian

def main():
    if len(sys.argv) != 4:
        print("Usage: python preprocess_rf.py iter evtid run_opt")
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
    pdict = read_params(f'solver/{mdir}/{evtname}/DATA/FWAT.PAR.yaml')['measure']['rf']
    
    # print log
    verbose = pdict['VERBOSE_MODE']
    CCODE = "." + pdict['CH_CODE']

    # read fkfile 
    #fktimes = np.loadtxt('solver/M%02d' %(iter) + f'.{evtname}/DATA/FKtimes',dtype=str,ndmin=2)

    # read station coordinates
    stationfile = 'solver/M%02d' %(iter) + f'.{evtname}/DATA/STATIONS_FILTERED'
    statxt = np.loadtxt(stationfile,dtype=str,ndmin=2)
    nsta = statxt.shape[0]
    
    # synthetic parameters
    syndir = f'solver/{mdir}/{evtname}/OUTPUT_FILES/'
    name = statxt[0,1] + "." + statxt[0,0] + CCODE + "Z.sac"
    syn_z_hd = SACTrace.read(syndir + name,headonly=True)
    npt_syn = syn_z_hd.npts
    dt_syn = syn_z_hd.delta
    t0_syn = syn_z_hd.b
    sac_head = syn_z_hd.copy()

    # get time window 
    win_tb,win_te = pdict['TIME_WINDOW']
    npt2 = int((win_te + win_tb) / dt_syn)
    if npt2 // 2 * 2 != npt2:
        npt2 += 1

    # get filters
    Flist = pdict['GAUSS_F0']
    tshift = float(pdict['TSHIFT'])

    # if run_opt == 1, save synthetic data to SYN
    if run_opt == 1:
        outdir = syndir + '/SYN/'
        os.makedirs(outdir,exist_ok=True)
        nsta = statxt.shape[0]

        for i in range(myrank,nsta,nprocs):
            # read data
            name = statxt[i,1] + "." + statxt[i,0] + CCODE + "Z.sac"
            syn_z_tr = SACTrace.read(syndir + "/" + name )
            name = statxt[i,1] + "." + statxt[i,0] + CCODE + "R.sac"
            syn_r_tr = SACTrace.read(syndir + "/" + name )

            # tr for saving 
            tr = syn_z_tr.copy()

            for ib in range(len(Flist)):
                # compute rf
                bandname = 'F%3.1f' %(Flist[ib])
                dt = syn_z_tr.delta 
                rf = deconit(syn_r_tr.data,syn_z_tr.data,dt,tshift,Flist[ib],0)

                # save 
                tr.data = rf 
                tr.b = -tshift 
                tr.user4 = Flist[ib]
                name = statxt[i,1] + "." + statxt[i,0] + CCODE + f"R.{bandname}.rf.sac"
                tr.write(outdir + '/' + name)
        
        return 0

    # now loop every band to estimate the stf 
    nbands = len(Flist)
    for ib in range(nbands):
        bandname = 'F%3.1f' %(Flist[ib])

        # gauss filter
        gauss = gauss_filter(npt_syn,dt_syn,Flist[ib])
        
        # mkdir to save all files
        os.makedirs(syndir + f"{bandname}",exist_ok=True)

        # global synthetics
        glob_syn_rf = np.zeros((nsta,npt_syn),'f4')
        glob_obs_rf = np.zeros((nsta,1,npt_syn),'f4')
        glob_syn = np.zeros((nsta,2,npt_syn),'f4')
        glob_tp = np.zeros((nsta),'f4')
        
        # loop every station
        nsta = statxt.shape[0]
        for i in range(myrank,nsta,nprocs):
            # read synthetic data
            name = statxt[i,1] + "." + statxt[i,0] + CCODE + "Z.sac"
            syn_z_tr = SACTrace.read(syndir + "/" + name )
            name = statxt[i,1] + "." + statxt[i,0] + CCODE +"R.sac"
            syn_r_tr = SACTrace.read(syndir + "/" + name )

            # read obs rf
            name = statxt[i,1] + "." + statxt[i,0] + CCODE +"R.rf.sac.obs"
            rf_tr = SACTrace.read(f'fwat_data/{evtid}/' + name)
            rf_tr.write(syndir + f"{bandname}/" + name)
            
            # resample obs rf
            rfdata = interpolate_syn(rf_tr.data,rf_tr.b,rf_tr.delta,rf_tr.npts,
                                    -tshift,dt_syn,npt_syn)

            # compute synthetic rf
            rfi = deconit(syn_r_tr.data,syn_z_tr.data,dt_syn,tshift,Flist[ib],0)
            tr = rf_tr.copy(); tr.data = rfi 
            tr.delta = dt_syn
            tr.b = -tshift 
            tr.user4 = Flist[ib]
            name = statxt[i,1] + "." + statxt[i,0] + CCODE + "R.rf.sac.syn"
            tr.write(syndir + f"{bandname}/" + name)

            # save data to global
            glob_obs_rf[i,0,:] = rfdata 
            glob_syn_rf[i,:] = rfi 
            glob_syn[i,0,:] = apply_gaussian(syn_z_tr.data,gauss,dt_syn)
            glob_syn[i,1,:] = apply_gaussian(syn_r_tr.data,gauss,dt_syn)

            # compute ttp 
            tp = np.argmax(np.abs(glob_syn[i,1,:])) * dt_syn + t0_syn 
            if tp - t0_syn < tshift:
                print('time length before P should greater than RF_TSHIFT')
                exit(1)
            glob_tp[i] = tp

        # all reduce
        tmp = comm.allreduce(glob_obs_rf); glob_obs_rf = tmp
        tmp = comm.allreduce(glob_syn_rf); glob_syn_rf = tmp
        tmp = comm.allreduce(glob_syn); glob_syn = tmp
        tmp = comm.allreduce(glob_tp); glob_tp = tmp

        # get average amplitude
        avgamp = get_average_amplitude(glob_obs_rf,0)
        if myrank == 0:
            print("average amplitude of data gather %g" %(avgamp))
        
        # compute adjoint source 
        os.makedirs(syndir + f"{bandname}/OUTPUT_FILES",exist_ok=True)
        for i in range(myrank,nsta,nprocs):
            synz_bp = apply_gaussian(glob_syn[i,0,:],gauss,dt_syn)
            synr_bp = apply_gaussian(glob_syn[i,0,:],gauss,dt_syn)
            syn_rf = glob_syn_rf[i,:]
            obs_rf = glob_obs_rf[i,:]

            # compute normalize factor
            zrf = deconit(synz_bp,synz_bp,dt_syn,10.,Flist[ib],0)
            nfac = np.max(zrf)
            syn_rf /= nfac; obs_rf /= nfac 

            # align synthetic waveforms with p
            nb = int((glob_tp[i] - tshift - t0_syn) / dt_syn)
            synr_shift = np.zeros_like(synr_bp)
            synz_shift = np.zeros(synz_bp)
            synr_shift[:npt_syn-nb] = synr_bp[nb:]
            synz_shift[:npt_syn-nb] = synz_bp[nb:] 

            # adjoint source in time domain
            diff_data = syn_rf - obs_rf 
            e = npt_syn * dt_syn - tshift
            r_rev = synr_shift[::-1]
            z_rev = synz_shift[::-1]
            adj_r_tmp = deconit(diff_data,z_rev,dt,e,Flist[ib],1)
            tmp_n = np.convolve(-diff_data,r_rev,'full')
            tmp_d = np.convolve(z_rev,z_rev,'full')
            adj_z_tmp = deconit(tmp_n,tmp_d,dt_syn,e,Flist[ib],1)

            # windowed rf
            nstart = (-win_tb + tshift) // dt 
            nend = (win_te + tshift) // dt + 1
            fac = 1. - np.cos(np.pi * np.linspace(0,1,nend-nstart)) ** 10
            syn_rf[nstart:nend] *= fac 
            obs_rf[nstart:nend] *= fac
            adj_r = np.zeros((npt_syn)); adj_z = adj_r * 1. 
            adj_r[nb:nb+nend-nstart] *= adj_r_tmp[nstart:nend] * fac 
            adj_z[nb:nb+nend-nstart] *= adj_z_tmp[nstart:nend] * fac

            # write adjoint source to semd
            d = np.zeros((npt_syn,2))
            d[:,0] = np.arange(npt_syn) * dt_syn - tshift
            d[:,1] = adj_r * 1.
            name = statxt[i,1] + "." + statxt[i,0] + CCODE +"R.adj"
            np.savetxt(syndir + f"{bandname}/OUTPUT_FILES/" + name,d)
            d[:,1] = adj_z * 1.
            name = statxt[i,1] + "." + statxt[i,0] + CCODE + "Z.adj"
            np.savetxt(syndir + f"{bandname}/OUTPUT_FILES/" + name,d)

        # write measure_adj window
        if myrank == 0:
            f = open(syndir + f"{bandname}" + "/MEASUREMENT.WINDOWS","w")
            f.write("%d\n" %(nsta))

            for i in range(nsta):
                # time window
                tstart = -tshift 
                tend = -tshift + (npt_syn - 1) * dt_syn
                name = statxt[i,1] + "." + statxt[i,0] + CCODE +"R" + ".sac"
                f.write("%s\n" %(name + '.obs'))
                f.write("%s\n" %(name + '.syn'))
                f.write("1\n")
                f.write("%f %f\n" %(tstart,tend))

            f.close()
    
    comm.Barrier()

if __name__ == "__main__":
    main()
