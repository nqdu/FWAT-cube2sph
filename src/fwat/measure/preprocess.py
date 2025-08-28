from mpi4py import MPI

def run(argv):

    if len(argv) != 4:
        print("Usage: fwat-main measure measure_type iter evtid run_opt")
        print("mpi should be enabled")
        exit(1)

    # get input parameter
    mtype = argv[0]
    iter0 = int(argv[1])
    evtid = argv[2]
    run_opt = int(argv[3])

    # init 
    if mtype == "tele":
        from .tele_preproc import Tele_PreOP
        op = Tele_PreOP(mtype,iter0,evtid,run_opt)

    elif mtype == "sks":
        from .sks_preproc import SKS_PreOP
        op = SKS_PreOP(mtype,iter0,evtid,run_opt)
        
    elif mtype == "noise":
        # from .noise_preproc import Noise_PreOP
        from .noise_mc_preproc import NoiseMC_PreOP
        op = NoiseMC_PreOP(mtype,iter0,evtid,run_opt)
    
    elif mtype == "rf":
        from .rf_preproc import RF_PreOP
        op = RF_PreOP(mtype,iter0,evtid,run_opt)
    else:
        print(f"{mtype} is not implemented!")
        exit(1)

    # run
    op.execute()

    # finalize
    MPI.Finalize()

def main():
    import sys 
    if len(sys.argv) != 5:
        print("Usage: python run_preprocess.py measure_type iter evtid run_opt")
        exit(1)

    # run 
    run(sys.argv[1:])

if __name__ == "__main__":
    main()