import sys 

def help_function():
    print("Usage fwat-main [module_name] [args]\n")


    print("fwat-main measure measure_type iter evtid run_opt")
    print("\trun measurements to generate adjoint source")
    print("\topt =1 forward, = 3 adjoint")
    print("\tmpi should be enabled")
    print("\texample: mpirun -np 31 fwat measure tele 0 XZ.FAF 3")

    print()
    print("fwat-main pack outname *list")
    print("\tpacking all seismograms to hdf5 file")
    print("\texample: fwat pack seismograms.h5 OUTPUT_FILES/all_seismograms.bin")

    print()
    print("fwat-main misfit MODEL simu_type")
    print("\tcompute misfit")
    print("\texample: fwat pack M00 noise")

    print()
    print("fwat-main bin2h5 WORK_DIR name nprocs combine=1,decompose=0")
    print("\tpack all kernels to hdf5 file")
    print("\texample: fwat bin2h5 solver/M00/P21/GRADIENT beta_kernel 160 1")

    print()
    print("fwat-main direc")
    print("\tget search direction, current iteration id will be read from fwat.yaml")
    print("\tmpi should be enabled")
    print("\texample:mpirun -np 4 fwat direc")

    print()
    print("fwat-main update MODEL OUT_DIR")
    print("\tupdate current model by using search direction")
    print("\tmpi should be enabled")
    print("\texample: mpirun -np 4 fwat update M06 OUT_DIR")

    print()
    print("fwat-main linesearch Model fcost_model fcost_ls")
    print("\tline search using Wolfe condition")
    print("\tmpi should be enabled")
    print("\texample: mpirun -np 4  fwat linesearch M06 10. 9.")

    print()
    print("fwat-main sum_kernel evtfile iter_cur MODEL ")
    print("\tsum all kernels from current iteration write to MODEL")
    print("\tmpi should be enabled")
    print("\texample: mpirun -np 4 fwat sum_kernel src_rec/sources.dat 0 M00.ls")

    print()
    print("fwat-main prepare forward/adjoint meatype iter evtid run_opt")
    print("\tprepare forward/adjoint files")
    print("\texample: fwat-main prepare forward tele 0 XZ.FAF 3")

def main():
    if len(sys.argv) < 2:
        help_function()
        sys.exit(1)
    
    # get cmd
    cmd = sys.argv[1]
    args = sys.argv[2:]

    # run command
    if cmd == "measure":
        from fwat.measure import preprocess  
        preprocess.run(args)
    elif cmd == "pack":
        from fwat.measure import pack_seismogram
        pack_seismogram.run(args)
    elif cmd == "misfit":
        from fwat.measure import cal_misfit
        cal_misfit.run(args)
    elif cmd == "bin2h5":
        from fwat.optimize import bin2h5
        bin2h5.run(args)
    elif cmd == "direc":
        from fwat.optimize import search_direction
        search_direction.run()
    elif cmd == "update":
        from fwat.optimize import model_update
        model_update.run(args)
    elif cmd == "linesearch":
        from fwat.optimize import std_linesearch
        std_linesearch.run(args)
    elif cmd == "sum_kernel":
        from fwat.optimize import sum_kernel
        sum_kernel.run(args)
    elif cmd == "prepare":
        from fwat.submit.submit import prepare_adj,prepare_fwd
        assert (args[0] in ['forward','adjoint'])
        if args[0] == "forward":
            prepare_fwd(args[1:])
        else:
            prepare_adj(args[1:])
    else:
        print(f"{cmd} is not a function in fwat!")
        sys.exit(1)