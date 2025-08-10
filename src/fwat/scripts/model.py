import sys 


def help_function():
    print("Usage fwat-model [module_name] [args]")
    
    print()
    print("fwat-model name type(grad/direc/model) ")
    print("\tprint names of current model parameterization, for gradient/search direction/model")
    print("\texample: fwat-model name grad")

    print()
    print("fwat-model visual  MODEL_DIR model_type kernel_type")
    print("\twrite out models for visualization")
    print("\tmodel_type and kernel_type can be find in class FwatModel")
    print("\texample: fwat-model visual M06 dtti 1")

    print()
    print("fwat-model kernel MODEL_DIR KERNEL_DIR model_type kernel_type")
    print("\trepack h5 kernel to binary")
    print("\tmpi should be enabled")
    print("\tusage: mpirun -np 4 fwat-model kernel M06 GRADIENT mdtype kltype")
    exit(1)

def main():
    if len(sys.argv) < 2:
        help_function()
        sys.exit(1)
    
    # get cmd
    cmd = sys.argv[1]
    args = sys.argv[2:]

    # run command
    if cmd == "name":
        from fwat.optimize import get_names
        get_names.run(args)
    elif cmd == "visual":
        from fwat.optimize import write_visual_models
        write_visual_models.run(args)
        pass 
    elif cmd == "kernel":
        pass
    else:
        print(f"{cmd} is not a function in fwat!")
        sys.exit(1)