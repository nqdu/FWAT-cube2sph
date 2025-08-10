import yaml 
from functools import reduce
from fwat.const import PARAM_FILE

def set_nested(d, path, value):
    keys = path.split('/')
    reduce(dict.__getitem__, keys[:-1], d)[keys[-1]] = value


def get_nested_value(d, path):
    keys = path.split('/')
    for key in keys:
        d = d.get(key, {})  # Default to an empty dict if key not found
    return d if isinstance(d, dict) else d  # Return the value if it's not a dictionary


def help_function():
    print("Usage fwat-utils [module_name] [args]\n")
    print()
    print("fwat clean MODEL evtid (deepclean) ")
    print("\tclean temporary diretories")
    print("\texample: fwat clean M00 XZ.FAF")

    print()
    print("fwat-utils setparam paramloc value file=[fwat_params/fwat.yaml] ")
    print("\tset parameters with value in the file")
    print("\texample: fwat-param set measure/tele/CH_CODE 'BH' ")

    print()
    print("fwat-utils getparam paramloc file=[fwat_params/fwat.yaml] ")
    print("\tget parameters in the file")
    print("\texample: fwat-param get measure/tele/CH_CODE ")

def set_param(argv):
    if len(argv) !=3 and len(argv) != 2:
        print("Usage: fwat-param set paramloc value [file=fwat_params/fwat.yaml]")
        exit(1)
    
    paramloc = argv[0]
    value = argv[1]
    if len(argv) == 3:
        filename = argv[2]
    else:
        filename = PARAM_FILE

    # load yaml 
    with open(filename,"r") as f:
        pdict = yaml.safe_load(f)

    # make sure value is str/int/float
    try:
        out = int(value)
    except:
        try:
            out = float(value)
        except:
            out = value

    # change parameters
    set_nested(pdict,paramloc,out)

    # write to file
    with open(filename,"w") as f:
        yaml.safe_dump(pdict,f)

def get_param(argv):
    if len(argv) !=1 and len(argv) != 2:
        print("Usage: fwat-param get paramloc [file=fwat_params/fwat.yaml]")
        exit(1)

    paramloc = argv[0]
    if len(argv) == 2:
        filename = argv[1]
    else:
        filename = PARAM_FILE
    
    # load yaml 
    with open(filename,"r") as f:
        pdict = yaml.safe_load(f)
    
    # get param
    value = get_nested_value(pdict,paramloc)
    print(value)

def main():
    import sys
    if len(sys.argv) < 2:
        help_function()
        sys.exit(1)

    # get cmd
    cmd = sys.argv[1]
    args = sys.argv[2:]

    # call modules
    if cmd == "setparam":
        set_param(args)
    elif cmd == "getparam":
        get_param(args)
    elif cmd == "clean":
        from fwat import clean
        clean.run(args)
    