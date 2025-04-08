import yaml 
import sys 

from functools import reduce

def set_nested(d, path, value):
    keys = path.split('/')
    reduce(dict.__getitem__, keys[:-1], d)[keys[-1]] = value

def main():
    if len(sys.argv) != 4 and len(sys.argv) != 3:
        print("Usage: set_param.py paramloc value file=[fwat_params/FWAT.PAR.yaml] ")
        print("example: set_param.py measure/tele/CH_CODE 'BH' ")
        exit(1)
    
    paramloc = sys.argv[1]
    value = sys.argv[2]
    if len(sys.argv) == 4:
        filename = sys.argv[3]
    else:
        filename = "fwat_params/FWAT.PAR.yaml"
    
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

if __name__ == "__main__":
    main()