import yaml 
import sys

def get_nested_value(d, path):
    keys = path.split('/')
    for key in keys:
        d = d.get(key, {})  # Default to an empty dict if key not found
    return d if isinstance(d, dict) else d  # Return the value if it's not a dictionary



def main():
    if len(sys.argv) != 3 and len(sys.argv) != 2:
        print("Usage: get_param.py paramloc file=[fwat_params/FWAT.PAR.yaml] ")
        print("example: get_param.py measure/tele/CH_CODE ")
        exit(1)
    
    paramloc = sys.argv[1]
    if len(sys.argv) == 3:
        filename = sys.argv[2]
    else:
        filename = "fwat_params/FWAT.PAR.yaml"
    
    # load yaml 
    with open(filename,"r") as f:
        pdict = yaml.safe_load(f)
    # get param
    value = get_nested_value(pdict,paramloc)
    print(value)

if __name__ == "__main__":
    main()