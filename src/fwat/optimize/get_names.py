from fwat import FwatModel
from fwat.const import PARAM_FILE


def run(argv):
    if len(argv) !=1 :
        print("usage: fwat model-name type (kernel/direc/model)")
        exit(1)
    
    type_ = argv[0]
    assert type_ in ['grad','direc','model']

    # init
    m = FwatModel(PARAM_FILE)

    if type_ == "grad":
        print(' '.join(m.grad_names()))
    elif type_ == "direc":
        print(' '.join(m.direc_names()))
    else:
        print(' '.join(m.model_names()))
