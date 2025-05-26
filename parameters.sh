# binary/python file location
fksem='/home/l/liuqy/nqdu//specfem3d-cube2sph/'    # solver location
FWATLIB=`pwd`/fwatlib    # fwatlib location
FWATPARAM=`pwd`/fwat_params

# nprocs for preprocessing
NPROC_PRE=4

# STOP HERE ###################


MEASURE_LIB=$FWATLIB/measure 
OPT_LIB=$FWATLIB/optimize

GET_GRAD_NAME(){
    local PFILE="$FWATPARAM/FWAT.PAR.yaml"
    local grad_list=`cd $OPT_LIB; python -c "from tools import FwatModel; m = FwatModel('$PFILE'); print(' '.join(m.get_grad_names()))"`
    echo $grad_list
}

GET_DIREC_NAME(){
    local PFILE="$FWATPARAM/FWAT.PAR.yaml"
    local grad_list=`cd $OPT_LIB; python -c "from tools import FwatModel; m = FwatModel('$PFILE'); print(' '.join(m.get_direc_names()))"`
    echo $grad_list
}

GET_MODEL_NAME(){
    local PFILE="$FWATPARAM/FWAT.PAR.yaml"
    local grad_list=`cd $OPT_LIB; python -c "from tools import FwatModel; m = FwatModel('$PFILE'); print(' '.join(m.get_model_names()))"`
    echo $grad_list
}