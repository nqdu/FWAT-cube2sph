# binary/python file location
fksem='/home/l/liuqy/nqdu//specfem3d-cube2sph/'    # solver location
FWATLIB=`pwd`/fwatlib    # fwatlib location
FWATPARAM=`pwd`/fwat_params


# STOP HERE ###################
MEASURE_LIB=$FWATLIB/measure 
OPT_LIB=$FWATLIB/optimize

GET_GRAD_NAME(){
    local PFILE="$FWATPARAM/FWAT.PAR.yaml"
   local grad_list=`cd $OPT_LIB; python -c "from tools import *; mod = get_gradname_list('$PFILE'); print(' '.join(mod))"`
    echo $grad_list
}

GET_DIREC_NAME(){
    local PFILE="$FWATPARAM/FWAT.PAR.yaml"
   local grad_list=`cd $OPT_LIB; python -c "from tools import *; mod = get_direc_name_list('$PFILE'); print(' '.join(mod))"`
    echo $grad_list
}

GET_MODEL_NAME(){
    local PFILE="$FWATPARAM/FWAT.PAR.yaml"
   local grad_list=`cd $OPT_LIB; python -c "from tools import *; mod = get_model_name_list('$PFILE'); print(' '.join(mod))"`
    echo $grad_list
}