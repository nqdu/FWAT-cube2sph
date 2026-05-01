#!/bin/bash
get_plot_latex_name() {
  local latex_name=""
  local param=$1
  local input_table=(vp vs rho vpv vph vsv vsh kappaa kappab G0 phi)
  local latex_table=("V_p" "V_s" "\rho" "\alpha_v" "\alpha_h" "\beta_v" "\beta_h" "\kappa_a" "\kappa_b" "G_0" "\phi")
  local units=("m/s" "m/s" "kg/m^3" "m/s" "m/s" "m/s" "m/s" "" "" "" "")
  local latex_dtable=("\delta \ln(V_p)" "\delta \ln(V_s)" "\delta \ln(\rho)" "\delta \ln(\alpha_v)" "\delta \ln(\alpha_h)" "\delta \ln(\beta_v)" "\delta \ln(\beta_h)" "\delta \kappa_a" "\delta \kappa_b" "\delta G_0" "\delta \phi")

  # check if "d" is in the $param 
  if [[ "$param" == d* ]]; then
    local base_param="${param:1}"  # remove the leading "d"
  else
    local base_param="$param"
  fi

  for i in "${!input_table[@]}"; do
    if [ "${input_table[$i]}" == "$base_param" ]; then
      if [[ "$param" == "$base_param" ]]; then
        latex_name="@[${latex_table[$i]}@[, ${units[$i]}"
      else
        latex_name="@[${latex_dtable[$i]}@[, \%"
      fi
      break
    fi
  done
  if [ "$latex_name" == "" ]; then
    latex_name="$1"
  fi
  echo $latex_name
}

specfem_dir="${HOME}/software/specfem3d-cube2sph/"
cube2sph_dir=$specfem_dir/utils/cube2sph/
NPROC=192

# horizontal slice
NSLICE_HORIZ=1
DEPTH_H=(225)  # add more if you want 
LON0_H=114; LON1_H=132;
LAT0_H=41; LAT1_H=46

# vertical slice
NSLICE_VERTI=1
MAX_DEP=400  # in km

# start points
LON0_V=(115 )
LAT0_V=(44)

#end points
LON1_V=(132 )
LAT1_V=(44)

# database_dir
DATABASE_DIR=../../DATABASES_MPI/
MODEL_DIR=../../optimize/
SOLVER_DIR=../../solver/

#parameters to interpolate
MODEL_SET="G0"
KERNEL_SET="dGsp dGcp"

# plot event kernels
PLOT_EVT_KERNEL=0
MDTYPE=dtti
KLTYPE=1

# interpolate options
INTP_KL=0 # if =1 interplate models
INTP_LS=0
run_indx=`seq 0 0`