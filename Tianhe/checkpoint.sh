#!/bin/bash

function CheckDebug
{
  local usr=${USER}
  local partition=$1

  if [[ $partition == "debug"  ]]; then
  	local njobs=`squeue -u $usr | grep debug| wc -l`
    while [ $njobs -eq 1  ]
    do
      echo "in debug partition, only one task could be submitted"
      squeue -u $usr | grep debug
      sleep 30
      njobs=`squeue -u $usr | grep debug | wc -l`
    done
  fi
}

function CheckForward
{
    if [ $# -ne 3 ]; then 
        echo "Usage: ./CheckForward model setb sete "
        echo "narg: " $#
        exit 1
    fi 
    # compute how many tasks
    local model=$1
    local setb=$2
    local sete=$3
    let local njobs=$sete-$setb+1

    # check for allocating resources
    echo "Waiting for forward simulation ... `date`"
    local noutfiles=0
    while [ $noutfiles -ne $njobs ]
    do 
        local ii 
        for ii in `seq $setb $sete`;
        do
            if [  -f output_fwat0_log_${model}.set${ii}.txt ] \
                && [ `grep 'Finished simulations here!!!' \
                output_fwat0_log_${model}.set${ii}.txt|wc -l` -eq 1 ]; 
            then
                noutfiles=`echo "$noutfiles+1" | bc`
            fi
        done 
        echo "$noutfiles jobs finished"
        if [ $noutfiles -eq $njobs ]; then
            break;
        else 
            noutfiles=0
            sleep 30 
        fi
    done 
}

function CheckFWDADJFinished 
{
    if [ $# -ne 3 ]; then 
        echo "Usage: ./CheckFWDADJFinished  model setb sete "
        echo "narg: " $#
        exit 1
    fi 
    # compute how many tasks
    local model=$1
    local setb=$2
    local sete=$3
    let local njobs=$sete-$setb+1

    # check for allocating resources
    echo "Waiting for forward/adjoint simulation ... `date`"
    local noutfiles=0
    while [ $noutfiles -ne $njobs ]
    do 
        local ii 
        for ii in `seq $setb $sete`;
        do
            if [  -f output_fwat1_log_${model}.set${ii}.txt ] \
                && [ `grep 'Finished simulations here!!!' \
                output_fwat1_log_${model}.set${ii}.txt|wc -l` -eq 1 ]; 
            then
                noutfiles=`echo "$noutfiles+1" | bc`
            fi
        done 
        echo "$noutfiles jobs finished"
        if [ $noutfiles -eq $njobs ]; then
            break;
        else 
            noutfiles=0
            sleep 30 
        fi
    done 
}

function CheckPostFinished
{
    if [ $# -ne 1 ]; then 
        echo "Usage: ./CheckPostFinished  model "
        echo "narg: " $#
        exit 1
    fi
    # compute how many steps
    local model=$1
    local temp=(`grep STEP_LENS fwat_params/FWAT.PAR |awk -F: '{print $2}'`)
    local nsteps=${#temp[@]}

    # check for allocating resources
    echo "Waiting for Post-processing ... `date`"
    local NPROC=`grep ^"NPROC" DATA/Par_file | cut -d'=' -f2`
    local noutfiles=0
    while [ $nsteps -ne $noutfiles ]
    do 
        local ii
        for ((ii=0;ii<$nsteps;ii++));do 
            if [ -d ./optimize/MODEL_${model}_step${temp[$ii]}  ] && \
                [ `ls ./optimize/MODEL_${model}_step${temp[$ii]} | grep external_mesh.bin | wc -l` -eq $NPROC ]; then 
                noutfiles=`echo "$noutfiles+1" | bc`
            fi 
        done 
        echo "$noutfiles jobs finished"
        if [ $noutfiles -eq $nsteps ];then 
            break;
        else 
            noutfiles=0
            sleep 30
        fi
    done 
}

function CheckLineSearch
{
    if [ $# -ne 1 ]; then 
        echo "Usage: ./CheckLineSearch model "
        echo "narg: " $#
        exit
    fi

    # compute how many steps
    local model=$1
    local temp=(`grep STEP_LENS fwat_params/FWAT.PAR |awk -F: '{print $2}'`)
    local nsteps=${#temp[@]}

    # check for allocating resources
    echo "Waiting for linear search ... `date`"
    local noutfiles=0
    while [ $noutfiles -ne $nsteps ]
    do 
        local ii
        for ((ii=0;ii<$nsteps;ii++));
        do
            if [  -f output_fwat3_log_${model}_step${temp[$ii]}.ls.txt ] \
                && [ `grep 'Finished simulations here!!!' \
                output_fwat3_log_${model}_step${temp[$ii]}.ls.txt|wc -l` -eq 1 ]; 
            then
                noutfiles=`echo "$noutfiles+1" | bc`
            fi
        done 
        echo "$noutfiles jobs finished"
        if [ $noutfiles -eq $nsteps ]; then
            break;
        else 
            noutfiles=0
            sleep 30 
        fi
    done 
}
