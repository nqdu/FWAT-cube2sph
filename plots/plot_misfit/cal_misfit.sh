#!/bin/bash
set -e 
#        ! KEY: write misfit function values to file (two for each window)
#        ! Here are the 20 columns of the vector window_chi
#        !  1: MT-TT chi,    2: MT-dlnA chi,    3: XC-TT chi,    4: XC-dlnA chi
#        !  5: MT-TT meas,   6: MT-dlnA meas,   7: XC-TT meas,   8: XC-dlnA meas
#        !  9: MT-TT error, 10: MT-dlnA error, 11: XC-TT error, 12: XC-dlnA error
#        ! WINDOW     : 13: data power, 14: syn power, 15: (data-syn) power, 16: window duration
#        ! FULL RECORD: 17: data power, 18: syn power, 19: (data-syn) power, 20: record duration
#        ! Example of a reduced file: awk '{print $2,$3,$4,$5,$6,$31,$32}' window_chi > window_chi_sub
 
if [ "$#" -eq "3" ]; then 
  step=$3
elif [ "$#" -eq "2" ]; then
  step=00
else 
  echo "Usage ./cal_misfit M00 simu_type (step_index)"
  exit 1
fi 

mod=$1
simu_type=$2
sumf=0.0
sumn=0

fwat_file=../fwat_params/FWAT.PAR.$simu_type

if [ $simu_type == "rf"  ];then
  f0=(`cat $fwat_file |grep 'F0:' |awk -F: '{print $2}'`)
  NUM_FILTER=`echo ${#f0[@]}`
  period_used=""
  for ((i=0;i<$NUM_FILTER;i++));
  do
    strt=`printf "F%2.1f " ${f0[$i]}`
    period_used="$period_used""$strt"
  done
else 
  SHORT_P=(`cat $fwat_file |grep 'SHORT_P:' |awk -F: '{print $2}'`)
  LONG_P=(`cat $fwat_file |grep 'LONG_P:' |awk -F: '{print $2}'`)
  NUM_FILTER=`echo ${#SHORT_P[@]}`
  period_used=""
  for ((i=0;i<$NUM_FILTER;i++));
  do
    strt=`printf "T%03d_T%03d " ${SHORT_P[$i]} ${LONG_P[$i]}`
    period_used="$period_used""$strt"
  done
fi 

for band in $period_used;
do
  if [ "$#" -eq "3" ];then
    cat ../src_rec/sources.dat.ls.$simu_type |awk '{printf"%s\n",$1}' >eid.dat
    files=../misfits/${mod}*.ls_${band}_${simu_type}_window_chi
    if (( $(echo "$step == 00" |bc -l) )); then
      files=`ls ../misfits/${mod}.*_${band}_${simu_type}_window_chi |grep -v '.ls'`
    fi

    for eid in `cat eid.dat`;do
      #eid=`echo $line |awk '{printf"%s.%s",$2,$1}'`
      misf=`cat $files |awk -v a=$eid '$1 == a {print $0}' | awk 'BEGIN{sum=0;} {if($29!=0) {sum=sum+$29}} END{print sum}' `
      nmeas=`cat $files |awk -v a=$eid '$1 == a {print $0}'| awk 'BEGIN{n=0;} {if($29!=0) {n=n+1}} END{print n}' `
      sumf=`echo "$sumf $misf" |awk '{printf "%f", $1+$2}'`
      sumn=`echo "$sumn + $nmeas" |bc -l`
    done
    \rm eid.dat
  else
    files=`ls ../misfits/${mod}.*_${band}_${simu_type}_window_chi |grep -v '.ls'`
    #echo $files
    misf=`cat $files |awk 'BEGIN{sum=0;} {if($29!=0) {sum=sum+$29}} END{print sum}' `
    nmeas=`cat $files |awk 'BEGIN{n=0;} {if($29!=0) {n=n+1}} END{print n}' `
    sumf=`echo "$sumf $misf" |awk '{printf "%f", $1+$2}'`
    sumn=`echo $sumn + $nmeas |bc -l`
  fi 
done 

echo $sumf $sumn

# for band in $period_used ;do
#   # calculate misfit of step 0.00
#   cat ../src_rec/sources_ls.dat.$simu_type |awk '{printf"%s\n",$1}' >eid.dat
#   files=../misfits/${mod}.set*_${band}_${simu_type}_window_chi
#   if [ "$#" -eq "3" ];then 
#     files=../misfits/${mod}_step${step}.ls_${band}_${simu_type}_window_chi
#   fi
#   #echo $files
#   for eid in `cat eid.dat`;do
#     #eid=`echo $line |awk '{printf"%s.%s",$2,$1}'`
#     misf=`cat $files |awk '{if($1==a) print $0}' a=$eid  |awk 'BEGIN{sum=0;} {if($29!=0) {sum=sum+$29}} END{print sum}' `
#     nmeas=`cat $files |awk '{if($1==a) print $0}' a=$eid |awk 'BEGIN{n=0;} {if($29!=0) {n=n+1}} END{print n}' `
#     sumf=`echo $sumf + $misf |bc -l`
#     sumn=`echo $sumn + $nmeas |bc -l`
#   done
# done 
