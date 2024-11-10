#!/bin/bash
module load gmt-6.0.0
set -e

if [ $# != 1 ];then 
    echo "./this Model(M00) "
    exit 1
fi 

mydir=../../solver/
proj=-JX8c/11c
sc=0.5c/-1
M0=M00
M1=$1
seisdir=seismograms_comp
mkdir -p $seisdir

cat ../../src_rec/sources.dat.tele |while read line;
do 
  # fktimes
  twb=5 #`head -1 $phasetab |awk '{print $4}'`
  twe=45 #`head -1 $phasetab |awk '{print $5}'`
  evid=`echo $line |awk '{print $1}'`

  # generate file
  file1=data.in.Z
  file2=syn.in.Z.init
  file3=syn.in.Z.$M1
  file4=data.in.R
  file5=syn.in.R.init
  file6=syn.in.R.$M1

  :>temp
  :> $file1 
  :> $file2 
  :> $file3 
  :> $file4 
  :> $file5 
  :> $file6

  cat ../../src_rec/STATIONS_${evid}_globe |while read line;
  do 
    net=`echo $line |awk  '{print $2}'` 
    stnm=`echo $line |awk  '{print $1}'`
    echo  $mydir/$M1/${evid}/OUTPUT_FILES/T005_T050/$net.$stnm.BXZ.sac.syn >> $file3
    echo  $mydir/$M1/${evid}/OUTPUT_FILES/T005_T050/$net.$stnm.BXR.sac.syn >> $file6
    echo  $mydir/$M0/${evid}/OUTPUT_FILES/T005_T050/$net.$stnm.BXZ.sac.syn >> $file2
    echo  $mydir/$M0/${evid}/OUTPUT_FILES/T005_T050/$net.$stnm.BXR.sac.syn >> $file5
    echo  $mydir/$M0/${evid}/OUTPUT_FILES/T005_T050/$net.$stnm.BXZ.sac.obs >> $file1
    echo  $mydir/$M0/${evid}/OUTPUT_FILES/T005_T050/$net.$stnm.BXR.sac.obs >> $file4
    saclst t0 f $mydir/$M0/${evid}/OUTPUT_FILES/T005_T050/$net.$stnm.BXR.sac.obs |awk '{print $2}' >> temp
  done  
  ymax=`wc -l $file1 |awk '{print $1}'`
  ty=`echo $ymax |awk '{print $1*1.05}'`
  bounds=-R0/130/-1/$ymax
  echo $evid $bounds

  gmt begin $seisdir/inv.$evid jpg
  gmt basemap $bounds $proj -Bxaf -Byaf -BWSen+t"BXZ before FWI" -X15c -Y20c
  cat $file1 | gmt sac -En -M$sc -W0.5p  
  cat $file2 | gmt sac -En -M$sc -W0.5p,red
  cat temp | awk '{printf"%f %f %f %f\n",$1-b,NR-1,0,0.2}' b=$twb  | gmt plot  -Sy0.2 -W0.5p,magenta
  cat temp | awk '{printf"%f %f %f %f\n",$1+e,NR-1,0,0.2}' e=$twe  | gmt plot  -Sy0.2 -W0.5p,magenta
#   gmt legend  -DjBR+w2.8c+o0.1c/0.1c -F+p1p << EOF
# S 0.25c - 0.5c - 0.25p 0.8c fwat
# S 0.25c - 0.5c - 0.25p,red,-- 0.8c my
# EOF

  gmt basemap $bounds $proj -Bxaf -Byaf -X9c -BWSen+t"BXZ after FWI"
  cat $file1 | gmt sac -En -M$sc -W0.5p
  cat $file3 | gmt sac -En -M$sc -W0.5p,red
  cat temp | awk '{printf"%f %f %f %f\n",$1-b,NR-1,0,0.2}' b=$twb  | gmt plot  -Sy0.2 -W0.5p,magenta
  cat temp | awk '{printf"%f %f %f %f\n",$1+e,NR-1,0,0.2}' e=$twe  | gmt plot  -Sy0.2 -W0.5p,magenta

  gmt basemap $bounds $proj -Bxaf+l"Time (s)" -Byaf -BWSen+t"BXR before FWI" -X-9c -Y-12c
  cat $file4 | gmt sac -En -M$sc -W0.5p  
  cat $file5 | gmt sac -En -M$sc -W0.5p,red
  cat temp | awk '{printf"%f %f %f %f\n",$1-b,NR-1,0,0.2}' b=$twb  | gmt plot  -Sy0.2 -W0.5p,magenta
  cat temp | awk '{printf"%f %f %f %f\n",$1+e,NR-1,0,0.2}' e=$twe  | gmt plot  -Sy0.2 -W0.5p,magenta
#   gmt legend  -DjBR+w2.8c+o0.1c/0.1c -F+p1p << EOF
# S 0.25c - 0.5c - 0.25p 0.8c fwat
# S 0.25c - 0.5c - 0.25p,red,-- 0.8c my
# EOF

  gmt basemap $bounds $proj -Bxaf+l"Time (s)" -Byaf -X9c -BWSen+t"BXR after FWI"
  cat $file4 | gmt sac -En -M$sc -W0.5p
  cat $file6 | gmt sac -En -M$sc -W0.5p,red
  cat temp | awk '{printf"%f %f %f %f\n",$1-b,NR-1,0,0.2}' b=$twb  | gmt plot  -Sy0.2 -W0.5p,magenta
  cat temp | awk '{printf"%f %f %f %f\n",$1+e,NR-1,0,0.2}' e=$twe  | gmt plot  -Sy0.2 -W0.5p,magenta

  gmt end 

  \rm -f $file1 $file2 $file3 $file4 $file5 $file6 
  #exit 
done