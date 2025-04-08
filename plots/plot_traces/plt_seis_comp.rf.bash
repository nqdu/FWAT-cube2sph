#!/bin/bash
set -e 
source activate pygmt

if [ $# != 1 ];then 
    echo "./this Model(M00) "
    exit 1
fi 

gmt gmtset ANNOT_FONT_SIZE_PRIMARY 12p
gmt gmtset ANNOT_OFFSET_PRIMARY 0.1c
gmt gmtset LABEL_FONT_SIZE 14p 
gmt gmtset LABEL_OFFSET 0.15c
gmt gmtset HEADER_FONT_SIZE 16p
gmt gmtset HEADER_OFFSET -0.5c 
gmt gmtset TICK_LENGTH -0.2c

seisdir=seismograms_comp
bt=-0.
sc=1.0
M0=M00
M1=$1
lags=F1.0
mkdir -p $seisdir
#IFS=$'\n'
for iset in `seq 21 26`;do
cat ../../src_rec/sources_set$iset.dat |
while read line;do
  ievt=`echo $line |awk '{printf $1}'`
  fwddir0=../../solver/$M0.set$iset/$ievt/OUTPUT_FILES #solver/M00.set1/5/OUTPUT_FILES
  fwddir1=../../solver/$M1.set$iset/$ievt/OUTPUT_FILES #solver/M00.set1/5/OUTPUT_FILES
  #saclst kstnm knetwk kcmpnm t0 f ../data/$ievt/*BXZ.sac |awk '{print $2,$3,$4,$5}' >$seisdir/$ievt/FKtimes.dat
  phasetab=../../src_rec/FKtimes_${ievt}
  twb=5 #`head -1 $phasetab |awk '{print $4}'`
  twe=30 #`head -1 $phasetab |awk '{print $5}'`
  echo $ievt $twb $twe
  #=======================================================

  ls $fwddir0/dat.*BXR*.sac.$lags |sort -V >file1.lst
  ls $fwddir1/dat.*BXR*.sac.$lags |sort -V >file3.lst
  cat /dev/null >file2.lst
  for dfile in ` cat file1.lst`;do
     net=`echo $dfile |awk -F$fwddir0/ '{print $2}' |awk -F. '{print $2}'` 
     stnm=`echo $dfile |awk -F$fwddir0/ '{print $2}' |awk -F. '{print $3}'` 
     echo $fwddir0/syn.$net.$stnm.BXR.rf.sac.$lags  >>file2.lst
  done
  cat /dev/null >file4.lst
  for dfile in ` cat file3.lst`;do
    net=`echo $dfile |awk -F$fwddir1/ '{print $2}' |awk -F. '{print $2}'` 
    stnm=`echo $dfile |awk -F$fwddir1/ '{print $2}' |awk -F. '{print $3}'` 
    echo $fwddir1/syn.$net.$stnm.BXR.rf.sac.$lags  >>file4.lst
  done

  nsta=`cat file1.lst |wc -l`
  len=45
  npart=`echo $nsta $len |awk '{print $1/$2}' | awk '{printf("%d\n",$0+=$0<0?0:0.9999)}'`
  echo "nsta,len,npart=" $nsta $len $npart
  for i in `seq $npart`;do
    ib=`echo $i $len |awk '{print ($1-1)*$2+1}'`
    ie=`echo $i $len |awk '{print $1*$2}'`
    if [ ${ie} -gt $nsta ];then
        ie=$nsta
    fi
    #ymax=`echo $len |awk '{print $1}'`
    ymax=`echo "$nsta + 1"|bc`
    echo $ib $ie $ymax
#    ##########################################
    if true;then
    out=$seisdir/fsismo_${ievt}_BX.seg${ib}-${ie}.$M1.ps
    ### plot reconstruted waveform
    ##(a)
    gmt psbasemap -R-5/55/-1/$ymax -JX8/12 -Ba20f10/a1weSn  -K -P -X3 -Y16 >$out
    cat file1.lst |sed -n "${ib},${ie}p" | gmt pssac -J -R -Entb -M$sc -W0.5p -O -K >> $out 
    cat file2.lst |sed -n "${ib},${ie}p" | gmt pssac -J -R -Entb -M$sc -W0.5p,red -O -K >> $out 
    ty=`echo $ymax |awk '{print $1*1.05}'`
    gmt pstext -J -R -N -O -K >>$out <<EOF
25 $ty 16 0 0 CM (a) BXR before FWI
EOF
    cat /dev/null >temp
    cat /dev/null >temp1
    for fname in `cat file1.lst`;do
      netwk=`echo $fname |awk -F$fwddir0/ '{print $2}' |awk -F. '{print $2}'`
      stnm=`echo $fname |awk -F$fwddir0/ '{print $2}' |awk -F. '{print $3}'`
      echo $netwk $stnm >>temp 
      grep $stnm $phasetab >>temp1 
    done
    cat temp |sed -n "${ib},${ie}p" |awk '{printf"%f %f 8 0 0 CM %s.%s\n",-5,NR-1,$1,$2}' |gmt pstext -J -R -W255/0/0 -N -K -O >>$out 
    cat temp1 |sed -n "${ib},${ie}p" |awk '{printf"%f %f %f %f\n",$3+t-b,NR-1,0,0.2}' t=$bt b=$twb | gmt psxy -J -R -Sy0.2 -W0.5p,magenta  -N -O -K >>$out
    cat temp1 |sed -n "${ib},${ie}p" |awk '{printf"%f %f %f %f\n",$3+t+e,NR-1,0,0.2}' t=$bt e=$twe | gmt psxy -J -R -Sy0.2 -W0.5p,magenta  -N -O -K  >>$out

    ##(b)
    gmt psbasemap -R-5/55/-1/$ymax -JX8/12 -Ba20f10/a1weSn  -O -K -X9.5 >>$out
    cat file1.lst |sed -n "${ib},${ie}p" | gmt pssac -J -R -Entb -M$sc -W0.5p -O -K >> $out 
    cat file4.lst |sed -n "${ib},${ie}p" | gmt pssac -J -R -Entb -M$sc -W0.5p,red -O -K >> $out 
    ty=`echo $ymax |awk '{print $1*1.05}'`
    gmt pstext -J -R -N -O -K >>$out <<EOF
25 $ty 16 0 0 CM (b) BXR after FWI
EOF

    cat temp |sed -n "${ib},${ie}p" |awk '{printf"%f %f 8 0 0 CM %s.%s\n",-5,NR-1,$1,$2}' |gmt pstext -J -R -W255/0/0 -N -K -O >>$out
    #cat temp1 |sed -n "${ib},${ie}p" |awk '{printf"%f %f %f %f\n",$3+t-b,NR-1,0,0.2}' t=$bt b=$twb | gmt psxy -J -R -Sy0.2 -W0.5p,magenta  -N -O -K >>$out
    #cat temp1 |sed -n "${ib},${ie}p" |awk '{printf"%f %f %f %f\n",$3+t+e,NR-1,0,0.2}' t=$bt e=$twe | gmt psxy -J -R -Sy0.2 -W0.5p,magenta  -N -O -K  >>$out

    gmt psxy -J -R -T -O >>$out

    gmt psconvert $out -A -Tg 
    rm $out
    #rm temp temp1
    fi ### end if plot stf
    done
    #########################################
    #rm file?.lst
done

#rm *.lst
done # end iset
