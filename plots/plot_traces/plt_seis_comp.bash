#!/bin/bash
set -e 
module load gmt-6.0.0

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
bt=-10.
sc=0.5/-1
M0=M00
M1=$1
mkdir -p $seisdir
#IFS=$'\n'
for iset in `seq 1 8`;do
band=T005_T050
cat ../../src_rec/sources.dat.tele |
while read line;do
  ievt=`echo $line |awk '{printf $1}'`
  fwddir0=../../solver/$M0.$ievt/OUTPUT_FILES/$band #solver/M00.set1/5/OUTPUT_FILES
  fwddir1=../../solver/$M1.$ievt/OUTPUT_FILES/$band #solver/M00.set1/5/OUTPUT_FILES
  #saclst kstnm knetwk kcmpnm t0 f ../data/$ievt/*BXZ.sac |awk '{print $2,$3,$4,$5}' >$seisdir/$ievt/FKtimes.dat
  phasetab=../../src_rec/FKtimes_${ievt}
  twb=10 #`head -1 $phasetab |awk '{print $4}'`
  twe=45 #`head -1 $phasetab |awk '{print $5}'`
  echo $ievt $twb $twe
  #=======================================================

  ls $fwddir0/*.BXZ*.sac.obs >file1.lst
  ls $fwddir0/*.BXR*.sac.obs >file3.lst
  ls $fwddir1/*.BXZ*.sac.obs >file5.lst
  ls $fwddir1/*.BXR*.sac.obs >file7.lst
  cat /dev/null >file2.lst
  cat /dev/null >file4.lst
  for dfile in ` cat file1.lst`;do
     net=`echo $dfile |awk -F$fwddir0/ '{print $2}' |awk -F. '{print $2}'` 
     stnm=`echo $dfile |awk -F$fwddir0/ '{print $2}' |awk -F. '{print $3}'` 
     echo $fwddir0/$net.$stnm.BXZ.sac.syn  >>file2.lst
     echo $fwddir0/$net.$stnm.BXR.sac.syn  >>file4.lst
  done
  cat /dev/null >file6.lst
  cat /dev/null >file8.lst
  for dfile in ` cat file5.lst`;do
     net=`echo $dfile |awk -F$fwddir1/ '{print $2}' |awk -F. '{print $2}'` 
     stnm=`echo $dfile |awk -F$fwddir1/ '{print $2}' |awk -F. '{print $3}'` 
     echo $fwddir1/$net.$stnm.BXZ.sac.syn  >>file6.lst
     echo $fwddir1/$net.$stnm.BXR.sac.syn  >>file8.lst
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
    ymax=$nsta
    echo $ib $ie $ymax
#    ##########################################
    if true;then
    out=$seisdir/fsismo_${ievt}_BX.seg${ib}-${ie}.$M1.ps
    ### plot reconstruted waveform
    ##(a)
    gmt psbasemap -R10/110/-1/$ymax -JX8/11 -Ba20f10/a1weSn  -K -P -X3 -Y16 >$out
    cat file1.lst |sed -n "${ib},${ie}p" | gmt pssac -J -R -En -M$sc -W0.5p -O -K >> $out 
    cat file2.lst |sed -n "${ib},${ie}p" | gmt pssac -J -R -En -M$sc -W0.5p,red -O -K >> $out 
    ty=`echo $ymax |awk '{print $1*1.05}'`
    gmt pstext -J -R -N -O -K >>$out <<EOF
55 $ty 16 0 0 CM (a) BXZ before FWI
EOF
    cat /dev/null >temp
    cat /dev/null >temp1
    for fname in `cat file1.lst`;do
      netwk=`echo $fname |awk -F$fwddir0/ '{print $2}' |awk -F. '{print $2}'`
      stnm=`echo $fname |awk -F$fwddir0/ '{print $2}' |awk -F. '{print $3}'`
      echo $netwk $stnm >>temp 
      grep $stnm $phasetab >>temp1 
    done
    cat temp |sed -n "${ib},${ie}p" |awk '{printf"%f %f 8 0 0 CM %s.%s\n",10,NR-1,$1,$2}' |gmt pstext -J -R -W255/0/0 -N -K -O >>$out 
    cat temp1 |sed -n "${ib},${ie}p" |awk '{printf"%f %f %f %f\n",$3+t-b,NR-1,0,0.2}' t=$bt b=$twb | gmt psxy -J -R -Sy0.2 -W0.5p,magenta  -N -O -K >>$out
    cat temp1 |sed -n "${ib},${ie}p" |awk '{printf"%f %f %f %f\n",$3+t+e,NR-1,0,0.2}' t=$bt e=$twe | gmt psxy -J -R -Sy0.2 -W0.5p,magenta  -N -O -K  >>$out

    ##(b)
    gmt psbasemap -R10/110/-1/$ymax -JX8/11 -Ba20f10/a1weSn  -O -K -X8.5 >>$out
    cat file5.lst |sed -n "${ib},${ie}p" | gmt pssac -J -R -En -M$sc -W0.5p -O -K >> $out 
    cat file6.lst |sed -n "${ib},${ie}p" | gmt pssac -J -R -En -M$sc -W0.5p,red -O -K >> $out 
    ty=`echo $ymax |awk '{print $1*1.05}'`
    gmt pstext -J -R -N -O -K >>$out <<EOF
55 $ty 16 0 0 CM (b) BXZ after FWI
EOF
    cat temp1 |sed -n "${ib},${ie}p" |awk '{printf"%f %f %f %f\n",$3+t-b,NR-1,0,0.2}' t=$bt b=$twb | gmt psxy -J -R -Sy0.2 -W0.5p,magenta  -N -O -K >>$out
    cat temp1 |sed -n "${ib},${ie}p" |awk '{printf"%f %f %f %f\n",$3+t+e,NR-1,0,0.2}' t=$bt e=$twe | gmt psxy -J -R -Sy0.2 -W0.5p,magenta  -N -O -K  >>$out
    
    #wb=`cat temp1 |sed -n "${ib},${ie}p" |awk '{printf"%f\n",$3+t-b}' t=$bt b=$twb`
    #we=`cat temp1 |sed -n "${ib},${ie}p" |awk '{printf"%f\n",$3+t+e}' t=$bt b=$twe`
    #echo $wb 
    #echo $we
    
    ##(c)
    gmt psbasemap -R10/110/-1/$ymax -JX8/11 -Ba20f10:"Time (s)":/a1weSn  -K -O -X-8.5 -Y-13 >>$out
    cat file3.lst |sed -n "${ib},${ie}p" | gmt pssac -J -R -En -M$sc -W0.5p -O -K >> $out 
    cat file4.lst |sed -n "${ib},${ie}p" | gmt pssac -J -R -En -M$sc -W0.5p,red -O -K >> $out 
    ty=`echo $ymax |awk '{print $1*1.05}'`
    gmt pstext -J -R -N -O -K >>$out <<EOF
55 $ty 16 0 0 CM (c) BXR before FWI
EOF
    cat temp |sed -n "${ib},${ie}p" |awk '{printf"%f %f 8 0 0 CM %s.%s\n",10,NR-1,$1,$2}' |gmt pstext -J -R -W255/0/0 -N -K -O >>$out 
    cat temp1 |sed -n "${ib},${ie}p" |awk '{printf"%f %f %f %f\n",$3+t-b,NR-1,0,0.2}' t=$bt b=$twb | gmt psxy -J -R -Sy0.2 -W0.5p,magenta  -N -O -K >>$out
    cat temp1 |sed -n "${ib},${ie}p" |awk '{printf"%f %f %f %f\n",$3+t+e,NR-1,0,0.2}' t=$bt e=$twe | gmt psxy -J -R -Sy0.2 -W0.5p,magenta  -N -O -K  >>$out
    
    ##(d)
    gmt psbasemap -R10/110/-1/$ymax -JX8/11 -Ba20f10:"Time (s)":/a1weSn  -K -O -X8.5 >>$out
    cat file7.lst |sed -n "${ib},${ie}p" | gmt pssac -J -R -En -M$sc -W0.5p -O -K >> $out 
    cat file8.lst |sed -n "${ib},${ie}p" | gmt pssac -J -R -En -M$sc -W0.5p,red -O -K >> $out 
    ty=`echo $ymax |awk '{print $1*1.05}'`
    gmt pstext -J -R -N -O -K >>$out <<EOF
55 $ty 16 0 0 CM (d) BXR after FWI
EOF
    cat temp1 |sed -n "${ib},${ie}p" |awk '{printf"%f %f %f %f\n",$3+t-b,NR-1,0,0.2}' t=$bt b=$twb | gmt psxy -J -R -Sy0.2 -W0.5p,magenta  -N -O -K >>$out
    cat temp1 |sed -n "${ib},${ie}p" |awk '{printf"%f %f %f %f\n",$3+t+e,NR-1,0,0.2}' t=$bt e=$twe | gmt psxy -J -R -Sy0.2 -W0.5p,magenta  -N -O -K  >>$out

    gmt psxy -J -R -T -O >>$out

    gmt psconvert $out -A -Tg 
    rm $out
    #rm temp temp1
    fi ### end if plot stf
    done
    #########################################
    #rm file?.lst
done

rm *.lst
done # end iset
