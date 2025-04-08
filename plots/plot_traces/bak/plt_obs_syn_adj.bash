#!/bin/bash

#gmt gmtset FONT_TITLE 24p
#gmt gmtset FONT_ANNOT_PRIMARY 14p
#gmt gmtset FONT_LABEL 16p
#gmt gmtset MAP_FRAME_TYPE plain
gmt gmtset ANNOT_FONT_SIZE_PRIMARY 12p
gmt gmtset ANNOT_OFFSET_PRIMARY 0.1c
gmt gmtset LABEL_FONT_SIZE 14p 
gmt gmtset LABEL_OFFSET 0.15c
gmt gmtset HEADER_FONT_SIZE 16p
gmt gmtset HEADER_OFFSET -0.5c 
gmt gmtset TICK_LENGTH -0.2c
gmt gmtset BASEMAP_TYPE plain 
export SAC_DISPLAY_COPYRIGHT=0

MOD=M02
tb=-20.
te=240.
sc=0.20
seisdir=$MOD
rm -rf $MOD
mkdir -p $seisdir
band1=T010_T025
band2=T010_T025
band3=T008_T020
num_filter_current=1
if (($num_filter_current==1));then
    band_all=$band1
elif (($num_filter_current==2));then
    band_all=`echo $band1 $band2 |awk '{print $1,$2}'`
elif (($num_filter_current==3));then
    band_all=`echo $band1 $band2 $band3|awk '{print $1,$2,$3}'`
fi
echo "bandall=$band_all"

for band in $band_all;do
for iset in 1 4 7;do
  cat ../../src_rec/sources_set$iset.dat |
  while read line;do
  ievt=`echo $line |awk '{printf $1}'`
  echo $ievt
  fwddir=../../solver/$MOD.set$iset/$ievt/OUTPUT_FILES #solver/M00.set1/5/OUTPUT_FILES
  mkdir -p $seisdir/$ievt
  hp=`echo $band |awk '{printf"%d",substr($1,2,3)}'`
  lp=`echo $band |awk '{printf"%d",substr($1,7,3)}'`
  R0=2.8
  R1=4.0
  saclst knetwk kstnm dist f ../../data/$ievt/*BXZ.sac |awk '{print $2,$3,$4}' >$seisdir/$ievt/dists.dat
  saclst dist f ../../data/$ievt/*BXZ.sac |awk '{printf"%f %f\n",$2/vmax-lp/2,$2/vmin+lp/2}' vmin=$R0 vmax=$R1 lp=$lp hp=$hp>$seisdir/$ievt/windows.dat
  paste $seisdir/$ievt/dists.dat $seisdir/$ievt/windows.dat >$seisdir/$ievt/MEASUREWINDOWS.dat
  #########get the time window##########
  for comp in Z ;do
    ls $fwddir/*BX${comp}.obs.sac.$band >file1.lst
    cat /dev/null >tmp_dist
    cat file1.lst|while read fileline;do
	    prefix=`echo ${fileline##*/}`
        tmp_stat=`echo $prefix |awk -F. '{print $1"."$2"."$3".sac"}'` 
        saclst dist f ../../data/$ievt/$tmp_stat |awk '{print $2}' >>tmp_dist
        
    done
    paste file1.lst tmp_dist>file1_dist
    sort -k2 -n file1_dist| awk '{print $1}'>file1.lst
    cat file1.lst |sed -e "s/obs/syn/g" >file2.lst
    cat /dev/null >file3.lst
    cat /dev/null >file4.lst
    #########cut data accoring to the time window##########
    for fname in `cat file1.lst`;do
      tracename=`echo $fname |awk -F$fwddir/ '{print $2}'`
      DELTA=`saclst delta f $fname |awk '{print $2}'`
      netwk=`echo $fname |awk -F$fwddir/ '{print $2}' |awk -F. '{print $1}'`
      stnm=`echo $fname |awk -F$fwddir/ '{print $2}' |awk -F. '{print $2}'`
      cutB=`grep $stnm $seisdir/$ievt/MEASUREWINDOWS.dat |awk '{print $4}'`
      cutB=`echo $cutB $tb|awk '{print $1-$2}'`
      cutE=`grep $stnm $seisdir/$ievt/MEASUREWINDOWS.dat |awk '{print $5}'`
      cutE=`echo $cutE $tb|awk '{print $1-$2}'`
      cutN=`echo $cutB $cutE $DELTA |awk '{print ($2-$1)/$3}' | awk '{printf("%d\n",$0+=$0<0?0:0.9999)}'`
      echo $seisdir/$ievt/${tracename}.cut >>file3.lst
      sac<<EOF
cut B $cutB N $cutN
r $fname
w $seisdir/$ievt/${tracename}.cut
q
EOF
    done
    for fname in `cat file2.lst`;do
      tracename=`echo $fname |awk -F$fwddir/ '{print $2}'`
      DELTA=`saclst delta f $fname |awk '{print $2}'`
      netwk=`echo $fname |awk -F$fwddir/ '{print $2}' |awk -F. '{print $1}'`
      stnm=`echo $fname |awk -F$fwddir/ '{print $2}' |awk -F. '{print $2}'`
      cutB=`grep $stnm $seisdir/$ievt/MEASUREWINDOWS.dat |awk '{print $4}'`
      cutB=`echo $cutB $tb|awk '{print $1-$2}'`
      cutE=`grep $stnm $seisdir/$ievt/MEASUREWINDOWS.dat |awk '{print $5}'`
      cutE=`echo $cutE $tb|awk '{print $1-$2}'`
      cutN=`echo $cutB $cutE $DELTA |awk '{print ($2-$1)/$3}' | awk '{printf("%d\n",$0+=$0<0?0:0.9999)}'`
      echo $seisdir/$ievt/${tracename}.cut >>file4.lst
      sac<<EOF
cut B $cutB N $cutN
r $fname
w $seisdir/$ievt/${tracename}.cut
q
EOF
    done
###adj
    adjdir=../../solver/$MOD.set$iset/$ievt/SEM 
    cat /dev/null >file7.lst
    cat file1.lst |sed -e "s/obs/adj/g" >file5.lst
    cat file5.lst |sed -e "s/OUTPUT_FILES/SEM/g" >file6.lst
    tbb=-1.0
    for fname in `cat file6.lst`;do
      tracename=`echo $fname |awk -F$adjdir/ '{print $2}'`
      DELTA=`saclst delta f $fname |awk '{print $2}'`
      netwk=`echo $fname |awk -F$adjdir/ '{print $2}' |awk -F. '{print $1}'`
      stnm=`echo $fname |awk -F$adjdir/ '{print $2}' |awk -F. '{print $2}'`
      cutB=`grep $stnm $seisdir/$ievt/MEASUREWINDOWS.dat |awk '{print $4}'`
      cutB=`echo $cutB $tbb|awk '{print $1-$2}'`
      cutE=`grep $stnm $seisdir/$ievt/MEASUREWINDOWS.dat |awk '{print $5}'`
      cutE=`echo $cutE $tbb|awk '{print $1-$2}'`
      cutN=`echo $cutB $cutE $DELTA |awk '{print ($2-$1)/$3}' | awk '{printf("%d\n",$0+=$0<0?0:0.9999)}'`
      echo $seisdir/$ievt/${tracename}.cut >>file7.lst
      sac<<EOF
cut B $cutB N $cutN
r $fname
w $seisdir/$ievt/${tracename}.cut
q
EOF
      depmax=`saclst depmax f ${seisdir}/${ievt}/${tracename}.cut |awk '{print $2}'`
      fsyn=$seisdir/$ievt/${netwk}.${stnm}.BX${comp}.syn.sac.${band}.cut
      fobs=$seisdir/$ievt/${netwk}.${stnm}.BX${comp}.obs.sac.${band}.cut
      if (( $(echo "$depmax<0.000001" |bc -l) ));then
        sac<<EOF
r $fsyn
mul 0.0
w $fsyn
q
EOF
        sac<<EOF
r $fobs
mul 0.0
w $fobs
q
EOF
      fi
    done


    nsta=`cat file2.lst |wc -l`
    len=56
    npart=`echo $nsta $len |awk '{print $1/$2}' | awk '{printf("%d\n",$0+=$0<0?0:0.9999)}'`
    echo "nsta,len,npart=" $nsta $len $npart
    for i in `seq $npart`;do
      ib=`echo $i $len |awk '{print ($1-1)*$2+1}'`
      ie=`echo $i $len |awk '{print $1*$2}'`
      if [ ${ie} -gt $nsta ];then
          ie=$nsta
      fi
      ymax=`echo $len |awk '{print $1+1}'`
      echo $ib $ie $ymax
    ##########################################
      if true;then
      out=$seisdir/$ievt/${MOD}_fsismo_${ievt}_BX${comp}_${band}_${ib}_${ie}.seg.ps
      xspace=8
      yspace=12
      gmt psbasemap -R$tb/$te/-1/$ymax -JX$xspace/$yspace -Ba40f20/a1:"Trace number":weSn -K -P > $out
      cat file3.lst |sed -n "${ib},${ie}p" | gmt pssac -J -R -Entb -M$sc -W0.5p -O -K>>$out
      cat file4.lst |sed -n "${ib},${ie}p" | gmt pssac -J -R -Entb -M$sc -W0.5p,255/0/0 -O -K>>$out
    cat /dev/null >temp
    cat /dev/null >temp1
    for fname in `cat file3.lst`;do
      netwk=`echo $fname |awk -F$seisdir/$ievt/ '{print $2}' |awk -F. '{print $1}'`
      stnm=`echo $fname |awk -F$seisdir/$ievt/ '{print $2}' |awk -F. '{print $2}'`
      echo $netwk $stnm >>temp 
      grep $stnm $seisdir/$ievt/MEASUREWINDOWS.dat >>temp1 
    done
    cat temp |sed -n "${ib},${ie}p" |awk '{printf"%f %f 8 0 0 CM %s.%s\n",-25,NR-1,$1,$2}' |gmt pstext -J -R -K -W255/0/0 -N  -O  >>$out 
    cat temp1 |sed -n "${ib},${ie}p" |awk '{printf"%f %f %f %f\n",$5,NR-1,0,0.2}' |gmt psxy -N -O  -J -R -K -W0.5p,blue >>$out
    cat temp1 |sed -n "${ib},${ie}p" |awk '{printf"%f %f %f %f\n",$4,NR-1,0,0.2}' |gmt psxy -N -O  -J -R -K -W0.5p,blue >>$out
#    #######################################
    if true;then
      dxspace=1.2
      dyspace=1.2
      XX=`echo $xspace $dxspace |awk '{print ($1+$2)}'`
      YY=`echo $yspace $dyspace |awk '{print ($1+$2)}'`
      echo $XX $YY
    gmt psbasemap -R-1/$te/-1/$ymax -JX$xspace/$yspace -Ba40f20/a1:"Trace number":weSn  -O -K -X${XX} >> $out
    cat file7.lst |sed -n "${ib},${ie}p" | gmt pssac -J -R -Entb -M$sc -W0.5p,255/0/0 -O -K >>$out
    ty=`echo $ymax |awk '{print $1*1.05}'`
    gmt pstext -J -R -N -K -O >>$out <<EOF
-150 $ty 14 0 0 CM (a) DATA and SYN at $band 
EOF
    cat /dev/null >temp2
    cat /dev/null >temp3
    for fname in `cat file7.lst`;do
      netwk=`echo $fname |awk -F$seisdir/$ievt/ '{print $2}' |awk -F. '{print $1}'`
      stnm=`echo $fname |awk -F$seisdir/$ievt/ '{print $2}' |awk -F. '{print $2}'`
      echo $netwk $stnm >>temp2
      grep $stnm $seisdir/$ievt/MEASUREWINDOWS.dat >>temp3 
    done
    cat temp2 |sed -n "${ib},${ie}p" |awk '{printf"%f %f 8 0 0 CM %s.%s\n",-16,NR-1,$1,$2}' |gmt pstext -J -R -W255/0/0 -N -K -O >>$out 
    cat temp3 |sed -n "${ib},${ie}p" |awk '{printf"%f %f %f %f\n",$5,NR-1,0,0.2}' |gmt psxy -N -O -K -J -R -W0.5p,blue >>$out
    cat temp1 |sed -n "${ib},${ie}p" |awk '{printf"%f %f %f %f\n",$4,NR-1,0,0.2}' |gmt psxy -N -O -K -J -R -W0.5p,blue >>$out
      ty=`echo $ymax |awk '{print $1*1.05}'`
      gmt pstext -J -R -N -K -O >>$out <<EOF
100 $ty 14 0 0 CM (b) Ajoint Source at $band 
EOF
    fi
    rm temp* 
#    #########################################
      tyy=`echo $ymax |awk '{print $1*1.15}'`
      gmt pstext -J -R -N -O >>$out <<EOF
0 $tyy 16 0 0 CM EventID: $ievt
EOF
    fi ### end if plot adj
    done
  done
  done
done # end iset
done # end band
