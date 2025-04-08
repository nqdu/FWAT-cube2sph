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

anno=("nul" "(a)" "(b)" "(c)" "(d)" "(e)" "(f)" "(g)" "(h)" "(i)" "(j)" "(k)" "(l)" "(m)" "(n)")
MOD=M20
tb=-20.
te=240.
sc=1.0
seisdir=$MOD
rm -rf $MOD
mkdir -p $seisdir
band1=T006_T015
band2=T010_T020
band3=T020_T035
num_filter_current=1
if (($num_filter_current==1));then
    band_all=$band1
elif (($num_filter_current==2));then
    band_all=`echo $band1 $band2 |awk '{print $1,$2}'`
elif (($num_filter_current==3));then
    band_all=`echo $band1 $band2 $band3|awk '{print $1,$2,$3}'`
fi
echo "bandall=$band_all"

cat source.lst | while read soline;do
evet=`echo $soline |awk '{print $2"."$1"_Z"}'`
evett=`echo $soline |awk '{print $2"."$1}'`
cat station.lst |sed -n '1,8p'| while read stline;do
  stat=`echo $stline |awk '{print $2"."$1}'`
  statt=`echo $stline |awk '{print $1}'`
  echo $evet $stat
  for band in $band_all;do
  fmist=$seisdir/misfit_${evett}_${stat}_$band 
  fwddir1=../../solver/M00.set1/$evet/OUTPUT_FILES 
  fwddir2=../../../mcr_west/solver/M00.set1/$evet/OUTPUT_FILES 
  fwddir3=../../../mcr_west/solver/M20.set1/$evet/OUTPUT_FILES
  echo $evet $statt
  cat ../../misfits/M00.set1_${band}_window_chi|grep $evet|grep $statt |awk '{print $13}' >$fmist
  cat ../../misfits/M00.set1_${band}_window_chi|grep $evet|grep $statt |awk '{print $13}'>>$fmist
  cat ../../misfits/M00.set1_${band}_window_chi|grep $evet|grep $statt |awk '{print $13}'>>$fmist
  hp=`echo $band |awk '{printf"%d",substr($1,2,3)}'`
  lp=`echo $band |awk '{printf"%d",substr($1,7,3)}'`
  R0=2.5
  R1=4.0
  saclst knetwk kstnm dist f ../../data/$evet/$stat.BXZ.sac |awk '{print $2,$3,$4}' >$seisdir/dists.dat
  saclst dist f ../../data/$evet/$stat.BXZ.sac |awk '{printf"%f %f\n",$2/vmax-lp/2,$2/vmin+lp/2}' vmin=$R0 vmax=$R1 lp=$lp hp=$hp>$seisdir/windows.${evett}_${stat}.$band.dat
  #########get the time window##########
  cat /dev/null >file1.$band.lst
  cat /dev/null >file.${evett}_${stat}.$band.lst
  ls $fwddir2/$stat.BXZ.obs.sac.$band >>file1.$band.lst
  ls $fwddir1/$stat.BXZ.syn.sac.$band >>file1.$band.lst
  ls $fwddir2/$stat.BXZ.syn.sac.$band >>file1.$band.lst
  ls $fwddir3/$stat.BXZ.syn.sac.$band >>file1.$band.lst
  #########cut data accoring to the time window##########
  itr=1
  for fname in `cat file1.$band.lst`;do
    tracename=${itr}.${evett}_${stat}.$band.sac
    echo $tracename
    DELTA=`saclst delta f $fname |awk '{print $2}'`
    netwk=`echo $fname |awk -F$fwddir/ '{print $2}' |awk -F. '{print $1}'`
    stnm=`echo $fname |awk -F$fwddir/ '{print $2}' |awk -F. '{print $2}'`
    cutB=`cat $seisdir/windows.${evett}_${stat}.$band.dat |awk '{print $1}'`
    cutE=`cat $seisdir/windows.${evett}_${stat}.$band.dat |awk '{print $2}'`
    cutB=`echo $cutB $tb|awk '{print $1-$2}'`
    cutE=`echo $cutE $tb|awk '{print $1-$2}'`
    cutN=`echo $cutB $cutE $DELTA |awk '{print ($2-$1)/$3}' | awk '{printf("%d\n",$0+=$0<0?0:0.9999)}'`
    echo $seisdir/${tracename}.cut >>file.${evett}_${stat}.$band.lst
    sac<<EOF
cut B $cutB N $cutN
r $fname
w $seisdir/${tracename}.cut
q
EOF
  ((itr++))
  done # end traces
  done # end band
done # end station
  ##########################################
  echo $evet
  if true;then
  out=$seisdir/fsismo_${evet}.ps
  ymax=0.8
  xspace=8
  yspace=`echo $ymax |awk '{print 2*$1}'`
  ist=1
  cat station.lst |sed -n "1,8p" | while read stline;do
  stat=`echo $stline |awk '{print $2"."$1}'`
  dx_move=`echo $xspace|awk '{print $1+0.5}'`
  dy_move=`echo $yspace |awk '{print $1}'`
  ii=$(((ist)%2))
  if [ $ii = 1 ];then
    XX=`echo $dx_move |awk '{print -1.0*$1}'`
    YY=`echo $dy_move $ymax |awk '{print -($1+$2)}'`
    echo "ii="$ii $XX $YY
  else
    XX=`echo $dx_move |awk '{print 1.0*$1}'`
    YY=`echo $dy_move |awk '{print 2*$1}'`
    echo "ii="$ii $XX $YY
  fi
  for band in $band_all;do
    cutB=`cat $seisdir/windows.${evett}_${stat}.$band.dat |awk '{print $1}'`
    cutE=`cat $seisdir/windows.${evett}_${stat}.$band.dat |awk '{print $2}'`
    fmist=$seisdir/misfit_${evett}_${stat}_$band 
    deltat1=`cat $fmist|sed -n '1,1p'|awk '{printf "%1.2f", $1}'`
    deltat2=`cat $fmist|sed -n '2,2p'|awk '{printf "%1.2f", $1}'`
    deltat3=`cat $fmist|sed -n '3,3p'|awk '{printf "%1.2f", $1}'`
    echo $evett $stat $fmist 'delta=' $deltat1 $deltat2 $deltat3
    echo $cutB $cutE
    if [ $band = "$band1" ];then
        if [ $ist = 1 ];then
          gmt psbasemap -R$cutB/$cutE/-${ymax}/$ymax -JX$xspace/$yspace -Ba20f10/a1:"Trace number":ws -K -P -Y25> $out
          cat file.${evett}_${stat}.$band.lst |sed -n "1,1p" | gmt pssac -J -R -Entb -M$sc -W0.5p -O -K>>$out
          cat file.${evett}_${stat}.$band.lst |sed -n "2,2p" | gmt pssac -J -R -Entb -M$sc -W0.5p,- -O -K>>$out
          cat file.${evett}_${stat}.$band.lst |sed -n "3,3p" | gmt pssac -J -R -Entb -M$sc -W0.5p,255/0/0 -O -K>>$out
          cat file.${evett}_${stat}.$band.lst |sed -n "4,4p" | gmt pssac -J -R -Entb -M$sc -W0.5p,0/255/0 -O -K>>$out
        else
          gmt psbasemap -R$cutB/$cutE/-${ymax}/$ymax -JX$xspace/$yspace -Ba20f10/a1:"Trace number":ws -K -O -Y${YY} -X${XX}>> $out
          cat file.${evett}_${stat}.$band.lst |sed -n "1,1p" | gmt pssac -J -R -Entb -M$sc -W0.5p -O -K>>$out
          cat file.${evett}_${stat}.$band.lst |sed -n "2,2p" | gmt pssac -J -R -Entb -M$sc -W0.5p,- -O -K>>$out
          cat file.${evett}_${stat}.$band.lst |sed -n "3,3p" | gmt pssac -J -R -Entb -M$sc -W0.5p,255/0/0 -O -K>>$out
          cat file.${evett}_${stat}.$band.lst |sed -n "4,4p" | gmt pssac -J -R -Entb -M$sc -W0.5p,0/255/0 -O -K>>$out
        fi
        tx=`echo $cutB $cutE |awk '{print ($1+$2)/2.0}'`
        ty=$ymax
    gmt pstext -J -R -N -K -O >>$out <<EOF
$tx $ty 14 0 0 CM ${anno[$ist]}  ${evett}-${stat}
EOF
    elif [ $band = "$band3" ];then
        if [ $ist -le 6 ];then
            gmt psbasemap -R$cutB/$cutE/-${ymax}/$ymax -JX$xspace/$yspace -Ba20f10/a1:"Trace number":wS -K -O -Y-${dy_move} -X0.0>> $out
        else
            gmt psbasemap -R$cutB/$cutE/-${ymax}/$ymax -JX$xspace/$yspace -Ba20f10:"Time (s)":wS -K -O -Y-${dy_move} -X0.0>> $out
        fi
        cat file.${evett}_${stat}.$band.lst |sed -n "1,1p" | gmt pssac -J -R -Entb -M$sc -W0.5p -O -K>>$out
        cat file.${evett}_${stat}.$band.lst |sed -n "2,2p" | gmt pssac -J -R -Entb -M$sc -W0.5p,- -O -K>>$out
        cat file.${evett}_${stat}.$band.lst |sed -n "3,3p" | gmt pssac -J -R -Entb -M$sc -W0.5p,255/0/0 -O -K>>$out
        cat file.${evett}_${stat}.$band.lst |sed -n "4,4p" | gmt pssac -J -R -Entb -M$sc -W0.5p,0/255/0 -O -K>>$out
    else
        gmt psbasemap -R$cutB/$cutE/-${ymax}/$ymax -JX$xspace/$yspace -Ba20f10/a1:"Trace number":ws -K -O -Y-${dy_move} -X0.0>> $out
        cat file.${evett}_${stat}.$band.lst |sed -n "1,1p" | gmt pssac -J -R -Entb -M$sc -W0.5p -O -K>>$out
        cat file.${evett}_${stat}.$band.lst |sed -n "2,2p" | gmt pssac -J -R -Entb -M$sc -W0.5p,- -O -K>>$out
        cat file.${evett}_${stat}.$band.lst |sed -n "3,3p" | gmt pssac -J -R -Entb -M$sc -W0.5p,255/0/0 -O -K>>$out
        cat file.${evett}_${stat}.$band.lst |sed -n "4,4p" | gmt pssac -J -R -Entb -M$sc -W0.5p,0/255/0 -O -K>>$out
    fi
    txx=`echo $cutB $cutE |awk '{print $1+($2-$1)/4}'`
    gmt pstext -R -J -W0.5p -G255/255/255 -O -K -N <<EOF >> $out
    $txx 0.4 8 0 1 ML $deltat1
EOF
    gmt pstext -R -J -W0.5p -G255/0/0 -O -K -N <<EOF >> $out
    $txx 0.0 8 0 1 ML $deltat2
EOF
    gmt pstext -R -J -W0.5p -G0/255/0 -O -K -N <<EOF >> $out
    $txx -0.4 8 0 1 ML $deltat3
EOF

  done # end band
  ((ist++))
  done # end station
  cat /dev/null |gmt psxy -J -R -O >>$out
  fi
##    #######################################
done # end source
