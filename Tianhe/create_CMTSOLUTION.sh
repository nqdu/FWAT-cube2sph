#/bin/bash

# loop around all teleset
for Set in src_rec/sources_set*.dat;
do
    # only keep tele source
    idx=`echo $Set| cut -d't' -f2 | cut -d'.' -f1`
    if [[ $idx -le 10  ]];then  continue; fi

    echo $Set

    # open file and read line by line
    cat $Set| while read line
    do
        srcidx=`echo $line | awk '{print $1}'`
        evla=`echo $line | awk '{print $2}'`
        evlo=`echo $line | awk '{print $3}'`
        evdp=`echo $line | awk '{print $4}'`
        cp DATA/CMTSOLUTION CMTSOLUTION_$srcidx
        sed -i "/latorUTM: /c\latorUTM:  $evla" CMTSOLUTION_$srcidx
        sed -i "/longorUTM: /c\longorUTM:  $evlo" CMTSOLUTION_$srcidx
        sed -i "/depth: /c\depth: $evdp" CMTSOLUTION_$srcidx
        mv CMTSOLUTION_$srcidx src_rec
        #echo $evla $evlo $evdp
    done
done