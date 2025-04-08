#!/bin/bash
set -e 

if [ $# -ne 3 ]; then 
    echo "Usage ./this setnum sourcename workdir"
    exit 1
fi 

setnum=$1
SOURCENAME=$2 
workdir=$3 

saclst dist f $workdir/data/$SOURCENAME/*.sac | while read line ;
do 
    stafile=$(basename `echo $line | awk '{print $1}'`)
    dist=`echo $line |awk '{print $2}'`
    netwk=`echo $stafile |cut -d'.' -f1`
    name=`echo $stafile |cut -d'.' -f2`
    sac << EOF 
    r sources_set$setnum/${name}.${netwk}.BXZ.*.sac*
    ch dist $dist 
    w over 
    q
EOF
done 


