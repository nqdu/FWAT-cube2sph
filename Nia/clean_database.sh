#!/bin/bash

if [ "$#" -ne "1" ]; then 
    echo "Usage ./clean_database.sh M00"
    exit 1
fi

modir=$1
for f in Database absorb_dsm inner normal;
do
    rm -f $modir/*${f}*
done
