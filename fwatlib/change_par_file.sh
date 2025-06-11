#!/bin/bash

# get input args
param=$1
value=$2
file=$3

# locate parameter
oldstr=`grep "^$param " $file`
newstr="$param           =     $value"

sed -i "s?$oldstr?$newstr?g" $file 
