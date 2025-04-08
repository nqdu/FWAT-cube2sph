#!/bin/bash
MODEL=M55

for f in solver/${MODEL}/*/*GRAD*;
do 
    rm -rf $f 
done 

for f in solver/${MODEL}.ls/*;
do 
    rm -rf $f 
done 