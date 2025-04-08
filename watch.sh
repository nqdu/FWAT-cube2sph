#!/bin/bash

for outf in `ls output_fwat*`;do
    echo $outf
    tail -2 $outf
done
