#!/bin/bash

for outf in `ls LOG/output_fwat*`;do
    echo $outf
    tail -2 $outf
done
