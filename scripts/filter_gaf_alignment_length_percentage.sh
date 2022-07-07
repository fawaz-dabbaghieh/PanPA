#!/bin/bash

in_file=$1
cutoff=$2

awk -F "\t" -v cutoff="$cutoff" '{if ((($4 - $3)*100)/$2 > cutoff) {print $0}}' $in_file
