#!/bin/bash

in_file=$1
cutoff=$2

# first we need to find which column has the NM tag
nm_column=`head -1 $in_file | awk '{for(i=1;i<=NF;i++){if($i~/^id:f/){print i}} }'`

awk -F "\t" -v cutoff="$cutoff" -v nm_column="$nm_column" '{split($nm_column,score,":"); if (score[3] >= cutoff) {print $0}}' $in_file

