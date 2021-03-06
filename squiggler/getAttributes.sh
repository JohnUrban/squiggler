#!/bin/bash

if [ $# -eq 0 ]; then
 echo "Usage: getAttributes.sh file.fast5"
 exit 0
fi

SEED=$RANDOM

h5dump -n 1 $1 | grep "attribute" | cut -f 4 -d " " > $SEED.tmp.txt

while read line; do
 value=`h5dump -a $line $1 | grep "(0)" | awk '{gsub(/\"/,""); print $2}'`
 echo -e $line "\t" $value
done < $SEED.tmp.txt

rm ${SEED}.tmp.txt
