#!/bin/bash
set -x
#sample="s1.fastq"
#sample="s3.fastq"
#sample="s5.fastq"
#sample="s11.fastq"
#opath="./output/log"

output='output'
opath='./'${output}'/log'


#for sample in 's1'
for sample in s11
do

mkdir ./${output}
mkdir ./${output}/log
mkdir ./${output}/preprocess
mkdir ./${output}/result
mkdir ./${output}/genome
mkdir ./${output}/other

bash ./preprocess.sh ${sample}.fastq  > ${opath}/preprocess.log 2>&1
bash ./genomemap.sh > ${opath}/genomemap.log 2>&1 
bash ./otherdbmap.sh > ${opath}/otherdbmap.log 2>&1 
#bash ./mergeresult.sh > ${opath}/mergeresult.log 2>&
mv ./output ./output${sample}

done
