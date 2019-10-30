# !/bin/bash
set -x

sample=$1
opath="./output/preprocess"

cutadapt -j 40 -a AGATCGGAAGAGCACACGTC -g GTTCAGAGTTCTACAGTCCGACGATC --trimmed-only --output=${opath}/cutadapt.fq --error-rate=0.1 --times=1 --overlap=3 --minimum-length=27 --maximum-length=55 ${sample}

umi_tools extract --extract-method=string --bc-pattern=NNNNNN --stdin=${opath}/cutadapt.fq --stdout ${opath}/5endumiextracted.fastq --quality-filter-threshold 20 --quality-encoding phred33 --log=${opath}/5endumiextracted.log

umi_tools extract --extract-method=string --bc-pattern=NNNNNNNNN --stdin=${opath}/5endumiextracted.fastq --stdout ${opath}/3endumiextracted.fastq --3prime --quality-filter-threshold 20 --quality-encoding phred33 --log=${opath}/3endumiextracted.log

python rename.py -i ${opath}/3endumiextracted.fastq -o ${opath}/rename.fastq

