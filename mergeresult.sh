#!bin/bash 
set -x

rpath="./output/result"
rname="human_rRNA_5.8S"

for ref in "human_rRNA_5S" "human_rRNA_12S" "human_rRNA_16S" "human_rRNA_18S" "human_rRNA_28S" "human_rRNA_45S" "human_rRNA_other" "hg19-tRNAs" "miRBase_21-hsa" "Rfam-12.3-human" "piR_human"
do

### merge by first col
sort -o ${rpath}/${rname}noheadsorted.txt ${rpath}/${rname}nohead.txt
sort -o ${rpath}/${ref}noheadsorted.txt  ${rpath}/${ref}nohead.txt
join -a1 -a2  ${rpath}/${rname}noheadsorted.txt ${rpath}/${ref}noheadsorted.txt  >  ${rpath}/${rname}_${ref}nohead.txt

rname=${rname}"_"${ref}

done

join -j 1 <(sort ${rpath}/genomenohead.txt) <(sort ${rpath}/${rname}nohead.txt)  >  ${rpath}/result.txt
