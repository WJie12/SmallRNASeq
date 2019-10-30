# !/bin/bash 
set -x
uppath="./output/preprocess"
opath="./output/other"
rpath="./output/result"

cp ${uppath}/rename.fastq ${opath}/unmap.fastq

for ref in 'hg19-tRNAs' 'human_rRNA_5.8S' 'human_rRNA_5S' 'human_rRNA_12S' 'human_rRNA_16S' 'human_rRNA_18S' 'human_rRNA_28S' 'human_rRNA_45S' 'human_rRNA_other' 'miRBase_21-hsa' 'piR_human' 'Rfam-12.3-human'
do
bowtie -q --threads 40 -v 1 -m 1 -a ${ref} ${opath}/unmap.fastq --un ${opath}/unmap.fq --sam > ${opath}/${ref}mapped.sam

samtools import ${ref}.fa ${opath}/${ref}mapped.sam ${opath}/${ref}mapped.bam

samtools sort ${opath}/${ref}mapped.bam -o ${opath}/${ref}sorted.bam

samtools index ${opath}/${ref}sorted.bam 

#umi_tools dedup -I ${opath}/${ref}sorted.bam --output-stats=${opath}/deduplicated -S ${opath}/${ref}deduplicated.bam
umi_tools dedup -I ${opath}/${ref}sorted.bam --method=unique --read-length -S ${opath}/${ref}deduplicated.bam

samtools view -h -o ${opath}/${ref}deduplicated.sam ${opath}/${ref}deduplicated.bam 

mv ${opath}/unmap.fq ${opath}/unmap.fastq

### 1 col $ 3 col
cat ${opath}/${ref}deduplicated.sam|sed -e '/^@/d'|awk '{print $1,$10,$3}' > ${rpath}/${ref}nohead.txt
### remove @ head

done

