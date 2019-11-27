# !/bin/bash
set -x

uppath="./output/preprocess"
opath="./output/genome"
rpath="./output/result"
ref="genome"

bowtie -q --threads 40 -v 1 -m 1 -a ${ref} ${uppath}/rename.fastq --un ${opath}/unmap.fastq --sam > ${opath}/${ref}mapped.sam

samtools import ${ref}.fa ${opath}/${ref}mapped.sam ${opath}/${ref}mapped.bam

samtools sort ${opath}/${ref}mapped.bam -o ${opath}/${ref}sorted.bam

samtools index ${opath}/${ref}sorted.bam 

#umi_tools dedup -I ${opath}/${ref}sorted.bam --output-stats=${opath}/${ref}deduplicated -S ${opath}/${ref}deduplicated.bam
umi_tools dedup -I ${opath}/${ref}sorted.bam --method=unique --read-length -S ${opath}/${ref}deduplicated.bam
#umi_tools group -I ${opath}/${ref}sorted.bam --read-length -S ${opath}/${ref}deduplicated.bam

samtools view -h -o ${opath}/${ref}deduplicated.sam ${opath}/${ref}deduplicated.bam 

### name, reads, ref, pos
cat ${opath}/${ref}deduplicated.sam|sed -e '/^@/d' | awk '{print $1,$10,$3,$4}' > ${rpath}/${ref}nohead.txt
### remove @ head
