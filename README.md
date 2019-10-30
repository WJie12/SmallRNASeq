# Step1: Mapping: shell scripts
```shell
$ bash analysis.sh
```
Required file folder framework: (put the shell scripts, fastq files and fa files at the same folder)
```
.
├── analysis.sh
├── preprocess.sh
├── rename.py
├── genomemap.sh
├── otherdbmap.sh
├── *.fastq (s1,s3,s5,s11)
├── *.fa (ref libraries)
```

After mapping, there should be output folders for sample files. 
Inside the output folder, such as 'outpus1' for s1, there should be 4 folders (log, preprocess, genome, other) as follow:

```
.
├── genome
│   ├── genomededuplicated.bam
│   ├── genomededuplicated.sam
│   ├── genomemapped.bam
│   ├── genomemapped.sam
│   ├── genomesorted.bam
│   └── genomesorted.bam.bai
├── log
│   ├── genomemap.log
│   ├── otherdbmap.log
│   └── preprocess.log
├── other
│   ├── hg19-tRNAsdeduplicated.bam
│   ├── hg19-tRNAsdeduplicated.sam
│   ├── hg19-tRNAsmapped.bam
│   ├── hg19-tRNAsmapped.sam
│   ├── hg19-tRNAssorted.bam
│   ├── hg19-tRNAssorted.bam.bai
│   ├── human_rRNA_12Sdeduplicated.bam
│   ├── human_rRNA_12Sdeduplicated.sam
│   ├── human_rRNA_12Smapped.bam
│   ├── human_rRNA_12Smapped.sam
│   ├── human_rRNA_12Ssorted.bam
│   ├── human_rRNA_12Ssorted.bam.bai
│   ├── human_rRNA_16Sdeduplicated.bam
│   ├── human_rRNA_16Sdeduplicated.sam
│   ├── human_rRNA_16Smapped.bam
│   ├── human_rRNA_16Smapped.sam
│   ├── human_rRNA_16Ssorted.bam
│   ├── human_rRNA_16Ssorted.bam.bai
│   ├── human_rRNA_18Sdeduplicated.bam
│   ├── human_rRNA_18Sdeduplicated.sam
│   ├── human_rRNA_18Smapped.bam
│   ├── human_rRNA_18Smapped.sam
│   ├── human_rRNA_18Ssorted.bam
│   ├── human_rRNA_18Ssorted.bam.bai
│   ├── human_rRNA_28Sdeduplicated.bam
│   ├── human_rRNA_28Sdeduplicated.sam
│   ├── human_rRNA_28Smapped.bam
│   ├── human_rRNA_28Smapped.sam
│   ├── human_rRNA_28Ssorted.bam
│   ├── human_rRNA_28Ssorted.bam.bai
│   ├── human_rRNA_45Sdeduplicated.bam
│   ├── human_rRNA_45Sdeduplicated.sam
│   ├── human_rRNA_45Smapped.bam
│   ├── human_rRNA_45Smapped.sam
│   ├── human_rRNA_45Ssorted.bam
│   ├── human_rRNA_45Ssorted.bam.bai
│   ├── human_rRNA_5.8Sdeduplicated.bam
│   ├── human_rRNA_5.8Sdeduplicated.sam
│   ├── human_rRNA_5.8Smapped.bam
│   ├── human_rRNA_5.8Smapped.sam
│   ├── human_rRNA_5.8Ssorted.bam
│   ├── human_rRNA_5.8Ssorted.bam.bai
│   ├── human_rRNA_5Sdeduplicated.bam
│   ├── human_rRNA_5Sdeduplicated.sam
│   ├── human_rRNA_5Smapped.bam
│   ├── human_rRNA_5Smapped.sam
│   ├── human_rRNA_5Ssorted.bam
│   ├── human_rRNA_5Ssorted.bam.bai
│   ├── human_rRNA_otherdeduplicated.bam
│   ├── human_rRNA_otherdeduplicated.sam
│   ├── human_rRNA_othermapped.bam
│   ├── human_rRNA_othermapped.sam
│   ├── human_rRNA_othersorted.bam
│   ├── human_rRNA_othersorted.bam.bai
│   ├── miRBase_21-hsadeduplicated.bam
│   ├── miRBase_21-hsadeduplicated.sam
│   ├── miRBase_21-hsamapped.bam
│   ├── miRBase_21-hsamapped.sam
│   ├── miRBase_21-hsasorted.bam
│   ├── miRBase_21-hsasorted.bam.bai
│   ├── piR_humandeduplicated.bam
│   ├── piR_humandeduplicated.sam
│   ├── piR_humanmapped.bam
│   ├── piR_humanmapped.sam
│   ├── piR_human_revdeduplicated.bam
│   ├── piR_human_revdeduplicated.sam
│   ├── piR_human_revmapped.sam
│   ├── piR_humansorted.bam
│   ├── piR_humansorted.bam.bai
│   ├── Rfam-12.3-humandeduplicated.bam
│   ├── Rfam-12.3-humandeduplicated.sam
│   ├── Rfam-12.3-humanmapped.bam
│   ├── Rfam-12.3-humanmapped.sam
│   ├── Rfam-12.3-humansorted.bam
│   ├── Rfam-12.3-humansorted.bam.bai
│   └── unmap.fastq
├── preprocess
│   ├── 3endumiextracted.fastq
│   ├── 3endumiextracted.log
│   ├── 5endumiextracted.fastq
│   ├── 5endumiextracted.log
│   ├── cutadapt.fq
│   └── rename.fastq
└── result
    ├── genomenohead.txt
    ├── hg19-tRNAsnohead.txt
    ├── human_rRNA_12Snohead.txt
    ├── human_rRNA_16Snohead.txt
    ├── human_rRNA_18Snohead.txt
    ├── human_rRNA_28Snohead.txt
    ├── human_rRNA_45Snohead.txt
    ├── human_rRNA_5.8Snohead.txt
    ├── human_rRNA_5Snohead.txt
    ├── human_rRNA_othernohead.txt
    ├── miRBase_21-hsanohead.txt
    ├── piR_humannohead.txt
    ├── piR_human_revnohead.txt
    └── Rfam-12.3-humannohead.txt


```


# Step2: Analysis: Python scripts

**resultanalysis.py:** Analysis mapping results. 

​	Input: *nohead.text files reduced by mapping. (in the 'result' folder)

​	Output: 'y_rnadb_y_genome_result.csv' and 'n_rnadb_y_genome_result.csv'

**unmapanalysis.py:** transfer fastq to fasta. Remove PCR duplication (deduplicate by UMI and reads). Count and remove dupication (deduplicate by reads). 

​	Input: unmapped reads in fastq format. 

​	output: unmapped reads in fastq format with idx and count number in name.

**graph.py:** draw pie graph for 'y_rnadb_y_genome_result.csv'.

Command Example:
```shell
python resultanalysis.py -i /zlab/data/Huang1/smallrnaseqTest/outputs3/result/ -o /zlab/data/Huang1/smallrnaseqTest/outputs3/result/
python graph.py -i /zlab/data/Huang1/smallrnaseqTest/outputs3/result/ -o /zlab/data/Huang1/smallrnaseqTest/outputs3/result/
python unmapanalysis.py -i /zlab/data/Huang1/smallrnaseqTest/outputs3/other/ -o /zlab/data/Huang1/smallrnaseqTest/outputs3/other/
python unmapanalysis.py -i /zlab/data/Huang1/smallrnaseqTest/outputs3/genome/ -o /zlab/data/Huang1/smallrnaseqTest/outputs3/genome/
```



# BLAST+ analysis

```shell
$ blastn -db nt -query unmap.fa -outfmt 11 -out 'rnadbunmap.blastn@nr.asn' -num_threads 10

#转M6格式
$ blast_formatter -archive "rnadbunmap.blastn@nr.asn" -outfmt "6 qseqid qlen sseqid sgi slen pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid ssciname" > "rnadbunmap.blastn@nr.blast"

#保留e value 最小的记录
$ awk '!a[$1]++{print}' rnadbunmap.blastn@nr.blast > uniqrnadbunmap.blastn@nr.blast
```

Blast -m 6格式就是 Blast -m 8格式。

-m 6格式共12列，每列意思如下：

1. Query id：查询序列ID标识

2. Subject id：比对上的目标序列ID标识

3. % identity：序列比对的一致性百分比

4. alignment length：符合比对的比对区域的长度

5. mismatches：比对区域的错配数

6. gap openings：比对区域的gap数目

7. q. start：比对区域在查询序列(Query id)上的起始位点

8. q. end：比对区域在查询序列(Query id)上的终止位点

9. s. start：比对区域在目标序列(Subject id)上的起始位点

10. s. end：比对区域在目标序列(Subject id)上的终止位点

11. e-value：比对结果的期望值，解释是大概多少次随即比对才能出现一次这个score, E-value越小，表明这种情况从概率上越不可能发生，那么发生了即说明这更有可能是真实的相似序列

12. bit score：比对结果的bit score值

    一般情况我们看第3、11、12两列，e值越小越可靠。