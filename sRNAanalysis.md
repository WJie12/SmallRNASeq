#Step1: Mapping 

#```$ bash analysis.sh```

##analysis.sh

```bash
#!/bin/bash
set -x
opath="./output/log"
# 分别对四个数据集s1,s3,s5,s11进行分析
for sample in "s1" "s3" "s5" "s11"
do

mkdir ./output
mkdir ./output/log
mkdir ./output/preprocess
mkdir ./output/result
mkdir ./output/genome
mkdir ./output/other
#预处理
bash ./preprocess.sh ${sample}.fastq  > ${opath}/preprocess.log 2>&1
#genome mapping
bash ./genomemap.sh > ${opath}/genomemap.log 2>&1
#其他small rna db mapping
bash ./otherdbmap.sh > ${opath}/otherdbmap.log 2>&1
#bash ./mergeresult.sh > ${opath}/mergeresult.log 2>&1
mv ./output ./output${sample}
done
```
##genomemap.sh

```bash
# !/bin/bash
set -x

uppath="./output/preprocess"
opath="./output/genome"
rpath="./output/result"
ref="genome"

#bowtie mapping
bowtie -q --threads 40 -v 1 -m 1 -a ${ref} ${uppath}/rename.fastq --sam > ${opath}/${ref}mapped.sam

#sam 2 bam
samtools import ${ref}.fa ${opath}/${ref}mapped.sam ${opath}/${ref}mapped.bam

#sort bam
samtools sort ${opath}/${ref}mapped.bam -o ${opath}/${ref}sorted.bam

#build index
samtools index ${opath}/${ref}sorted.bam

#umi_tool deduplicate based on UMI and read sequence, considering read-length, method=unique
umi_tools dedup -I ${opath}/${ref}sorted.bam --method=unique --read-length -S ${opath}/${ref}deduplicated.bam

#bam 2 sam
samtools view -h -o ${opath}/${ref}deduplicated.sam ${opath}/${ref}deduplicated.bam

### remove @ head, extract col 1,10,3 (name, sequence, ref)
cat ${opath}/${ref}deduplicated.sam|sed -e '/^@/d' | awk '{print $1, $10,$3}' > ${rpath}/${ref}nohead.txt
```
##otherdbmap.sh

```bash
# !/bin/bash 
set -x
uppath="./output/preprocess"
opath="./output/other"
rpath="./output/result"

cp ${uppath}/rename.fastq ${opath}/unmap.fastq

for ref in 'hg19-tRNAs' 'human_rRNA_5.8S' 'human_rRNA_5S' 'human_rRNA_12S' 'human_rRNA_16S' 'human_rRNA_18S' 'human_rRNA_28S' 'human_rRNA_45S' 'human_rRNA_other' 'miRBase_21-hsa' 'piR_human' 'piR_human_rev' 'Rfam-12.3-human'
do
bowtie -q --threads 40 -v 1 -m 1 -a ${ref} ${opath}/unmap.fastq --un ${opath}/unmap.fq --sam > ${opath}/${ref}mapped.sam

samtools import ${ref}.fa ${opath}/${ref}mapped.sam ${opath}/${ref}mapped.bam

samtools sort ${opath}/${ref}mapped.bam -o ${opath}/${ref}sorted.bam

samtools index ${opath}/${ref}sorted.bam

umi_tools dedup -I ${opath}/${ref}sorted.bam --method=unique --read-length -S ${opath}/${ref}deduplicated.bam

samtools view -h -o ${opath}/${ref}deduplicated.sam ${opath}/${ref}deduplicated.bam

mv ${opath}/unmap.fq ${opath}/unmap.fastq

### remove @ head, extract col 1,10,3 (name, sequence, ref)
cat ${opath}/${ref}deduplicated.sam|sed -e '/^@/d'|awk '{print $1,$10,$3}' > ${rpath}/${ref}nohead.txt

done
```

# Step1 Result:

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



#Step2: Result Analysis

#```python resultanalysis.py```

resultanalysis.py

```python
import pandas as pd

#读入txt（小rna数据库结果)
rna_dbs=['hg19-tRNAs','human_rRNA_5.8S','human_rRNA_5S','human_rRNA_12S','human_rRNA_16S','human_rRNA_18S','human_rRNA_28S','human_rRNA_45S','human_rRNA_other','miRBase_21-hsa','piR_human','piR_human_rev','Rfam-12.3-human']
genome='genome'
for idx in range(0, len(rna_dbs)):
    if idx == 0:
        rs_rnadb = pd.read_table(rna_dbs[idx] + 'nohead.txt',header=None,sep=' ',names=['name','reads', rna_dbs[idx]])
        print(rna_dbs[idx],rs_rnadb.shape)
    else:
        df = pd.read_table(rna_dbs[idx] + 'nohead.txt',header=None,sep=' ', names=['name','reads', rna_dbs[idx]])
        print(rna_dbs[idx],df.shape)
        rs_rnadb = pd.merge(rs_rnadb, df, how = 'outer', on = ['name','reads'], suffixes=[rna_dbs[idx-1], rna_dbs[idx]]).fillna('*')
        print(rs_rnadb.shape)

#读入txt（genome mapping结果)
rs_genome = pd.read_table(genome + 'nohead.txt',header=None,sep=' ', names=['name','reads', 'genome'])
print(rs_genome.shape)

#merge genome and small RNA DB mapping results
rs = pd.merge(rs_genome, rs_rnadb, how = 'outer', on = ['name']).fillna('NA')
print(rs.shape)

#group by small rna ref
rs_group=rs.groupby(rna_dbs).describe().reset_index()

#write to csv
rs_group.to_csv('./finalresult.csv',index=False)

```

#Step2 result: finalresult.csv

![1571588668590](C:\Users\admin\AppData\Roaming\Typora\typora-user-images\1571588668590.png)

```
hg19-tRNAs,human_rRNA_5.8S,human_rRNA_5S,human_rRNA_12S,human_rRNA_16S,human_rRNA_18S,human_rRNA_28S,human_rRNA_45S,human_rRNA_other,miRBase_21-hsa,piR_human,piR_human_rev,Rfam-12.3-human,name,name,name,name,reads_x,reads_x,reads_x,reads_x,genome,genome,genome,genome,reads_y,reads_y,reads_y,reads_y
,,,,,,,,,,,,,count,unique,top,freq,count,unique,top,freq,count,unique,top,freq,count,unique,top,freq
*,*,*,*,*,*,*,*,*,*,*,*,AADB02000969.1/1180739-1181019,7,7,15404473_ACCCGGGCGTTCTGT,1,7,1,NA,7,7,1,NA,7,7,7,GGGAGGTGGAGGTGGAGT,1
*,*,*,*,*,*,*,*,*,*,*,*,AADB02001193.1/31306-31018,1,1,13017009_GAAACACTTTTCCCT,1,1,1,NA,1,1,1,NA,1,1,1,AATCTGCAAGTGT,1
*,*,*,*,*,*,*,*,*,*,*,*,AADB02013261.1/2344330-2344071,4,4,1219071_TCCAGGCCTGCCACT,1,4,1,NA,4,4,1,NA,4,4,1,TGGTACTCTTTTT,4
*,*,*,*,*,*,*,*,*,*,*,*,AADB02013937.1/195-495,1,1,10507606_GAGGAGCTCTCAAAT,1,1,1,NA,1,1,1,NA,1,1,1,CTGCACTCCAGCCTGGGCGACGGAGCGAGACTCCC,1
```

# Step3 : Unmapped reads analysis

没有比对到small RNA DB 的reads: unmap.fastq

```
collapse: unmap.fastq -> unmap.fa
```

# Blast+

```shell
$ blastn -db nt -query unmap.fa -outfmt 11 -out 'rnadbunmap.blastn@nr.asn' -num_threads 10

#转M6格式
$ blast_formatter -archive "rnadbunmap.blastn@nr.asn" -outfmt "6 qseqid qlen sseqid sgi slen pident length mismatch gapopen qstart qend sstart send evalue bitscore staxid ssciname" > "rnadbunmap.blastn@nr.blast"

#保留e value 最小的记录
$ awk '!a[$1]++{print}' rnadbunmap.blastn@nr.blast > uniqrnadbunmap.blastn@nr.blast
```

结果目录：/zlab/data/Huang1/smallrnaseqTest/nogenomeresults/temp/

![1571644824063](C:\Users\admin\AppData\Roaming\Typora\typora-user-images\1571644824063.png)

从图中可以看出共12列，下面来列举一下这12列的意思
1、Query id：查询序列ID标识
2、Subject id：比对上的目标序列ID标识
3、% identity：序列比对的一致性百分比
4、alignment length：符合比对的比对区域的长度
5、mismatches：比对区域的错配数
6、gap openings：比对区域的gap数目
7、q. start：比对区域在查询序列(Query id)上的起始位点
8、q. end：比对区域在查询序列(Query id)上的终止位点
9、s. start：比对区域在目标序列(Subject id)上的起始位点
10、s. end：比对区域在目标序列(Subject id)上的终止位点
11、e-value：比对结果的期望值，解释是大概多少次随即比对才能出现一次这个score, E-value越小，表明这种情况从概率上越不可能发生，那么发生了即说明这更有可能是真实的相似序列
12、bit score：比对结果的bit score值
一般情况我们看第3、11、12两列，e值越小越可靠。

python resultanalysis.py -i /zlab/data/Huang1/smallrnaseqTest/outputs3/result/ -o /zlab/data/Huang1/smallrnaseqTest/outputs3/result/
python graph.py -i /zlab/data/Huang1/smallrnaseqTest/outputs3/result/ -o /zlab/data/Huang1/smallrnaseqTest/outputs3/result/
python unmapanalysis.py -i /zlab/data/Huang1/smallrnaseqTest/outputs3/other/ -o /zlab/data/Huang1/smallrnaseqTest/outputs3/other/
python unmapanalysis.py -i /zlab/data/Huang1/smallrnaseqTest/outputs3/other/ -o /zlab/data/genome/smallrnaseqTest/outputs3/genome/
