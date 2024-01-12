# cancer_genomics_workshop

## Align raw sequencing data to human reference genome
```
mkdir bams
bwa mem -t 1 -R '@RG\tID:1\tSM:PT1\tPL:PacBio' GRCh38/genome_chr17_0_10Mb.fa fastq/PT1_R1.fastq.gz fastq/PT1_R2.fastq.gz > bams/PT1.sam
```
