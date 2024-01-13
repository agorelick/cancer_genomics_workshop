# cancer_genomics_workshop

## Clone this GitHub repo and install prerequisites

## Examine the raw sequencing data (.fastq files)

1. What does raw Next Generation Sequencing data look like? 
2. What do we actually **want to know** from this data?
3. What do we need to do to make this data _interpretable_?

## Examine the reference genome

## Align raw sequencing data to human reference genome
```
mkdir bams
bwa mem -t 1 -R '@RG\tID:9\tSM:N\tPL:PacBio' GRCh38/genome_chr17_0_10Mb fastq/N_R1.fastq.gz fastq/N_R2.fastq.gz > bams/N.sam
bwa mem -t 1 -R '@RG\tID:1\tSM:PT1\tPL:PacBio' GRCh38/genome_chr17_0_10Mb fastq/PT1_R1.fastq.gz fastq/PT1_R2.fastq.gz > bams/PT1.sam
bwa mem -t 1 -R '@RG\tID:2\tSM:PT2\tPL:PacBio' GRCh38/genome_chr17_0_10Mb fastq/PT2_R1.fastq.gz fastq/PT2_R2.fastq.gz > bams/PT2.sam
bwa mem -t 1 -R '@RG\tID:3\tSM:PT3\tPL:PacBio' GRCh38/genome_chr17_0_10Mb fastq/PT3_R1.fastq.gz fastq/PT3_R2.fastq.gz > bams/PT3.sam
bwa mem -t 1 -R '@RG\tID:4\tSM:LN1\tPL:PacBio' GRCh38/genome_chr17_0_10Mb fastq/LN1_R1.fastq.gz fastq/LN1_R2.fastq.gz > bams/LN1.sam
bwa mem -t 1 -R '@RG\tID:5\tSM:Liv1\tPL:PacBio' GRCh38/genome_chr17_0_10Mb fastq/Liv1_R1.fastq.gz fastq/Liv1_R2.fastq.gz > bams/Liv1.sam
bwa mem -t 1 -R '@RG\tID:6\tSM:Lun1\tPL:PacBio' GRCh38/genome_chr17_0_10Mb fastq/Lun1_R1.fastq.gz fastq/Lun1_R2.fastq.gz > bams/Lun1.sam
bwa mem -t 1 -R '@RG\tID:7\tSM:Lun2\tPL:PacBio' GRCh38/genome_chr17_0_10Mb fastq/Lun2_R1.fastq.gz fastq/Lun2_R2.fastq.gz > bams/Lun2.sam
bwa mem -t 1 -R '@RG\tID:8\tSM:Lun3\tPL:PacBio' GRCh38/genome_chr17_0_10Mb fastq/Lun3_R1.fastq.gz fastq/Lun3_R2.fastq.gz > bams/Lun3.sam
```

## Format aligned reads so that variant callers can use them
```
samtools view -hb bams/N.sam > bams/N.bam; samtools sort bams/N.bam > bams/N_sorted.bam; samtools index bams/N_sorted.bam
samtools view -hb bams/PT1.sam > bams/PT1.bam; samtools sort bams/PT1.bam > bams/PT1_sorted.bam; samtools index bams/PT1_sorted.bam
samtools view -hb bams/PT2.sam > bams/PT2.bam; samtools sort bams/PT2.bam > bams/PT2_sorted.bam; samtools index bams/PT2_sorted.bam
samtools view -hb bams/PT3.sam > bams/PT3.bam; samtools sort bams/PT3.bam > bams/PT3_sorted.bam; samtools index bams/PT3_sorted.bam
samtools view -hb bams/LN1.sam > bams/LN1.bam; samtools sort bams/LN1.bam > bams/LN1_sorted.bam; samtools index bams/LN1_sorted.bam
samtools view -hb bams/Liv1.sam > bams/Liv1.bam; samtools sort bams/Liv1.bam > bams/Liv1_sorted.bam; samtools index bams/Liv1_sorted.bam
samtools view -hb bams/Lun1.sam > bams/Lun1.bam; samtools sort bams/Lun1.bam > bams/Lun1_sorted.bam; samtools index bams/Lun1_sorted.bam
samtools view -hb bams/Lun2.sam > bams/Lun2.bam; samtools sort bams/Lun2.bam > bams/Lun2_sorted.bam; samtools index bams/Lun2_sorted.bam
samtools view -hb bams/Lun3.sam > bams/Lun3.bam; samtools sort bams/Lun3.bam > bams/Lun3_sorted.bam; samtools index bams/Lun3_sorted.bam
```






