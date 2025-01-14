# Welcome to Naxerova Lab's Cancer Genomics Workshop!

## 1. Install miniconda on your computer

(See Google Slides presentation slide for instructions)

## 2. Install prerequisite software into a new conda environment

Open a Terminal (Go -> Utilities -> Terminal), then run the following commands:

```
# After conda is installed, create a new conda environment for this workshop:
conda create --name naxerova_lab
conda activate naxerova_lab
conda config --add channels bioconda

# if you have an older intel chip mac:
conda install gatk4 igv samtools git bwa

# if you have a newer M1/M2 mac:
conda install gatk4 igv git
conda install nanoporetech::samtools nanoporetech::bwa
```

## 3. Clone this GitHub repo to access raw data and scripts

```
# clone the repo to your computer and change directories to this location
git clone https://github.com/agorelick/cancer_genomics_workshop.git; cd cancer_genomics_workshop
```

## 4. Examine the raw sequencing data (.fastq files)

1. What does raw Next Generation Sequencing data look like? 

```
# view the first 10 lines of the compressed data
head -10 fastq/N_R1.fastq.gz

# view the first 10 lines of uncompressed data
cat fastq/N_R1.fastq.gz | gunzip | head -10
```

2. What do we actually **want to know** from this data?
3. What do we need to do to make this data _interpretable_?

## 5. Examine the reference genome

1. What is a reference genome?
2. Why do we need it?
3. Are there potential downsides to using a reference genome?
4. Do the downsides affect us if we are calling _somatic_ mutations?


## 6. Align raw sequencing data to human reference genome
```
mkdir bams
bwa mem -t 1 -R '@RG\tID:9\tSM:N' GRCh38/genome_chr17_0_10Mb fastq/N_R1.fastq.gz fastq/N_R2.fastq.gz > bams/N.sam
bwa mem -t 1 -R '@RG\tID:1\tSM:PT1' GRCh38/genome_chr17_0_10Mb fastq/PT1_R1.fastq.gz fastq/PT1_R2.fastq.gz > bams/PT1.sam
bwa mem -t 1 -R '@RG\tID:2\tSM:PT2' GRCh38/genome_chr17_0_10Mb fastq/PT2_R1.fastq.gz fastq/PT2_R2.fastq.gz > bams/PT2.sam
bwa mem -t 1 -R '@RG\tID:3\tSM:PT3' GRCh38/genome_chr17_0_10Mb fastq/PT3_R1.fastq.gz fastq/PT3_R2.fastq.gz > bams/PT3.sam
bwa mem -t 1 -R '@RG\tID:4\tSM:LN1' GRCh38/genome_chr17_0_10Mb fastq/LN1_R1.fastq.gz fastq/LN1_R2.fastq.gz > bams/LN1.sam
bwa mem -t 1 -R '@RG\tID:5\tSM:Liv1' GRCh38/genome_chr17_0_10Mb fastq/Liv1_R1.fastq.gz fastq/Liv1_R2.fastq.gz > bams/Liv1.sam
bwa mem -t 1 -R '@RG\tID:6\tSM:Lun1' GRCh38/genome_chr17_0_10Mb fastq/Lun1_R1.fastq.gz fastq/Lun1_R2.fastq.gz > bams/Lun1.sam
bwa mem -t 1 -R '@RG\tID:7\tSM:Lun2' GRCh38/genome_chr17_0_10Mb fastq/Lun2_R1.fastq.gz fastq/Lun2_R2.fastq.gz > bams/Lun2.sam
bwa mem -t 1 -R '@RG\tID:8\tSM:Lun3' GRCh38/genome_chr17_0_10Mb fastq/Lun3_R1.fastq.gz fastq/Lun3_R2.fastq.gz > bams/Lun3.sam
```

1. What is the difference between a SAM file and a FASTQ file? (see: https://en.wikipedia.org/wiki/SAM_(file_format))
2. Are SAM files compressed?
3. Are the reads sorted?
4. What does the CIGAR field tell us? (https://en.wikipedia.org/wiki/Sequence_alignment#CIGAR_Format) 
5. Can you find any reads that differ from the reference genome? Can we tell if these are biological differences (e.g. mutations) or technical (e.g. noise)?


## 7. Format aligned reads so that variant callers can use them
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


## 8. Look at aligned sequencing data on IGV

```
# open IGV (Integrative Genomics Viewer) from the command line
igv
```

### Instructions
        1. Select "Human hg38" as the reference genome
        2. Open the patient's normal bam file (N_sorted.bam) on IGV (File -> Load From File -> Navigate to bams/N_sorted.bam)
        3. Go to gene TP53 (Type "TP53" in the search bar and click "Go"), then zoom in to view the reads.

1. Can you find any differences in the patient's genome from the reference genome (e.g. single-nucleotide polymorphisms)? Are they heterozygous or homozygous?
2. Can you find any somatic mutations? Can we tell if they likely occured either early or late in the patient's cancer?
3. Based on IGV (if you had to guess) do either the patient's lung or liver metastases seem more closely related to the locoregional lymph node (LN) metastasis?


## 9. Call somatic mutations in paired tumor/normal mode (~15min)
```
gatk Mutect2 -R GRCh38/genome_chr17_0_10Mb.fa \
        -I bams/PT1_sorted.bam \
        -I bams/PT2_sorted.bam \
        -I bams/PT3_sorted.bam \
        -I bams/LN1_sorted.bam \
        -I bams/Liv1_sorted.bam \
        -I bams/Lun1_sorted.bam \
        -I bams/Lun2_sorted.bam \
        -I bams/Lun3_sorted.bam \
        -I bams/N_sorted.bam \
        -L chr17:7000000-8000000 \
        -normal "N" \
        -O unfiltered.vcf

# quick look at the output
less -RNS unfiltered.vcf
```
1. What does a VCF file show?
2. Where can you find the number of reads supporting each mutation's reference and alternate allele in each sample?
3. Can you see any mutations that look like they are real (true positives) or artifacts? What information might help?


## 10. Apply filters to try to remove false-positive mutations
```
gatk FilterMutectCalls -R GRCh38/genome_chr17_0_10Mb.fa -V unfiltered.vcf -O filtered.vcf

# quick look at the output
less -RNS filtered.vcf
```
1. What are some of the filters to flag potential artifacts/false positive mutations?


## 11. Make heatmap and phylogenetic tree

Start R-studio, then let's step through the script `make_heatmap_and_tree.R`


## Optional: uninstall Conda once we're done

conda install anaconda-clean
anaconda-clean --yes
rm -rf ~/anaconda3







