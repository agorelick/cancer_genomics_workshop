# Welcome to Naxerova Lab's Cancer Genomics Workshop!

## Install Prequisite software

### Instructions for Mac

1. Open a terminal on your computer: `Finder > Go > Utilities > Terminal`
2. clone this GitHub repo to your computer
   * In terminal, type `git clone https://github.com/agorelick/cancer_genomics_workshop.git; cd cancer_genomics_workshop` and hit enter.
3. Install Conda on your computer following the instructions here: https://docs.anaconda.com/miniconda/install/#quick-command-line-install
   * Older Intel-chip Macs: Select `MacOS > Intel`
   * Newer Apple Silicon Macs: Select `MacOS > Apple Silicon`
   * Copy and paste the commands into your terminal and hit enter.
   * Close the terminal and open a new terminal. You should see the text `(base)` in the prompt at the bottom.
4. Use Conda to install all required software into a contained environment called *naxerova_workshop:
   ```
   conda env create -f environment.yml
   ```
6. Activate the naxerova_workshop Conda environment: `conda activate naxerova_workshop`. You should see the text `(naxerova_workshop)` at the bottom.
7. Download IGV (Integrative Genomics Viewer): https://igv.org/doc/desktop/#DownloadPage/
   * Select either `IGV For MacOS (Apple Chip) - Java Included` or `IGV For MacOS (Intel Chip) - Java Included`
   

### Instructions for Windows

You will need to install a linux terminal within your windows computer. Windows 10+ has a built-in utility for this called WSL (Windows Subsystem for Linux). To use enable WSL, follow the instructions here. *this may take a few minutes and will require restarting.* 
1. Follow these instructions to install a linux command-line on your windows computer: https://www.geeksforgeeks.org/how-to-install-wsl2-windows-subsystem-for-linux-2-on-windows-10/
2. Restart your computer
3. Search for WSL to start your linux terminal. Choose a username and password (something easy to remember).
4. Install Conda on your **linux terminal** following the instructions here: https://docs.anaconda.com/miniconda/install/#quick-command-line-install
   * Select `Linux > 64-bit`
   * Copy and paste the commands into your terminal and hit enter.
   * Close the terminal and open a new terminal. You should see the text `(base)` in the prompt at the bottom.
5. Use Conda to install all required software into a contained environment called *naxerova_workshop:
   ```
   conda env create -f environment.yml
   ```
7. Activate the naxerova_workshop Conda environment: `conda activate naxerova_workshop`. You should see the text `(naxerova_workshop)` at the bottom.
8. Download IGV (Integrative Genomics Viewer): https://igv.org/doc/desktop/#DownloadPage/
   * Select `IGV For Windows - Java Included`

## 1. Examine the raw sequencing data (.fastq files)

1. What does raw Next Generation Sequencing data look like? 

```
# view the first 10 lines of compressed data
head -10 fastq/N_R1.fastq.gz

# view the first 10 lines of uncompressed data
cat fastq/N_R1.fastq.gz | gunzip | head -10
```

2. What do we actually **want to know** from this data?
3. What do we need to do to make this data _interpretable_?

## 2. Align raw sequencing data to human reference genome

1. What is a reference genome?
2. Why do we need it?
3. Are there potential downsides to using a reference genome?
4. Do the downsides affect us if we are calling _somatic_ mutations?
   
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

5. What is the difference between a SAM file and a FASTQ file? (see: https://en.wikipedia.org/wiki/SAM_(file_format))
6. Are SAM files compressed?
7. Are the reads sorted?
8. What does the CIGAR field tell us? (https://en.wikipedia.org/wiki/Sequence_alignment#CIGAR_Format) 
9. Can you find any reads that differ from the reference genome? Can we tell if these are biological differences (e.g. mutations) or technical (e.g. noise)?


## 4. Format aligned reads so that variant callers can use them
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


## 5. Look at aligned sequencing data on IGV

### Instructions
1. Open IGV
    * Mac: Extract the IGV application .zip file and run it.
    * Windows: open it from the installation
3. Select "Human (GRCh38/hg38)" as the reference genome
4. Load the patient's normal bam (N_sorted.bam) onto IGV (Go to: File -> Load From File -> Navigate to bams/N_sorted.bam)
    * Windows: Home > cancer_genomics_workshop > bams > N_sorted.bam
6. Go to gene TP53 (Type "TP53" in the search bar and click "Go"), then zoom in to view the reads.
    * Can you find any differences in the patient's genome from the reference genome (e.g. single-nucleotide polymorphisms, SNPs)? Are they heterozygous or homozygous?
7. Without removing the N1 bam file, add a primary tumor bam file (PT1_sorted.bam) onto IGV (Go to: File -> Load From File -> Navigate to bams/N_sorted.bam)
    * Can you find any somatic mutations?
8. Without removing N1 or PT1, add *sorted* bam files for LN1, Liv1, Lung1 onto IGV. (Right click on reads and select "squished" to see all reads.)
    * Can you find somatic mutations that arose early in the patient's cancer? Late in the cancer? How do we know?
    * Based on IGV (if you had to guess) are either the patient's lung or liver metastases seem more closely related to the locoregional lymph node (LN) metastasis?
    * Can you find any mutations present in multiple metastasis samples and absent in the primary tumor sample? What scenarios in cancer evolution can explain this (multiple answers!)


## 6. Call somatic mutations in paired multi-sample tumor/normal mode (~15min)
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

# quick look at the output (scroll with keyboard arrows; q to exit)
less -RNS unfiltered.vcf
```
1. What does a VCF file show?
2. Where can you find the number of reads supporting each mutation's reference and alternate allele in each sample?
3. Can you see any mutations that look like they are real (true positives)? What about false-positives/artifacts? What information might help?


## 7. Apply filters to flag potential false-positive mutations
```
gatk FilterMutectCalls -R GRCh38/genome_chr17_0_10Mb.fa -V unfiltered.vcf -O filtered.vcf

# quick look at the output (scroll with keyboard arrows; q to exit)
less -RNS filtered.vcf
```
1. What are some filters to flag potential artifacts/false positive mutations? Why might they be a common source of false positives in mutation calling?


## 8. Make heatmap and phylogenetic tree

1. Let's look through the included R script `make_heatmap_and_tree.R`:
```
# quick look at the R script (scroll with keyboard arrows; q to exit)
less make_heatmap_and_tree.R
```
2. Run the R script to create a heatmap and phylogenetic tree for the tumor samples, then open the PDFs on your computer. 
```
Rscript make_heatmap_and_tree.R
```



## Optional: Revert your computer

### Remove the Conda environment and all newly installed software
```
conda deactivate
conda remove --name naxerova_workshop --all -y
```
### Remove Conda from your computer
```
~/miniconda3/uninstall.sh
```

### (Mac) Remove IGV
Move the IGV .zip file and extracted application to the Trash.

### (Windows) Remove IGV
Uninstall IGV: Start > Search for IGV > Open File Location > Run uninstaller

### (Windows) Remove WSL terminal from your computer
Completely remove WSL with following these instructions (select your Windows version): 
* Windows 10: https://medium.com/@bonguides25/how-to-completely-uninstall-the-subsystem-for-linux-on-windows-10-20c5c1377117
* Windows 11: https://www.elevenforum.com/t/uninstall-windows-subsystem-for-linux-wsl-distro-in-windows-11.12250/




