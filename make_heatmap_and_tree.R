library(pheatmap) # for heatmap 
library(ape) # for neighbor-joining tree
library(phytools) # for rooting tree at the germline/normal sample
library(ggplot2) # plotting
library(ggtree) # plotting

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# load mutation data from the VCF file
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## load data
d <- read.csv('filtered.vcf', comment.char='#', sep='\t', header=F)

## extract field names from the header
txt <- read.csv('filtered.vcf', sep='\n', header=F)[[1]]
header <- strsplit(grep('#CHROM', txt, value=T),'\t')[[1]]
header[1] <- 'CHROM'
names(d) <- header

## subset for mutations that passed all QC filters
d <- d[d$FILTER=='PASS',]


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# extract REF/ALT read counts for each variant in each 
# sample from the VCF
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d$ID <- paste0(d$CHROM,'.',d$REF,d$POS,d$ALT)
res <- d[,c('ID','PT1','PT2','PT3','LN1','Liv1','Lun1','Lun2','Lun3')]

## extract the number of ALT/REF reads in each sample
extract_read_counts <- function(string) {
    AD <- sapply(string, function(x) strsplit(x, '[:]')[[1]][2])
    ref_count <- sapply(AD, function(x) strsplit(x, '[,]')[[1]][1])
    alt_count <- sapply(AD, function(x) strsplit(x, '[,]')[[1]][2])
    ref_count <- as.integer(ref_count)
    alt_count <- as.integer(alt_count)
    data.frame(ref=ref_count, alt=alt_count)
}

AD_PT1 <- extract_read_counts(res$PT1)
AD_PT2 <- extract_read_counts(res$PT2)
AD_PT3 <- extract_read_counts(res$PT3)
AD_LN1 <- extract_read_counts(res$LN1)
AD_Liv1 <- extract_read_counts(res$Liv1)
AD_Lun1 <- extract_read_counts(res$Lun1)
AD_Lun2 <- extract_read_counts(res$Lun2)
AD_Lun3 <- extract_read_counts(res$Lun3)

## create a matrix of VAFs for all samples/mutations
alt_counts<-cbind(PT1=AD_PT1$alt, PT2=AD_PT2$alt, PT3=AD_PT3$alt, LN1=AD_LN1$alt, Liv1=AD_Liv1$alt, Lun1=AD_Lun1$alt, Lun2=AD_Lun2$alt, Lun3=AD_Lun3$alt)
ref_counts<-cbind(PT1=AD_PT1$ref, PT2=AD_PT2$ref, PT3=AD_PT3$ref, LN1=AD_LN1$ref, Liv1=AD_Liv1$ref, Lun1=AD_Lun1$ref, Lun2=AD_Lun2$ref, Lun3=AD_Lun3$ref)
total_counts <- alt_counts + ref_counts
vaf <- alt_counts / total_counts
rownames(vaf) <- res$ID

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# remove any variants with insufficient read support 
# (possible false positives due to seq. error)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# first remove any mutations without at least 20 read coverage in all samples
number_of_samples_with_good_coverage <- rowSums(total_counts >= 20) 
bad_variants <- number_of_samples_with_good_coverage < ncol(vaf)
vaf <- vaf[-bad_variants,]
total_counts <- total_counts[-bad_variants,]
alt_counts <- alt_counts[-bad_variants,]
ref_counts <- ref_counts[-bad_variants,]

# next remove any mutations not mutated with at least 5% VAF and 3 ALT reads in at least 1 sample
samples_with_evidence <- rowSums(alt_counts >= 3 & vaf >= 0.05)
bad_variants <- samples_with_evidence < 1
vaf <- vaf[-bad_variants,]
total_counts <- total_counts[-bad_variants,]
alt_counts <- alt_counts[-bad_variants,]
ref_counts <- ref_counts[-bad_variants,]

# add germline to the VAF data (with VAF=0 for all mutations)
vaf <- cbind(vaf, germline=0)

## save the VAF matrix
write.table(vaf, file = 'vaf_matrix.txt', sep = "\t", quote = FALSE, col.names = NA)


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot a heatmap of the clustered samples and mutations
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## make a heatmap from the VAFs for each sample/variant
rownames(vaf) <- rep('', nrow(vaf))
pheatmap(t(vaf),file='heatmap.png',width=8, height=5.5)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# generate cancer evolution tree
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## make a NJ tree based on the mutations' VAFs
dm <- dist(t(vaf), method='euclidean')
tree <- nj(dm)
tree <- phytools::reroot(tree, which(tree$tip.label=='germline'))

groups <- rbind(data.frame(label=grep('^germline',tree$tip.label,value=T), group='Normal'),
                data.frame(label=grep('^PT',tree$tip.label,value=T), group='Primary'),
                data.frame(label=grep('^Liv',tree$tip.label,value=T), group='Liver'),
                data.frame(label=grep('^Lun',tree$tip.label,value=T), group='Lung'),
                data.frame(label=grep('^LN',tree$tip.label,value=T), group='Lymph node'))

cols <- c('black','forestgreen','steelblue','magenta','orange')
names(cols) <- c('Normal','Primary','Liver','Lung','Lymph node')
p <- ggtree(tree, layout='ape') 
p <- p %<+% groups
p <- p + geom_tiplab(aes(color=group), angle=0) + theme(legend.position='none')
p <- p + xlim(min(p$data$x-1.0), max(p$data$x+3.0)) + 
    scale_color_manual(values=cols, name='Sample type') 
ggsave('tree.png')



