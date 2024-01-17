## install ape package if it is not installed
if(!'ape' %in%  rownames(installed.packages())) {
    options(install.packages='always')
    install.packages('ape')
}

library(ape)

## load data
d <- read.csv('filtered.vcf', comment.char='#', sep='\t', header=F)

## extract field names from the header
txt <- read.csv('filtered.vcf', sep='\n', header=F)[[1]]
header <- strsplit(grep('#CHROM', txt, value=T),'\t')[[1]]
header[1] <- 'CHROM'
names(d) <- header

## subset for PASS mutations
d <- d[d$FILTER=='PASS',]
d$ID <- paste0(d$CHROM,'.',d$REF,d$POS,d$ALT)
res <- d[,c('ID','PT1','PT2','PT3','LN1','Liv1','Lun1','Lun2','Lun3')]

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
vaf <- alt_counts / (alt_counts + ref_counts)
rownames(vaf) <- res$ID

## remove any variants with insufficient read support (possible artifacts)
samples_with_evidence <- rowSums(alt_counts >= 3)
bad_variants <- samples_with_evidence==0
vaf <- vaf[-bad_variants,]
vaf <- cbind(vaf, N=0)

## save the VAF matrix
write.table(vaf, file = 'vaf_matrix.txt', sep = "\t", quote = FALSE, col.names = NA)

## make a heatmap from the VAFs for each sample/variant
rownames(vaf) <- rep('', nrow(vaf))
pdf('heatmap.pdf',width=9, height=7)
heatmap(t(vaf))
dev.off()

## make a NJ tree based on the true CCF values
dm <- dist(t(vaf), method='euclidean')
tree <- nj(dm)
tree <- root(tree, outgroup='N')
pdf('tree.pdf',width=10, height=8)
plot(tree, type='u')
dev.off()



