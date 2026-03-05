### added 'row.names=1' to each time we read in file
me <- read.table("Raw_Data/AncestryDNA_LI.txt", sep="\t", header=TRUE,row.names=1)
mom <- read.table("Raw_Data/AncestryDNA_BM.txt", sep="\t", header=TRUE,row.names=1)
dad <- read.table("Raw_Data/AncestryDNA_DI.txt", sep="\t", header=TRUE,row.names=1)

pgf <- read.table("Raw_Data/AncestryDNA_FI.txt", sep="\t", header=TRUE,row.names=1)
mgm <- read.table("Raw_Data/AncestryDNA_MH.txt", sep="\t", header=TRUE,row.names=1)
mgf <- read.table("Raw_Data/AncestryDNA_RM.txt", sep="\t", header=TRUE,row.names=1)


gwas1 <- read.csv('Raw_Data/gwas1.csv', header = TRUE)
gwas2 <- read.csv('Raw_Data/gwas2.csv', header = TRUE)
gwas <- rbind(gwas1, gwas2)
rm(gwas1, gwas2)



### Look into this
chr21 <- get_phased(chrom = 21, mom = mom, dad = dad, mgm = mgm, mgf = mgf, me = me)

chr21_traits <- gwas[which(rownames(chr21) %in% gwas$SNPS),]

## 
table(chr21_traits$DISEASE.TRAIT)


snp <- rownames(chr21)[28]

### This shows you the my genotype
chr21[snp,]

### This shows you the traits associated with that snp
gwas[ gwas$SNPS == snp,]