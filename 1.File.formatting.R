## ################ ##
## File formatting  ##
## Narinder Singh   ##
## ################ ##

# load functions and required libraries
source('functions.R')
req.libraries <- c('data.table', 'ape', 'phyclust', 'rrBLUP', 'scatterplot3d', 'rgl', 'ggplot2', 
                   'corrplot', 'plyr', 'multtest', 'gplots', 'LDheatmap', 'genetics', 'compiler',
                   'mapdata', 'maps', 'ggmap', 'snpStats')
if (!require("pacman")) install.packages("pacman")
pacman::p_load(req.libraries)

if(!file.exists('project.log')) file.create('project.log')

# read and format the data file
hap <- fread("data/Ae_tauschiiDiversity.hmp.txt", header = T, check.names = F, data.table = F)
table(hap$chrom)
hap <- hap[grep('NWV', hap$chrom, invert = T), ]
table(hap$chrom)

# update chromosome designation
chrom.designation <- read.table(file = 'data/pseudoChrom2chrom.txt', header = T)
for (i in 1:nrow(chrom.designation)) {
  hap$chrom[hap$chrom == chrom.designation$pseudoChrom[i]] = chrom.designation$chrom[i]
}
table(hap$chrom)

# retain only biallelic SNP sites
biallelic = sapply(strsplit(as.character(hap$alleles), '/'), length) == 2
sum(!biallelic)
hap = hap[biallelic, ]

# replace het IUPAC codes with H and N
missAmbiguous <- c('0', '+', '-')
hetCodes <- c('R','Y','S','W','K','M','B','D','H','V')

   hapgeno <- as.matrix(hap[, 12:ncol(hap)])
   hapgeno[hapgeno %in% missAmbiguous] = 'N'
   hapgeno[hapgeno %in% hetCodes] = 'H'
   sort(unique(c(hapgeno)))

hap=cbind(hap[,1:11], hapgeno)
rm(biallelic, hapgeno, hetCodes, missAmbiguous)

# using empty columns
colnames(hap)[5:11]=c('alleleA', 'alleleB', 'het', 'het_prop', 'missing', 'propmiss', 'maf')

# removing known wrong sample
hap <- hap[, !colnames(hap) %in% remove.samples]

# writing hap file in text file format
write.table(hap, file="data/hapFile.txt" , col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

# write sample names to a file
if(!file.exists('output')) dir.create('output')
initial.sample.names <- cbind('TA'=sort(colnames(hap)[-c(1:11)]), 'Sample'=NA)
write.table(initial.sample.names, file = "output/initial.sample.names.txt", row.names = F, quote = F, sep = '\t')

#########
## END ##
#########