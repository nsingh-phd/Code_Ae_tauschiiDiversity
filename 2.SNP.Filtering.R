
# read the data file
hap <- fread("data/hapFile.txt", header = T, check.names = F, data.table = F)

# compute basic stats - missing data, heterozygosity, MAF and allele count
hap <- basic.stats(hap = hap)

# some statistics
cat('Missing data range =', round(range(hap$propmiss), digits = 2))
cat('# SNPs having less than 20% missing data =', sum(hap$propmiss < 0.2))
cat('Heterozygosity range =', round(range(hap$het_prop), digits = 2))

# remove bad snps
hap <- hap[hap$propmiss < 0.2 & hap$het_prop < 0.05 & hap$maf >= 0.01, ]

# fisher exact test
p.values <- c(apply(cbind(hap$het, hap$alleleB, hap$alleleA, hap$missing), 1, fisher.p.values))

# filtering snps based on 
hap <- hap[p.values < 0.001/length(p.values), ]

# computing per sample heterozygosity and missing data
per.sample.het <- colSums(hap == 'H', na.rm = T) / nrow(hap)
per.sample.miss <- colSums(hap == 'N', na.rm = T) / nrow(hap)
cat('Per sample missing data range =', round(range(per.sample.miss), digits = 2))

  # per sample heterozygosity for
  cat('Per sample het range (wheat) =', round(range(per.sample.het[names(per.sample.het) %in% tauschii.lines]), digits = 2))  # tauschii
  cat('Per sample het range (wheat) =', round(range(per.sample.het[names(per.sample.het) %in% wheat.lines]), digits = 2))  # wheat
  
# removing wheat specific SNPs
snp.levels <- check.snp.levels(hap = hap[, !colnames(hap) %in% wheat.lines])
table(snp.levels)   # number of monomorphic and polymorphic SNPs
hap <- hap[snp.levels == 2, ] # removing wheat specific SNPs

# removing bad samples
hap <- hap[, per.sample.het <= 0.05 & per.sample.miss <= 0.8]

# write hapfile
write.table(hap, file="data/hapFile.filtered.txt", col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
# hap <- fread("data/hapFile.filtered.txt", header = T, check.names = F, data.table = F)

# write out sample names
write.table(colnames(hap)[-c(1:11)], file = "output/samples_after_filtering_all.txt", col.names = F, row.names = F, quote = F)

## convert major allele to -1 and minor allele to 1, hets to 0
hap01 = hap2numeric(hap = hap)

## writing hap01 file into a text format
write.table(hap01, file="data/hap01.filtered.txt", col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
# hap01 <- fread("data/hap01.filtered.txt", header = T, check.names = F, data.table = F)

## STRUCTURE input for all samples ##
write.table(file4faststructure(hap = hap), 
            file = "data/str_all.f.str", sep = '\t', quote = F, row.names = F, col.names = F)

## STRUCTURE input for tauschii samples ##
write.table(file4faststructure(hap = hap[, !colnames(hap) %in% wheat.lines]), 
            file = "data/str_tauschii.f.str", sep = '\t', quote = F, row.names = F, col.names = F)

#################