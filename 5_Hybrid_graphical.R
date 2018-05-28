
# datasets for two lineages
lin1 <- basic.stats(hap[, !colnames(hap) %in% lininfo$TA[lininfo$Lineage %in% c(2,3,4)]])
lin2 <- basic.stats(hap[, !colnames(hap) %in% lininfo$TA[lininfo$Lineage %in% c(1,3,4)]])

############################
### finding private snps ###
############################
## lin2 private alleles based on lin1 MAF=0
lin1.fix.allele <- data.frame('snp'=which(lin1$maf == 0), 
                              'lin1.fixed.allele'=find.fixed.allele(hap = lin1[lin1$maf == 0, ]), 
                              'tag'=lin1$`rs#`[lin1$maf == 0], stringsAsFactors = F)
  # lin2 private alleles
  lin2.private.allele <- find.private.allele(hap = lin2[lin1.fix.allele$snp, ], fixed = lin1.fix.allele$lin1.fixed.allele)
  # lin1 fixed and lin2 private alleles
  lin1.fix.allele <- data.frame(lin1.fix.allele, 
                                'lin2.private.allele'=lin2.private.allele, stringsAsFactors = F)
  lin1.fix.allele <- lin1.fix.allele[lin1.fix.allele$lin2.private.allele != 'N', ]

## lin1 private alleles based on lin2 MAF=0
lin2.fix.allele <- data.frame('snp'=which(lin2$maf == 0), 
                              'lin2.fixed.allele'=find.fixed.allele(hap = lin2[lin2$maf == 0, ]), 
                              'tag'=lin2$`rs#`[lin2$maf == 0], stringsAsFactors = F)
  # lin1 private alleles
  lin1.private.allele <- find.private.allele(hap = lin1[lin2.fix.allele$snp, ], fixed = lin2.fix.allele$lin2.fixed.allele)
  # lin2 fixed and lin1 private alleles
  lin2.fix.allele <- data.frame(lin2.fix.allele, 
                                'lin1.private.allele'=lin1.private.allele, stringsAsFactors = F)
  lin2.fix.allele <- lin2.fix.allele[lin2.fix.allele$lin1.private.allele != 'N', ]
  
# both datasets together
  lin1.lin2.private <- merge(lin1.fix.allele, lin2.fix.allele, by = 'snp', all = T)
  
## COLORING hybrids and wheat
  # getting x positions for bars in barplot
    chr.length <- read.table('data/chromLength.txt', header = T, as.is = T)
    chr.length <- cbind(chr.length, 'barplot.x'=as.numeric(barplot(chr.length$chrom, plot = F)))
    chr.length <- merge(chr.length, read.table('data/chrCentromeres.txt', header = T, as.is = T), by = 'chrom')

  # hybrids and wheat hap (a subset of hapfile)
    rows.to.keep <- lin1.lin2.private$snp
    cols.to.keep <- which(colnames(hap) %in% lininfo$TA[lininfo$Lineage %in% c(3,4)])
    hybrids <- as.character(lininfo$TA[lininfo$Lineage == 3])
    set.seed(1234)
    lin1.acc <- which(colnames(hap) %in% colnames(lin1)[sample(12:ncol(lin1), 1)])
    lin2.acc <- which(colnames(hap) %in% colnames(lin2)[sample(12:ncol(lin2), 1)])
    color.hybrids <- as.matrix(hap[rows.to.keep, c(1:11, lin1.acc, lin2.acc, cols.to.keep)])
    
  # colors
    color.hybrids[color.hybrids == lin1.lin2.private$lin1.private.allele] = 'red'
    color.hybrids[color.hybrids == lin1.lin2.private$lin2.private.allele] = 'blue'
    color.hybrids[color.hybrids %in% c('A','C','G','T','N','H')] = NA
    color.hybrids = as.data.frame(color.hybrids, stringsAsFactors = F)
    color.hybrids$pos = as.numeric(color.hybrids$pos)
    # remove rows where all NA genotypes
    color.hybrids = color.hybrids[rowSums(!is.na(color.hybrids[, -c(1:11)]), na.rm = T) > 0, ]
    # merge with chrom info
    color.hybrids = merge(chr.length, color.hybrids, by = 'chrom')
    
  # get wheat consensus
    wheat.cols <- which(colnames(color.hybrids) %in% wheat.lines)
    hybrid.cols <- which(colnames(color.hybrids) %in% hybrids)
    color.hybrids <- data.frame('lineage1' = rowSums(color.hybrids[, wheat.cols]=='red', na.rm = T), 
                                'lineage2' = rowSums(color.hybrids[, wheat.cols]=='blue', na.rm = T),
                                'Wheat.D.genome.consensus' = 'blue',
                                color.hybrids, stringsAsFactors = F, check.names = F)
    # update wheat consensus color column
    color.hybrids$Wheat.D.genome.consensus[color.hybrids$lineage1 > 0] = 'red'
    color.hybrids$Wheat.D.genome.consensus[(color.hybrids$lineage1 + color.hybrids$lineage2) == 0] = NA

    color.hybrids <- color.hybrids[, c(3:6, 10, 18:19, which(colnames(color.hybrids) %in% hybrids))]

# plotting perfect hybrid and wheat consensus
alpha.val = 0.75
pdf('output/Fig.6_TA3429_Wheat.pdf', width = 11, height = 4.25)
par(mfrow=c(1,2))
  for (i in 1:ncol(color.hybrids)) { # loop to plot chromosomes
    if(colnames(color.hybrids)[i] %in% c('TA3429', 'Wheat.D.genome.consensus')) plot.hybrid(i = i)
  }
dev.off()

# paint chromosomes in hybrids
pdf('output/S7_Hybrids.pdf', width = 11, height = 11)
par(mfrow=c(3,2))
for (i in 6:ncol(color.hybrids)) { # loop to plot chromosomes
  if(colnames(color.hybrids)[i] != 'TA3429') plot.hybrid(i = i, cex = 1.35)
}
dev.off()

# distribution of lineage-specific allele on wheat chromosomes
table(color.hybrids[, 1:2])

## identity matrix with private allele SNPs
source('R/allele.match.R')

lin1.id <- allele.match(hap = hap.read(hap[rows.to.keep, !colnames(hap) %in% lininfo$TA[lininfo$Lineage %in% c(2,4)]], 
                                       data.col = 12))
lin2.id <- allele.match(hap = hap.read(hap[rows.to.keep, !colnames(hap) %in% lininfo$TA[lininfo$Lineage %in% c(1,4)]], 
                                       data.col = 12))

save(lin1.id, file = 'data/lin1.id.RData')
save(lin2.id, file = 'data/lin2.id.RData')
load('data/lin1.id.RData'); load('data/lin2.id.RData')

lin1.id[lin1.id > 1] = NA; diag(lin1.id) = NA
lin1.id[lower.tri(lin1.id)] = t(lin1.id)[lower.tri(lin1.id)]
lin1.id <- lin1.id[-which(rownames(lin1.id) %in% hybrids), colnames(lin1.id) %in% hybrids]

lin2.id[lin2.id > 1] = NA; diag(lin2.id) = NA
lin2.id[lower.tri(lin2.id)] = t(lin2.id)[lower.tri(lin2.id)]
lin2.id <- lin2.id[-which(rownames(lin2.id) %in% hybrids), colnames(lin2.id) %in% hybrids]

# hybrids parents
hybrid.parents <- data.frame('Putative hybrid' = hybrids, 
                             'Putative L1 parent'=NA, 
                             'Putative L2 parent'=NA)
all(colnames(lin1.id) == colnames(lin2.id))

for (i in 1:nrow(hybrid.parents)) {
  hybrid = as.character(hybrid.parents[i,1])
  lin1.id.col <- as.data.frame(lin1.id[, colnames(lin1.id) %in% hybrid])
  l1.parent = which(lin1.id.col == max(lin1.id.col, na.rm = T))
  hybrid.parents[i,2] = paste0(rownames(lin1.id)[l1.parent], ' (', round(lin1.id.col[l1.parent,1]*100, 2), '%)')
  lin2.id.col <- as.data.frame(lin2.id[, colnames(lin2.id) %in% hybrid])
  l2.parent = which(lin2.id.col == max(lin2.id.col, na.rm = T))
  hybrid.parents[i,3] = paste0(rownames(lin2.id)[l2.parent], ' (', round(lin2.id.col[l2.parent,1]*100, 2), '%)')
}

write.csv(hybrid.parents, file = 'output/hybrid_parents.csv', quote = F, row.names = F)

# hybrids related to each other
id.mat <- allele.match(hap = hap.read(hap.obj = hap[, !colnames(hap) %in% lininfo$TA[lininfo$Lineage %in% c(1,2,4)]], 
                                      data.col = 12))
identical.accessions(t(id.mat), threshold = 0.98)
##############
