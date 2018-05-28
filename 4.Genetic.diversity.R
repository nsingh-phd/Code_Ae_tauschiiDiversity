
lin1 <- basic.stats(hap[, !colnames(hap) %in% lininfo$TA[lininfo$Lineage %in% c(2,3,4)]])
lin2 <- basic.stats(hap[, !colnames(hap) %in% lininfo$TA[lininfo$Lineage %in% c(1,3,4)]])
mc <- read.table('data/mc.txt', header = T, as.is = T)
mc.hap <- basic.stats(hap[, c(1:11, which(colnames(hap) %in% mc$TA))])
# L1
nei.lin1 <- nei.index(lin1)
# L2
nei.lin2 <- nei.index(lin2)
# possible hybrids
nei.hybrids <- nei.index(hap = basic.stats(hap[, !colnames(hap) %in% lininfo$TA[lininfo$Lineage %in% c(1,2,4)]]))
# Wheat
nei.wheat <- nei.index(hap = basic.stats(hap[, !colnames(hap) %in% lininfo$TA[lininfo$Lineage %in% c(1,2,3)]]))
# all tauschii
nei.tauschii <- nei.index(hap = basic.stats(hap[, !colnames(hap) %in% lininfo$TA[lininfo$Lineage %in% c(4)]]))
# whole dataset including wheat
nei.all.samples <- nei.index(hap = basic.stats(hap))
# minicore
nei.mc <- nei.index(hap = mc.hap)
  # segregating snps
  sum(mc.hap$maf > 0)

# MAF L1 and L2, and joint
pdf('output/S8_MAF-spectrum.pdf', width = 12, height = 4)
par(mfrow=c(1,3))
hist((basic.stats(hap[, !colnames(hap) %in% lininfo$TA[lininfo$Lineage %in% c(2,3,4)]]))$maf,
     breaks = 30, col = 'black', xlab = 'Minor allele frequency', ylab = '# SNPs', main = '(A). Lineage1 (L1)', 
     cex.main = 2, cex.lab = 1.5)
hist((basic.stats(hap[, !colnames(hap) %in% lininfo$TA[lininfo$Lineage %in% c(1,3,4)]]))$maf, 
     breaks = 30, col = 'black', xlab = 'Minor allele frequency', ylab = '# SNPs', main = '(B). Lineage2 (L2)',
     cex.main = 2, cex.lab = 1.5)
plot((basic.stats(hap[, !colnames(hap) %in% lininfo$TA[lininfo$Lineage %in% c(2,3,4)]]))$maf, 
     (basic.stats(hap[, !colnames(hap) %in% lininfo$TA[lininfo$Lineage %in% c(1,3,4)]]))$maf, 
     pch=20, cex=0.75, xlab='L1 MAF', ylab='L2 MAF', main = '(C). L1 vs. L2',
     cex.main = 2, cex.lab = 1.5)
dev.off()

# chromosome wise MAF
pdf('output/S9_chrom-MAF-spectrum.pdf', width = 8.5, height = 11)
par(mfrow=c(7,1), mar=c(4,5,2,1))
for (i in 1:7) {
  if (i %in% 1:6) {
    window.maf=sliding.window.maf(chrom = i, window.size = 1)
    plot(window.maf$window, window.maf$lin1.maf, type = 'h', ylim = c(-0.5, 0.5), xlab = '', ylab = 'MAF', axes = F, col = 'red')
    axis(side = 2, at = c(-0.5,0,0.5), labels = c(0.5,0,0.5))
    axis(side = 1, at = c(0,100,200,300,400,500,600,700), labels = c(0,100,200,300,400,500,600,700))
    lines(window.maf$window, 0 - window.maf$lin2.maf, col = 'blue', type = 'h')
    points(chr.length$center[i]/1000000, 0, pch = 18, cex = 1.5)
    mtext(text = paste0('Chr ', i), cex = 0.75, font = 2)
  } else {
    window.maf=sliding.window.maf(chrom = i, window.size = 1)
    plot(window.maf$window, window.maf$lin1.maf, type = 'h', ylim = c(-0.5, 0.5), xlab = 'Megabase', ylab = 'MAF', axes = F, col = 'red')
    axis(side = 2, at = c(-0.5,0,0.5), labels = c(0.5,0,0.5))
    axis(side = 1, at = c(0,100,200,300,400,500,600,700), labels = c(0,100,200,300,400,500,600,700))
    lines(window.maf$window, 0 - window.maf$lin2.maf, col = 'blue', type = 'h')
    points(chr.length$center[i]/1000000, 0, pch = 18, cex = 1.5)
    mtext(text = paste0('Chr ', i), cex = 0.75, font = 2)
  }
}
dev.off()

## Fst
source('R/hap.read.R')
source('R/gbs.diversity.R')
hap.converted <- hap.read(hap.obj = hap, data.col = 12)
diversity.stats <- gbs.diversity(hap = hap.converted, popGroups = lininfo)
save(diversity.stats, file = 'output/diversity.stats.RData')
