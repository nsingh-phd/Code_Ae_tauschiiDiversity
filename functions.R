## ################### ##
## General R functions ##
## Narinder Singh      ##
## ################### ##

# samples info
lininfo <- read.table('data/initial.sample.names.txt', header = T, as.is = T)
wheat.lines <- lininfo$TA[lininfo$Sample == 'wheat']
tauschii.lines <- lininfo$TA[lininfo$Sample == 'tauschii']
remove.samples <- c('TA10148') # these samples are not tauschii

# Compute basic stats
basic.stats <- function(hap = x) {
  # subsetting alleles
    a = substring(hap$alleles,1,1)   
    b = substring(hap$alleles,3,3)
    sum(a == b)
      # compute number of homozygous individuals
        cat('Counting alleles...')
        hap$alleleA = alleleA = rowSums(hap[, 12:ncol(hap)] == a, na.rm = T)
        hap$alleleB = alleleB = rowSums(hap[, 12:ncol(hap)] == b, na.rm = T)
        cat('Done. \n')
      # updating major minor alleles
        a[alleleA<alleleB] = substring(hap$alleles,3,3)[alleleA<alleleB]
        hap$alleleA[alleleA<alleleB] = alleleB[alleleA<alleleB]
        b[alleleA<alleleB] = substring(hap$alleles,1,1)[alleleA<alleleB]
        hap$alleleB[alleleA<alleleB] = alleleA[alleleA<alleleB]
        if (sum(a == b) == 0) cat('All SNPs are bi-allelic. \n')
        if (all(hap$alleleA >= hap$alleleB))  {
            cat('Alleles updated successfully. \n')
            hap$alleles = paste(a, b, sep = '/')
        } else cat ('Allele counts are irrespective of major/minor state. \n')
  # compute per SNP heterozygosity
    cat('Counting het calls...')
    hap$het = rowSums(hap[, 12:ncol(hap)] == 'H', na.rm = T) 
    hap$het_prop = (hap$het / (ncol(hap) - 11))
    cat('Done. \n')
  # compute MAF
    hap$maf = (2*hap$alleleB + hap$het) / (2*(hap$alleleA + hap$alleleB + hap$het))
    cat('MAF computation done. \n')
  # compute missing data
    cat('Computing missing data...')
    hap$missing = rowSums(hap[, 12:ncol(hap)] == 'N', na.rm = T) # compute per SNP missing data
    hap$propmiss = hap$missing/(ncol(hap)-11)
    cat('Done.')
  return(hap)
}

# compute fisher exact test p values
fisher.p.values <- function(x) {
  allele.testing <- matrix(c(x[1], x[2], x[3], x[4]),
                           nrow = 2,
                           dimnames = list(alleleA = c("Y", "N"),
                                           alleleB = c("Y", "N")))
  return(fisher.test(allele.testing)$p.value)
}

# convert hap to numeric
hap2numeric <- function(hap = x) {
  cat('Numericalizing the alleles...')
  hap01 = hap
  hap01[, 12:ncol(hap01)] = NA
  hap01[1:5, 1:50]
  
  a = substring(hap$alleles,1,1)   # subsetting alleles
  b = substring(hap$alleles,3,3)
  sum(a == b)
  
  cat ('almost done...Numericalizing...')
  
  hap01[hap == a] = -1
  hap01[hap == b] = 1
  hap01[hap == "H"] = 0
  hap01[1:10,1:20]
  cat('Done.')
  return(hap01)
}

# input file for faststructure
file4faststructure <- function(hap = x) {
  file4struct <- as.data.frame(t(hap[, 12:ncol(hap)]))
  file4struct <- cbind(rownames(file4struct), file4struct)
  file4struct <- file4struct[sort(rep(seq(nrow(file4struct)), 2)), ]
  cat('Getting alleles at each SNP...')
  alleles <- apply(X = file4struct[, -1], MARGIN = 2, function(x) intersect(unique(x), c('A', 'C', 'G', 'T')))
  
  # replace hets with alleles
  cat('Replacing IUPAC het codes with H...')
  for(i in 1:ncol(alleles)) {
    file4struct[which(file4struct[, i+1]=='H'), i+1] <- alleles[, i]
  }
  file4struct <- droplevels(file4struct)
  
  # convert to matrix for faster computation
  cat('Numericalizing...')
  file4struct <- as.matrix(file4struct)
  unique(c(file4struct[, -1]))
  file4struct[file4struct=="A"]=1
  file4struct[file4struct=="C"]=2
  file4struct[file4struct=="G"]=3
  file4struct[file4struct=="T"]=4
  file4struct[file4struct=="H"]=-9
  file4struct[file4struct=="N"]=-9
  file4struct <- gsub(' ', '', file4struct)
  cat('Done. \n')

  # format file for fastStructure
  file4struct <- as.data.frame(cbind(file4struct[, 1], NA, NA, NA, NA, NA, file4struct[, -1]))
  cat('There are', nrow(file4struct)/2, 'individuals and', ncol(file4struct)-6, 'SNPs')
  return(file4struct)
}

# check if snp is segregating or not
check.snp.levels <- function(hap = y) {
  snp.levels <- apply(hap[, 12:ncol(hap)], 1, function(x) {
    lvls <- levels(as.factor(as.character(x)))
    return(length(unique(lvls[!lvls %in% c('N', 'H')])))
  })
}

# nei's diversity index
nei.index <- function(hap = NULL) {
  maf = hap$maf
  nei = mean(1 - ((maf^2) + ((1-maf)^2)), na.rm = T)
  return(nei)
}

# plot MAF histogram
maf.plot <- function(hap = x) {
  hist(hap$maf, col = 'black', breaks = 20, xlab = 'Minor allele frequency', ylab = '#SNPs')
}

# find fixed allele
find.fixed.allele <- function(hap = x) {
  fixed.allele = apply(hap[, 12:ncol(hap)], 1, function(snp) unique(snp)[!unique(snp) %in% 'N'])
  if (all(nchar(fixed.allele) == 1)) {
    cat(length(fixed.allele), 'fixed alleles found successfully.')
    return(fixed.allele)
  } else cat('One or more SNP is segregating or totally missing!')
}

# find private allele
find.private.allele <- function(hap = x, fixed = y) {
  alleles <- apply(hap[, 12:ncol(hap)], 1, function(snp) unique(snp)[!unique(snp) %in% c('H', 'N')])
  private.allele <- sapply(seq_along(fixed), function(i) setdiff(alleles[[i]], fixed[i]))
  private.allele[lengths(private.allele) == 0] <- 'N'
  return(unlist(private.allele))
}
  
# plot hybrids
plot.hybrid <- function(i = x, cex=0.85) {
  # paint tags on chromosomes
  spread = 0.45
  barplot(height = chr.length$length/1000000, col = 'azure', xlab = 'Mb', ylab = 'Chromosome',
          xlim = c(0,750), ylim = c(0,10), horiz = T, offset = 1)
  axis(side = 2, at = chr.length$barplot.x, labels = chr.length$chrom, las = 2, cex=2, line = NA)
  title(main = colnames(color.hybrids)[i], line = 2.5)
  rect(color.hybrids$pos/1000000, color.hybrids$barplot.x - spread,  
       color.hybrids$pos/1000000, color.hybrids$barplot.x + spread,  
       border = alpha(as.character(color.hybrids[,i]), alpha.val))
  # place centromeres
  rect(chr.length$center/1000000 - 2, chr.length$barplot.x - spread, 
       chr.length$center/1000000 + 2, chr.length$barplot.x + spread, col = 'black')
  legend('top', pch = '|', lty = NA, lwd = 1, pt.cex = 1.5, bty = 'n', horiz = T, cex = cex,
         col = c('black', 'red', 'blue'),
         legend = c('Centromere',
                    paste0('L1 alleles (', sum(color.hybrids[, i] == 'red', na.rm = T),')'),
                    paste0('L2 alleles (', sum(color.hybrids[, i] == 'blue', na.rm = T),')')))
}

# sliding window MAF
sliding.window.maf <- function(chrom = x, window.size = y) {
  window.size = window.size * 1000000                           # convert Mb to bp
  window.num <- 1:ceiling(chr.length$length[i] / window.size)   # list of window numbers
  lower.bound <- (window.num*window.size) - window.size + 1     # lower bound of a window
  upper.bound = window.num * window.size                        # upper bound of a window
  window.maf <- data.frame('window' = window.num, 'lin1.maf' = NA, 'lin2.maf' = NA)
  
  for (i in 1:nrow(window.maf)) {
    chrom = as.integer(chrom)
    window.maf$lin1.maf[i] <- mean(lin1$maf[lin1$chrom == chrom & lin1$pos >= lower.bound[i] & lin1$pos <= upper.bound[i]])
    window.maf$lin2.maf[i] <- mean(lin2$maf[lin2$chrom == chrom & lin2$pos >= lower.bound[i] & lin2$pos <= upper.bound[i]])
  }
  return(window.maf[complete.cases(window.maf), ])
}

# Summary of missing SNPs data
missing.summary <- function(missing.SNP.data.col = hap$propmiss) {
   threshold.values <- seq(from = 0.1, to = 0.95, by = 0.05)
   cat('\n')
   cat('Number of SNPs with percent missing data', '\n')
   cat('========================================', '\n')
   cat('\n')
   for (i in 1:length(threshold.values)) {
      cat(' ', threshold.values[i]*100, '% -', sum(as.numeric(missing.SNP.data.col) < threshold.values[i]))
      cat('\n')
   }
   cat('\n')
}

# Summary of blank wells
check.blank.wells <- function(data.obj = data.obj, blank.name = 'blank') {
   blank.wells <- data.obj[, grep(blank.name, colnames(data.obj), ignore.case = T)]
   n.cols <- ncol(blank.wells)

   if (is.null(n.cols)) {
      blank.wells <- as.character(blank.wells)
      if(is.null(n.cols)) cat('Blank well has', sum(blank.wells != "N"), 'sequence reads')
   } else if (n.cols == 0) {
      'There are no blank wells'
   } else {
      blankSum <- apply(blank.wells!="N", 2, sum)
      blankSum <- data.frame(reads.present = blankSum, percent.total = paste(round((blankSum/nrow(blank.wells))*100, 3), '%'),
                             row.names = names(blankSum))
      if (exists('blankSum')) print(blankSum)
   }
}

# Summary of non-blank wells
check.non.blank.wells <- function(data.obj = hap, 
                                  cols.to.exclude = c('blank', 'void', 'empty'), 
                                  data.col = 12, miss.threshold = 1) {
   data.obj <- data.obj[, grep(paste(cols.to.exclude, collapse = '|'), colnames(data.obj), ignore.case = T, invert = T)]
   data.obj <- data.obj[, data.col:ncol(data.obj)]
      reads.present <- apply(data.obj!="N", 2, sum)
      reads.present <- data.frame(reads.present, paste(round(((nrow(data.obj)-reads.present)/nrow(data.obj))*100, 3), '%'))
      colnames(reads.present) <- c('SNPs.Present', 'Percent.SNPs.Missing')
      cat('Samples with >', miss.threshold*100, '% missing data', '\n')
      cat('===============================', '\n')
      print(reads.present[nrow(data.obj)-reads.present$SNPs.Present > miss.threshold*(nrow(data.obj)), ])
}

# convert to IUPAC code
convert2iupac <- function(hap = x) {
  ## substitute all values of H for IUPAC symbol
  hap=ddply(hap,.(alleles), function(y) {
    if(y$alleles[1]=="A/G" | y$alleles[1]=="G/A") y[y=="H"]="R"
    if(y$alleles[1]=="C/T" | y$alleles[1]=="T/C") y[y=="H"]="Y"
    if(y$alleles[1]=="G/C" | y$alleles[1]=="C/G") y[y=="H"]="S"
    if(y$alleles[1]=="A/T" | y$alleles[1]=="T/A") y[y=="H"]="W"
    if(y$alleles[1]=="G/T" | y$alleles[1]=="T/G") y[y=="H"]="K"
    if(y$alleles[1]=="A/C" | y$alleles[1]=="C/A") y[y=="H"]="M"
    return(y)
    }
  )
  hap[order(hap$chrom, hap$pos), ]
}

# plot MAF
plot.maf <- function(hap = NULL, type = 'h') {
  for (i in 1:7) {
    chrName = paste(i, "D", sep = '')
    chr <- tau.only[grep(chrName, hap$chrom, ignore.case = T), c(3,4,10)]
    plot(chr$pos/1000000, chr$maf, 
         type = type, xlab = 'Mb', ylab = 'MAF', main = paste('Chromosome ', i, 'D', sep = ''),
         xlim = c(0, (max(chr$pos)/1000000) + 1))
    cen.start = cen.positions[grep(chrName, cen.positions$agp_chr, ignore.case = T), 2] / 1000000
    abline(v = cen.start, lwd = 3, col = 'red')
  }
}

# identify identical accessions
identical.accessions <- function(id, rows = 1:nrow(id), cols = 1:ncol(id), threshold = 0.99) {
  
  for (i in 1:ncol(id)) { # loop to set replicates to 0
    id[which(rownames(id) %in% colnames(id)[i]), i] = 0
  }
  
  id[id > 1] = NA
  tmp.csv.out <- matrix(NA, nrow = length(rows), ncol = 2)
  tmpIdMat <- t(id[rows, cols])
  for (i in 1:nrow(tmpIdMat)) {
    tmp.csv.out[i, 1] <- rownames(tmpIdMat)[i]
    for (j in i:ncol(tmpIdMat)) {
      if (tmpIdMat[i, j] >= threshold)
        tmp.csv.out[i, 2] <- ifelse(is.na(tmp.csv.out[i, 2]),
                                    colnames(tmpIdMat)[j],
                                    paste(tmp.csv.out[i, 2], colnames(tmpIdMat)[j], sep=", "))
    }
  }
  return(tmp.csv.out)
}

# unique alleles in lin1 and lin2 (to be used only with fixed alleles)
unique.alleles <- function(x = NULL) {
  levels.x <- unique(x)
  levels.x <- levels.x[levels.x %in% c('A', 'C', 'G', 'T')]
  return(levels.x)
}