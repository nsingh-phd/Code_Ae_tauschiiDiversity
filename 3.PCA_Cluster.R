
source('functions.R')

###############
## all samples
###############

# read hap01 file
hap01 <- fread("data/hap01.filtered.txt", header = T, check.names = F, data.table = F)
lininfo = read.table(file="data/lininfo.txt", header=T)
  geno = t( as.matrix( hap01[, 12:ncol(hap01)] ) ) # only genotypes
  geno[1:15,1:5]
  
####################
## CLUSTER ANALYSIS
####################
  # compute genetic distance
    ifelse(test = file.exists('data/distMat.RData'),
             yes = load('data/distMat.RData'), no = distMat <- dist(geno))
  # saving distance matrix because it's a computationally intensive step
    if (!file.exists('data/distMat.RData')) save(distMat, file = "data/distMat.RData")
    
  # cluster accessions and convert to phylo object for ape
    hc2 <- as.phylo(hclust(distMat))
    
  # cluster coloring
    edgecols <- cbind('TA'=NA, 1:nrow(hc2$edge), color='black') # create data frame
    edgecols[, 1] <- hc2$tip.label[hc2$edge[, 2]] # get labels
    edgecols <- as.matrix(merge(edgecols, lininfo, by = 'TA', all.x = T)) # get samples info
    edgecols <- edgecols[order(as.numeric(edgecols[, 2])), ] # get samples in original order
    edgecols[,4] <- as.numeric(edgecols[,4])
    
    # coloring diff lineages with different colors
    edgecols[, 3][edgecols[, 4] == 1] = "red"
    edgecols[, 3][edgecols[, 4] == 2] = "blue"
    edgecols[, 3][edgecols[, 4] == 3] = "gold"
    edgecols[, 3][edgecols[, 4] == 4] = "green"
    
    # tips coloring
    tipcols <- as.matrix(merge(hc2$tip.label, edgecols, by = 1, sort = F))
    
    # plotting tree
    pdf("output/Fig.3_Cluster.pdf", width = 8.5, height = 8.5)
      plotnj(unrooted.tree = hc2, type = 'u',
             show.tip.label = T, lab4ut = "axial", label.offset = 1, cex = 0.25,
             edge.color = edgecols[, 3], edge.width = 1, tip.color = tipcols[, 3], rotate.tree = -50)
      legend(5, 110, lty = 1, lwd = 2, cex = 0.75,
             legend = c('Lineage 1 (L1)', 'Lineage 2 (L2)', 'Possible L1 & L2 Hybrids', 'Wheat'),
             text.col = 'black', col = c("red", "blue", "gold", "green"))
    dev.off()
    
#######
## PCA
#######
  # compute A matrix for eigen values
    A = A.mat(geno, impute.method="mean", return.imputed = T)
  # PCA label names
    names = as.data.frame(rownames(A$A)); colnames(names)=c("TA")
    names = merge(names, edgecols, by='TA', sort = F, all.x = T)
  # give shapes to points
    names$color = as.character(names$color)
    names$Lineage = as.numeric(names$Lineage)

  # eigenvectors
  e = eigen(A$A)
    ## PCA plot
    pdf(file = 'output/S4_PCA.pdf', width = 8.25, height = 11)
    par(mfrow=c(2,1))
      plot(e$vectors[,1], e$vectors[,2], pch=names$Lineage, cex=1.25, col=names$color, lwd =1.5,
           xlab=paste("PC 1 (", round(e$values[1]/sum(e$values), digits=2)*100, "%)", sep = ''),
           ylab=paste("PC 2 (", round(e$values[2]/sum(e$values), digits=2)*100, "%)", sep = ''),
           main = substitute(paste("PCA of WGRC ", italic("Ae. tauschii")," Collection and Wheat")),
           font.main=1, cex.main=1.5, cex.lab=1.25)
      legend("topleft", pch=c(1,2,3,4), cex=1.25, pt.cex = 2, pt.lwd = 2,
             col=c("red", "blue", "gold", "green"), 
             legend=c("Lineage1 (L1)","Lineage2 (L2)","Possible L1 & L2 hybrids", "Wheat"))

#################
## tauschii only
#################
  geno = geno[rownames(geno) %in% tauschii.lines, ]
  geno[1:15,1:5]
  # compute A mat
    A = A.mat(geno, impute.method="mean", return.imputed = T)
    # PCA label names
    names = as.data.frame(rownames(A$A)); colnames(names)=c("TA")
    names = merge(names, edgecols, by='TA', sort = F, all.x = T)
    # give shapes to points
    names$color = as.character(names$color)
    names$Lineage = as.numeric(names$Lineage)

  # eigenvectors
  e = eigen(A$A)
    ## PCA plot
      plot(e$vectors[,1], e$vectors[,2], pch=names$Lineage, cex=1.25, col=names$color, lwd =1.5,
           xlab=paste("PC 1 (", round(e$values[1]/sum(e$values), digits=2)*100, "%)", sep = ''),
           ylab=paste("PC 2 (", round(e$values[2]/sum(e$values), digits=2)*100, "%)", sep = ''),
           main = substitute(paste("PCA of WGRC ", italic("Ae. tauschii")," Collection")),
           font.main=1, cex.main=1.5, cex.lab=1.25)
      legend("topright", pch=c(1,2,3), cex=1.25, pt.cex = 2, pt.lwd = 2,
             col=c("red", "blue", "gold"), 
             legend=c("Lineage1 (L1)","Lineage2 (L2)","Possible L1 & L2 hybrids"))
    dev.off()

##########################
### 3D PCA, tauschii only 
##########################
passport <- read.csv("data/long_lat_alt.csv", header = T, as.is = T, fill = T)
passport <- merge(lininfo, passport, by = 'TA', all.x = T)

## color vector for PCA
passport <- merge(data.frame('TA'=rownames(A$A), 'col'="black", 'ord'=1:nrow(A$A), stringsAsFactors = F), 
                  passport, by = 'TA', sort = F, all.x = T)

## function for creating lin1 long gradient colors
gradient.lin1.long <- colorRampPalette(colors = c("yellow", "orange", "red"))
lin1.long <- sort(unique(passport$longitude[passport$Lineage == 1]))
cols.lin1.long <- data.frame('longitude'=lin1.long, 'cols.long'=gradient.lin1.long(length(lin1.long)), stringsAsFactors = F)

## function for creating lin2 altitude gradient colors
gradient.lin2.alt <- colorRampPalette(colors = c("blue", "lightgreen", "green"))
lin2.alt <- sort(unique(passport$altitude[passport$Lineage == 2]))
cols.lin2.alt <- data.frame('altitude'=lin2.alt, 'cols.alt'=gradient.lin2.alt(length(lin2.alt)), stringsAsFactors = F)

# merge colors with labels
passport <- merge(passport, cols.lin1.long, by = 'longitude', all.x=T)
passport <- merge(passport, cols.lin2.alt, by = 'altitude', all.x=T)

## fill color column with gradient colors
passport$col[passport$Lineage == 1] <- passport$cols.long[passport$Lineage == 1]
passport$col[passport$Lineage == 2] <- passport$cols.alt[passport$Lineage == 2]
passport$col[passport$Lineage == 3]="gold"
passport <- passport[order(passport$ord),]

# ## plot and color PCA plot

pca1v <- paste("PC 1 (", round(e$values[1]/sum(e$values), digits=2)*100, "%)", sep = '')
pca2v <- paste("PC 2 (", round(e$values[2]/sum(e$values), digits=2)*100, "%)", sep = '')
pca3v <- paste("PC 3 (", round(e$values[3]/sum(e$values), digits=2)*100, "%)", sep = '')

pdf('output/Fig.4_PCA.pdf', width = 11, height = 6.38)
layout(matrix(1:2, ncol=2), width = c(3, 1), height = c(1, 1))
scatterplot3d(e$vectors[,1], e$vectors[,3], e$vectors[,2],
              pch = passport$Lineage, color = passport$col, cex.symbols = 1.5, lwd = 1.5,
              xlab = pca1v, ylab = pca3v, zlab = pca2v,
              zlim = c(-0.1, 0.1), ylim = c(-0.1, 0.1))
legend(x = 0.3, y = 3.77, pch=c(1,2,3), pt.cex = 1, pt.lwd = 1.5, 
       box.lwd = 1, cex = 0.7, bg = 'white', bty = 'n',
       legend = c("Lineage1 (L1)", "Lineage2 (L2)", "Possible L1 & L2 hybrids"), 
       col = c(1, 1, 'gold'))
x.for.legend <- c(rep(1, nrow(cols.lin1.long)), rep(1.7, nrow(cols.lin1.long)))
y.for.legend <- c(1:nrow(cols.lin1.long), 1:nrow(cols.lin1.long))
plot(x = x.for.legend, y = y.for.legend, 
     pch = 15, cex = 2.5, 
     col = c(gradient.lin1.long(length(x.for.legend)/2), gradient.lin2.alt(length(x.for.legend)/2)), 
     ann = F, axes = F, xlim = c(1, 2))
axis(side = 2, at = c(1, 44, 87, 130, 173, 216), 
     labels = c(26, 43, 60, 77, 94, 111),
     line = 0.15)
axis(side = 2, at = seq(1, 216, length.out = 9), 
     labels = round(seq(-30, 2250, length.out = 9), 0),
     line = -4.55)
mtext(text = 'Longitude gradient', side = 2, line = -1.7)
mtext(text = 'Altitude gradient (m)', side = 2, line = -6.5)
mtext(text = "L1", side = 3, adj = 0)
mtext(text = "L2", side = 3, adj = 0.7)

dev.off()

par(mfrow=c(1,1))

##############
## BiPlots ##
#############
passport <- merge(data.frame('TA'=rownames(A$A), 'PC2'=e$vectors[, 2], 'PC3'=e$vectors[, 3], stringsAsFactors = F),
                  passport, by = 'TA', sort = F, all.x = T)
non.outliers.l1 <- passport$Lineage==1 & passport$longitude < 100 & passport$longitude > 38 # removing extreme values of longitude

pdf('output/S5_PC_correlations.pdf', width = 11, height = 6)
par(mfrow=c(1,2))
corrplot(cor(passport[non.outliers.l1, c(2,4,5,9)], use = "complete.obs"),
         type = 'upper', order = 'original', diag = T, addCoef.col = 1, tl.srt = 45)
corrplot(cor(passport[passport$Lineage==2, c(3,4,5,9)], use = "complete.obs"),
         type = 'upper', order = 'original', diag = T, addCoef.col = 1, tl.srt = 45)
dev.off()

# pc2-long
pdf(file = 'output/Fig.5_pc2_longitude.pdf', height = 8.5, width = 11)
plot(passport$longitude[non.outliers.l1], passport$PC2[non.outliers.l1], pch = 1,
     cex = 1.25, cex.axis = 1.5, cex.lab = 1.5, xlab = 'Longitude', ylab = 'PC2', main = 'Longitudinal distribution along PC2')
segments(51.38, 0, 51.38, -0.08, col = 'red', lty = 3, lwd = 2)
text(52.5, -0.02, labels = 'Tehran', srt = 90, cex=1.5)
legend('topright', cex = 1.5, bg = 'azure',
       legend = paste('r =', 
                      round(cor(passport$longitude[non.outliers.l1], passport$PC2[non.outliers.l1], use = "complete.obs"),2)))
dev.off()

# pc3-alt
pdf(file = 'output/S6_pc3_altitude.pdf', height = 8.5, width = 11)
plot(passport$altitude[passport$Lineage==2], passport$PC3[passport$Lineage==2], pch = 2,
     cex = 1.25, cex.axis = 1.5, cex.lab = 1.5,
     xlab = 'Altitude', ylab = 'PC3', main = 'Altitudinal distribution along PC3')
segments(150, -0.15, 150, -0.06, col = 'red', lty = 3, lwd = 2)
legend('topright', cex = 1.5, bg = 'azure',
       legend = paste('r =', 
                      round(cor(passport$altitude[passport$Lineage==2], passport$PC3[passport$Lineage==2], use = "complete.obs"),2)))
dev.off()

#################
### violin plots
#################
passport$Lineage = as.factor(passport$Lineage)

pdf(file = 'output/violin_altitude.pdf')
ggplot(data=subset(passport[passport$Lineage != 3, ], !is.na(altitude)), aes(x=Lineage, y=altitude)) +
  geom_violin(trim = F) + coord_flip() + ylab('Altitude (m)') + xlab(NULL) +
  stat_summary(fun.y=median, geom="point", size=2, color="red") +
  scale_x_discrete(labels=c("1" = "L1", "2" = "L2")) + labs(title = "(A)") +
  theme(axis.text=element_text(size=17), 
        axis.title=element_text(size=20,face="bold"),
        plot.title=element_text(size=22))
dev.off()

pdf(file = 'output/violin_longitude.pdf')
ggplot(data=subset(passport[passport$Lineage != 3 & passport$longitude < 100 & passport$longitude > 38, ], !is.na(longitude)), aes(x=Lineage, y=longitude)) +
  geom_violin(trim = F) + coord_flip() + ylab('Longitude') + xlab(NULL) +
  stat_summary(fun.y=median, geom="point", size=2, color="red") +
  scale_x_discrete(labels=c("1" = "", "2" = "")) + labs(title = "(B)") +
  theme(axis.text=element_text(size=17), 
        axis.title=element_text(size=20,face="bold"),
        plot.title=element_text(size=22))
dev.off()

pdf(file = 'output/violin_latitude.pdf')
ggplot(data=subset(passport[passport$Lineage != 3, ], !is.na(latitude)), aes(x=Lineage, y=latitude)) +
  geom_violin(trim = F) + coord_flip() + ylab('Latitude') + xlab(NULL) +
  stat_summary(fun.y=median, geom="point", size=2, color="red") +
  scale_x_discrete(labels=c("1" = "", "2" = "")) + labs(title = "(C)") +
  theme(axis.text=element_text(size=17), 
        axis.title=element_text(size=20,face="bold"),
        plot.title=element_text(size=22))
dev.off()

###
