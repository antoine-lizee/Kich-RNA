
library(dplyr)
library(ggplot2)

# Read files --------------------------------------------------------------

# Read raw file
allDf <- readr::read_csv("Input/RNAseq-PEMF_reads_count_data/counts_raw_all.csv")
# Read normalized file
allNormDf <- readr::read_csv("Input/RNAseq-PEMF_reads_count_data/counts.norm_DEseq_all.csv", skip = 3) 
filtNormDf <- readr::read_csv("Input/RNAseq-PEMF_reads_count_data/counts.norm_DEseq_filt.csv", skip = 3) 
# Read column names seperatly
allNormDfCN <- readr::read_csv("Input/RNAseq-PEMF_reads_count_data/counts.norm_DEseq_all.csv", col_names = FALSE, n_max = 1) 
stopifnot(colnames(allDf[-1]) == allNormDfCN[-1])
colnames(allNormDf) <- colnames(allDf)

samps <- colnames(allDf)[-1]
sampsDf <- data_frame(name = samps, 
                      day = as.numeric(sub("D([[:digit:]]+)\\.([[:alpha:]]?)([[:digit:]])([SL])", replacement = "\\1", samps)),
                      type = sub("D([[:digit:]]+)\\.([[:alpha:]]?)([[:digit:]])([SL])", replacement = "\\2", samps),
                      mouse_num = as.numeric(sub("D([[:digit:]]+)\\.([[:alpha:]]?)([[:digit:]])([SL])", replacement = "\\3", samps)),
                      tissue = factor(sub("D([[:digit:]]+)\\.([[:alpha:]]?)([[:digit:]])([SL])", replacement = "\\4", samps))) %>% 
  mutate(mouseName = paste0(ifelse(type == "", "Z", type), mouse_num),
         type = factor(type, labels = c("Zero", "Control", "PEMF")),
         tissue = factor(tissue, labels = c("Lymph", "Spinal")))

## Create the matrices
# Raw data
allMatFull <- as.matrix(allDf[-1])
rownames(allMatFull) <- allDf[[1]]
noExpr <- apply(allMatFull, 1, function(row) all(row == 0))
allMat <- allMatFull[!noExpr,]
# Normlaized expression
allNormMat <- as.matrix(allNormDf[-1])
rownames(allNormMat) <- allNormDf[[1]]


# PCAs for non-normalized (raw) data --------------------------------------------------------------------

pcaGenes <- prcomp(allMat)
plot(pcaGenes)
plot(pcaGenes$x[, 1:2])
pcaGenes$rotation[,1]

pdf("Output/PCA_Raw.pdf")

pca <- prcomp(t(allMat))#, scale. = TRUE)
barplot(pcaNorm$sdev^2 / sum(pcaNorm$sdev^2) * 100,
        main = "Variance plot, scaled PCA")
plot(cumsum(pcaNorm$sdev^2 / sum(pcaNorm$sdev^2)) * 100,
     main = "Cumulative Variance plot, scaled PCA", ylim = c(0,100))
plot(pca$x[, 1:2], col = RColorBrewer::brewer.pal(3, "Dark2")[sampsDf$type], pch = 20,
     main = "PCs plot, non-scaled PCA, per type")
plot(pca$x[, 1:2], col = RColorBrewer::brewer.pal(3, "Dark2")[sampsDf$tissue], pch = 20,
     main = "PCs plot, non-scaled PCA, per tissue")
plot(pca$rotation[,1:2],
     main = "Rotation plot (genes), non-scaled PCA")
plot(pca$x[, 3:4], col = RColorBrewer::brewer.pal(3, "Dark2")[sampsDf$type], pch = 20,
     main = "PCs plot, non-scaled PCA, per type")
plot(pca$rotation[,3:4],
     main = "Rotation plot (genes), non-scaled PCA")

# No outlier sample on PC1, here is the main (outlier) gene on PC1:
hist(allMat[pca$rotation[,1] < -0.8, ], 100,
     main = "feature plot for gene outlier on PC1")
# Outlier sample on PC2, driven by the outlier gene:
hist(allMat[pca$rotation[,2] > 0.8, ], 100,
     main = "feature plot for gene outlier on PC2")

dev.off()

# PCAs for normalized data --------------------------------------------------------------------

pcaNormGenes <- prcomp(allNormMat)
plot(pcaNormGenes)
plot(pcaNormGenes$x[, 1:2])
barplot(pcaNormGenes$rotation[,1])

pdf("Output/PCA_Norm.pdf")

## No scaling:
pcaNorm <- prcomp(t(allNormMat))
barplot(pcaNorm$sdev^2 / sum(pcaNorm$sdev^2) * 100,
        main = "Variance plot, scaled PCA")
plot(cumsum(pcaNorm$sdev^2 / sum(pcaNorm$sdev^2)) * 100,
     main = "Cumulative Variance plot, scaled PCA", ylim = c(0,100))
plot(pcaNorm$x[, 1:2], col = RColorBrewer::brewer.pal(3, "Dark2")[sampsDf$type], pch = 20,
     main = "PCs plot, non-scaled PCA, per sample")
plot(pcaNorm$x[, 1:2], col = RColorBrewer::brewer.pal(3, "Dark2")[sampsDf$type], pch = 20,
     main = "PCs plot, non-scaled PCA, per type")
plot(pcaNorm$x[, 1:2], col = RColorBrewer::brewer.pal(3, "Dark2")[sampsDf$tissue], pch = 20,
     main = "PCs plot, non-scaled PCA, per tissue")
plot(pcaNorm$rotation[,1:2],
     main = "Rotation plot (genes), non-scaled PCA")
plot(pcaNorm$x[, 3:4], col = RColorBrewer::brewer.pal(3, "Dark2")[sampsDf$type], pch = 20,
     main = "PCs plot, non-scaled PCA, per type")
plot(pcaNorm$x[, 3:4], col = RColorBrewer::brewer.pal(3, "Dark2")[sampsDf$tissue], pch = 20,
     main = "PCs plot, non-scaled PCA, per tissue")
plot(pcaNorm$rotation[,3:4],
     main = "Rotation plot (genes), non-scaled PCA")

# No outlier sample on PC1, here is the main (outlier) gene on PC1:
hist(allNormMat[pcaNorm$rotation[,1] < -0.8, ], 100,
     main = "feature plot for gene outlier on PC1")
# Outlier sample on PC2, driven by the outlier gene:
hist(allNormMat[pcaNorm$rotation[,2] > 0.8, ], 100,
     main = "feature plot for gene outlier on PC2")
# No one driver to PC3:
hist(allNormMat[pcaNorm$rotation[,3] > 0.8, ], 100,
     main = "feature plot for gene outlier 1/2 on PC3")
hist(allNormMat[pcaNorm$rotation[,3] < -0.2, ], 100,
     main = "feature plot for gene outlier 2/2 on PC3")

## PCA with scaling each gene separately:
pcaNorm <- prcomp(t(allNormMat), scale. = TRUE)
barplot(pcaNorm$sdev^2 / sum(pcaNorm$sdev^2) * 100,
        main = "Variance plot, scaled PCA")
plot(cumsum(pcaNorm$sdev^2 / sum(pcaNorm$sdev^2)) * 100,
     main = "Cumulative Variance plot, scaled PCA", ylim = c(0,100))
plot(pcaNorm$x[, 1:2], col = RColorBrewer::brewer.pal(3, "Dark2")[sampsDf$type], pch = 20,
     main = "PCs plot, scaled PCA, per type")
plot(pcaNorm$x[, 1:2], col = RColorBrewer::brewer.pal(3, "Dark2")[sampsDf$tissue], pch = 20,
     main = "PCs plot, scaled PCA, per tissue")
text(pcaNorm$x[, 1:2], col = RColorBrewer::brewer.pal(3, "Dark2")[sampsDf$tissue], label = as.character(1:40),
     main = "PCs plot, scaled PCA, per tissue")
plot(pcaNorm$rotation[,1:2], pch = 20, col = adjustcolor("black", a = 0.1),
     main = "Rotation plot (genes), scaled PCA"); abline(v = 0, col = "red"); abline(h = 0, col = "red")
plot(pcaNorm$x[, 3:4], col = RColorBrewer::brewer.pal(3, "Dark2")[sampsDf$type], pch = 20,
     main = "PCs plot, scaled PCA, per type")
text(pcaNorm$x[, 3:4], col = RColorBrewer::brewer.pal(3, "Dark2")[sampsDf$type], label = as.character(1:40),
     main = "PCs plot, scaled PCA, per tissue")
plot(pcaNorm$x[, 3:4], col = RColorBrewer::brewer.pal(3, "Dark2")[sampsDf$tissue], pch = 20,
     main = "PCs plot, scaled PCA, per tissue")
plot(pcaNorm$rotation[,3:4], pch = 20, col = adjustcolor("black", a = 0.1),
     main = "Rotation plot (genes), scaled PCA"); abline(v = 0, col = "red"); abline(h = 0, col = "red")
plot(pcaNorm$x[, 5:6], col = RColorBrewer::brewer.pal(3, "Dark2")[sampsDf$type], pch = 20,
     main = "PCs plot, scaled PCA, per type")
text(pcaNorm$x[, 5:6], col = RColorBrewer::brewer.pal(3, "Dark2")[sampsDf$type], label = as.character(1:40),
     main = "PCs plot, scaled PCA, per tissue")

# outliers:
print(sampsDf[pcaNorm$x[,4] < -60, ])

dev.off()


# Parking lot -------------------------------------------------------------

# pcaDf <- as.data.frame(ttt$x)
# qplot(data = as.data.frame(ttt$x[,1:2], type = ), x = PC1, y = PC2, 
#       col = type, alpha = 0.2, shape = I(20)) +
#   scale_color_brewer("Dark2")
