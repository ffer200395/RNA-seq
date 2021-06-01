## --- Paquetes necesarios para el análisis
if (!require(DESeq2))BiocManager::install("DESeq2")
if(!(require(dplyr))) install.packages("dplyr")

## --- Cargamos los datos
counts <- read.csv('data/counts.csv', sep = ';', header = T, row.names = 1)
targets <- read.csv('data/targets.csv', sep = ',', header = T)

# --- Creación del subconjunto del análisis (usar datos de expresión RNA-Seq) y preparación de los datos
library(dplyr)
set.seed(1234)
# Filtramos para tener sólo datos RNA Seq (NGS) y elegimos 10 muestras de cada grupo
my.targets <- targets %>% group_by(Group) %>% sample_n(10)

# Modificamos los nombres de las muestras para que sea igual que las columnas de la tabla counts
my.targets$Sample_Name <- gsub('-', '.', my.targets$Sample_Name)

# Filtramos la tabla counts
my.counts <- counts[, my.targets$Sample_Name]

# Los nombres de las filas de Targets deben coincidir con las columnas de Counts
rownames(my.targets) <- my.targets$Sample_Name
sum(colnames(my.counts) == rownames(my.targets))

# Convertimos la variable Group en factor (El grupo no tratado es el de referencia)
my.targets$Group %>% as.factor() %>% relevel('NIT') -> my.targets$Group

# Construimos el objeto DESeqDataSet a partir de Targets y Counts
library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData = my.counts, colData = my.targets, design = ~ Group)

# --- Preprocesado de los datos: filtraje, normalización y estabilización de la varianza

# Nos quedamos con aquellos genes que tengan al menos 10 reads en total
nrow(dds)
dds <- dds[ rowSums(counts(dds)) >= 10, ]
nrow(dds)

# Número de reads por muestra expresado en millones
sort(colSums(my.counts))/1e6

# log2 normalized counts 
log.norm.counts <- log2(counts(estimateSizeFactors(dds), normalized=TRUE) + 1)

# Regularized log transformation
rld <- rlog(dds)

# Variance stabilizing transformation
vsd <- vst(dds)

# --- EDA: Scatter plots, boxplots, cluster y PCA

# Scatter plots
par(mfrow=c(2,2),mar=c(1,1,1,1))
plot(log2(counts(dds)[,1:2] + 1), cex=.1, main = 'Sin normalizar')
plot(log.norm.counts[,1:2], cex=.1, main = 'Log norm')
plot(assay(rld)[,1:2], cex=.1, main = 'Regularized Log')
plot(assay(vsd)[,1:2], cex=.1, main = 'VST')

# Boxplots
par(mfrow=c(2,2))
boxplot(log2(counts(dds) + 1), names=my.targets$SRA_Sample, las=2, main = 'Sin normalizar')
boxplot(log.norm.counts, names=my.targets$SRA_Sample, las=2, main = 'Log norm')
boxplot(assay(rld), names=my.targets$SRA_Sample, las=2, main = 'Regularized Log')
boxplot(assay(vsd), names=my.targets$SRA_Sample, las=2, main = 'VST')

# Cluster
par(mfrow=c(3,1),mar=c(1,1,1,1))
plot(hclust(dist(t(log.norm.counts))), labels=colData(dds)$Group, main = 'Log norm')
plot(hclust(dist(t(assay(rld)))), labels=colData(rld)$Group, main = 'Regularized Log')
plot(hclust(dist(t(assay(vsd)))), labels=colData(vsd)$Group, main = 'VST')

# PCA
my.pca <- function(datos, grupos){
  pca_res <- prcomp(t(assay(datos)), scale = FALSE)
  pca_df <- as.data.frame(pca_res$x)
  pca_df$Group <- grupos
  pcts <- round(pca_res$sdev^2/sum(pca_res$sdev^2)*100,1)
  plot(x = pca_df$PC1, y = pca_df$PC2, col = pca_df$Group, xlab = paste0('PC1 ',pcts[1],'%'), ylab = paste0('PC2 ',pcts[2],'%'))
  legend("bottomleft", legend = unique(pca_df$Group), col = unique(pca_df$Group), pch =19, horiz = T,y.intersp=0.5,x.intersp=0.5,text.width=0.1,bty = 'n')
}
my.pca(rld, my.targets$Group)
my.pca(vsd, my.targets$Group)

# --- Identificación de genes diferencialmente expresados

# Ejecutamos el modelo DESeq2
dds <- DESeq(dds)

# Hacemos los contrastes entre grupos
res_sfi_nit <- results(dds, contrast=c("Group",'SFI','NIT'))
res_eli_nit <- results(dds, contrast=c("Group",'ELI','NIT'))
res_eli_sfi <- results(dds, contrast=c("Group",'ELI','SFI'))

# Número de genes con p-val ajustado menor de 0.05
# If we consider a fraction of 10% false positives acceptable, we can consider all genes 
# with an adjusted p value below 10% = 0.1 as significant. How many such genes are there?
sum(res_sfi_nit$padj < 0.05, na.rm=TRUE)
sum(res_eli_nit$padj < 0.05, na.rm=TRUE)
sum(res_eli_sfi$padj < 0.05, na.rm=TRUE)

# Genes más significativos
res_sfi_nit_Sig <- subset(res_sfi_nit, padj < 0.05)
res_eli_nit_Sig <- subset(res_eli_nit, padj < 0.05)
res_eli_sfi_Sig <- subset(res_eli_sfi, padj < 0.05)

# Counts plot
#A quick way to visualize the counts for a particular gene is to use the plotCounts
plotCounts(dds, gene = rownames(res_sfi_nit_Sig)[which.min(res_sfi_nit_Sig$padj)], intgroup=c("Group"))
plotCounts(dds, gene = rownames(res_eli_nit_Sig)[which.min(res_eli_nit_Sig$padj)], intgroup=c("Group"))
plotCounts(dds, gene = rownames(res_eli_sfi_Sig)[which.min(res_eli_sfi_Sig$padj)], intgroup=c("Group"))

# p-val hist
hist(res_sfi_nit$pvalue[res_sfi_nit$baseMean > 1], xlab = 'p-valores', main = 'Histograma para SFI-NIT')
hist(res_eli_nit$pvalue[res_eli_nit$baseMean > 1], xlab = 'p-valores', main = 'Histograma para ELI-NIT')
hist(res_eli_sfi$pvalue[res_eli_sfi$baseMean > 1], xlab = 'p-valores', main = 'Histograma para ELI-SFI')

# MA-plot
plotMA(res_sfi_nit, ylim=c(-4,4))
plotMA(res_eli_nit, ylim=c(-4,4))
plotMA(res_eli_sfi, ylim=c(-4,4))

# --- Anotación
library("AnnotationDbi")
library(org.Hs.eg.db)

my.annotation <- function(data){
  claves <- substr(rownames(data), 1, 15)
  data$symbol <- mapIds(org.Hs.eg.db, keys=claves, column="SYMBOL", keytype="ENSEMBL", multiVals="first")
  data$entrez <- mapIds(org.Hs.eg.db, keys=claves, column="ENTREZID", keytype="ENSEMBL", multiVals="first")
  return(data[complete.cases(data), ])
}

# Anotamos los genes más significativos
res_sfi_nit_Sig <- my.annotation(res_sfi_nit_Sig)
res_eli_nit_Sig <- my.annotation(res_eli_nit_Sig)
res_eli_sfi_Sig <- my.annotation(res_eli_sfi_Sig)

# Cuántos quedan
res_sfi_nit_Sig@nrows
res_eli_nit_Sig@nrows
res_eli_sfi_Sig@nrows









# Criterio más estricto, incrementando el umbral de log2 fold change
#res_sfi_nit_LFC1 <- results(dds, contrast=c("Group",'SFI','NIT'), lfcThreshold=1)
#res_eli_nit_LFC1 <- results(dds, contrast=c("Group",'ELI','NIT'), lfcThreshold=1)
#res_eli_sfi_LFC1 <- results(dds, contrast=c("Group",'ELI','SFI'), lfcThreshold=1)

# Genes up-regulated
#head(res_sfi_nit_Sig[ order(res_sfi_nit_Sig$log2FoldChange, decreasing = TRUE), ])