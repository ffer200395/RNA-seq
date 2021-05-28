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
vsd <- varianceStabilizingTransformation(dds)

# Scatter plots
par(mfrow=c(2,2))
plot(log2(counts(dds)[,1:2] + 1), cex=.1)
plot(log.norm.counts[,1:2], cex=.1)
plot(assay(rld)[,1:2], cex=.1)
plot(assay(vsd)[,1:2], cex=.1)

# Boxplots
par(mfrow=c(2,2))
boxplot(log2(counts(dds) + 1), names=my.targets$SRA_Sample, las=2)
boxplot(log.norm.counts, names=my.targets$SRA_Sample, las=2)
boxplot(assay(rld), names=my.targets$SRA_Sample, las=2)
boxplot(assay(vsd), names=my.targets$SRA_Sample, las=2)

# PCA

# Cluster

# Heatmap

# --- Identificación de genes diferencialmente expresados

# Ejecutamos el modelo DESeq2
dds <- DESeq(dds)

# Hacemos los contrastes entre grupos
res_sfi_nit <- results(dds, contrast=c("Group",'SFI','NIT'))
res_eli_nit <- results(dds, contrast=c("Group",'ELI','NIT'))
res_eli_sfi <- results(dds, contrast=c("Group",'ELI','SFI'))

# Número de genes con p-val ajustado menor de 0.1
table(res_sfi_nit$padj < 0.1)
table(res_eli_nit$padj < 0.1)
table(res_eli_sfi$padj < 0.1)

# Criterio más estricto, incrementando el umbral de log2 fold change
res_sfi_nit_LFC1 <- results(dds, contrast=c("Group",'SFI','NIT'), lfcThreshold=1)
res_eli_nit_LFC1 <- results(dds, contrast=c("Group",'ELI','NIT'), lfcThreshold=1)
res_eli_sfi_LFC1 <- results(dds, contrast=c("Group",'ELI','SFI'), lfcThreshold=1)
table(res_sfi_nit_LFC1$padj < 0.1)
table(res_eli_nit_LFC1$padj < 0.1)
table(res_eli_sfi_LFC1$padj < 0.1)

# If we consider a fraction of 10% false positives acceptable, we can consider all genes 
# with an adjusted p value below 10% = 0.1 as significant. How many such genes are there?
sum(res_sfi_nit$padj < 0.1, na.rm=TRUE)
sum(res_eli_nit$padj < 0.1, na.rm=TRUE)
sum(res_eli_sfi$padj < 0.1, na.rm=TRUE)

# Genes más significativos
res_sfi_nit_Sig <- subset(res_sfi_nit, padj < 0.1)
res_eli_nit_Sig <- subset(res_eli_nit, padj < 0.1)
res_eli_sfi_Sig <- subset(res_eli_sfi, padj < 0.1)

# Genes up-regulated
head(res_sfi_nit_Sig[ order(res_sfi_nit_Sig$log2FoldChange, decreasing = TRUE), ])

# MA-plot

# p-val hist

# Plot of counts

# Heatmap of top genes

# Anotación
library("AnnotationDbi")
columns(org.Hs.eg.db)

res_sfi_nit_Sig$symbol <- mapIds(org.Hs.eg.db,
                     keys=row.names(res_sfi_nit_Sig),
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")

res_sfi_nit_Sig$entrez <- mapIds(org.Hs.eg.db,
                    keys=row.names(res_sfi_nit_Sig),
                    column="ENTREZID",
                    keytype="ENSEMBL",
                    multiVals="first")













