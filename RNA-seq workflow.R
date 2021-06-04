## --- Paquetes necesarios para el análisis
if (!require(DESeq2))BiocManager::install("DESeq2")
if (!require(GOstats))BiocManager::install("GOstats")
if (!(require(dplyr))) install.packages("dplyr")
if (!require(VennDiagram)) install.packages("VennDiagram")
if (!require(pheatmap)) install.packages("pheatmap")

## --- Cargamos los datos
counts <- read.csv('data/counts.csv', sep = ';', header = T, row.names = 1)
targets <- read.csv('data/targets.csv', sep = ',', header = T)

# --- Creación del subconjunto del análisis (usar datos de expresión RNA-Seq) y preparación de los datos
library(dplyr)
set.seed(1234)
# Elegimos 10 muestras de cada grupo
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
par(mfrow=c(2,2),mar=c(2,2,2,2))
plot(log2(counts(dds)[,1:2] + 1), cex=.1, main = 'Sin normalizar')
plot(log.norm.counts[,1:2], cex=.1, main = 'Log norm')
plot(assay(rld)[,1:2], cex=.1, main = 'Regularized Log')
plot(assay(vsd)[,1:2], cex=.1, main = 'VST')

# Boxplots
par(mfrow=c(2,2),mar=c(2,2,2,2))
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
library(ggplot2)
plotPCA(DESeqTransform(SummarizedExperiment(log.norm.counts, colData=colData(dds))), intgroup="Group")+ggtitle('PCA para log2 norm')+theme(plot.title = element_text(hjust = 0.5))
plotPCA(rld, intgroup="Group")+ggtitle('PCA para Rlog')+theme(plot.title = element_text(hjust = 0.5))
plotPCA(vsd, intgroup="Group")+ggtitle('PCA con VST')+theme(plot.title = element_text(hjust = 0.5))

# --- Identificación de genes diferencialmente expresados

# Ejecutamos el modelo DESeq2
dds <- DESeq(dds)

# Hacemos los contrastes entre grupos
res_sfi_nit <- results(dds, contrast=c("Group",'SFI','NIT'))
res_eli_nit <- results(dds, contrast=c("Group",'ELI','NIT'))
res_eli_sfi <- results(dds, contrast=c("Group",'ELI','SFI'))

# Número de genes con p-val ajustado menor de 0.1 
#(Si aceptamos un 10% de falsos positivos entonces los genes con p-val<0.1 son significativos)
sum(res_sfi_nit$padj < 0.1, na.rm=TRUE)
sum(res_eli_nit$padj < 0.1, na.rm=TRUE)
sum(res_eli_sfi$padj < 0.1, na.rm=TRUE)

# Genes más significativos
res_sfi_nit_Sig <- subset(res_sfi_nit, padj < 0.1)
res_eli_nit_Sig <- subset(res_eli_nit, padj < 0.1)
res_eli_sfi_Sig <- subset(res_eli_sfi, padj < 0.1)

# Counts plot para visualizar el conteo para un gen en concreto
plotCounts(dds, gene = rownames(res_sfi_nit_Sig)[which.min(res_sfi_nit_Sig$padj)], intgroup=c("Group"))

# p-val hist
par(mfrow=c(1,3))
hist(res_sfi_nit_Sig$padj, xlab = 'p-val adj.', main = 'Histograma para SFI-NIT', freq=F, ylim = c(0,50))
hist(res_eli_nit_Sig$padj, xlab = 'p-val adj.', main = 'Histograma para ELI-NIT', freq=F, ylim = c(0,50))
hist(res_eli_sfi_Sig$padj, xlab = 'p-val adj.', main = 'Histograma para ELI-SFI', freq=F, ylim = c(0,50))
 
# MA-plot
par(mfrow=c(1,3),mar=c(2,2,2,2))
plotMA(res_sfi_nit, ylim=c(-4,4), main = 'MA plot SFI-NIT')
plotMA(res_eli_nit, ylim=c(-4,4), main = 'MA plot ELI-NIT')
plotMA(res_eli_sfi, ylim=c(-4,4), main = 'MA plot ELI-SFI')

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

# --- Patrones de expresión y agrupación de las muestras

# Volcano plot
get.volcano <- function(data, pval, fc, title){
  df = as.data.frame(data@listData)
  # Inicialmente asumimos que no se expresan
  df$diffexpressed <- "NO"
  # Si el pvalor es inferior a pval y el foldchange es mayor de fc se considera UP
  df$diffexpressed[df$log2FoldChange > fc & df$pvalue < pval] <- "UP"
  # Si el foldchange es menor que -fc entonces se considera DOWN
  df$diffexpressed[df$log2FoldChange < -fc & df$pvalue < pval] <- "DOWN"
  # Creamos el volcano plot
  ggplot(data=df, aes(x=log2FoldChange, y=-log10(padj), col=diffexpressed)) + geom_point() + geom_text(aes(label=symbol),check_overlap=T)+ ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
}

get.volcano(res_sfi_nit_Sig, 0.05, 2.5, 'Volcano plot SFI_NIT')
get.volcano(res_eli_nit_Sig, 0.05, 2.5, 'Volcano plot ELI_NIT')
get.volcano(res_eli_sfi_Sig, 0.05, 2.5, 'Volcano plot ELI_SFI')

# Diagrama de Venn para visualizar el número de genes en común
library(VennDiagram)
get.num.inter <- function(data1,data2){return(length(intersect(as.array(row.names(data1)),as.array(row.names(data2)))))}
dev.off()
draw.triple.venn(area1 = nrow(res_eli_nit_Sig), area2 = nrow(res_eli_sfi_Sig), area3 = nrow(res_sfi_nit_Sig), 
                 n12 = get.num.inter(res_eli_nit_Sig,res_eli_sfi_Sig), 
                 n23 = get.num.inter(res_eli_sfi_Sig,res_sfi_nit_Sig), 
                 n13 = get.num.inter(res_eli_nit_Sig,res_sfi_nit_Sig), 
                 n123 = length(intersect(intersect(as.array(row.names(res_eli_nit_Sig)),as.array(row.names(res_sfi_nit_Sig))),as.array(row.names(res_eli_sfi_Sig)))), 
                 category =  c("eli_nit" , "eli_sfi" , "sfi_nit"), lty = "blank", 
                 fill = c("skyblue", "pink1", "mediumorchid"))

# Heatmap
library(pheatmap)
# Ordenamos los datos de mayor a menor conteo normalizado y seleccionamos los primeros 20
resSort <- order(rowVars(assay(vsd)), decreasing = TRUE)[c(1:20)]
# Obtenemos la id de esos genes
topgenes <-  row.names(assay(vsd))[resSort]
# Seleccionamos los conteos normalizados de esos genes
mat <- assay(vsd)[topgenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(vsd)[, c("Group")])
rownames(df) <- colnames(mat)
pheatmap(mat, annotation_col = df)

# --- Análisis de significación biológica
library(GOstats)

# Tabla a estudiar
topTable <- res_eli_nit_Sig

# Almacenamos todo el conjunto de genes y los genes más diferenciados según un umbral del p-valor
entrezUni <- unique(my.annotation(counts)$entrez) # Usamos la tabla counts que contiene todos los genes del análisis
whichgenes <- which(topTable$padj<0.01)
geneIds <- unique(topTable[whichgenes,]$entrez)

# Definimos los hiperparámetros para realizar las búsquedas
paramsGO <- new("GOHyperGParams", geneIds =geneIds, universeGeneIds = entrezUni, annotation = "org.Hs.eg.db", ontology = "BP", testDirection = "over", pvalueCutoff=0.001) 

# Ejecutamos el análisis
hypGO <- hyperGTest(paramsGO)

# Almacenamos los resultados del análisis en un informe html
compare <- 'ELIvsNIT'
htmlReport(hypGO, file = paste0(compare,'_GO.html'))

# --- Términos
library(tm)
# Palabras que no queremos incluir
my.stopWords <- c('of', 'the', 'to', 'in', 'and', 'terms')

# Obtenemos los términos obtenidos en el análisis y los filtramos
raw_terms <- removeWords(paste(tolower(summary(hypGO)$Term),collapse=""),my.stopWords) 
terms <- as.array(strsplit(raw_terms, "\\W"))

# Ordenamos los terminos por su aparición
topTerms <- as.data.frame(sort(table(terms), decreasing=TRUE)[c(2:11)])

# Ploteamos el Top 10
barplot(height=topTerms$Freq, names=topTerms$terms, col="#CCFF66",las=2,horiz = F,main = 'Términos Top 10',ylim = c(0,200))