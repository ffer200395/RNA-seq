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
my.targets <- targets %>% filter(molecular_data_type == 'RNA Seq (NGS)') %>% group_by(Group) %>% sample_n(8)

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
ddsMat <- DESeqDataSetFromMatrix(countData = my.counts, colData = my.targets, design = ~ Group)

# Ejecutamos el modelo DESeq2
dds <- DESeq(ddsMat)

# --- Preprocesado de los datos: filtraje y normalización

# Nos quedamos con aquellos genes que tengan al menos 10 reads en total
nrow(dds)
dds <- dds[ rowSums(counts(dds)) >= 10, ]
nrow(dds)

# Número de reads por muestra expresado en millones
sort(colSums(my.counts))/1e6

# Normalizamos
my.counts.norm <- log2(counts(dds, normalized=TRUE)+1)

# Mostramos los datos normalizados versus sin normalizar
par(mfrow=c(1,2))
plot(my.counts[,1:2], main='Muestra 1 vs muestra 2 sin normalizar')
plot(my.counts.norm[,1:2],  main='Muestra 1 vs muestra 2 tras normalizar')
par(mfrow=c(1,2))
boxplot(log2(counts(dds)+1), main='Muestras sin normalizar', names=my.targets$SRA_Sample, las=2)
boxplot(my.counts.norm, main='Muestras normalizadas', names=my.targets$SRA_Sample, las=2)











