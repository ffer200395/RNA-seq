## --- Paquetes necesarios para el análisis
if (!require(DESeq2))BiocManager::install("DESeq2")

## --- Cargamos los datos
counts <- read.csv('data/counts.csv', sep = ';', header = T, row.names = 1)
targets <- read.csv('data/targets.csv', sep = ',', header = T)

# --- Creación del subconjunto del análisis (usar datos de expresión RNA-Seq) y preparación de los datos
library(dyplr)
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

# Eliminamos filas del objeto DESeqDataSet que tienen un conteo de 0 o 1 en las muestras
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 1, ]
nrow(dds)

# Normalizamos...........








