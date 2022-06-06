
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biovizBase")

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")
#
#
#
#
#
# create a gene activity matrix
activity.matrix <- CreateGeneActivityMatrix(
  peak.matrix = peaks,
  annotation.file = "/Users/aaaaaa/Downloads/Homo_sapiens.GRCh37.82.gtf",
  seq.levels = c(1:22, "X", "Y"),
  upstream = 2000,
  verbose = TRUE
) 


# before installing signac: 
# Install bioconductor
# Install bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()

# To automatically install Bioconductor dependencies
setRepositories(ind=1:2)

install.packages("Signac")
#uninstall.packages

# the packages needed:  
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)
set.seed(1234)


#------------------------------------------------------------------------> Gene activity matrix: 

# source: https://satijalab.org/signac/articles/pbmc_vignette.html#create-a-gene-activity-matrix-1 
#-----------------------------------> seurat object:


counts <- Read10X_h5(filename = "/Users/aaaaaa/Downloads/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5")
metadata <- read.csv(
  file = "/Users/aaaaaa/Downloads/atac_v1_pbmc_10k_singlecell.csv",
  header = TRUE,
  row.names = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg19',
  fragments = '/Users/aaaaaa/Downloads/atac_v1_pbmc_10k_fragments.tsv.gz',
  min.cells = 10,
  min.features = 200
)

pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)


pbmc

pbmc[['peaks']]


granges(pbmc)


# 2nd part of pre processing: 

# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'

# add the gene information to the object
Annotation(pbmc) <- annotations

#-----------------------------------------------------

gene.activities <- GeneActivity(pbmc)


# add the gene activity matrix to the Seurat object as a new assay and normalize it
pbmc[['RNA']] <- CreateAssayObject(counts = gene.activities)
pbmc <- NormalizeData(
  object = pbmc,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(pbmc$nCount_RNA)
)

DefaultAssay(pbmc) <- 'RNA'
# filter
# pbmc <- subset(pbmc, subset = nCount_ATAC > 5000)

pbmc$tech <- "atac"

pbmc.atac

pbmc

#



# converting to file for python processing:

library(SeuratDisk)

SaveH5Seurat(pbmc, filename = "atacrnapbmc.h5Seurat")
Convert("atacrnapbmc.h5Seurat", dest = "h5ad")

# ?SaveH5Seurat


annotations

pbmc

gene.activities

pbmc[['RNA']]
pbmc$gene_id <- pbmc[['RNA']]

pbmc$gene_id


# stuff to try: 
#meta <- meta[colnames(pbmc), ]
#pbmc <- AddMetaData(pbmc, metadata = meta)


# how to load the seurat object: 
atpbmc <- LoadH5Seurat("atacrnapbmc.h5Seurat")

atpbmc

onlyrna <- atpbmc[['RNA']]

onlyrna


if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("mojaveazure/seurat-disk")

# loading the data: 
bcell1_dataset <- readRDS("/Users/aaaaaa/Downloads/bcell1data.rds")

bcell1_dataset[[1]] # this works!

library(SeuratDisk)

SaveH5Seurat(bcell1_dataset[[1]], filename = "bcells1dataset.h5Seurat")
Convert("/Users/aaaaaa/Desktop/bcells1dataset.h5Seurat", dest = "h5ad")

# adata succesfully created. 

# how to load:
brain2 <- LoadH5Seurat("anterior1.h5Seurat")

summary(bcell1_dataset)

rm(list = ls())


cvdata <- readRDS("/Users/aaaaaa/Downloads/corr_10xlymphoma.rds")
cvdata

cvdata@assays$ATAC

cvdata@assays$chromvar

cvdata@assays$chromvar

head(cvdata@assays[["chromvar"]]@data, 1)

# summary(cvdata@assays[["chromvar"]]@data)

rownames(cvdata[["chromvar"]])

colnames(cvdata[["chromvar"]]) 

cvdata@assays$chromvar@counts

library(stringr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.NCBI.GRCh38)

library(ggplot2)



cvdata@assays[["chromvar"]]@data[5, ] 

cvdata@assays[["chromvar"]]@data[, 2] 

cvdata@assays[["chromvar"]]@data[1:5, 1:5] 


GetAssayData(cvdata@assays[["chromvar"]], slot=c("data"))

mtrx <- cvdata@assays[["chromvar"]]@data


write.table(mtrx, file = "/Users/aaaaaa/Downloads/10xlymphomamtrx.csv", sep = ",")


h5data <- Read10X_h5("/Users/aaaaaa/Downloads/filtered_feature_bc_matrix.h5")

h5data

library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)
library(dplyr)
library(EnsDb.Hsapiens.v86)

# extract RNA and ATAC data
rna_counts <- h5data$`Gene Expression`
atac_counts <- h5data$Peaks


# Create Seurat object
bcell <- CreateSeuratObject(counts = rna_counts)
bcell[["percent.mt"]] <- PercentageFeatureSet(bcell, pattern = "^MT-")

# Now add in the ATAC-seq data
# we'll only use peaks in standard chromosomes
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

frag.file <- "/Users/aaaaaa/Downloads/atac_fragments.tsv.gz"
chrom_assay <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = frag.file,
  min.cells = 10,
  annotation = annotations
)
bcell[["ATAC"]] <- chrom_assay


# make save file here in case: 
saveRDS(bcell, file = "/Users/aaaaaa/Downloads/bcellsbeforepwm.rds")
# works so far!

DefaultAssay(bcell) <- "ATAC"

#pbmc

pwm_set <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

id2sym <- data.frame(
  id=sapply(1:length(pwm_set),  function(x) ID(pwm_set[[x]]) ),
  symbol=sapply(1:length(pwm_set),  function(x) name(pwm_set[[x]]) )
)
# Scan the DNA sequence of each peak for the presence of each motif, and create a Motif object
#DefaultAssay(h5data) <- "ATAC"
# pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
motif.matrix <- CreateMotifMatrix(features = granges(bcell), pwm = pwm_set, genome = 'hg38', use.counts = FALSE)
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)
bcell <- SetAssayData(bcell, assay = 'ATAC', slot = 'motifs', new.data = motif.object)

# save point: 
saveRDS(bcell, file = "/Users/aaaaaa/Downloads/alltfbeforechromvarbcells.rds")
# save successful

# Note that this step can take 30-60 minutes 
bcell <- RunChromVAR(
  object = bcell,
  genome = BSgenome.Hsapiens.UCSC.hg38
)



bcell@assays[["chromvar"]]@data[2, 3]

#cvdata@assays[["chromvar"]]@data


# save point: 
saveRDS(bcell, file = "/Users/aaaaaa/Downloads/afterchromvarbcells.rds")



mtxx <- bcell@assays[["chromvar"]]@data


write.table(mtxx, file = "/Users/aaaaaa/Downloads/bcell10ktf.csv", sep = ",")


bcell@assays[["chromvar"]]@meta.features


#rownames(bcell@assays[["chromvar"]]@meta.features)

#motif.name <- ConvertMotifID(bcell, rownames(bcell@assays[["chromvar"]]@meta.features))

#rownames(bcell@assays[["chromvar"]]@meta.features)


#motifs <- getJasparMotifs()

pwm_set 

annotations@elementMetadata@listData$gene_name

id2sym

write.csv(id2sym, "/Users/aaaaaa/Downloads/all746tfs.csv")






