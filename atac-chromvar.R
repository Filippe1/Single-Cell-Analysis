
#
#
#
#
#
# the packages needed:  
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)
set.seed(1234)


#------------------------------------------------------------------------> 

library(SeuratDisk)
library(stringr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Hsapiens.UCSC.hg38)
library(BSgenome.Hsapiens.NCBI.GRCh38)
library(ggplot2)

# Data from Cell Ranger: 

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

# prepare to run chromvar:
DefaultAssay(bcell) <- "ATAC"


pwm_set <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

id2sym <- data.frame(
  id=sapply(1:length(pwm_set),  function(x) ID(pwm_set[[x]]) ),
  symbol=sapply(1:length(pwm_set),  function(x) name(pwm_set[[x]]) )
)
# Scan the DNA sequence of each peak for the presence of each motif, and create a Motif object
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


# save point: 
saveRDS(bcell, file = "/Users/aaaaaa/Downloads/afterchromvarbcells.rds")

# converting tf activities to csv file: 

mtxx <- bcell@assays[["chromvar"]]@data


write.table(mtxx, file = "/Users/aaaaaa/Downloads/bcell10ktf.csv", sep = ",")
# can convert csv to python dataframe and analyze in Scanpy!
            
            

            

bcell@assays[["chromvar"]]@meta.features


# might be useful to have: 
write.csv(id2sym, "/Users/aaaaaa/Downloads/all746tfs.csv")






