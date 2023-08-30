if (!require(ggplot2)) install.packages('ggplot2')
library(ggplot2)

if (!require(dplyr)) install.packages('dplyr')
library(dplyr)

if (!require(Seurat)) install.packages('Seurat')
library(Seurat)

if (!require(tidyverse)) install.packages('tidyverse')
library(tidyverse)

if (!require(data.table)) install.packages('data.table')
library(data.table)

if (!require(sctransform)) install.packages('sctransform')
library(sctransform)

library(magrittr)
library(clustree)
library(futile.logger)


in.path = "../../BCexh_scRNAseq/data/doublet_output/"
out.path = "../../BCexh_scRNAseq/data/out/"


flog.info("Starting SCTransform")
all.merged <- readRDS(file = paste0(out.path, "complete_merged_postFilter.rds"))

future::plan("multisession", workers = 2)
options(future.globals.maxSize = 5* 1000 * 1024^2)

# run sctransform
all.merged <- SCTransform(all.merged, verbose = TRUE)

saveRDS(all.merged, file = paste0(out.path, "complete_merged_sctransformed.rds"))

#read
all.merged <- readRDS (file = paste0(out.path, "complete_merged_sctransformed.rds"))


flog.info("Starting DimReduc analysis")
all.merged <- RunPCA(object = all.merged, features = VariableFeatures(object = all.merged), verbose = FALSE)
all.merged <- FindNeighbors(object = all.merged, dims = 1:27)
all.merged<- FindClusters(object = all.merged, resolution = 2)
all.merged <- RunUMAP(object = all.merged, dims = 1:27)

saveRDS(all.merged, file = paste0(out.path, "merged_complete_umap.rds"))
all.merged <- readRDS(file = paste0(out.path, "merged_complete_umap.rds"))

flog.info("Starting cell cluster analysis")
Idents(all.merged) <-all.merged$SCT_snn_res.2
cluster.ids <- read.csv(file = "../03_Additional_files/cluster_celltypes_res2_v2.csv")

new.cluster.ids <- cluster.ids$cell.type
new.cluster.ids <- as.character(new.cluster.ids)
names(x = new.cluster.ids) <- levels(x = all.merged)
all.merged <- RenameIdents(object = all.merged, new.cluster.ids)

saveRDS(all.merged, file = paste0(out.path, "merged_complete_inclCelltype.rds"))




flog.info("Starting immune cluster subsetting")
run1.immune <- subset(x = all.merged, idents = c("T/NK cell", "myeloid", "B cell", "granulocyte", "cDC", "pDC", "plasma cell"))
#Test
DimPlot(object = run1.immune, reduction = 'umap', label = TRUE, pt.size = 0.5)
#Save Seurat object
saveRDS(run1.immune, file = paste0(out.path, "immune.rds"))
rm(run1.immune)

## T/NK cell subset
run1.TNK <- subset(x = all.merged, idents = "T/NK cell")
saveRDS(run1.TNK, file = paste0(out.path, "T_NK_cells.rds"))
rm(run1.TNK)

## Myeloid subset (incl. DCs)
run1.myeloid <- subset(x = all.merged, idents = c("myeloid", "cDC", "pDC"))
saveRDS(run1.myeloid, file = paste0(out.path, "myeloid_inclDC.rds"))
rm(run1.myeloid)

## B cell subset
run1.B <- subset(x = all.merged, idents = "B cell")
saveRDS(run1.B, file = paste0(out.path, "B_cells.rds"))
rm(run1.B)

## Plasma cell subset
run1.PC <- subset(x = all.merged, idents = "plasma cell")
saveRDS(run1.PC, file = paste0(out.path, "plasma_cells.rds"))
rm(run1.PC)

## Granulocyte subset
run1.gran <- subset(x = all.merged, idents = "granulocyte")
saveRDS(run1.gran, file = paste0(out.path, "granulocytes.rds"))
rm(run1.gran)

## tumor cell subset
run1.ep <- subset(x = all.merged, idents = "epithelial")
saveRDS(run1.ep, file = paste0(out.path, "epithelial.rds"))
rm(run1.ep)

## fibroblast subset
run1.fib <- subset(x = all.merged, idents = "fibroblast")
saveRDS(run1.fib, file = paste0(out.path, "fibroblast.rds"))
rm(run1.fib)

## endothelial cell subset
run1.endo <- subset(x = all.merged, idents = "endothelial")
saveRDS(run1.endo, file = paste0(out.path, "endothelial.rds"))
rm(run1.endo)

flog.info("Finished")