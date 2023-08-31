suppressMessages({
  library(Seurat)
  library(dplyr)
  library(data.table)
  library(tools)
  library(stringr)
  library(futile.logger)
  library(devtools)
  library(ouija)
})


### Prepare ###
all.Tcell <- readRDS("../data/out/custom_merge_tcell.rds")

## Subset CD8_exhausted, cytotoxic and naive T cells
Idents(all.Tcell) <- all.Tcell$Tcell_cluster
CD8 <- subset(all.Tcell, idents = c("T-naive", "T-cytotoxic-1", "T-cytotoxic-2","T-cytotoxic-3","T-cytotoxic-4", "T-CD8-exhausted"))

# Subsample max 800 cells per sample
Idents(CD8) <- CD8$orig.ident
set.seed(11)
CD8.sub.800 <- subset(x = CD8, downsample=800)
CD8.sub.800 <- NormalizeData(CD8.sub.800)

# CUSTOM: Recalculate ouija pseudotime
## Load genes from publication
ouija_genes <- read.table(text = "CCR7	switch
IL7R	switch
SELL	switch
CD69	switch
PDCD1	switch
CXCL13	switch
LAG3	switch
HAVCR2	switch
CD27	switch
CD38	switch
TIGIT	switch
CTLA4	switch
ENTPD1	switch
GZMB	switch
FASLG	switch
TCF7	transient
KLRG1	transient
CX3CR1	transient
FCGR3A	transient
PRF1	transient
TNF	transient
IFNG	transient
GZMK	transient")
colnames(ouija_genes) <- c("gene",	"putative_response_type")

## Get relevant matrix from Seurat
count_matrix <- LayerData(CD8.sub.800, layer = "data")
ouija_matrix <- count_matrix[rownames(count_matrix) %in% ouija_genes$gene, ]

## Reorder gene table to get correct gene types to ouija
ouija_genes <- ouija_genes[match(rownames(ouija_matrix), ouija_genes$gene), ]

# Calculate ouija
options(mc.cores = parallel::detectCores())
oui <- ouija(t(as.matrix(ouija_matrix)), ouija_genes$putative_response_type, iter = 6000) 

saveRDS(oui, "../data/out/ouija.rds")