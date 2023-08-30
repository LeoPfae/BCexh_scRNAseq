library(Seurat)
library(futile.logger)

flog.info("Starting merge")
pilot_1 <- readRDS("../data/doublet_output/TBB330_singlet.rds")
pilot_2 <- readRDS("../data/doublet_output/TBB338_singlet.rds")
flog.info("Loaded pilot files")
run1_1 <- readRDS("../data/doublet_output/TBB075_singlet.rds")
run1_2 <- readRDS("../data/doublet_output/TBB102_singlet.rds")
run1_3 <- readRDS("../data/doublet_output/TBB111_singlet.rds")
run1_4 <- readRDS("../data/doublet_output/TBB129_singlet.rds")
run1_5 <- readRDS("../data/doublet_output/TBB165_singlet.rds")
run1_6 <- readRDS("../data/doublet_output/TBB171_singlet.rds")
run1_7 <- readRDS("../data/doublet_output/TBB214_singlet.rds")
run1_8 <- readRDS("../data/doublet_output/TBB226_singlet.rds")
flog.info("Loaded run1 files")
flog.info("Merging pilot and run1 files")
run1.pilot <- merge(x = pilot_1, y = c(pilot_2, run1_1, run1_2, run1_3, run1_4, run1_5, run1_6, run1_7, run1_8), add.cell.ids = c("TBB330", "TBB338", "TBB075", "TBB102", "TBB111", "TBB129", "TBB165", "TBB171", "TBB214", "TBB226"))
flog.info("Finished merging pilot and run1 files")
flog.info("Saving merged files")
saveRDS(run1.pilot, "../data/doublet_output/run1_pilot_merge.RDS")
flog.info("Saved merged files")

rm(pilot_1)
rm(pilot_2)
rm(run1_1)
rm(run1_2)
rm(run1_3)
rm(run1_4)
rm(run1_5)
rm(run1_6)
rm(run1_7)
rm(run1_8)

in.path = "../../BCexh_scRNAseq/data/doublet_output/"
out.path = "../../BCexh_scRNAseq/data/out/"

TBB035 <- readRDS(file = paste0(in.path, 'TBB035_singlet.rds'))
TBB184 <- readRDS(file = paste0(in.path, 'TBB184_singlet.rds'))
TBB011 <- readRDS(file = paste0(in.path, 'TBB011_singlet.rds'))
TBB212 <- readRDS(file = paste0(in.path, 'TBB212_singlet.rds'))

merge1 <- merge(x=run1.pilot, y = TBB011, add.cell.ids = c("", "TBB011"))
merge2 <- merge(x=merge1, y = TBB035, add.cell.ids = c("", "TBB035"))
merge3 <- merge(x=merge2, y = TBB184, add.cell.ids = c("", "TBB184"))
all.merged <- merge(x=merge3, y = TBB212, add.cell.ids = c("", "TBB212"))

#save
saveRDS(all.merged, file = paste0(out.path, "complete_merged_preFilter.rds"))


all.merged$run <- "Run1"
all.merged$run[all.merged$batch == "pilot"] <- "pilot"
all.merged$run[all.merged$batch == "B5" | all.merged$batch == "B6"] <- "Run2"
table(all.merged$run)

# store mitochondrial percentage in object meta data
all.merged <- PercentageFeatureSet(all.merged, pattern = "^MT-", col.name = "percent.mt")
all.merged$percent.mito <- NULL
all.merged$DF_hi.lo <- NULL

all.merged <- subset(x = all.merged, subset = nFeature_RNA > 200 & nFeature_RNA < 7500 & percent.mt < 20  & nCount_RNA < 75000)

saveRDS(all.merged, file = paste0(out.path, "complete_merged_postFilter.rds"))
