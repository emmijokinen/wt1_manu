# AML-BM-SC data processing

library(Seurat)
library(harmony)
library(ggplot2)
library(stringr)
library(dplyr)
library(Azimuth)


# Read in data and subset CD8 T cells of baseline samples
# preprocessed object is available from Zenodo: https://zenodo.org/doi/10.5281/zenodo.10947495
sc_aml <- readRDS("AML_CD8.rds")

sc_aml$pubclusters <- Idents(sc_aml)
sc_aml$cohort <- !grepl("_dg",sc_aml$orig.ident)
sc_aml <- subset(sc_aml, cohort)

# Name patients consistently
sc_aml$ID <- ""
u.o.i <- unique(sc_aml$orig.ident)
for (i in 1:length(u.o.i)){
  sc_aml$ID[sc_aml$orig.ident==u.o.i[i]] <- paste0("AML-BM-SC-",as.character(i))
}

# normalize data
sc_aml <- NormalizeData(sc_aml, normalization.method = "LogNormalize", scale.factor = 10000)

# 2000 most highly variable genes not including TCR clonality genes
sc_aml <- FindVariableFeatures(sc_aml, selection.method = "vst", nfeatures = 2500)
hvg <- VariableFeatures(sc_aml)
clonality_genes <- grep("^(TRAV|TRAJ|TRBV|TRBJ|TRBD|TRGV|TRGJ|TRDV|TRDJ)", 
                        rownames(sc_aml), value = T)
hvg  <- hvg[!(hvg %in% clonality_genes)][1:2000]
VariableFeatures(sc_aml) <- hvg

# Phase correction, scaling, batch correction
sc_aml <- ScaleData(sc_aml, features = VariableFeatures(sc_aml), verbose=F)

sc_aml <- CellCycleScoring(sc_aml, s.features = cc.genes$s.genes, 
                           g2m.features = cc.genes$g2m.genes, set.ident = T)
sc_aml <- ScaleData(sc_aml, features = VariableFeatures(sc_aml),verbose=F,
                    vars.to.regress = c("S.Score", "G2M.Score")) 

sc_aml <- RunPCA(sc_aml, features=VariableFeatures(sc_aml), nPCs = 50,verbose = F)
sc_aml <- RunUMAP(sc_aml, dims = 1:20) 
sc_aml <- RunHarmony(sc_aml, "patient", nPCs = 20)

# UMAP and clustering
sc_aml <- RunUMAP(sc_aml,reduction = "harmony",reduction.name = "HUMAP20",dims=1:20)
sc_aml <- FindNeighbors(sc_aml, reduction = "harmony", dims = 1:20,verbose = F)
sc_aml <- FindClusters(object = sc_aml,resolution = 0.5, verbose=F)

# Name clusters, combine clusters 4 and 5 as one Tn cluster
sc_aml$Clusters <- factor(recode(sc_aml$RNA_snn_res.0.5, 
                                 `0`="Temra", `1`="Teff", `2`="Tem/rm", `3`= "NK-like Temra", 
                                 `4`="Tn",`5`="Tn",`6`="IFN CTL"),
                          levels=c("Tn","Tem/rm","Teff","Temra","NK-like Temra","IFN CTL"))


##############################################################################
# AML-BM-SC: WT1 expression in HSPC, Mono, and DC cells
# Preprocessed object can be downloaded from Zenodo
seurat_object <- readRDS("sc_aml_all-cells.RDS")

# Subsetting based on Azimuth cell types and preprocessing
seurat_object <- RunAzimuth(query=seurat_object,reference="bonemarrowref")
seurat_bm <- subset(seurat_object,predicted.celltype.l1 %in% c("HSPC","Mono","DC"))
seurat_bm <- subset(seurat_bm, predicted.celltype.l2 %in% c("CD4 Memory", "CD8 Naive","CD8 Effector_2","Late Eryth"),invert=T)
seurat_bm <- ScaleData(seurat_bm, features = VariableFeatures(seurat_bm),verbose=FALSE) 
seurat_bm <- RunHarmony(seurat_bm, "ID")
seurat_bm <- RunUMAP(seurat_bm,reduction = "harmony",reduction.name = "HUMAP20",dims=1:20)



