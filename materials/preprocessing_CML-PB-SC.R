library(Seurat)
library(harmony)
library(ggplot2)
library(stringr)
library(scico)
library(dplyr)

# Read in data and subset CD8 T cells of baseline samples
# object is available from Zenodo: https://zenodo.org/doi/10.5281/zenodo.10947495
# sc_cml <- readRDS("sc_cml.RDS")
sc_cml$pubclusters <- Idents(sc_cml)
sc_cml <- subset(sc_cml, (timepoint == "baseline")  &
                   (CD3E > 1 | CD8A > 1 | CD8B > 1) &
                   (pubclusters %in% c("CD8 Tem/emra", "CD8 Tem/rm", 
                                       "CD8 Tcm/n", "Unconv: NKT", 
                                       "Unconv: MAIT", "Unconv: Tgd")))

# Name patients consistently
sc_cml$ID <- sc_cml$orig.ident
u.o.i <- unique(sc_cml$orig.ident)
for (i in 1:length(u.o.i)){
  sc_cml$ID[sc_cml$ID==u.o.i[i]] <- paste0("CML-PB-SC-",as.character(i))
}

# normalize data
sc_cml <- NormalizeData(sc_cml, normalization.method = "LogNormalize", scale.factor = 10000)

# 2000 most highly variable genes not including TCR clonality genes
sc_cml <- FindVariableFeatures(sc_cml, selection.method = "vst", nfeatures = 2500)
hvg <- VariableFeatures(sc_cml)
clonality_genes <- grep("^(TRAV|TRAJ|TRBV|TRBJ|TRBD|TRGV|TRGJ|TRDV|TRDJ)", 
                        rownames(sc_cml), value = T)
hvg  <- hvg[!(hvg %in% clonality_genes)][1:2000]
VariableFeatures(sc_cml) <- hvg

# Phase correction, scaling, batch correction
sc_cml <- ScaleData(sc_cml, features = VariableFeatures(sc_cml), verbose=F)

sc_cml <- CellCycleScoring(sc_cml, s.features = cc.genes$s.genes, 
                           g2m.features = cc.genes$g2m.genes, set.ident = T)
sc_cml <- ScaleData(sc_cml, features = VariableFeatures(sc_cml),verbose=F,
                    vars.to.regress = c("S.Score", "G2M.Score")) 

sc_cml <- RunPCA(sc_cml, features=VariableFeatures(sc_cml), nPCs = 50,verbose = F)
sc_cml <- RunUMAP(sc_cml, dims = 1:20) 
sc_cml <- RunHarmony(sc_cml, "patient", nPCs = 20)

# UMAP and clustering
sc_cml <- RunUMAP(sc_cml,reduction = "harmony",reduction.name = "HUMAP20",dims=1:20)
sc_cml <- FindNeighbors(sc_cml, reduction = "harmony", dims = 1:20,verbose = F)
sc_cml <- FindClusters(object = sc_cml,resolution = 0.2, verbose=F)

# Remove cluster 5 and name remaining clusters
sc_cml <- subset(sc_cml, RNA_snn_res.0.2 != 5)
sc_cml$Clusters <- factor(recode(sc_cml$RNA_snn_res.0.2, 
                                 `0`="NK-like Temra", `1`="Tem/rm", `2`="Tn", `3`= "MAIT-like", `4`="IFN CTL"),
                          levels=c("Tn","Tem/rm","NK-like Temra","IFN CTL", "MAIT-like"))
