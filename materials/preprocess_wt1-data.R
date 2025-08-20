library(Seurat)
library(dplyr)
library(ggplot2)
library(viridis)
library(harmony)
library(scico)

"%ni%" <- Negate("%in%")

addTCRinfo <- function(seurat_object, filename){
  
  # Read TCR table
  df.TCR <- read.csv(filename)
  
  # A and B chains separately
  df.TRA <- df.TCR[df.TCR$chain=="TRA",]
  df.TRB <- df.TCR[df.TCR$chain=="TRB",]
  
  # which cells have duplicate TRA or TRB ?
  # Count occurrences and find cells with two TRA or TRB
  n.occur <- data.frame(table(df.TRA$barcode))
  duplicate.A <- n.occur[n.occur$Freq > 1,]
  n.occur <- data.frame(table(df.TRB$barcode))
  duplicate.B <- n.occur[n.occur$Freq > 1,]
  
  # Remove duplicates from each
  df.TRA <- subset(df.TRA, (barcode %ni% duplicate.A$Var1) )
  df.TRB <- subset(df.TRB, (barcode %ni% duplicate.B$Var1) )
  
  # Merge into one table
  df.TRAB <- merge(df.TRA,df.TRB,by="barcode",suffixes=c(".A",".B"),all=TRUE)
  df.TRAB$chains <- paste(df.TRAB$chain.A,df.TRAB$chain.B,sep="_")
  
  ## Add TCR information to seurat object
  #Add CDR3 A+B and TCR A+B (CDR3+V+J) and mark available chains
  
  seurat_object$CDR3AB <- ""
  seurat_object$CDR3AB[df.TRAB$barcode] <- paste(df.TRAB$cdr3.A,df.TRAB$cdr3.B,sep="_")
  
  seurat_object$CDR3A <- ""
  seurat_object$CDR3A[df.TRAB$barcode] <- df.TRAB$cdr3.A
  
  seurat_object$CDR3B <- ""
  seurat_object$CDR3B[df.TRAB$barcode] <- df.TRAB$cdr3.B
  
  seurat_object$TCRA <- ""
  seurat_object$TCRA[df.TRAB$barcode] <- paste(df.TRAB$cdr3.A,df.TRAB$v_gene.A,df.TRAB$j_gene.A,sep="_")
  
  seurat_object$TCRB <- ""
  seurat_object$TCRB[df.TRAB$barcode] <- paste(df.TRAB$cdr3.B,df.TRAB$v_gene.B,df.TRAB$j_gene.B,sep="_")
  
  seurat_object$availableTCRchains <- ""
  seurat_object$availableTCRchains[df.TRAB$barcode] <- gsub("(NA_|_NA)","",df.TRAB$chains)
  
  return(seurat_object)
}


create_WT1object_with_TCR_data <- function(dir.data.matrix,file.vdj,project.name,donor_group,epitope){
  # Read in sample data
  seurat_object <- CreateSeuratObject(counts=Read10X(data.dir = dir.data.matrix), project = project.name, min.cells = 3, min.features = 200)
  
  #make labels
  seurat_object$designation = c(donor_group)
  seurat_object$epitope = c(epitope)
  
  # Add TCR data
  # As this data will be used for training a TCRGP model, 
  # we will only add TCRs for T cells that have a unique alpha and beta chain, 
  # as we want that all the TCRs in the training data are WT1 specific
  seurat_object <- addTCRinfo(seurat_object,file.vdj)
  
  seurat_object <- subset(seurat_object, availableTCRchains!="")
  
  return(seurat_object)
}


# Define correct locations of the raw data here

dir.VLDh <- "VLD_healthy/filtered_feature_bc_matrix"
vdj.file.VLDh <- "VLD_healthy/filtered_contig_annotations.csv"

dir.VLDaml <- "VLD_AML/filtered_feature_bc_matrix/"
vdj.file.VLDaml <- "VLD_AML/filtered_contig_annotations.csv"

dir.RMFh <- "RMF_healthy/filtered_feature_bc_matrix"
vdj.file.RMFh <- "RMF_healthy/filtered_contig_annotations.csv"

dir.RMFaml <- "RMF_AML/filtered_feature_bc_matrix"
vdj.file.RMFaml <- "RMF_AML/filtered_contig_annotations.csv"

gliph2.tcrs.VLD <- "selected_TCRs_VLD.csv"
gliph2.tcrs.RMF <- "selected_TCRs_RMF.csv"


# Load data and create objects for VLD from healthy and AML patients
scVLDh <- create_WT1object_with_TCR_data(dir.VLDh,vdj.file.VLDh,"scVLDh","Healthy","VLD")
scVLDaml <- create_WT1object_with_TCR_data(dir.VLDaml,vdj.file.VLDaml,"scVLDaml","AML","VLD")
# combine
scVLDcomb <- merge(scVLDh, scVLDaml,add.cell.ids = c("VLDh", "VLDaml"),project="VLDcomb")
scVLDcomb <- JoinLayers(scVLDcomb)
# mark GLIPH2 selected TCRs
df.GLIPH2_VLD <- read.csv(gliph2.tcrs.VLD)
scVLDcomb$GLIPH2_TCR <- ifelse(test = scVLDcomb$CDR3B %in% df.GLIPH2_VLD$CDR3B, "GLIPH2","NA")

# Load data and create objects for RMF from healthy and AML patients
scRMFh <- create_WT1object_with_TCR_data(dir.RMFh,vdj.file.RMFh,"scRMFh","Healthy","RMF")
scRMFaml <- create_WT1object_with_TCR_data(dir.RMFaml,vdj.file.RMFaml,"scRMFaml","AML","RMF")
# combine
scRMFcomb <- merge(scRMFh, scRMFaml,add.cell.ids = c("RMFh", "RMFaml"),project="RMFcomb")
scRMFcomb <- JoinLayers(scRMFcomb)
# mark GLIPH2 selected TCRs
df.GLIPH2_RMF <- read.csv(gliph2.tcrs.RMF)
scRMFcomb$GLIPH2_TCR <- ifelse(test = scRMFcomb$CDR3B %in% df.GLIPH2_RMF$CDR3B, "GLIPH2","NA")

# Merge VLD and RMF
sc_wt1 <- merge(scRMFcomb, scVLDcomb)
sc_wt1 <- JoinLayers(sc_wt1)

# normalize data
sc_wt1 <- NormalizeData(sc_wt1, normalization.method = "LogNormalize", 
                           scale.factor = 10000,verbose = F)

# 2000 Variable features without TCR clonality genes
sc_wt1 <- FindVariableFeatures(sc_wt1, selection.method = "vst", 
                                  nfeatures = 2500,verbose = F)
clonality_genes <- grep("^(TRAV|TRAJ|TRBV|TRBJ|TRBD|TRGV|TRGJ|TRDV|TRDJ)", 
                        rownames(sc_wt1), value = T)
hvg_WT1 <- VariableFeatures(sc_wt1)
hvg_WT1  <- hvg_WT1[!(hvg_WT1 %in% clonality_genes)][1:2000]
VariableFeatures(sc_wt1) <- hvg_WT1

# Cell cycle scoring and data scaling
sc_wt1 <- CellCycleScoring(sc_wt1, s.features = cc.genes$s.genes, 
                              g2m.features = cc.genes$g2m.genes, set.ident = F)
sc_wt1 <- ScaleData(sc_wt1, features = VariableFeatures(sc_wt1),
                       vars.to.regress = c("S.Score", "G2M.Score"),verbose=F) 
## PCA
sc_wt1 <- RunPCA(sc_wt1, features=VariableFeatures(sc_wt1), npcs = 20,verbose = F)
sc_wt1 <- RunUMAP(sc_wt1, dims = 1:20)  # UMAP without batch correction

## Harmony batch correction
sc_wt1 <- RunHarmony(sc_wt1, "orig.ident")
sc_wt1 <- RunUMAP(sc_wt1, dims = 1:20,reduction = "harmony",reduction.name="HUMAP20")

# Clustering
sc_wt1 <- FindNeighbors(sc_wt1, reduction = "harmony", dims = 1:20,verbose = F)
sc_wt1 <- FindClusters(object = sc_wt1,resolution = 0.3, verbose=F)
sc_wt1$Clusters <- factor(recode(sc_wt1$RNA_snn_res.0.3, 
                                    `0`="Activated Tem/Temra", `1`="Transitional Tem", `2`="Tn/Tcm_1", 
                                    `3`= "Tn/Tcm_2", `4`="Activated Teff", `5`="Activated Tem/Temra", 
                                    `6`="Tpex", `7`= "Transitional Tem"),
                             levels=c("Tn/Tcm_1","Tn/Tcm_2","Transitional Tem",
                                      "Activated Teff", "Activated Tem/Temra","Tpex"))

# Cell origins / conditions
sc_wt1$Origin <- factor(recode(sc_wt1$orig.ident,"scRMFh"="RMF Healthy",
                                  "scRMFaml"="RMF AML","scVLDh"="VLD Healthy","scVLDaml"="VLD AML"),
                           levels=c("RMF AML","RMF Healthy","VLD AML","VLD Healthy"))

# GLIPH2 selected T cells
sc_wt1$GLIPH2.RMF.VLD <- sc_wt1$epitope
sc_wt1$GLIPH2.RMF.VLD[grepl("NA",sc_wt1$GLIPH2_TCR)] <- "unselected"

