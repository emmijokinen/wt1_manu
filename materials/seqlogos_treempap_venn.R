
# Preliminaries
#setwd("Desktop/projects/WT1/")
## Libraries  
library(stringr)
library(ggseqlogo)
library(ggplot2)
library(readr)
library(treemap)
library(VennDiagram)
library(pheatmap)
library(gridExtra)
library(viridis)

## Functions

simpson_clonality <- function(t){
  tab_t = table(t)
  N = sum(unname(tab_t))
  p = unname(tab_t)/N
  R = length(names(tab_t))
  return(sqrt(sum(p**2)))
}

get_CDR3_list <- function(fname){
  tcr_string <- read_file(fname)
  tcrs <- str_split(tcr_string,",")[[1]]
  n <- length(tcrs)
  tcrs[n] <- str_replace(tcrs[n],"\n","") # remove possible line change from last one
  return(tcrs)
}

cdr3_treemap <- function(tcrs,lens,l,cohort,epi,lb=0.3,figpath=""){
  t <- tcrs[lens==l]
  clo = simpson_clonality(t)
  cdr3counts <- sort(table(t),decreasing = T)
  data = data.frame(cdr3counts)
  title = paste0(cohort,": ",epi,", clonality =  ",format(clo,digits=3))
  
  pdf(paste0(figpath,"treemap_",cohort,"_",epi,"_cdr3_",as.character(l),".pdf"),width=9,height = 3)
  treemap(data,index="t", vSize="Freq", type="index", lowerbound.cex.labels = lb,fontface.labels = 1,
          border.lwds=0.2, palette = viridis(6,option="viridis")[2:6],title = title,aspRatio = 3)
  dev.off()
}

## Load data

### Set correct paths
cdr3path <- "/Users/jokiemmi/Desktop/projects/TCRGP/wt1_manu/data/CDR3_lists/"
figpath <- "figures_test/"
filename_tcrgp_rmf_train <- "/Users/jokiemmi/Desktop/projects/TCRGP/training_data/complete_1to10_RMFPNAPYL-ghruh+gtcrex.tsv"
filename_tcrgp_vld_train <- "/Users/jokiemmi/Desktop/projects/TCRGP/training_data/complete_1to10_VLDFAPPGA-ghruh+gtcrex.tsv"

cohorts <- c("AML-BM","MDS-BM","CML-BM","CML-PB","HC-PB","training data")


### Read all tcrs in one named list
tcrRMF_list <- list()
tcrVLD_list <- list()
for (cohort in cohorts[1:5]){
  tcrs_RMF <- get_CDR3_list(paste0(cdr3path,cohort,"_predRMF-CDR3B.txt"))
  tcrs_VLD <- get_CDR3_list(paste0(cdr3path,cohort,"_predVLD-CDR3B.txt"))
  tcrRMF_list <- c(tcrRMF_list, list(tcrs_RMF))
  tcrVLD_list <- c(tcrVLD_list, list(tcrs_VLD))
} 

df <- read.csv(filename_tcrgp_rmf_train,sep="\t")
tcrs_RMF_train <- df$CDR3B[df$Epitope=="RMFPNAPYL"]
tcrRMF_list <- c(tcrRMF_list, list(tcrs_RMF_train))

df <- read.csv(filename_tcrgp_vld_train,sep="\t")
tcrs_VLD_train <- df$CDR3B[df$Epitope=="VLDFAPPGA"]
tcrVLD_list <- c(tcrVLD_list, list(tcrs_VLD_train))

names(tcrRMF_list) <- cohorts
names(tcrVLD_list) <- cohorts




# Seqlogos and treemaps

l=15 # CDR3B length to be used in seqlogos

for (cohort in cohorts){
  for (epi in c("RMF","VLD")){
    tcrs <- ifelse(epi=="RMF", tcrRMF_list[cohort], tcrVLD_list[cohort])[[1]]
    lens = sapply(tcrs,nchar)
    Il <- lens==l
    nc <- sum(Il)
    
    loglist <- list(tcrs[Il])
    names(loglist) <- c(paste0(cohort,": ",epi,", ", "length: ",as.character(l),", count: ",as.character(nc)))
    ggseqlogo(loglist,ncol=1, seq_type='aa',method="prob")
    ggsave(paste0(figpath,"seqlogo_",cohort,"_cdr3_",as.character(l),".pdf"),device=pdf,
           width = 7,height=3)
    
    cdr3_treemap(tcrs,lens,l,cohort,epi,lb=0.2,figpath=figpath)
  }
}

# Venn Diagrams

palvenn <- c("#b53319","#f59551","#217196","#5fb0c2","#f5f0cb")

## RMF

### Venn diagram of shared predicted RMF-specific CDR3Bs between cohorts
venn.diagram(
  x = tcrRMF_list[1:5],
  filename = paste0(figpath,'venn_RMF.png'),
  output=TRUE,fill=palvenn[1:5],
  lwd=1, col= palvenn[1:5],
  width=3000, height=2350,
  fontfamily = "sans", cat.fontfamily = "sans",cex=0.8, cat.cex=0.9,margin=0.1,
  disable.logging = T,force.unique = T
)


### RMF-specific TCRs in all cohorts

message("Predicted RMF-specific CDR3Bs shared between all cohorts and their numbers in the cohorts:")
RMFshared <- intersect(tcrRMF_list[[1]],intersect(tcrRMF_list[[2]],intersect(tcrRMF_list[[3]],intersect(tcrRMF_list[[4]],tcrRMF_list[[5]]))))
message(paste(RMFshared,collapse=", "))
for (i in 1:5){
  message(cohorts[i])
  for (v in RMFshared){
    message(v,": ",sum(str_detect(tcrRMF_list[[i]],v)))
  }
}


## VLD

### Venn diagram of shared predicted VLD-specific CDR3Bs between cohorts
venn.diagram(
  x = tcrVLD_list[1:5],
  filename = paste0(figpath,'venn_VLD.png'),
  output=TRUE,fill=palvenn[1:5],
  lwd=1, col= palvenn[1:5],
  width=3000, height=2350,
  fontfamily = "sans", cat.fontfamily = "sans",cex=0.8, cat.cex=0.9,margin=0.1,
  disable.logging = T,force.unique = T
)


### Shared VLD specific TCRs in all cohorts

message("Predicted VLD-specific CDR3Bs shared between all cohorts and their numbers in the cohorts:")
VLDshared <- intersect(tcrVLD_list[[1]],intersect(tcrVLD_list[[2]],intersect(tcrVLD_list[[3]],intersect(tcrVLD_list[[4]],tcrVLD_list[[5]]))))
message(paste(VLDshared,collapse=", "))
for (i in 1:5){
  message(cohorts[i])
  for (v in VLDshared){
    message(v,": ",sum(str_detect(tcrVLD_list[[i]],v)))
  }
}



