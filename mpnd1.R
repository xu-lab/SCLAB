# This is the script to reproduce the PND1 mouse lung cell types.
# Please note that Seurat 1.4 objects and functions were used in this scripts.
# Author: Minzhe Guo (minzhe.guo@cchmc.org)

library(dplyr)
library(Matrix)
library(reshape2)
library(Seurat)
library(outliers)
library(gmodels)
library(RANN) 
library(Biobase)
library(gridExtra)
library(pheatmap)

source("https://raw.githubusercontent.com/xu-lab/SINCERA/umi/R/functions.R")

project.name <- "Mouse.Lung.PND1"

# We will update the path to access the raw data file
mpnd1 <- readRDS(file="mpnd1.rds")

min.cells <- 2

##########################################################################
######################## major cell type #################################
##########################################################################

# select highly variable genes
batch1.cells <- rownames(subset(mpnd1@data.info, batch=="PND1.Batch1"))
batch2.cells <- rownames(subset(mpnd1@data.info, batch=="PND1.Batch2"))
mpnd1@var.genes <- pc.genes.selection(mpnd1@data, rownames(mpnd1@data), batch1.cells, batch2.cells)$pc.genes
cat("Reducing dimensions using ", length(mpnd1@var.genes), " highly variable genes\n", sep="")

# pca using highly variable genes
mpnd1 <- DoPCA(mpnd1, genes=mpnd1@var.genes, fast=F)
PCElbowPlot(mpnd1, num.pc = 50)

r <- 5 # ret$r
pc.use <- c(1:r)

# tSNE plot
mpnd1 <- RunTSNE(mpnd1, dims.use = pc.use, do.fast = T, dim_embed = 2, k.seed = 1)
plot1 <- TSNEPlot(mpnd1, pt.size=0.1)

# clustering
mpnd1 <- graph.clustering(mpnd1, pcs.use=pc.use, do.jaccard=T, method = "Louvain", num.nn=30)
mpnd1 <- BuildClusterTree(mpnd1, do.reorder = T, reorder.numeric = T, pcs.use = pc.use, do.plot=T)

# assign major cell types
mpnd1@data.info$major.cluster <- mpnd1@data.info$tree.ident
mpnd1@data.info$major.type <- NA
mpnd1@data.info$major.type[which(mpnd1@data.info$tree.ident %in% c(2:5))] <- "Endothelial"
mpnd1@data.info$major.type[which(mpnd1@data.info$tree.ident %in% c(16:20))] <- "Epithelial"
mpnd1@data.info$major.type[which(mpnd1@data.info$tree.ident %in% c(9:10,12:15))] <- "Mesenchyme"
mpnd1@data.info$major.type[which(mpnd1@data.info$tree.ident %in% c(6:8))] <- "Immune"
mpnd1@data.info$major.type[which(mpnd1@data.info$tree.ident %in% c(1))] <- "Contam.Endo+Epi+"
mpnd1@data.info$major.type[which(mpnd1@data.info$tree.ident %in% c(11))] <- "Contam.Endo+FB+"

mpnd1 <- SetIdent(mpnd1, ident.use=mpnd1@data.info$major.type)
TSNEPlot(mpnd1, do.label = F, pt.size=0.1)


# compare with paper results
cmp <- table(mpnd1@data.info$tree.ident, cellmeta$major.cluster)
cmp <- cmp/rowSums(cmp)
pheatmap(cmp, cluster_cols = FALSE, cluster_rows = FALSE)
cmp <- table(mpnd1@data.info$tree.ident, cellmeta$major.type)
cmp <- cmp/rowSums(cmp)
pheatmap(cmp, cluster_cols = FALSE, cluster_rows = FALSE)
cmp <- table(mpnd1@data.info$major.type, cellmeta$major.type)
cmp <- cmp/rowSums(cmp)
pheatmap(cmp, cluster_cols = FALSE, cluster_rows = FALSE)



##########################################################################
################### mesenchymal cell subtype #############################
##########################################################################

mesen.cells <- colnames(mpnd1@data)[mpnd1@ident == "Mesenchyme"]
batch1.mesen.cells <- mesen.cells[which(mesen.cells %in% batch1.cells)]
batch2.mesen.cells <- mesen.cells[which(mesen.cells %in% batch2.cells)]
mesen.cells <- c(batch1.mesen.cells, batch2.mesen.cells)
print(length(mesen.cells))

batch1.mesen.genes <- prefiltering.genes(mpnd1@raw.data[, batch1.mesen.cells], 
                                         min.exp=0, min.cells=min.cells, 
                                         min.counts=0)
batch2.mesen.genes <- prefiltering.genes(mpnd1@raw.data[, batch2.mesen.cells], 
                                         min.exp=0, min.cells=min.cells, 
                                         min.counts=0)
mesen.genes <- intersect(batch1.mesen.genes, batch2.mesen.genes)
print(length(mesen.genes))
mesen.count <- mpnd1@raw.data[mesen.genes, mesen.cells]
mesen.lognorm <- mpnd1@data[mesen.genes, mesen.cells]
  

mpnd1.mesen <- new("seurat", raw.data = mesen.lognorm)
mpnd1.mesen <- Setup(mpnd1.mesen, min.cells = 0, min.genes = 0, 
                     is.expr=0, do.logNormalize = F, total.expr = 1e4, 
                     do.scale=TRUE, do.center=TRUE, 
                     project = paste(project.name, ".mesen", sep=""))
mpnd1.mesen@raw.data <- mesen.count

mpnd1.mesen@data.info$orig.ident <- mpnd1@data.info[rownames(mpnd1.mesen@data.info), "orig.ident"]
mpnd1.mesen@data.info$major.cluster <- mpnd1@data.info[rownames(mpnd1.mesen@data.info), "major.cluster"]
mpnd1.mesen@ident <- mpnd1.mesen@data.info$orig.ident
names(mpnd1.mesen@ident) <- rownames(mpnd1.mesen@data.info)

mesen.data.combat <- pc.batch.combat(t(mpnd1.mesen@data), 
                                     sample=mpnd1.mesen@data.info$orig.ident, 
                                     db=paste(project.name, ".mesen", sep=""))  
mpnd1.mesen@data <- t(data.frame(mesen.data.combat, check.names = F))
mpnd1.mesen <- ScaleData(mpnd1.mesen,do.scale = T, do.center = T)

mpnd1.mesen@var.genes <- pc.genes.selection(mpnd1.mesen@data, 
                                            mesen.genes, batch1.mesen.cells, 
                                            batch2.mesen.cells, 
                                            x.low.cutoff = 0.2, 
                                            x.high.cutoff=4, 
                                            y.cutoff=0)$pc.genes
cat("Reducing dimensions using ", length(mpnd1.mesen@var.genes), " highly variable genes\n", sep="")

mpnd1.mesen <- DoPCA(mpnd1.mesen, genes=mpnd1.mesen@var.genes, fast=F)
PCElbowPlot(mpnd1.mesen)
    
    
pc.use <- 1:5
mpnd1.mesen <- RunTSNE(mpnd1.mesen, dims.use = pc.use, 
                       do.fast = T, dim_embed = 2)
plot1 <- TSNEPlot(mpnd1.mesen, pt.size=0.1)


mpnd1.mesen <- graph.clustering(mpnd1.mesen, pcs.use=pc.use, 
                                do.jaccard=T, method = "Louvain", num.nn=20)
mpnd1.mesen <- BuildClusterTree(mpnd1.mesen, do.reorder = T, 
                                reorder.numeric = T, pcs.use = pc.use, 
                                do.plot=T)
TSNEPlot(mpnd1.mesen, pt.size=0.1, do.label = T, label.size = 6)


# compare with paper results
cmp <- table(mpnd1.mesen@data.info$tree.ident, cellmeta[rownames(mpnd1.mesen@data.info), "SubType"])
cmp <- cmp/rowSums(cmp)
pheatmap(cmp, cluster_cols = FALSE, cluster_rows = FALSE)

# define mesenchymal subpopulations
mpnd1.mesen@data.info$mesen.cluster <- mpnd1.mesen@data.info$tree.ident

mpnd1.mesen@data.info$mesen.type.1 <- NA
mpnd1.mesen@data.info$mesen.type.1[which(mpnd1.mesen@data.info$mesen.cluster %in% 13)] <- "Mesen.D"
mpnd1.mesen@data.info$mesen.type.1[which(mpnd1.mesen@data.info$mesen.cluster %in% 9)] <- "Mesen.B"
mpnd1.mesen@data.info$mesen.type.1[which(mpnd1.mesen@data.info$mesen.cluster %in% 10)] <- "Mesen.C"
mpnd1.mesen@data.info$mesen.type.1[which(mpnd1.mesen@data.info$mesen.cluster %in% 1:8)] <- "Mesen.A"
mpnd1.mesen@data.info$mesen.type.1[which(mpnd1.mesen@data.info$mesen.cluster %in% c(11))] <- "Mesen.E"
mpnd1.mesen@data.info$mesen.type.1[which(mpnd1.mesen@data.info$mesen.cluster %in% c(12))] <- "Mesen.F"
mpnd1.mesen@data.info$mesen.type.1[which(mpnd1.mesen@data.info$mesen.cluster %in% c(14))] <- "Mesen.G"

mpnd1.mesen@data.info$mesen.type.2 <- NA
mpnd1.mesen@data.info$mesen.type.2[which(mpnd1.mesen@data.info$mesen.cluster %in% 13)] <- "Pericyte.2"
mpnd1.mesen@data.info$mesen.type.2[which(mpnd1.mesen@data.info$mesen.cluster %in% 9)] <- "MatrixFB.2"
mpnd1.mesen@data.info$mesen.type.2[which(mpnd1.mesen@data.info$mesen.cluster %in% 10)] <- "Pericyte.1"
mpnd1.mesen@data.info$mesen.type.2[which(mpnd1.mesen@data.info$mesen.cluster %in% 1:8)] <- "MatrixFB.1"
mpnd1.mesen@data.info$mesen.type.2[which(mpnd1.mesen@data.info$mesen.cluster %in% c(11))] <- "MyoFB.1"
mpnd1.mesen@data.info$mesen.type.2[which(mpnd1.mesen@data.info$mesen.cluster %in% c(12))] <- "MyoFB.2"
mpnd1.mesen@data.info$mesen.type.2[which(mpnd1.mesen@data.info$mesen.cluster %in% c(14))] <- "Smooth Muscle"

# compare with paper results
cmp <- table(mpnd1.mesen@data.info$mesen.type.2, cellmeta[rownames(mpnd1.mesen@data.info), "SubType"])
cmp <- cmp/rowSums(cmp)
pheatmap(cmp, cluster_cols = FALSE, cluster_rows = FALSE)


mpnd1.mesen <- SetIdent(mpnd1.mesen, ident.use=mpnd1.mesen@data.info$mesen.type.1)
TSNEPlot(mpnd1.mesen, do.label = F, pt.size=0.1)

saveRDS(mpnd1.mesen, file="mpnd1.mesen.rds")


##########################################################################
#################### epithelial cell subtype #############################
##########################################################################
  
epi.cells <- colnames(mpnd1@data)[mpnd1@ident == "Epithelial"]
batch1.epi.cells <- epi.cells[which(epi.cells %in% batch1.cells)]
batch2.epi.cells <- epi.cells[which(epi.cells %in% batch2.cells)]
epi.cells <- c(batch1.epi.cells, batch2.epi.cells)
print(length(epi.cells))

batch1.epi.genes <- prefiltering.genes(mpnd1@raw.data[, batch1.epi.cells], 
                                       min.exp=0, min.cells=min.cells, 
                                       min.counts=0)
batch2.epi.genes <- prefiltering.genes(mpnd1@raw.data[, batch2.epi.cells], 
                                       min.exp=0, min.cells=min.cells, 
                                       min.counts=0)
epi.genes <- intersect(batch1.epi.genes, batch2.epi.genes)
print(length(epi.genes))
epi.count <- mpnd1@raw.data[epi.genes, epi.cells]
epi.lognorm <- mpnd1@data[epi.genes, epi.cells]

db <- "mpnd1.epi"

mpnd1.epi <- new("seurat", raw.data = epi.lognorm)
mpnd1.epi <- Setup(mpnd1.epi, min.cells = 0, min.genes = 0, 
                   is.expr=0, do.logNormalize = F, total.expr = 1e4, 
                   do.scale=TRUE, do.center=TRUE, project = db)
mpnd1.epi@raw.data <- epi.count
mpnd1.epi@data.info$orig.ident <- mpnd1@data.info[rownames(mpnd1.epi@data.info), "orig.ident"]
mpnd1.epi@data.info$major.cluster <- mpnd1@data.info[rownames(mpnd1.epi@data.info), "major.cluster"]
mpnd1.epi@ident <- mpnd1.epi@data.info$orig.ident
names(mpnd1.epi@ident) <- rownames(mpnd1.epi@data.info)
  
epi.data.combat <- pc.batch.combat(t(mpnd1.epi@data), 
                                   sample=mpnd1.epi@data.info$orig.ident, 
                                   db="mpnd1.epi")  
mpnd1.epi@data <- t(data.frame(epi.data.combat, check.names = F))
mpnd1.epi <- ScaleData(mpnd1.epi,do.scale = T, do.center = T)
mpnd1.epi@var.genes <- pc.genes.selection(mpnd1.epi@data, 
                                          epi.genes, batch1.epi.cells, 
                                          batch2.epi.cells, 
                                          x.low.cutoff = 0.2, 
                                          x.high.cutoff=4, 
                                          y.cutoff=0.25)$pc.genes
cat("Reducing dimensions using ", length(mpnd1.epi@var.genes), " highly variable genes\n", sep="")
mpnd1.epi <- DoPCA(mpnd1.epi, genes=mpnd1.epi@var.genes, fast=F)
PCAPlot(mpnd1.epi)
PCElbowPlot(mpnd1.epi)

pc.use <- 1:5
mpnd1.epi <- RunTSNE(mpnd1.epi, dims.use = pc.use, do.fast = T, dim_embed = 2)
plot1 <- TSNEPlot(mpnd1.epi,pt.size=0.1)


mpnd1.epi <- graph.clustering(mpnd1.epi, pcs.use=pc.use, 
                              do.jaccard=T, method = "Louvain", num.nn=10)
mpnd1.epi <- BuildClusterTree(mpnd1.epi, do.reorder = T, reorder.numeric = T, pcs.use = pc.use, do.plot=T)
TSNEPlot(mpnd1.epi, do.label = T, label.size = 6, pt.size=0.1)

# compare with paper results
cmp <- table(mpnd1.epi@data.info$tree.ident, 
             cellmeta[rownames(mpnd1.epi@data.info), "SubType"])
cmp <- cmp/rowSums(cmp)
pheatmap(cmp, cluster_cols = FALSE, cluster_rows = FALSE)

# define epithelial subpopulations
mpnd1.epi@data.info$epi.cluster <- mpnd1.epi@data.info$tree.ident
mpnd1.epi@data.info$epi.type.1 <- NA
mpnd1.epi@data.info$epi.type.1[which(mpnd1.epi@data.info$epi.cluster %in% c(1,4:6))] <- "Epi.A"
mpnd1.epi@data.info$epi.type.1[which(mpnd1.epi@data.info$epi.cluster %in% c(2,3,7:9))] <- "Epi.B"
mpnd1.epi@data.info$epi.type.1[which(mpnd1.epi@data.info$epi.cluster %in% c(10:17))] <- "Epi.C"
mpnd1.epi@data.info$epi.type.1[which(mpnd1.epi@data.info$epi.cluster %in% 18)] <- "Epi.D" 
mpnd1.epi@data.info$epi.type.1[which(mpnd1.epi@data.info$epi.cluster %in% 19)] <- "Epi.E" 
mpnd1.epi@data.info$epi.type.1[which(mpnd1.epi@data.info$epi.cluster %in% 20)] <- "Epi.F"

mpnd1.epi@data.info$epi.type.2 <- NA
mpnd1.epi@data.info$epi.type.2[which(mpnd1.epi@data.info$epi.cluster %in% c(1,4:6))] <- "AT1)"
mpnd1.epi@data.info$epi.type.2[which(mpnd1.epi@data.info$epi.cluster %in% c(2,3,7:9))] <- "AT1+AT2+"
mpnd1.epi@data.info$epi.type.2[which(mpnd1.epi@data.info$epi.cluster %in% c(10:17))] <- "AT2"
mpnd1.epi@data.info$epi.type.2[which(mpnd1.epi@data.info$epi.cluster %in% 18)] <- "Club" 
mpnd1.epi@data.info$epi.type.2[which(mpnd1.epi@data.info$epi.cluster %in% 19)] <- "Sox2_hi" 
mpnd1.epi@data.info$epi.type.2[which(mpnd1.epi@data.info$epi.cluster %in% 20)] <- "Ciliated"

# compare with paper results
cmp <- table(mpnd1.epi@data.info$epi.type.2, 
             cellmeta[rownames(mpnd1.epi@data.info), "SubType"])
cmp <- cmp/rowSums(cmp)
pheatmap(cmp, cluster_cols = FALSE, cluster_rows = FALSE)

mpnd1.epi <- SetIdent(mpnd1.epi, ident.use=mpnd1.epi@data.info$epi.type.1)
TSNEPlot(mpnd1.epi, do.label = F, pt.size=0.1)

saveRDS(mpnd1.epi, file="mpnd1.epi.rds")
  

##########################################################################
#################### epithelial cell subtype #############################
##########################################################################

endo.cells <- colnames(mpnd1@data)[mpnd1@ident == "Endothelial"]
batch1.endo.cells <- endo.cells[which(endo.cells %in% batch1.cells)]
batch2.endo.cells <- endo.cells[which(endo.cells %in% batch2.cells)]
endo.cells <- c(batch1.endo.cells, batch2.endo.cells)

batch1.endo.genes <- prefiltering.genes(mpnd1@raw.data[, batch1.endo.cells], 
                                       min.exp=0, min.cells=min.cells, 
                                       min.counts=0)
batch2.endo.genes <- prefiltering.genes(mpnd1@raw.data[, batch2.endo.cells], 
                                        min.exp=0, min.cells=min.cells, 
                                        min.counts=0)
endo.genes <- intersect(batch1.endo.genes, batch2.endo.genes)
print(length(endo.genes))
endo.count <- mpnd1@raw.data[endo.genes, endo.cells]
endo.lognorm <- mpnd1@data[endo.genes, endo.cells]

db <- "mpnd1.endo"

mpnd1.endo <- new("seurat", raw.data = endo.lognorm)
mpnd1.endo <- Setup(mpnd1.endo, min.cells = 0, min.genes = 0, is.expr=0, do.logNormalize = F, total.expr = 1e4, do.scale=TRUE, do.center=TRUE, project = db)
mpnd1.endo@raw.data <- endo.count
mpnd1.endo@data.info$orig.ident <- mpnd1@data.info[rownames(mpnd1.endo@data.info), "orig.ident"]
mpnd1.endo@data.info$major.cluster <- mpnd1@data.info[rownames(mpnd1.endo@data.info), "major.cluster"]
mpnd1.endo@ident <- mpnd1.endo@data.info$orig.ident
names(mpnd1.endo@ident) <- rownames(mpnd1.endo@data.info)

endo.data.combat <- pc.batch.combat(t(mpnd1.endo@data), 
                                    sample=mpnd1.endo@data.info$orig.ident, 
                                    db="mpdn1.endo")  
mpnd1.endo@data <- t(data.frame(endo.data.combat, check.names = F))
mpnd1.endo <- ScaleData(mpnd1.endo, do.scale = T, do.center = T)
mpnd1.endo@var.genes <- pc.genes.selection(mpnd1.endo@data, endo.genes, 
                                           batch1.endo.cells, 
                                           batch2.endo.cells, 
                                           x.low.cutoff = 0.2, 
                                           x.high.cutoff=4, 
                                           y.cutoff=0.25)$pc.genes
cat("Reducing dimensions using ", length(mpnd1.endo@var.genes), " highly variable genes\n", sep="")
mpnd1.endo <- DoPCA(mpnd1.endo, genes=mpnd1.endo@var.genes, fast=F)
PCElbowPlot(mpnd1.endo)

pc.use <- 1:5
mpnd1.endo <- RunTSNE(mpnd1.endo, dims.use = pc.use, do.fast = T, dim_embed = 2)
plot1 <- TSNEPlot(mpnd1.endo, pt.size=0.1)

mpnd1.endo <- graph.clustering(mpnd1.endo, pcs.use=pc.use, 
                             do.jaccard=T, method = "Louvain", num.nn=20)
mpnd1.endo <- BuildClusterTree(mpnd1.endo, do.reorder = T, 
                             reorder.numeric = T, pcs.use = pc.use, do.plot=T)
TSNEPlot(mpnd1.endo, do.label = T, label.size = 5)

cluster.id <- as.numeric(mpnd1.endo@ident)
cluster.id[which(cluster.id==1)] <- 0
cluster.id[which(cluster.id==2)] <- 1
cluster.id[which(cluster.id==0)] <- 2
names(cluster.id) <- rownames(mpnd1.endo@data.info)
mpnd1.endo@data.info$tree.ident <- factor(cluster.id)
mpnd1.endo <- SetIdent(mpnd1.endo, ident.use=factor(cluster.id))
TSNEPlot(mpnd1.endo, do.label = T, label.size = 6, pt.size=0.1)


# compare with paper results
cmp <- table(mpnd1.endo@data.info$tree.ident, 
             cellmeta[rownames(mpnd1.endo@data.info), "SubType"])
cmp <- cmp/rowSums(cmp)
pheatmap(cmp, cluster_cols = FALSE, cluster_rows = FALSE)

# define endothelial subpopulations
mpnd1.endo@data.info$endo.cluster <- mpnd1.endo@data.info$tree.ident
mpnd1.endo@data.info$endo.type.1 <- NA
mpnd1.endo@data.info$endo.type.1[which(mpnd1.endo@data.info$endo.cluster %in% c(2:17))] <- "Endo.A"
mpnd1.endo@data.info$endo.type.1[which(mpnd1.endo@data.info$endo.cluster %in% c(1))] <- "Endo.B"

mpnd1.endo@data.info$endo.type.2 <- NA
mpnd1.endo@data.info$endo.type.2[which(mpnd1.endo@data.info$endo.cluster %in% c(2:17))] <- "VasEndo"
mpnd1.endo@data.info$endo.type.2[which(mpnd1.endo@data.info$endo.cluster %in% c(1))] <- "LymEndo"

# compare with paper results
cmp <- table(mpnd1.endo@data.info$endo.type.1, 
             cellmeta[rownames(mpnd1.endo@data.info), "SubType"])
cmp <- cmp/rowSums(cmp)
pheatmap(cmp, cluster_cols = FALSE, cluster_rows = FALSE)

mpnd1.endo <- SetIdent(mpnd1.endo, ident.use=mpnd1.endo@data.info$endo.type.1)
TSNEPlot(mpnd1.endo, do.label = F, pt.size=0.1)

saveRDS(mpnd1.endo, file="mpnd1.endo.rds")

  
##########################################################################
######################## immune cell subtype #############################
##########################################################################

  
immune.cells <- colnames(mpnd1@data)[mpnd1@ident == "Immune"]
batch1.immune.cells <- immune.cells[which(immune.cells %in% batch1.cells)]
batch2.immune.cells <- immune.cells[which(immune.cells %in% batch2.cells)]
immune.cells <- c(batch1.immune.cells,  batch2.immune.cells)
print(length(immune.cells))

batch1.immune.genes <- prefiltering.genes(mpnd1@raw.data[, batch1.immune.cells], 
                                          min.exp=0, min.cells=min.cells,
                                          min.counts=0)
batch2.immune.genes <- prefiltering.genes(mpnd1@raw.data[, batch2.immune.cells],
                                           min.exp=0, min.cells=min.cells, 
                                           min.counts=0)
immune.genes <- intersect(batch1.immune.genes, batch2.immune.genes)
print(length(immune.genes))
immune.count <- mpnd1@raw.data[immune.genes, immune.cells]
immune.lognorm <- mpnd1@data[immune.genes, immune.cells]

db <- "mpnd1.immune"

mpnd1.immune <- new("seurat", raw.data = immune.lognorm)
mpnd1.immune <- Setup(mpnd1.immune, min.cells = 0, min.genes = 0, is.expr=0, do.logNormalize = F, total.expr = 1e4, do.scale=TRUE, do.center=TRUE, project = db)
mpnd1.immune@raw.data <- immune.count
mpnd1.immune@data.info$orig.ident <- mpnd1@data.info[rownames(mpnd1.immune@data.info), "orig.ident"]
mpnd1.immune@data.info$major.cluster <- mpnd1@data.info[rownames(mpnd1.immune@data.info), "major.cluster"]
mpnd1.immune@ident <- mpnd1.immune@data.info$orig.ident
names(mpnd1.immune@ident) <- rownames(mpnd1.immune@data.info)
  
immune.data.combat <- pc.batch.combat(t(mpnd1.immune@data), 
                                      sample=mpnd1.immune@data.info$orig.ident,
                                      db="mpnd1.immune")  
mpnd1.immune@data <- t(data.frame(immune.data.combat, check.names = F))
mpnd1.immune <- ScaleData(mpnd1.immune, do.scale = T, do.center = T)
mpnd1.immune@var.genes <- pc.genes.selection(mpnd1.immune@data, immune.genes, 
                                             batch1.immune.cells, 
                                             batch2.immune.cells, 
                                             x.low.cutoff = 0.15, 
                                             x.high.cutoff = 4, 
                                             y.cutoff = 0.3)$pc.genes
cat("Reducing dimensions using ", length(mpnd1.immune@var.genes), " highly variable genes\n", sep="")
mpnd1.immune <- DoPCA(mpnd1.immune, genes=mpnd1.immune@var.genes, fast=F)
PCElbowPlot(mpnd1.immune, num.pc = 50)

pc.use <- 1:14
mpnd1.immune <- RunTSNE(mpnd1.immune, dims.use = pc.use, do.fast = T, dim_embed = 2)
plot1 <- TSNEPlot(mpnd1.immune, pt.size=0.1)


mpnd1.immune <- graph.clustering(mpnd1.immune, pcs.use=pc.use, 
                                 do.jaccard=T, method = "Louvain", num.nn=5)
mpnd1.immune <- BuildClusterTree(mpnd1.immune, do.reorder = T, 
                                 reorder.numeric = T, pcs.use = pc.use, 
                                 do.plot=T)
TSNEPlot(mpnd1.immune, do.label = T, label.size = 4, pt.size=0.1)


# compare with paper results
cmp <- table(mpnd1.immune@data.info$tree.ident, 
             cellmeta[rownames(mpnd1.immune@data.info), "SubType"])
cmp <- cmp/rowSums(cmp)
pheatmap(cmp, cluster_cols = FALSE, cluster_rows = FALSE)


# define immune subpopulations
mpnd1.immune@data.info$immune.cluster <- mpnd1.immune@data.info$tree.ident
mpnd1.immune@data.info$immune.type.1 <- NA
mpnd1.immune@data.info$immune.type.1[which(mpnd1.immune@data.info$immune.cluster %in% c(1))] <- "Contam.Platelet+"
mpnd1.immune@data.info$immune.type.1[which(mpnd1.immune@data.info$immune.cluster %in% c(8))] <- "Immune.D"
mpnd1.immune@data.info$immune.type.1[which(mpnd1.immune@data.info$immune.cluster %in% c(19:20))] <- "Contam.Immune+FB+"
mpnd1.immune@data.info$immune.type.1[which(mpnd1.immune@data.info$immune.cluster %in% c(2))] <- "Contam.RBC+"
mpnd1.immune@data.info$immune.type.1[which(mpnd1.immune@data.info$immune.cluster %in% c(3,4))] <- "Immune.C"
mpnd1.immune@data.info$immune.type.1[which(mpnd1.immune@data.info$immune.cluster %in% c(5))] <- "Immune.A"
mpnd1.immune@data.info$immune.type.1[which(mpnd1.immune@data.info$immune.cluster %in% c(6,7))] <- "Immune.B"
mpnd1.immune@data.info$immune.type.1[which(mpnd1.immune@data.info$immune.cluster %in% c(9:18,21:26))] <- "Immune.E"

mpnd1.immune@data.info$immune.type.2 <- NA
mpnd1.immune@data.info$immune.type.2[which(mpnd1.immune@data.info$immune.cluster %in% c(1))] <- "Platelet+"
mpnd1.immune@data.info$immune.type.2[which(mpnd1.immune@data.info$immune.cluster %in% c(8))] <- "Basophils"
mpnd1.immune@data.info$immune.type.2[which(mpnd1.immune@data.info$immune.cluster %in% c(19:20))] <- "Immune+FB+"
mpnd1.immune@data.info$immune.type.2[which(mpnd1.immune@data.info$immune.cluster %in% c(2))] <- "RBC+"
mpnd1.immune@data.info$immune.type.2[which(mpnd1.immune@data.info$immune.cluster %in% c(3,4))] <- "T Lymphocyte"
mpnd1.immune@data.info$immune.type.2[which(mpnd1.immune@data.info$immune.cluster %in% c(5))] <- "Proliferative Lymphocyte"
mpnd1.immune@data.info$immune.type.2[which(mpnd1.immune@data.info$immune.cluster %in% c(6,7))] <- "B Lymphocyte"
mpnd1.immune@data.info$immune.type.2[which(mpnd1.immune@data.info$immune.cluster %in% c(9:18,21:26))] <- "Macrophage"

# compare with paper results
cmp <- table(mpnd1.immune@data.info$immune.type.2, 
             cellmeta[rownames(mpnd1.immune@data.info), "SubType"])
cmp <- cmp/rowSums(cmp)
pheatmap(cmp, cluster_cols = FALSE, cluster_rows = FALSE)



mpnd1.immune <- SetIdent(mpnd1.immune, ident.use=mpnd1.immune@data.info$immune.type.1)
TSNEPlot(mpnd1.immune, do.label = F, pt.size=0.1)

saveRDS(mpnd1.immune, file="mpnd1.immune.rds")
  
  
##########################################################################
######################## collecting subtype info #############################
##########################################################################

# adding subclusters
mpnd1@data.info$cluster <- mpnd1@data.info$major.type
mpnd1@data.info[rownames(mpnd1.mesen@data.info), "cluster"] <- mpnd1.mesen@data.info$mesen.cluster
mpnd1@data.info[rownames(mpnd1.epi@data.info), "cluster"] <- mpnd1.epi@data.info$epi.cluster
mpnd1@data.info[rownames(mpnd1.endo@data.info), "cluster"] <- mpnd1.endo@data.info$endo.cluster
mpnd1@data.info[rownames(mpnd1.immune@data.info), "cluster"] <- mpnd1.immune@data.info$immune.cluster

# update subtype 
mpnd1@data.info$celltype.2<- mpnd1@data.info$major.type
mpnd1@data.info[rownames(mpnd1.mesen@data.info), "celltype.2"] <- as.character(mpnd1.mesen@data.info$mesen.type.1)
mpnd1@data.info[rownames(mpnd1.epi@data.info), "celltype.2"] <- as.character(mpnd1.epi@data.info$epi.type.1)
mpnd1@data.info[rownames(mpnd1.endo@data.info), "celltype.2"] <- as.character(mpnd1.endo@data.info$endo.type.1)
mpnd1@data.info[rownames(mpnd1.immune@data.info), "celltype.2"] <- as.character(mpnd1.immune@data.info$immune.type.1)

mpnd1@data.info$celltype.3<- mpnd1@data.info$major.type
mpnd1@data.info[rownames(mpnd1.mesen@data.info), "celltype.3"] <- as.character(mpnd1.mesen@data.info$mesen.type.2)
mpnd1@data.info[rownames(mpnd1.epi@data.info), "celltype.3"] <- as.character(mpnd1.epi@data.info$epi.type.2)
mpnd1@data.info[rownames(mpnd1.endo@data.info), "celltype.3"] <- as.character(mpnd1.endo@data.info$endo.type.2)
mpnd1@data.info[rownames(mpnd1.immune@data.info), "celltype.3"] <- as.character(mpnd1.immune@data.info$immune.type.2)

# updating major types
MajorType <- as.character(mpnd1@data.info$major.type)
MajorType[which(substr(MajorType, 1, 6) == "Contam")] <- "Contaminants"
MajorType[which(substr(as.character(mpnd1@data.info$celltype.2), 1, 6) == "Contam")] <- "Contaminants"
mpnd1@data.info$celltype.1 <- MajorType


mpnd1@data.info <- mpnd1@data.info[, c("nGene","nUMI","orig.ident","batch","cluster", "celltype.1","celltype.2","celltype.3")]

cmp <- table(mpnd1@data.info$celltype.3, pnd1$All@data.info$SubType)
cmp <- cmp/rowSums(cmp)
pheatmap(cmp, cluster_cols = FALSE, cluster_rows = FALSE)


mpnd1 <- SetIdent(mpnd1, ident.use=mpnd1@data.info$celltype.2)
TSNEPlot(mpnd1)


ret <- list(All=mpnd1, 
              Mesenchyme=mpnd1.mesen, 
              Immune=mpnd1.immune, 
              Epithelial=mpnd1.epi, 
              Endothelial=mpnd1.endo)

save(ret, file="mpnd1.final.rds")

  
  
  
  



