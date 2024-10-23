library(phytools)
library(Momocs)
library(readxl)
#remember to change directory to wherever you have the files
setwd("~/Uncarina")

#read in outlines from the path (update this as well)
jpg.list<-list.files(path="/outlines",pattern = ".jpg",full.names = T)
import_jpg(jpg.list)->returns_leaves
Out(returns_leaves)-> leaves

#read in groups and set
groups <- read_excel("Uncarina_leaf_images_metadata.xls")
groups<-as.data.frame(groups)
groups$species<-as.factor(groups$species)
groups$clades<-as.factor(groups$clade)
leaves$fac=data.frame(groups$species)

images.groups<-gsub("..jpg","",groups$File_Name)
images.names<-names(returns_leaves)
groups.verbatum<-cbind.data.frame(images.names,images.groups,groups$species,groups$clade)

#landmark outlines to ensure consistent orientation
par(mar=c(1,1,1,1))
ldks1 <- vector("list", length=length(leaves$coo))
names(ldks1) <- names(leaves$coo)
coo_centpos(leaves) -> centroids
#4 ldks
par(mar=c(1,1,1,1))
ldks1 <- vector("list", length=length(leaves$coo))
names(ldks1) <- names(leaves$coo)
coo_centpos(leaves) -> centroids
for (i in 1:length(leaves)) {
  ldk <- numeric(4)
  leaves$coo[[i]] -> x
  which.min(x[,1]) -> min.x
  which.max(x[,1]) -> max.x
  which.min(x[,2]) -> min.y
  which.max(x[,2]) -> max.y
  p <- c(centroids[i,1], min(x[,2]))
  l <- apply(x, 1, function(y) sqrt(sum((p - y)^2)))
  ldk[1] <- which.min(l)
  p <- c(centroids[i,1], max(x[,2]))
  l <- apply(x, 1, function(y) sqrt(sum((p - y)^2)))
  ldk[2] <- which.min(l)
  p <- c(min(x[,1]), centroids[i,2])
  l <- apply(x, 1, function(y) sqrt(sum((p - y)^2)))
  ldk[3] <- which.min(l)
  p <- c(max(x[,1]), centroids[i,2])
  l <- apply(x, 1, function(y) sqrt(sum((p - y)^2)))
  ldk[4] <- which.min(l)
  ldk -> ldks1[[i]]
}

#Procrustes alignment
leaves$ldk <- ldks1
coo_slide(leaves, ldk=2) -> leaves
fgProcrustes(leaves) -> leaves2
stack(leaves2, ldk.col="red", ldk.cex=1)

#Elliptic fourier analysis

#run to figure out numer of harmonics
#calibrate_harmonicpower_efourier(x=leaves2)

efou_leaves<-efourier(leaves2, norm=F, nb.h=12, smooth.it=10)
efou_leaves$fac<-tibble::as_tibble(groups[,1:4])

#generate mean shapes and figures
UNC.ms.sp <- MSHAPES(efou_leaves, ~species)

pdf("mean.shapes.panel.pdf")
Out(UNC.ms.sp$shp) %>% panel(cols="black",borders="grey",names=T,cex.names=0.1)
dev.off()

pdf("mean.shapes.names.pdf")
Out(UNC.ms.sp$shp) %>% panel(cols="white",borders="grey",names=T,cex.names=0.7)
dev.off()

UNC.ms.cl <- MSHAPES(efou_leaves, ~clade)
pdf("mean.shapes.clade.panel.pdf")
Out(UNC.ms.cl$shp) %>% panel(cols="black",borders="grey",names=T,cex.names=0.1)
dev.off()

pdf("mean.shapes.clade.names.pdf")
Out(UNC.ms.cl$shp) %>% panel(cols="white",borders="grey",names=T,cex.names=0.7)
dev.off()

#PCA
pca_efou_leaves <- PCA(efou_leaves)
summary(pca_efou_leaves)

#look at variation explained by each component
pdf("PC_contrib.pdf", height=5, width=10)
PCcontrib(pca_efou_leaves, nax=1:4, sd.r = c(-2, -1, 0, 1, 2))
dev.off()

pca_efou_leaves$fac <- tibble::as_tibble(groups[,1:4])
pca_efou_leaves$fac$species<-as.factor(pca_efou_leaves$fac$species)
pca_efou_leaves$fac$clade<-as.factor(pca_efou_leaves$fac$clade)

scores_coos<-groups
scores_pscs<-pca_efou_leaves$x[,1:2]
scores_coos$PC1<-scores_pscs[,1]
scores_coos$PC2<-scores_pscs[,2]
scores_decaryi<-scores_coos[scores_coos$species=="U. decaryi",]

#Plot#
pdf("Uncarina_leaf_morphospace_groups.pdf", height=12, width=10)
plot_PCA(pca_efou_leaves, ~species,  palette = col_qual,morphospace=T,points=T,legend=F,points_transp =0,chull=T,labelpoints = F, center_origin = F,axes=c(1,2),eigen=F,box=T,axesvar =F) %>%
  layer_chullfilled(alpha=.9) %>% layer_points(cex=3,pch=20,transp=0) %>% layer_axesvar(cex=1) %>% layer_morphospace_PCA(pos="xy",size=.15,col="black") %>%   layer_legend(cex=5)
dev.off()

pdf("Uncarina_leaf_morphospace_groups_legend.pdf", height=12, width=10)
plot_PCA(pca_efou_leaves, ~species,  palette = col_qual,morphospace=T,points=T,legend=T,points_transp =0,chull=T,labelpoints = T, center_origin = F,axes=c(1,2),eigen=F,box=T,axesvar =F) %>%
  layer_chullfilled(alpha=.9) %>% layer_points(cex=3,pch=20,transp=0) %>% layer_axesvar(cex=1) %>% layer_morphospace_PCA(pos="xy",size=.15,col="black") %>% 
  layer_legend(cex=5)
dev.off()

pdf("Uncarina_leaf_morphospace_clade_legend.pdf", height=12, width=10)
plot_PCA(pca_efou_leaves, ~clade,  palette = col_qual,morphospace=T,points=T,legend=T,points_transp = 0,chull=T,labelpoints = F, center_origin = F,axes=c(1,2),eigen=F,box=T,axesvar =F) %>%
  layer_chull() %>% layer_points(cex=2.5,pch=20,transp=0) %>% layer_axesvar(cex=1) %>% layer_morphospace_PCA(pos="xy",size=.2,col="black") %>%   layer_legend(cex=10)
dev.off()

pdf("Uncarina_leaf_morphospace_centroids.pdf", height=12, width=10)
plot_PCA(pca_efou_leaves, ~species,  palette = col_black,morphospace=T,points=F,points_transp =0,chull=F,labelpoints = F, legend=F,center_origin = F,axes=c(1,2),eigen=T,box=T,axesvar =F) %>%
  layer_ellipsesaxes(conf = 0.05, lwd = 1, alpha = 0) %>% layer_labelgroups(cex = .8,  font = 4,rect = F,alpha = 0,abbreviate = FALSE) 
dev.off()