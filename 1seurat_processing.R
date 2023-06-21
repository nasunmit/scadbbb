###########################
id='ROSMAP.VascularCells.seurat.harmony.final'
##########################
library(ggrepel)
library(scales)
library(plyr)
library(Seurat)
library(dplyr)
library(harmony)
library(pheatmap)
library(RColorBrewer)
wr <- colorRampPalette(colors = c( "white", "red"))(100)
rwb <- colorRampPalette(colors = c("blue", "white", "red"))(100)
ryb=colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100)
prgn=colorRampPalette(rev(brewer.pal(n = 10, name ="PRGn")))(100)
wb <- colorRampPalette(colors = c( "white", "blue"))(100)
#####################################
data=readRDS('ROSMAP.VascularCells.counts_full,rds')
meta=readRDS('ROSMAP.VascularCells.meta_full.rds')
brain=CreateSeuratObject(counts = data, project = "vascular",meta.data = meta,min.cells = 50, min.features = 500)
brain

brain[["percent.mt"]] <- PercentageFeatureSet(brain, pattern = "^MT-")
brain[["percent.rp"]] <- PercentageFeatureSet(brain, pattern = "^RP")
VlnPlot(brain, features = c("nFeature_RNA", "percent.rp", "percent.mt"), ncol = 3)

brain <- subset(brain, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 5 & percent.rp<5)
brain

meta=brain@meta.data
meta=meta[meta$projid!='50106280',] ## 50106280 has so many fibroblast in hipp
meta=droplevels(meta)

brain=subset(brain,cells=rownames(meta))

brain <- NormalizeData(brain, normalization.method = "LogNormalize", scale.factor = 10000)
brain <- FindVariableFeatures(brain, selection.method = "vst", nfeatures = 2000)

## scaling the data
brain <- ScaleData(brain, features = rownames(brain))
## Perform linear dimensional reduction
brain <- RunPCA(brain, features = VariableFeatures(object = brain))
DimPlot(brain, reduction = "pca")
ElbowPlot(brain,ndims = 50)
k=1:30

## cluster the cells
brain <- FindNeighbors(brain, dims = k)
brain <- FindClusters(brain, resolution = 0.5)
head(Idents(brain), 5)
brain <- RunUMAP(brain, dims = k)

brain <- brain %>% RunHarmony("batch", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(brain, 'harmony')

brain <- brain %>% 
  RunUMAP(reduction = "harmony", dims = k) %>% 
  FindNeighbors(reduction = "harmony", dims = k) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

mycol=col.all[1:length(table(Idents(brain)))]
table(Idents(brain))

pdf(file=paste(id,'.umap.pdf',sep=''),width = 6,height = 5)
DimPlot(brain, reduction = "umap",label=T,cols=mycol)
DimPlot(brain, reduction = "umap",label=F,group.by = 'brain_region',cols=col.dark2[c(1:5,7)])
DimPlot(brain, reduction = "umap",label=F,group.by = 'batch',cols=c(col.set1,col.set3))
DimPlot(brain, reduction = "umap",label=F,group.by = 'projid',cols=hue_pal()(length(table(brain$projid))))+NoLegend()
FeaturePlot(brain,features = 'nFeature_RNA')
FeaturePlot(brain,features = 'percent.mt')
FeaturePlot(brain,features = 'percent.rp')
dev.off()

### recolor 

col1=c('#e69138','#3c78d8','#6aa84f','#c27ba0','#9900ff')
#col1=adjustcolor(col1, alpha.f = 0.5)
names(col1)=c('Endo','Per','Fib','SMC','Ependymal')
DimPlot(brain, reduction = "umap",label=T,group.by='celltype',cols = col1)

##
endo=c('#b45f06','#ffd966','#f6b26b')
names(endo)=c('capEndo','aEndo','vEndo')
per=c('#6d9eeb','#1155cc')
names(per)=c('Per1','Per2')
fib=c('#b6d7a8','#38761d','#00ff00')
names(fib)=c('Fib1','Fib2','Fib3')
smc=c('#a64d79','#d5a6bd')
names(smc)=c('vSMC','aSMC')
epd='#9900ff'
names(epd)='Ependymal'
col2=c(endo,per,fib,smc,epd)
DimPlot(brain, reduction = "umap",label=F,group.by='cellsubtype',cols = col2)

pdf(file=paste(id,'.umap.celltype_annotation.recolored.pdf',sep=''),width = 6,height = 5)
DimPlot(brain, reduction = "umap",label=T,group.by='celltype',cols = col1)
DimPlot(brain, reduction = "umap",label=F,group.by='celltype',cols = col1)
DimPlot(brain, reduction = "umap",label=T,group.by='cellsubtype',cols = col2)
DimPlot(brain, reduction = "umap",label=F,group.by='cellsubtype',cols = col2)
dev.off()



pdf(file=paste(id,'.umap.celltype_annotation.pdf',sep=''),width = 6,height = 5)
DimPlot(brain, reduction = "umap",label=T,group.by='celltype',cols=col.all[1:length(table(brain$celltype))])
DimPlot(brain, reduction = "umap",label=F,group.by='celltype',cols=col.all[1:length(table(brain$celltype))])
DimPlot(brain, reduction = "umap",label=T,group.by='celltype')
DimPlot(brain, reduction = "umap",label=T,group.by='cellsubtype',cols=col.all[1:length(table(brain$cellsubtype))])
DimPlot(brain, reduction = "umap",label=F,group.by='cellsubtype',cols=col.all[1:length(table(brain$cellsubtype))])
DimPlot(brain, reduction = "umap",label=T,group.by='cellsubtype')
dev.off()

saveRDS(brain,file=paste(id,'.rds',sep=''))
write.table(brain@meta.data,file=paste(id,'.metadata.txt',sep = ''),sep = '\t',quote = F)


####
##########
# find markers for every cluster compared to all remaining cells, report only the positive ones
Idents(brain)=brain$celltype
brain.markers <- FindAllMarkers(brain, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
table(brain.markers$cluster)
write.table(brain.markers,file=paste(id,'.celltypeDEGs.txt',sep=''),sep = '\t',quote = F)

clusters=names(table(brain.markers$cluster))
degset=list()
for(c in clusters){
  genes=as.character(brain.markers[brain.markers$cluster==c,7])
  degset[[c]]=genes
}
str(degset)
save(degset,file=paste(id,'.celltypeDEGs.list.RData',sep=''))

# find markers for every cluster compared to all remaining cells, report only the positive ones
Idents(brain)=brain$cellsubtype
brain.markers <- FindAllMarkers(brain, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
table(brain.markers$cluster)
write.table(brain.markers,file=paste(id,'.cellsubtypeDEGs.txt',sep=''),sep = '\t',quote = F)

clusters=names(table(brain.markers$cluster))
degset=list()
for(c in clusters){
  genes=as.character(brain.markers[brain.markers$cluster==c,7])
  degset[[c]]=genes
}
str(degset)
save(degset,file=paste(id,'.cellsubtypeDEGs.list.RData',sep=''))

## pseudo-bulk
newmeta=brain@meta.data
data=brain@assays$RNA@data

## celltype
bulkdata=c()
myclusters=names(table(brain$celltype))
for (c in myclusters){
  cells=rownames(newmeta[newmeta$celltype==c,])
  tmp=data[,cells]
  bulkdata=cbind(bulkdata,rowMeans(as.matrix(tmp)))
}
colnames(bulkdata)=myclusters
write.table(bulkdata,file=paste(id,'.celltype_pseudobulk.txt',sep=''),sep = '\t',quote = F)
## cellsubtype
bulkdata=c()
myclusters=names(table(brain$cellsubtype))
for (c in myclusters){
  cells=rownames(newmeta[newmeta$cellsubtype==c,])
  tmp=data[,cells]
  bulkdata=cbind(bulkdata,rowMeans(as.matrix(tmp)))
}
colnames(bulkdata)=myclusters
write.table(bulkdata,file=paste(id,'.cellsubtype_pseudobulk.txt',sep=''),sep = '\t',quote = F)

##########
regions=names(table(brain$brain_region))
colors=col.dark2[c(1:5,7)]
names(colors)=regions

pdf(file=paste(id,'.umap.eachBrainRegion.pdf',sep=''),width = 4,height = 4)
for (r in regions){
  cells=rownames(newmeta[newmeta$brain_region==r,])
  print(DimPlot(brain,cells=cells, reduction = "umap",label=F,group.by = 'brain_region',cols=colors[r],pt.size = 0.1)+NoLegend()+labs(title=r))
  print(DimPlot(brain,cells=cells, reduction = "umap",label=F,group.by = 'celltype',cols=col.all[1:length(table(brain$celltype))],pt.size = 0.1)+NoLegend()+labs(title=r))
  print(DimPlot(brain,cells=cells, reduction = "umap",label=F,group.by = 'cellsubtype',cols=col.all[1:length(table(brain$cellsubtype))],pt.size = 0.1)+NoLegend()+labs(title=r))
}
dev.off()

pdf(file=paste(id,'.umap.eachBrainRegion_overlay.pdf',sep=''),width = 4,height = 4)
for (r in regions){
  myorder=c(r,regions[regions!=r])
  print(DimPlot(brain,reduction = "umap",label=F,group.by = 'brain_region',order=myorder,cols=c(rep('grey90',5),colors[r]),pt.size = 0.1)+NoLegend()+labs(title=r))
}
dev.off()

##


