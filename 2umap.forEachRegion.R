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
col1=c('#e69138','#3c78d8','#6aa84f','#c27ba0','#9900ff')
names(col1)=c('Endo','Per','Fib','SMC','Ependymal')

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

###
brain=readRDS('ROSMAP.VascularCells.seurat.harmony.final.rds')
meta=brain@meta.data
regions=names(table(brain$brain_region))
colors=col.dark2[c(1:5,7)]
names(colors)=regions

ctps=names(table(brain$celltype))
cols=c(unname(col1),'grey90')
names(cols)=c(names(col1),'other')
cols

pdf(file=paste(id,'.umap.eachBrainRegion_overlay.coloredbyCelltype.pdf',sep=''),width = 4,height = 4)
for (r in regions){
  tmp1=meta[meta$brain_region==r,]
  tmp1$mycol=tmp1$celltype
  tmp2=meta[meta$brain_region!=r,]
  tmp2$mycol=rep('other',nrow(tmp2))
  tmp=rbind(tmp1,tmp2)
  tmp=tmp[rownames(meta),]
  brain.tmp=brain
  brain.tmp@meta.data=tmp
  myorder=names(cols)
  print(DimPlot(brain.tmp,reduction = 'umap',label = F,group.by = 'mycol',order=myorder,cols=cols,pt.size = 0.1)+NoLegend()+labs(title=r))
}
dev.off()


### subtype

ctps=names(table(brain$cellsubtype))
cols=c(unname(col2),'grey90')
names(cols)=c(names(col2),'other')
cols

pdf(file=paste(id,'.umap.eachBrainRegion_overlay.coloredbyCellsubtype.pdf',sep=''),width = 4,height = 4)
for (r in regions){
  tmp1=meta[meta$brain_region==r,]
  tmp1$mycol=tmp1$cellsubtype
  tmp2=meta[meta$brain_region!=r,]
  tmp2$mycol=rep('other',nrow(tmp2))
  tmp=rbind(tmp1,tmp2)
  tmp=tmp[rownames(meta),]
  brain.tmp=brain
  brain.tmp@meta.data=tmp
  myorder=names(cols)
  print(DimPlot(brain.tmp,reduction = 'umap',label = F,group.by = 'mycol',order=myorder,cols=cols,pt.size = 0.1)+NoLegend()+labs(title=r))
}
dev.off()


