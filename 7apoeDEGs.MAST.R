##########################
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
###########################

## Identify DEGs
#######
library(edgeR)
library(limma)
library(MAST)
library(tidyr)
######
files='ROSMAP.VascularCells.seurat.harmony.final.rds'
####
for (f in files){
  fid=strsplit(rev(strsplit(f,'/')[[1]])[1],'[.]rds')[[1]][1]
  brain=readRDS(f)
  meta=brain@meta.data
  
  ## in the context of AD
  #meta=meta[meta$ADdiag2types=='AD',]
  tmp33=meta[meta$apoe_genotype=='33',]
  tmp33$apoe_genotype.binary=rep('e3',nrow(tmp33))
  tmp4=meta[meta$apoe_genotype %in% c('34','44'),]
  tmp4$apoe_genotype.binary=rep('e4',nrow(tmp4))
  meta=rbind(tmp33,tmp4)
  meta$apoe_genotype.binary=as.factor(meta$apoe_genotype.binary)
  table(meta$apoe_genotype.binary)
  brain=subset(brain,cells=rownames(meta))
  
  data=brain@assays$RNA@counts
  lib_size <- colSums(data)
  norm <- t(t(data)/lib_size * median(lib_size))
  ## cell type
  tb=table(meta$cellsubtype)
  tb=tb[tb>50]
  ctps=names(sort(tb,decreasing = T))
  
  for (c in ctps){
    tmpmeta=meta[meta$cellsubtype==c,]
    tmpmeta=droplevels(tmpmeta)
    tmpdata=data[,rownames(tmpmeta)]
    tmpnorm=norm[,rownames(tmpmeta)]
    group=as.character(tmpmeta$apoe_genotype.binary)
    table(group)
    if (min(table(group))>20){
      ad=tmpmeta$ADdiag2types
      nfeat=tmpmeta$nFeature_RNA
      batch=tmpmeta$batch
      dxpark=tmpmeta$dxpark
      dlb=tmpmeta$clin_dlb
      age=tmpmeta$age_death
      sex=tmpmeta$msex
      pmi=tmpmeta$pmi
      mt=tmpmeta$percent.mt
      rp=tmpmeta$percent.rp
      #### MAST   age+sex+pmi+tdp
      log_counts <- log2(tmpdata + 1)
      fData <- data.frame(names = rownames(log_counts))
      rownames(fData) <- rownames(log_counts)
      cData <- data.frame(cond = group,ad=ad,batch=batch,dxpark=dxpark,dlb=dlb,age=age,sex=sex,pmi=pmi,mt=mt,rp=rp)
      cData$cond=relevel(factor(cData$cond),"e3")
      rownames(cData) <- colnames(log_counts)
      obj <- FromMatrix(as.matrix(log_counts), cData, fData)
      colData(obj)$cngeneson <- scale(colSums(assay(obj) > 0))
      cond <- factor(colData(obj)$cond)
      if (length(table(batch))>1){
        zlmCond <- zlm(~ cond + cngeneson+ad+batch+dxpark+dlb+age+sex+pmi+mt+rp, obj) 
      }else{zlmCond <- zlm(~ cond + cngeneson+ad+dxpark+dlb+age+sex+pmi+mt+rp, obj) }
      summaryCond <- summary(zlmCond, doLRT = "conde4")
      
      summaryDt <- summaryCond$datatable
      res <- merge(summaryDt[contrast=='conde4' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                   summaryDt[contrast=='conde4' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
      
      res[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
      res=as.data.frame(res)
      
      write.table(res,file=paste(paste(fid,c,sep = '.'),'.APOE.MASTres.txt',sep=''),sep = '\t',row.names = F,quote = F)
    }
  }
 
}







