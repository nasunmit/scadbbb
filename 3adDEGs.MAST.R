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
    group=as.character(tmpmeta$ADdiag2types)
    table(group)
    
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
    cData <- data.frame(cond = group,batch=batch,dxpark=dxpark,dlb=dlb,age=age,sex=sex,pmi=pmi,mt=mt,rp=rp)
    cData$cond=relevel(factor(cData$cond),"nonAD")
    rownames(cData) <- colnames(log_counts)
    obj <- FromMatrix(as.matrix(log_counts), cData, fData)
    colData(obj)$cngeneson <- scale(colSums(assay(obj) > 0))
    cond <- factor(colData(obj)$cond)
    if (length(table(batch))>1){
      zlmCond <- zlm(~ cond + cngeneson+batch+dxpark+dlb+age+sex+pmi+mt+rp, obj)
    }else{zlmCond <- zlm(~ cond + cngeneson+dxpark+dlb+age+sex+pmi+mt+rp, obj) }
    summaryCond <- summary(zlmCond, doLRT = "condAD")
    
    summaryDt <- summaryCond$datatable
    res <- merge(summaryDt[contrast=='condAD' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                 summaryDt[contrast=='condAD' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
    
    res[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
    res=as.data.frame(res)
    
    write.table(res,file=paste(paste(fid,c,sep = '.'),'.MASTres.txt',sep=''),sep = '\t',row.names = F,quote = F)
  }
  
 
}
