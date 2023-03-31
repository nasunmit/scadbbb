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
## APOE33
for (f in files){
  fid=strsplit(rev(strsplit(f,'/')[[1]])[1],'[.]rds')[[1]][1]
  brain=readRDS(f)
  meta=brain@meta.data
 # brain=subset(brain,cells=rownames(meta[meta$celltype != 'Ependymal',]))
  brain=subset(brain,cells=rownames(meta[meta$apoe_genotype=='33',])) ## apoe 33 
  meta=brain@meta.data
  meta=droplevels(meta)
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
    cog=tmpmeta$cogn_global_lv
    
    group=tmpmeta$ADdiag2types
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
    cData <- data.frame(cond = cog,group=group,batch=batch,dxpark=dxpark,dlb=dlb,age=age,sex=sex,pmi=pmi,mt=mt,rp=rp)
    #cData$cond=relevel(factor(cData$cond),"nonAD")
    rownames(cData) <- colnames(log_counts)
    obj <- FromMatrix(as.matrix(log_counts), cData, fData)
    colData(obj)$cngeneson <- scale(colSums(assay(obj) > 0))
    cond <- factor(colData(obj)$cond)
    if (length(table(batch))>1){
      zlmCond <- zlm(~ cond + group+cngeneson+batch+dxpark+dlb+age+sex+pmi+mt+rp, obj)
    }else{zlmCond <- zlm(~ cond + group+cngeneson+dxpark+dlb+age+sex+pmi+mt+rp, obj) }
    summaryCond <- summary(zlmCond, doLRT = "cond")
    
    summaryDt <- summaryCond$datatable
    res <- merge(summaryDt[contrast=='cond' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                 summaryDt[contrast=='cond' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
    
    res[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
    res=as.data.frame(res)
    
    write.table(res,file=paste(paste(fid,c,sep = '.'),'.APOE3.MASTres.txt',sep=''),sep = '\t',row.names = F,quote = F)
  }
  
 
}

## APOE34 + APOE44
for (f in files){
  fid=strsplit(rev(strsplit(f,'/')[[1]])[1],'[.]rds')[[1]][1]
  brain=readRDS(f)
  meta=brain@meta.data
 # brain=subset(brain,cells=rownames(meta[meta$celltype != 'Ependymal',]))
  brain=subset(brain,cells=rownames(meta[meta$apoe_genotype %in% c('34','44'),])) ## apoe 4
  meta=brain@meta.data
  meta=droplevels(meta)
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
    cog=tmpmeta$cogn_global_lv
    
    group=tmpmeta$ADdiag2types
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
    cData <- data.frame(cond = cog,group=group,batch=batch,dxpark=dxpark,dlb=dlb,age=age,sex=sex,pmi=pmi,mt=mt,rp=rp)
    #cData$cond=relevel(factor(cData$cond),"nonAD")
    rownames(cData) <- colnames(log_counts)
    obj <- FromMatrix(as.matrix(log_counts), cData, fData)
    colData(obj)$cngeneson <- scale(colSums(assay(obj) > 0))
    cond <- factor(colData(obj)$cond)
    if (length(table(batch))>1){
      zlmCond <- zlm(~ cond + group+cngeneson+batch+dxpark+dlb+age+sex+pmi+mt+rp, obj)
    }else{zlmCond <- zlm(~ cond + group+cngeneson+dxpark+dlb+age+sex+pmi+mt+rp, obj) }
    summaryCond <- summary(zlmCond, doLRT = "cond")
    
    summaryDt <- summaryCond$datatable
    res <- merge(summaryDt[contrast=='cond' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                 summaryDt[contrast=='cond' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
    
    res[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
    res=as.data.frame(res)
    
    write.table(res,file=paste(paste(fid,c,sep = '.'),'.APOE4.MASTres.txt',sep=''),sep = '\t',row.names = F,quote = F)
  }
  
  
}


