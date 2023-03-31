###########3
load('vas.adDEGs.pval01_coef02.genelist.RData')
str(degset)
### 
library(enrichR)
options(base_address = "http://amp.pharm.mssm.edu/Enrichr/")
dbs <- listEnrichrDbs()

cutoff=0.01
mydbs <- c("TRANSFAC_and_JASPAR_PWMs","ChEA_2016",'ENCODE_TF_ChIP-seq_2015')
for (c in names(degset)){
  genes=as.character(degset[[c]])
  enriched <- enrichr(genes, mydbs)
  
  out=enriched$TRANSFAC_and_JASPAR_PWMs
  sel=out[out$P.value<cutoff,]
  
  out=enriched$ChEA_2016
  sel=rbind(sel,out[out$P.value<cutoff,])
  
  out=enriched[['ENCODE_TF_ChIP-seq_2015']]
  sel=rbind(sel,out[out$P.value<cutoff,])
  if (nrow(sel)>0){
    write.table(sel,file=paste(paste('vas.adDEGs',c,sep='.'),'.enrichedTFs.txt',sep=''),sep='\t',quote = F)
  }
  
  
}

