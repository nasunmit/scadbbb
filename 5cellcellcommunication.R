##########################
library(scales)
library(plyr)
library(dplyr)
library(ComplexHeatmap)
library(RColorBrewer)
wr <- colorRampPalette(colors = c( "white", "red"))(100)
rwb <- colorRampPalette(colors = c("blue", "white", "red"))(100)
ryb=colorRampPalette(rev(brewer.pal(n = 10, name ="RdYlBu")))(100)
prgn=colorRampPalette(rev(brewer.pal(n = 10, name ="PRGn")))(100)
wb <- colorRampPalette(colors = c( "white", "blue"))(100)
###########
###  STEP 1: 
##### vas vs all
### co-variation
allfiles=list.files(pattern = 'module.expMat.txt') ### Matrix of average expression of module across all individuals: row (each module), column (each individual)
vasfiles=allfiles[grep('vas',allfiles)]

res=c()
for (f1 in vasfiles){
  vasdata=read.table(f1,header = T,sep = '\t')
  for (f2 in allfiles){
    if (f1 != f2){
      otherdata=read.table(f2,header = T,sep = '\t')
      ov=intersect(colnames(vasdata),colnames(otherdata))
      if (length(ov)>30){
        vasdatasel=vasdata[,ov]
        otherdatasel=otherdata[,ov]
        dim(vasdatasel)
        dim(otherdatasel)
        
        for (i in c(1:nrow(vasdatasel))){
          for (j in c(1:nrow(otherdatasel))){
            test=cor.test(as.numeric(vasdatasel[i,]),as.numeric(otherdatasel[j,]))
            pval=test$p.value
            pc=test$estimate
            res=rbind(res,c(rownames(vasdatasel)[i],rownames(otherdatasel)[j],pval,pc))
          }
        }
      }
      
    }
  }
  
}
dim(res)
res=as.data.frame(res)
colnames(res)=c('celltype1','celltype2','pvalue','pcc')
res$fdr=p.adjust(res$pvalue,method = 'fdr',n=nrow(res))
write.table(res,file='vas_all.modules_correlation_test.result.txt',sep = '\t',quote = F)

#### Step2: filtering 
ligrec=read.table('CellChat_Phone_Talk_SingleCellSignalR.combined.ligrec.txt',header = T,sep = '\t') ## ligand-receptor pairs: Lig_Rec as id, Ligand, Receptor
dim(ligrec)
rownames(ligrec)=ligrec$LigRec

load('adDEGs.modules.allcelltype.RData') ## rename it as module if not


res=res[res$fdr<0.01,]

output=c()
for (i in c(1:nrow(res))){
  gset1=module[[res[i,1]]]
  gset2=module[[res[i,2]]]
  ligrec.sel1=ligrec[ligrec$Lig %in% gset1 & ligrec$Rec %in% gset2,]
  if(nrow(ligrec.sel1)>0){
    print(ligrec.sel1)
    ligrec.sel1=cbind(ligrec.sel1,rep(res[i,1],nrow(ligrec.sel1)),rep(res[i,2],nrow(ligrec.sel1)))
    ligrec.sel1=cbind(ligrec.sel1,rep(res[i,4],nrow(ligrec.sel1)),rep(res[i,5],nrow(ligrec.sel1)))
    colnames(ligrec.sel1)=c(colnames(ligrec),'sender','receiver','pcc','fdr')
    output=rbind(output,ligrec.sel1)
  }
  ligrec.sel2=ligrec[ligrec$Lig %in% gset2 & ligrec$Rec %in% gset1,]
  if(nrow(ligrec.sel2)>0){
    print(ligrec.sel2)
    ligrec.sel2=cbind(ligrec.sel2,rep(res[i,2],nrow(ligrec.sel2)),rep(res[i,1],nrow(ligrec.sel2)))
    ligrec.sel2=cbind(ligrec.sel2,rep(res[i,4],nrow(ligrec.sel2)),rep(res[i,5],nrow(ligrec.sel2)))
    colnames(ligrec.sel2)=c(colnames(ligrec),'sender','receiver','pcc','fdr')
    output=rbind(output,ligrec.sel2)
  }
}
output=as.data.frame(output)
dim(output)
head(output)

load('adDEGs.pval01_coef02.genelist.RData') ## all DEGs in other cell types
otherdegset=degset
load('vas.adDEGs.pval01_coef02.genelist.RData') ## DEGs in vascular cell types
vasdegset=degset
degset=c(vasdegset,otherdegset)

tmp=c()
for (i in c(1:nrow(output))){
  c1=strsplit(as.character(output[i,6]),'[.]M')[[1]][1]
  c2=strsplit(as.character(output[i,7]),'[.]M')[[1]][1]
  g1=as.character(output[i,1])
  g2=as.character(output[i,2])
  val=c()
  if (g1 %in% degset[[paste(c1,'up',sep = '.')]]){val=c(val,'up')}else{val=c(val,'down')}
  if (g2 %in% degset[[paste(c2,'up',sep = '.')]]){val=c(val,'up')}else{val=c(val,'down')}
  tmp=rbind(tmp,val)
}
output=cbind(output,tmp)
colnames(output)[10:11]=c('sender.AD','receiver.AD')

outputnew=c()
for (i in names(output)){
  outputnew=cbind(outputnew,as.character(output[[i]]))
}
outputnew=as.data.frame(outputnew)
colnames(outputnew)=names(output)

dim(outputnew)
head(outputnew)

write.table(outputnew,file='vas_all.modules_interaction_lig_rec.res.txt',row.names = F,sep = '\t',quote = F)

#### Step3: filtering by GO terms and plotting
data=read.table('vas_all.modules_interaction_lig_rec.res.txt',header = T,sep = '\t')

data$type=paste(data$sender.AD,data$receiver.AD,sep = '.')

gofiles=list.files(pattern = 'enrichedGO',full.names = T) ### enriched GO terms of adDEGs in each module (output of enrichr in R). This step can be ignored if don't want to filter pairs based on functions of module.

datasel=c()
for (i in c(1:nrow(data))){
  s=data[i,6]
  r=data[i,7]
  sfile=gofiles[grep(paste(s,'.enriched',sep=''),gofiles)]
  n=0
  if (length(sfile)>0){
    for (f in sfile){
      df=read.table(f,header = T,sep = '\t')
      genes=unique(strsplit(paste(df$Genes,collapse = ';'),';')[[1]])
      if (data[i,1] %in% genes){n=n+1}
    }
  }
  
  rfile=gofiles[grep(paste(r,'.enriched',sep=''),gofiles)]
  m=0
  if (length(rfile)>0){
    for (f in rfile){
      df=read.table(f,header = T,sep = '\t')
      genes=unique(strsplit(paste(df$Genes,collapse = ';'),';')[[1]])
      if (data[i,2] %in% genes){m=m+1}
    }
  }
  
  if(n>0 & m>0){datasel=rbind(datasel,data[i,])}
}

dim(data)
dim(datasel)
datasel=datasel[!duplicated(datasel),]
dim(datasel)
datasel=datasel[-(intersect(grep('Endo',datasel$sender),grep('Endo',datasel$receiver))),]
dim(datasel)
c1=c()
c2=c()
for (i in c(1:nrow(datasel))){c1=c(c1,strsplit(datasel[i,6],'[.]')[[1]][1]);c2=c(c2,strsplit(datasel[i,7],'[.]')[[1]][1])}
datasel$sender.celltype=c1
datasel$receiver.celltype=c2
head(datasel)
table(datasel$type)
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

ctps=unique(c(datasel$sender.celltype,datasel$receiver.celltype))
ov=intersect(names(col2),ctps)
mycolor <- col.all[1:(length(ctps)-length(ov))]
names(mycolor)=setdiff(ctps,ov)
mycolor=c(mycolor,col2[ov])

##########################  gain of interactions
write.table(datasel,file='all.interation.txt',sep = '\t',quote = F,row.names = F)

data.up=datasel[datasel$type=='up.up',]

write.table(data.up,file='gain_interation.txt',sep = '\t',quote = F,row.names = F)


df=data.up[,c(13,14)]
df$num=rep(1,nrow(df))
library(dplyr)
dfsum=df %>% 
  group_by(sender.celltype,receiver.celltype) %>% 
  summarise(across(everything(), sum))

data_long=as.data.frame(dfsum)
colnames(data_long)=c('rowname','key','value')

n=length(unique(c(data_long$rowname,data_long$key)))
# parameters
circos.clear()
circos.par(start.degree = 90, gap.degree = 4, track.margin = c(-0.1, 0.1), points.overflow.warning = FALSE)
par(mar = rep(0, 4))



pdf('gain_interation.chordDiagram.pdf',width = 8,height = 8)
# Base plot
chordDiagram(
  x = data_long, 
  grid.col = mycolor,
  transparency = 0.25,
  directional = 1,
  direction.type = c("arrows", "diffHeight"), 
  diffHeight  = -0.04,
  annotationTrack = "grid", 
  annotationTrackHeight = c(0.05, 0.1),
  link.arr.type = "big.arrow", 
  link.sort = TRUE, 
  link.largest.ontop = TRUE)

# Add text and axis
circos.trackPlotRegion(
  track.index = 1, 
  bg.border = NA, 
  panel.fun = function(x, y) {
    
    xlim = get.cell.meta.data("xlim")
    sector.index = get.cell.meta.data("sector.index")
    
    # Add names to the sector. 
    circos.text(
      x = mean(xlim), 
      y = 2.5, 
      labels = sector.index, 
      facing = "bending", 
      cex = 1
    )

  }
)

dev.off()

##########################  loss of interactions

data.down=datasel[!(datasel$type=='up.up'),]

write.table(data.down,file='loss_interation.txt',sep = '\t',quote = F,row.names = F)


df=data.down[,c(13,14)]
df$num=rep(1,nrow(df))
library(dplyr)
dfsum=df %>% 
  group_by(sender.celltype,receiver.celltype) %>% 
  summarise(across(everything(), sum))

data_long=as.data.frame(dfsum)
colnames(data_long)=c('rowname','key','value')

n=length(unique(c(data_long$rowname,data_long$key)))
# parameters
circos.clear()
circos.par(start.degree = 90, gap.degree = 4, track.margin = c(-0.1, 0.1), points.overflow.warning = FALSE)
par(mar = rep(0, 4))


pdf('loss_interation.chordDiagram.pdf',width = 8,height = 8)
# Base plot
chordDiagram(
  x = data_long, 
  grid.col = mycolor,
  transparency = 0.25,
  directional = 1,
  direction.type = c("arrows", "diffHeight"), 
  diffHeight  = -0.04,
  annotationTrack = "grid", 
  annotationTrackHeight = c(0.05, 0.1),
  link.arr.type = "big.arrow", 
  link.sort = TRUE, 
  link.largest.ontop = TRUE)

# Add text and axis
circos.trackPlotRegion(
  track.index = 1, 
  bg.border = NA, 
  panel.fun = function(x, y) {
    
    xlim = get.cell.meta.data("xlim")
    sector.index = get.cell.meta.data("sector.index")
    
    # Add names to the sector. 
    circos.text(
      x = mean(xlim), 
      y = 2.5, 
      labels = sector.index, 
      facing = "bending", 
      cex = 1
    )
  }
)

dev.off()


datasel=datasel[,c(1,2,3,12:14)]
datasel=datasel[!duplicated(datasel),]
write.table(datasel,file = 'all.interation.short.txt',sep = '\t',quote = F,row.names = F)








