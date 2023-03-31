
##########
library(qqman)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(qlcMatrix)
library(tidyr)

dataselnew=readRDS('data/GWASdata_20.forManhattanPlot.rds')
#####
mat=read.table('data/adDEGs.GWASgenes.matrix.txt',header = T,sep = '\t')
genes=rownames(mat)


# Prepare the dataset
don <- dataselnew %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(dataselnew, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate( BPcum=BP+tot) %>%
  
  # Add highlight and annotation information
  mutate( is_annotate=ifelse(Gene %in% genes, "yes", "no")) 

axisdf = don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

test=don[don$is_annotate=='yes',]
table(test$Gene)

# Make the plot
options(ggrepel.max.overlaps = Inf)

p1=ggplot(don, aes(x=BPcum, y=-log10(Pfinal))) +
  # Show all points
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=0.5) +
  #geom_point( data=subset(don, is_annotate=="yes"),aes(color=as.factor(CHR)), alpha=0.8, size=0.5) +
  scale_color_manual(values = rep(c("grey", "black"), 22 )) +
  ylim(0,20)+
  
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
 # scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=Gene), size=2) +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

p2=ggplot(don, aes(x=BPcum, y=-log10(Pfinal))) +
  # Show all points
  #geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=0.5) +
  geom_point( data=subset(don, is_annotate=="yes"),aes(color=as.factor(CHR)), alpha=0.8, size=0.5) +
  scale_color_manual(values = rep(c("grey", "black"), 22 )) +
  ylim(0,20)+
  
  # custom X axis:
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  # scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  
  # Add label using ggrepel to avoid overlapping
  geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=Gene), size=2) +
  
  # Custom the theme:
  theme_bw() +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank()
  )

pdf('AD.GWAS_20.ManhattanPlot.withDEGslabled.pdf',height = 6,width = 16) ## max p-value 50
p1+ geom_hline(yintercept = 5,col='red')
p2+ geom_hline(yintercept = 5,col='red')
dev.off()







