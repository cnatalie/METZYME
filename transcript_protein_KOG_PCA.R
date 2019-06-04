######PCA plots
#TRANSCRIPTS
http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/113-ca-correspondence-analysis-in-r-essentials/#graph-of-row-variables
  
library(factoextra)
library(FactoMineR)
library(vegan)
library(labdsv)
library(dplyr)
library(tibble)
library(ggplot2)

a<-read.csv('TPM_TRANSCRIPTS_Dino.lpi0.8_only_annotations_KOGpost.csv')
head(a)
rownames(a)<-a$X
a<-a[,-1]
data<-a
meta<-read.csv('CCA_meta.csv')
rownames(meta)<-meta$ID
meta<-meta[,-1]
stand<-decostand(data, method = "hellinger") #hellinger transformed PCA
pca<- rda(stand, scale = FALSE)
ef <- envfit(pca, meta, permu = 999)
scores<-scores(pca)
scores <- data.frame(scores$sites)
uscores <- inner_join(rownames_to_column(data), rownames_to_column(data.frame(scores)), type = "right", by = "rowname")
vscores <- data.frame(ef$vectors$arrows)
vscores$env<-rownames(vscores)

eig<- eigenvals(pca)
eig<- eig / sum(eig)

ggplot(uscores) + 
  geom_point(aes(x = PC1, y = PC2, col = meta$Depth,
                 shape = meta$Site, size=4)) +
  scale_colour_gradient(low="#A2FEFF",high="black") +
  scale_fill_manual(values=c("#A2FEFF",'black')) +
  scale_shape_manual(values = c(21:25,8,9)) +
  geom_segment(data = vscores, aes(x = 0, y = 0, xend = vscores$PC1*2, yend = vscores$PC2*2), arrow=arrow(length=unit(0.2,"cm")),
               alpha = 1, color = 'black')+ geom_label(data=vscores, aes(x=PC1, y=PC2, label = env)) +
  theme_bw() +
  theme(strip.text.y = element_text(angle = 0))
#add in eig values in adobe illustrator