#Principal Coordinate Analysis (PCA) on dinoflagellate TPM-normalized transcript counts 

'''
Helpful resources:
http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/113-ca-correspondence-analysis-in-r-essentials/#graph-of-row-variables
http://cc.oulu.fi/~jarioksa/opetus/metodi/sessio2.pdf
http://cc.oulu.fi/~jarioksa/opetus/metodi/ordination101.html#94
'''

library(vegan)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(tibble)

a<-read.csv('TPM_TRANSCRIPTS_Dino.lpi0.8_only_annotations_orf_allcontigs.csv')
head(a)
rownames(a)<-a$X
rownames(a) <- gsub("[.y]", "", rownames(a))
a<-a[,-1]
a<-t(a)
data<-a
meta<-read.csv('CCA_meta.csv')
rownames(meta)<-meta$ID
meta<-meta[,-1]

data<-log(data+1)
pca <- rda(data, distance = "euclidean")
plot(pca)


meta1<-meta #Clean metadata file and remove highly co-linear parameters (e.g., nitrate/phosphate, oxygen/temperature)
meta<-meta[,-1]
meta<-meta[,-1]
meta<-meta[,-1]
meta<-meta[,-2]
ef <- envfit(ca, meta, permu = 999)
meta<-meta1 #Revert back to original table with site and depth information
plot(ef, col = 'red')
ef #Results of linear regression vector fitting of environmental parameters to ordination
meta<-meta1 #Revert back to original table with site and depth information

ef.score <- scores(ef, "vectors", choices = 1:2)
efvec <- ef.score * ordiArrowMul(ef)

scores<-scores(pca)
scores <- data.frame(scores$sites)
uscores <- inner_join(rownames_to_column(data), rownames_to_column(data.frame(scores)), type = "right", by = "rowname")
vscores <- as.data.frame(efvec)
vscores <- cbind(vscores, Species = rownames(vscores))
vscores$env<-rownames(vscores)

eig<- eigenvals(pca)
eig<- eig / sum(eig) #Percent variation explained by each component/axis

#Graph ordination with ggplot2
ggplot(uscores) + theme_bw() + scale_shape_manual(values = c(21:25,8,9)) + theme(strip.text.y = element_text(angle = 0)) + geom_segment(data = vscores, aes(x = 0, y = 0, xend = vscores$PC1, yend = vscores$PC2),
                                                                                                                                        alpha = 1, color = 'black')+ geom_label(data=vscores, aes(x=PC1, y=PC2, label = env)) +
  geom_point(aes(x = PC1, y = PC2, col = meta$Depth,
                 shape = meta$Site, size=4)) +
  scale_colour_gradient(low="#A2FEFF",high="black") +
  scale_fill_manual(values=c("#A2FEFF",'black') +
                      geom_segment(data = vscores, aes(x = 0, y = 0, xend = vscores$PC1, yend = vscores$PC2), 
                                   alpha = 1, color = 'black')+ geom_label(data=vscores, aes(x=PC1, y=PC2, label = env)))

