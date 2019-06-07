#CCA plots
'''
Example for dinoflagellate TPM-normalized transcript counts annotated at the KOG level
Helpful resources:
http://www.sthda.com/english/articles/31-principal-component-methods-in-r-practical-guide/113-ca-correspondence-analysis-in-r-essentials/#graph-of-row-variables
http://cc.oulu.fi/~jarioksa/opetus/metodi/sessio2.pdf
http://cc.oulu.fi/~jarioksa/opetus/metodi/ordination101.html#94
'''
library(factoextra)
library(FactoMineR)
library(vegan)
library(labdsv)
library(dplyr)
library(tibble)
library(ggplot2)

a<-read.csv('TPM_TRANSCRIPTS_Dino.lpi0.8_only_annotations_KOGpost.csv') #Read in count data
rownames(a)<-a$X
a<-a[,-1]
data<-a
meta<-read.csv('CCA_meta.csv') #Read in metadata
rownames(meta)<-meta$ID
meta<-meta[,-1]
cca<-cca(data ~ Co + Temperature + Fe + NO3 + NH4, data=meta)
vif.cca(cca) #variance inflation factor to address collinearity of environmental parameters
anova(cca)
anova(cca, by="term") #permutation test (sequential test)
anova(cca, by="margin") #permutation test (marginal test)
summary(cca) #Proportion explained for x and y axes (CCA1 and CCA2)

#Ploting the ordination
scores<-scores(cca)
scores <- data.frame(cca$CCA$u)
uscores <- inner_join(rownames_to_column(data), rownames_to_column(data.frame(scores)), type = "right", by = "rowname")
vscores <- data.frame(cca$CCA$biplot)
vscores$env<-rownames(vscores)

ggplot(uscores) + theme_bw() + scale_shape_manual(values = c(21:25,8,9)) + theme(strip.text.y = element_text(angle = 0)) + geom_segment(data = vscores, aes(x = 0, y = 0, xend = vscores$CCA1*2, yend = vscores$CCA2*2), arrow=arrow(length=unit(0.2,"cm")),
                               alpha = 1, color = 'black')+ geom_label(data=vscores, aes(x=CCA1, y=CCA2, label = env)) +
  geom_point(aes(x = CCA1, y = CCA2, col = meta$Depth,
                 shape = meta$Site, size=4)) +
  scale_colour_gradient(low="#A2FEFF",high="black") +
  scale_fill_manual(values=c("#A2FEFF",'black') +
  geom_segment(data = vscores, aes(x = 0, y = 0, xend = vscores$CCA1*2, yend = vscores$CCA2*2), arrow=arrow(length=unit(0.2,"cm")),
               alpha = 1, color = 'black')+ geom_label(data=vscores, aes(x=CCA1, y=CCA2, label = env)))
