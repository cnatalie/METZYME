#Correspondence analysis (CA) on dinoflagellate TPM-normalized transcript counts annotated at the KOG level

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

a<-read.csv('TPM_TRANSCRIPTS_Dino.lpi0.8_only_annotations_orf_KOGpost.csv')
head(a)
rownames(a)<-a$X
a<-a[,-1]
data<-a
meta<-read.csv('CCA_meta.csv')
rownames(meta)<-meta$ID
meta<-meta[,-1]
stand<-decostand(data, method = "hellinger") #hellinger transformed
ca<-cca(stand)
summary(ca) #Summary includes proportion of variation explained by each component

meta1<-meta #Clean metadata file and remove highly co-linear parameters (e.g., nitrate/phosphate, oxygen/temperature)
meta<-meta[,-1]
meta<-meta[,-1]
meta<-meta[,-1]
meta<-meta[,-2]
ef <- envfit(ca, meta, permu = 999)
meta<-meta1 #Revert back to original table with site and depth information
plot(ef)
ef #Results of linear regression vector fitting of environmental parameters to ordination

scores<-scores(ca)
scores <- data.frame(ca$CA$u)
uscores <- inner_join(rownames_to_column(data), rownames_to_column(data.frame(scores)), type = "right", by = "rowname")
vscores <- as.data.frame(scores(ef, display = "vectors"))
vscores <- cbind(vscores, Species = rownames(vscores))
vscores$env<-rownames(vscores)

#Graph ordination with ggplot2
ggplot(uscores) + theme_bw() + scale_shape_manual(values = c(21:25,8,9)) + theme(strip.text.y = element_text(angle = 0)) + geom_segment(data = vscores, aes(x = 0, y = 0, xend = vscores$CA1, yend = vscores$CA2), arrow=arrow(length=unit(0.2,"cm")),
                                                                                                                                        alpha = 1, color = 'black')+ geom_label(data=vscores, aes(x=CA1, y=CA2, label = env)) +
  geom_point(aes(x = CA1, y = CA2, col = meta$Depth,
                 shape = meta$Site, size=4)) +
  scale_colour_gradient(low="#A2FEFF",high="black") +
  scale_fill_manual(values=c("#A2FEFF",'black') +
                      geom_segment(data = vscores, aes(x = 0, y = 0, xend = vscores$CA1, yend = vscores$CA2), arrow=arrow(length=unit(0.2,"cm")),
                                   alpha = 1, color = 'black')+ geom_label(data=vscores, aes(x=CA1, y=CA2, label = env)))


