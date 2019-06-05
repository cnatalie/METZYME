library("phyloseq")
library("ggplot2")
library("vegan")
library("ape")

#Followed Joey's phyloseq guidelines: https://joey711.github.io/phyloseq/import-data.html
#http://joey711.github.io/phyloseq-demo/phyloseq-demo.html

#Read in OTU counts table
a<-read.csv('OTU.csv')
rownames(a)<-a$id
a<-a[,-1]
a<-a[,-42]
otu<-as.matrix(sapply(a, as.numeric))

#Hellinger-normalized OTU counts for PCA
rownames(otu)<-rownames(a)
otumat<-otu
otumat<-decostand(otumat, method = "hellinger")
OTU = otu_table(otumat, taxa_are_rows = TRUE)

#Read in OTU taxonomy 
b<-read.csv('TAXA_newPR2.csv')
rownames(b)<-b$id
b<-b[,-1]
taxa<-as.matrix(sapply(b, as.character))
rownames(taxa)<-rownames(b)
taxmat<-taxa


TAX = tax_table(taxmat)
physeq = phyloseq(OTU, TAX)
tree = read.tree("18S.tre") #Alignment with MUSCLE on OTU fasta file, Maximum-likelihood tree with partial deletions and 100 bootstraps

#PCA
library(dplyr)
library(tibble)
b<-t(otu)
b<-data.frame(b)
scores<-scores(pca)
scores <- data.frame(scores$sites)
uscores <- inner_join(rownames_to_column(b), rownames_to_column(data.frame(scores)), type = "right", by = "rowname")
vscores <- data.frame(ef$vectors$arrows)
vscores$env<-rownames(vscores)
eig<- eigenvals(pca)
eig<- eig / sum(eig)
ggplot(uscores) + 
  geom_point(aes(x = PC1, y = PC2, col = x$Depth,
                 shape = x$Station, size=4)) + scale_shape_manual(values = c(21:25,8,9)) + 
  geom_segment(data = vscores, aes(x = 0, y = 0, xend = vscores$PC1*2, yend = vscores$PC2*2), arrow=arrow(length=unit(0.2,"cm")),
               alpha = 1, color = 'black')+ geom_label(data=vscores, aes(x=PC1, y=PC2, label = env)) + theme_bw() + theme(strip.text.y = element_text(angle = 0)) +
  scale_colour_gradient(low="#A2FEFF",high="black") +
  scale_fill_manual(values=c("#A2FEFF","black")) 

#HEATMAP on non-transformed counts. Helpful: https://joey711.github.io/phyloseq/plot_heatmap-examples.html#subset_a_smaller_dataset_based_on_an_archaeal_phylum
physeq1 = merge_phyloseq(physeq, sampledata, tree)
gpt <- subset_taxa(physeq1, Division=="Dinoflagellata")
gpt <- prune_taxa(names(sort(taxa_sums(gpt),TRUE)[1:50]), gpt)

plot_heatmap(gpt, taxa.label = "Taxa", low = "#000033", method = "NMDS", distance = "bray",
             high = "#66CCFF", na.value = "black", title = NULL, sample.order = NULL, taxa.order = NULL,
             first.sample = NULL, first.taxa = NULL) 



