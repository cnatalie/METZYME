library("phyloseq")
library("ggplot2")
library("vegan")
library("ape")

#Followed Joey's phyloseq guidelines: https://joey711.github.io/phyloseq/import-data.html
#http://joey711.github.io/phyloseq-demo/phyloseq-demo.html

#Read in OTU counts table
a<-read.csv('OTU.csv')
rownames(a)<-a$id
a<-a[,-1] #Remove ID
a<-a[,-42] #Remove 1900m sample
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

#Read in environmental data
sampledata<-read.csv('sampledata.csv', stringsAsFactors=FALSE)
x<-sampledata
rownames(sampledata)<-sampledata$X
sampledata = sample_data(data.frame(
Station = x$Station,
Depth = x$Depth,
NO3 = x$NO3,
Co = x$Co,
NH4 = x$NH4,
O2 = x$O2,
PO4=x$PO4,
Fe=x$Fe,
temp=x$Temperature,
row.names=sample_names(physeq),
stringsAsFactors=FALSE
))

physeq1 = merge_phyloseq(physeq, sampledata)
physeq1.dino = subset_taxa(physeq1, Division == "Dinoflagellata")
physeq1.dino = prune_samples(names(which(sample_sums(physeq1.dino) >= 5)), physeq1.dino) #trim OTUs with less than 5 reads

#CA plot
b<-t(otu)
b<-data.frame(b)

ca <- ordinate(physeq = physeq1.dino, method = "CCA")

#Delete all columns except parameters of interest: Co, Fe, NO3, temp and NH4
rownames(x)<-x$X
xx<-x
x<-x[,-1]
x<-x[,-1]
x<-x[,-1]
x<-x[,-1]
x<-x[,-2]
ef <- envfit(ca, x, permu = 999)
plot(ca)
plot(ef)
#After fitting vectors, revert back to original table with station/depth info
x<-xx

#Plot with ggplot2, color code by depth
scores<-scores(ca)
scores <- data.frame(ca$CA$u)
uscores <- inner_join(rownames_to_column(b), rownames_to_column(data.frame(scores)), type = "right", by = "rowname")
vscores <- data.frame(ef$vectors$arrows)
vscores$env<-rownames(vscores)

ggplot(uscores) + theme_bw() + scale_shape_manual(values = c(21:25,8,9)) +
  theme(strip.text.y = element_text(angle = 0)) + 
  geom_segment(data = vscores, aes(x = 0, y = 0, xend = vscores$CA1, yend = vscores$CA2), arrow=arrow(length=unit(0.2,"cm")), alpha = 1, color = 'black') +
  geom_label(data=vscores, aes(x=CA1, y=CA2, label = env)) +
  geom_point(aes(x = CA1, y = CA2, col = x$Depth, shape = x$Station, size=4)) +
  scale_colour_gradient(low="#A2FEFF",high="black") +
  scale_fill_manual(values=c("#A2FEFF",'black') +
                      geom_segment(data = vscores, aes(x = 0, y = 0, xend = CA1, yend = CA2), arrow=arrow(length=unit(0.2,"cm")), alpha = 1, color = 'black')+
                      geom_label(data=vscores, aes(x=CA1, y=CA2, label = env)))


#Heatmap on non-transformed counts. Helpful: https://joey711.github.io/phyloseq/plot_heatmap-examples.html#subset_a_smaller_dataset_based_on_an_archaeal_phylum
a<-read.csv('OTU.csv')
rownames(a)<-a$id
a<-a[,-1] #Remove ID
a<-a[,-42] #Remove 1900m sample
otu<-as.matrix(sapply(a, as.numeric))
rownames(otu)<-rownames(a)
otumat<-otu
OTU = otu_table(otumat, taxa_are_rows = TRUE)
b<-read.csv('TAXA_newPR2.csv')
rownames(b)<-b$id
b<-b[,-1]
taxa<-as.matrix(sapply(b, as.character))
rownames(taxa)<-rownames(b)
taxmat<-taxa
TAX = tax_table(taxmat)
physeq = phyloseq(OTU, TAX)
physeq1 = merge_phyloseq(physeq, sampledata, tree)
gpt <- subset_taxa(physeq1, Division=="Dinoflagellata")
gpt <- prune_taxa(names(sort(taxa_sums(gpt),TRUE)[1:50]), gpt)
plot_heatmap(gpt, taxa.label = "Taxa", low = "#000033", method = "NMDS", distance = "bray",
             high = "#66CCFF", na.value = "black", title = NULL, sample.order = NULL, taxa.order = NULL,
             first.sample = NULL, first.taxa = NULL) 



