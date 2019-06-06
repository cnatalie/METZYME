###Calculating normalized transcripts as TPM, caculated following Micheal Love's blog post: https://support.bioconductor.org/p/91218/

#Read in dinoflagellate only raw transcripts, annotations and gene length 
a<-read.csv('annotation_all.filtered.grps.go_TRANSCRIPTS_0.8lpi_Dino_only_sansAnnotations.csv')
rownames(a)<-a$orf_id
x <- a[3:44] / a$contig_length
tpm.mat <- t( t(x) * 1e6 / colSums(x) )
write.csv(tpm.mat, 'TPM_TRANSCRIPTS_Dino_only.csv')
b<-read.csv('TPM_TRANSCRIPTS_Dino_only.csv')
colnames(b)[colnames(b) == 'X'] <- 'orf_id' #rownames = orf_id
#Read in file with all annotations
a<-read.csv('annotation_all.filtered.grps.go_TRANSCRIPTS_0.8pli_Dino_only.csv')
c<-merge(a, b, by='orf_id')
write.csv(c, 'TPM_TRANSCRIPTS_Dino.lpi0.8_only_annotations.csv')

library(pheatmap)
library(genefilter)
library(viridis)
library(RColorBrewer)

### TPM-normalized KEGG heatmap (three different taxa)
a<-read.csv('TPM_TRANSCRIPTS_Dino_diatom_pelagophyte.lpi0.8_annotations_KOpresum.csv')
dia<-c('Diatoms')
b<-a[a$group %in% dia,]
head(b)
library(caroline)
diatoms<-groupBy(b, by='KO',clmns=(3:43),aggregation='sum')
colnames(diatoms)<-colnames(b[3:43])
colnames(diatoms) <- paste("Diatoms", colnames(diatoms), sep = "_")
diatoms$KO<-rownames(diatoms)
dino<-c('Dinophyta')
b<-a[a$group %in% dino,]
dino<-groupBy(b, by='KO',clmns=(3:43),aggregation='sum')
colnames(dino)<-colnames(b[3:43])
colnames(dino) <- paste("Dinoflagellates", colnames(dino), sep = "_")
dino$KO<-rownames(dino)
dia<-c('Other Stramenopiles')
b<-a[a$group %in% dia,]
head(b)
pelago<-groupBy(b, by='KO',clmns=(3:43),aggregation='sum')
colnames(pelago)<-colnames(b[3:43])
colnames(pelago) <- paste("Pelagophytes", colnames(pelago), sep = "_")
pelago$KO<-rownames(pelago)
c<-merge(diatoms, dino, by='KO')
d<-merge(c, pelago, by='KO')
library(viridis)
library(pheatmap)
paletteLength=100
myColor <- rev(viridis_pal(option = "B")(paletteLength))
annotation <- data.frame(Var1 = factor(1:124, labels = c('1')))
rownames(annotation)<-colnames(d)
annotation$Var1<-rownames(annotation)
fix(annotation) #Fill in dataframe with Group, Station and Depth
rownames(annotation) <- colnames(d)
colnames(annotation)<-c('Group','Station','Depth')
annotation<-annotation[-1,]
tail(annotation)
rownames(d)<-d$KO
d<-d[,-1]
z<-log2(d+1)
pheatmap(z, color=myColor, cluster_cols=T, fontsize_row=8, cluster_rows=T, show_colnames=F, show_rownames = F, annotation = annotation)

#PFams heatmaps with proteins
setwd("~/METZYME/metaproteome 3um")
a<-read.csv('exclusive_counts_annotations_dino_lpi0.8_pseudoTPM_PFAMS_post.csv')
rownames(a)<-a$X
b<-a[,-1]
b<-b[,-40]
head(b)
c<-log2(b+1)
library(pheatmap)
library(genefilter)
library(viridis)
library(RColorBrewer)
paletteLength <- 50
myColor <- colorRampPalette(brewer.pal(9, "YlGnBu"))(100)
rv <- apply(b, 1, var)
idx <- names(V[order(V, decreasing = T)][1:40])
nrow(b)
head(c)
annotation <- data.frame(Var1 = factor(1:41, labels = c('1')))
rownames(annotation)<-colnames(b)
annotation$Var1<-rownames(annotation)
fix(annotation)
colnames(annotation)<-c('Depth',"Station")

#To flip a dendrogram branch
callback = function(hc, mat){
         sv = svd(t(mat))$v[,1]
         dend = reorder(as.dendrogram(hc), wts = sv)
         as.hclust(dend)
}

pheatmap(c[idx,], scale="row", color=myColor, 
         cluster_cols=T, fontsize_row=8, 
         cluster_rows=T, cellwidth=9, 
         cellheight=9, show_colnames=T, 
         show_rownames = T, annotation = annotation, clustering_callback = callback)
}

### transcripts KEGG heatmap
library(caroline)
a<-read.csv('TPM_TRANSCRIPTS_Dino.lpi0.8_KO_postsum.csv')
x<-read.tab('kodef.tab') #Read in a table of KEGG IDs and definitions
x$KO_def<-paste(x$KO, x$def, sep='_')
a$KO<-rownames(a)
rownames(a)<-a$X
a$KO<-rownames(a)
head(a)
z<-merge(a, x, by='KO')
rownames(z)<-z$KO_def
head(z)
z<-z[,-1] #Clean up by deleting text columns, leave only data matrix
z<-z[,-45]
z<-z[,-44]
head(z)
z<-z[,-43]
d<-log2(z+1)
head(z)
z<-z[,-42] #Delete 1900m sample
rv<- rowVars(z) #Sort by variance
idx<- order(-rv)[1:50] #Show top 50 genes with highest variances
d<-log2(z+1)
annotation <- data.frame(Var1 = factor(1:41, labels = c('1')))
rownames(annotation)<-colnames(d)
annotation$Var1<-rownames(annotation)
fix(annotation) #Fill in annotation table
colnames(annotation)<-c('Depth',"Station")
pheatmap(d[idx,], scale="row", color=myColor,
         cluster_cols=T, fontsize_row=8,
         cluster_rows=T, cellwidth=9,
         cellheight=9, show_colnames=T,
         show_rownames = T, annotation = annotation, clustering_callback = callback, cutree_rows=2)

## Heatmap to explore differences in gene expression across the biogeochemical gradient - surface samples only (<100 m)
library(pheatmap)
library(genefilter)
library(RColorBrewer)
a<-read.csv('merged.csv', stringsAsFactors=F)
rownames(a)<-a$X
a<-a[,-1]
d<-log2(a+1)
rv<- rowVars(a)
idx<- order(-rv)[1:50]
annotation <- data.frame(Var1 = factor(1:15, labels = c('1')))
rownames(annotation)<-colnames(d)
annotation$Var1<-rownames(annotation)
fix(annotation)
colnames(annotation)<-c('Depth',"Station")
myColor <- colorRampPalette(brewer.pal(9, "YlGnBu"))(100)
pheatmap(d[idx,], cluster_cols=T, fontsize_row=8, cluster_rows=T, annotation = annotation, cellheight=9,cellwidth=9, show_rownames = T, color = myColor, cutree_rows=2, cutree_cols=2, scale="row")
