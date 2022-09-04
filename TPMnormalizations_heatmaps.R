###Calculating normalized transcripts as TPM, caculated following Micheal Love's blog post: https://support.bioconductor.org/p/91218/

#Transcripts normalization - Read in dinoflagellate raw transcripts, annotations and open reading frame length in base pairs and TPM normalize
a<-read.csv('annotation_all.filtered.grps.go_TRANSCRIPTS_0.8lpi_Dino_only_sansAnnotations.csv') #Subset dinoflagellate contigs and raw counts with LPI score > 0.8
b<-read.csv('orflength.csv') #ORF lengths from PhyloDB output
c<-merge(a, b, by='orf_id')
rownames(c)<-c$orf_id
x <- c[3:44] / c$orf_length #TPM normalize: https://support.bioconductor.org/p/91218/
tpm.mat <- t( t(x) * 1e6 / colSums(x) )
rownames(tpm.mat)<-c$orf_id
write.csv(tpm.mat, 'TPM_TRANSCRIPTS_Dino_only_orf.csv')
b<-read.csv('TPM_TRANSCRIPTS_Dino_only_orf.csv')
colnames(b)[colnames(b) == 'X'] <- 'orf_id' #rownames = orf_id
#Read in file with all annotations - Supplemental Dataset 2 (the dataset is split in 2 parts due to the large size and upload limits - merge to one big file!)
a<-read.csv('dataset2.csv')
c<-merge(a, b, by='orf_id') #Should be left with dinoflagellate TPM counts and functional annotations
write.csv(c, 'TPM_TRANSCRIPTS_Dino.lpi0.8_only_annotations_orf.csv')

#Protein normalization - NSAF normalization using amino acid residues as length and multiplying by scaler of 1000 
a<-read.csv('exclusive_counts_annotated_withscaffoldgroups_dinolpi0.8_NSAFpre.csv')
rownames(a)<-a$X
x <- a[5:44] / a$orf_length_aa
xx <- t( t(x) * 1e3 / colSums(x) )
y<-colSums(xx) #Should all equal 1,000
write.csv(xx, 'exclusive_counts_annotations_dino_lpi0.8_post_NSAF_orf.csv')
a<-read.csv('exclusive_counts_annotations_dino_lpi0.8_post_NSAF_orf.csv')
b<-read.csv('exclusive_counts_annotations_dino_lpi0.8_pre.csv') #import original annotation dataframe
c<-merge(a, b, by="X") #merge new counts with annotations

library(pheatmap)
library(genefilter)
library(viridis)
library(RColorBrewer)
library(stringr)

### TPM-normalized KEGG heatmap (three different taxa). Sum together same KO annotations within each taxonomic group, and then merge into one dataframe
a<-read.csv('annotation_all.filtered.grps.go.lpi_0.8_dino_diatom_hapto_TPM_KOpre.csv')
dia<-c('Diatom')
b<-a[a$GROUP %in% dia,]
head(b)
library(caroline)
diatoms<-groupBy(b, by='KO',clmns=(4:44),aggregation='sum')
colnames(diatoms)<-colnames(b[4:44])
colnames(diatoms) <- paste("Diatoms", colnames(diatoms), sep = "_")
diatoms$KO<-rownames(diatoms)
dino<-c('Dinophyta')
b<-a[a$GROUP %in% dino,]
dino<-groupBy(b, by='KO',clmns=(4:44),aggregation='sum')
colnames(dino)<-colnames(b[4:44])
colnames(dino) <- paste("Dinoflagellates", colnames(dino), sep = "_")
dino$KO<-rownames(dino)
dia<-c('Haptophyta')
b<-a[a$GROUP %in% dia,]
head(b)
hapto<-groupBy(b, by='KO',clmns=(4:44),aggregation='sum')
colnames(hapto)<-colnames(b[4:44])
colnames(hapto) <- paste("Haptophytes", colnames(hapto), sep = "_")
hapto$KO<-rownames(hapto)
c<-merge(diatoms, dino, by='KO')
d<-merge(c, hapto, by='KO')
rownames(d)<-d$KO
d<-d[,-1]
library(viridis)
library(pheatmap)
paletteLength=100
myColor <- rev(viridis_pal(option = "B")(paletteLength))
annotation <- data.frame(Var1 = factor(1:123, labels = c('1')))
rownames(annotation)<-colnames(d)
annotation$Var1<-rownames(annotation)
annotation$Depth<-str_extract(annotation$Var1, "(?<=\\_)\\d+")
annotation$Depth<- as.numeric(as.character(annotation$Depth))
annotation$Station<-str_extract(annotation$Var1, "(?=\\_?)\\d+")
annotation$Group<-str_extract(annotation$Var1, "\\w+(?=_X)")
annotation<-annotation[,-1]
rownames(d)<-d$KO
z<-log2(d+1)
pheatmap(z, color=myColor, cluster_cols=T, fontsize_row=8, cluster_rows=T, show_colnames=F, show_rownames = F, annotation = annotation)

#PFams heatmaps with proteins
setwd("~/METZYME/metaproteome 3um")
a<-read.csv('exclusive_counts_annotations_dino_lpi0.8_post_NSAF_orf_062520_Pfampost_annotated.csv')
rownames(a)<-a$PFams
b<-a[,-1]
b<-b[,-40]
head(b)
c<-log2(b+1)
library(pheatmap)
library(genefilter)
library(viridis)
library(RColorBrewer)
library(stringr)
paletteLength <- 50
myColor <- colorRampPalette(brewer.pal(9, "YlGnBu"))(100)
#NOTE on rowVars: Because of the mean-variance relationship in RNA count data, more highly expressed genes will have higher variance. Log-transforming will have the reverse effect, and log-transformation of highly expressed genes will have lower variance. 
#Alternative options include variance stabilization approaches, or using fold-change cutoffs to evaluate differentially expressed genes 
#Popular RNASeq differential expression tools were not used here as they have assumptions about most genes not being differentially expressed, which is likely violated in mixed natural communities that could be taxonomically and/or functionally changing with ocean region
#MA plots can be helpful to visualize fold-changes in expression alongside transcript abundance
#Very useful summary of the effects of log (or not log) transforming to calculate variance: https://github.com/friedue/Notes/blob/master/RNA_heteroskedasticity.md
rv <- apply(b, 1, var)
idx <- names(V[order(V, decreasing = T)][1:50])
nrow(b)
head(c)
annotation <- data.frame(Var1 = factor(1:39, labels = c('1')))
rownames(annotation)<-colnames(c)
annotation$Var1<-rownames(annotation)
annotation$Depth<-str_extract(annotation$Var1, "(?<=\\_)\\d+")
annotation$Station<-str_extract(annotation$Var1, "(?=\\_?)\\d+")
annotation<-annotation[,-1]

#To flip a dendrogram branch
callback = function(hc, mat){
         sv = svd(t(mat))$v[,1]
         dend = reorder(as.dendrogram(hc), wts = sv)
         as.hclust(dend)
}

pheatmap(c[idx,], color=myColor, 
         cluster_cols=T, fontsize_row=8, 
         cluster_rows=T, cellwidth=9, 
         cellheight=9, show_colnames=T, 
         show_rownames = T, annotation = annotation, clustering_callback = callback)


### transcripts KEGG heatmap
library(caroline)

a<-read.csv('TPM_TRANSCRIPTS_Dino.lpi0.8_only_annotations_orf_KOpost.csv') #TPM counts from the same KO IDs were summed together.
x<-read.tab('kodef.tab')
x$KO_def<-paste(x$KO, x$def, sep='_')
rownames(a)<-a$KO
head(a)
#a<-a[,-1]
z<-merge(a, x, by='KO')
rownames(z)<-z$KO_def
head(z)
z<-z[,-1]
z<-z[,-45]
z<-z[,-44]
z<-z[,-43]
z<-z[,-42] #Delete non-numeric columns
o<-a[1660:1663,] #carry over kegg genes without descriptions 
o<-o[,-44] #Delete non-numeric columns
o<-o[,-43]
o<-o[,-1]
x<-rbind(o,z)
z<-x
d<-log2(z+1)

#NOTE on rowVars: Because of the mean-variance relationship in RNA count data, more highly expressed genes will have higher variance. Log-transforming will have the reverse effect, and log-transformation of highly expressed genes will have lower variance. 
#Alternative options include variance stabilization approaches, or using fold-change cutoffs to evaluate differentially expressed genes 
#Popular RNASeq differential expression tools were not used here as they have assumptions about most genes not being differentially expressed, which is likely violated in mixed natural communities that could be taxonomically and/or functionally changing with ocean region
#MA plots can be helpful to visualize fold-changes in expression alongside transcript abundance
#Very useful summary of the effects of log (or not log) transforming to calculate variance: https://github.com/friedue/Notes/blob/master/RNA_heteroskedasticity.md

rv<- rowVars(z)
idx<- order(-rv)[1:50]
paletteLength <- 50
myColor <- rev(viridis_pal(option = "C")(paletteLength))
myColor <- colorRampPalette(brewer.pal(9, "YlGnBu"))(100)
annotation <- data.frame(Var1 = factor(1:41, labels = c('1')))
rownames(annotation)<-colnames(d)
annotation$Var1<-rownames(annotation)
annotation$Depth<-str_extract(annotation$Var1, "(?<=\\_)\\d+")
annotation$Station<-str_extract(annotation$Var1, "(?=\\_?)\\d+")
annotation<-annotation[,-1]
pheatmap(d[idx,], scale="row", color=myColor,
         cluster_cols=T, fontsize_row=8,
         cluster_rows=T, cellwidth=9,
         cellheight=9, show_colnames=T,
         show_rownames = T, annotation = annotation, clustering_callback = callback, cutree_rows=2)

## Heatmap to explore differences in gene expression across the biogeochemical gradient - surface samples only (<100 m). 
library(pheatmap)
library(genefilter)
library(RColorBrewer)
a<-read.csv('TPM_transcripts_Dino.lpi0.8_KOdef_top100m.csv', stringsAsFactors=F) #TPM counts from the same KO IDs were summed together.
rownames(a)<-a$X
a<-a[,-1]
d<-log2(a+1)
rv<- rowVars(a)
idx<- order(-rv)[1:45]
annotation <- data.frame(Var1 = factor(1:15, labels = c('1')))
rownames(annotation)<-colnames(d)
annotation$Var1<-rownames(annotation)
annotation$Depth<-str_extract(annotation$Var1, "(?<=\\_)\\d+")
annotation$Station<-str_extract(annotation$Var1, "(?=\\_?)\\d+")
annotation<-annotation[,-1]
pheatmap(d[idx,], cluster_cols=T, fontsize_row=8, cluster_rows=T, annotation = annotation, cellheight=9,cellwidth=9, show_rownames = T, color = myColor, cutree_rows=2, cutree_cols=2, scale="row")
