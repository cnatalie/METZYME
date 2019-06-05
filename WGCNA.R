#WGCNA pipeline adapted from Sarah Davies (http://sites.bu.edu/davieslab/)


#Data input and cleaning
library(WGCNA)
options(stringsAsFactors=FALSE)
allowWGCNAThreads()

#Read in normalized, gene-annotated transcript counts
a<-read.csv('TPM_TRANSCRIPTS_Dino.lpi0.8_KO_postsum.csv')
a<-a[,-43] #Remove outlier sample from 1900m depth
datExpr0 = as.data.frame(t(a[, -c(1)]))
names(datExpr0) = a$X
rownames(datExpr0) = names(a)[-c(1)]
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK #if TRUE, no outlier genes

#log2 transform
datExpr0<-log2(datExpr0+1)

#Read in trait data
traitData<-read.csv('CCA_meta_WCGNA_short.csv', stringsAsFactors=F)
rownames(traitData)=rownames(datExpr0)
 traitData$Sample= NULL 
datTraits=traitData
table(rownames(datTraits)==rownames(datExpr0)) #should return TRUE if datasets align correctly, otherwise your names are out of order
datTraits<-datTraits[,-1]

#Dendrogram and trait heat map showing outliers
A=adjacency(t(datExpr0),type="signed") #calculate network adjacency
k=as.numeric(apply(A,2,sum))-1 #Sum columns, not sure why (-1) but doesn't seem to make a difference
Z.k=scale(k)
thresholdZ.k=-2.5 # often -2.5
outlierColor=ifelse(Z.k<thresholdZ.k,"red","black") #Outliers in red
library(flashClust)
sampleTree = flashClust(as.dist(1-A), method = "average")
#datTraits<-datTraits[,-1] #Get rid of text ID column
traitColors=data.frame(numbers2colors(datTraits,signed=FALSE))
dimnames(traitColors)[[2]]=paste(names(datTraits))
datColors=data.frame(outlierC=outlierColor,traitColors)
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree,groupLabels=names(datColors), colors=datColors,main="Sample dendrogram and trait heatmap")

#Remove outlying samples from expression and trait data
remove.samples= Z.k<thresholdZ.k | is.na(Z.k)
datExpr=datExpr0[!remove.samples,]
datTraits=datTraits[!remove.samples,]
datExpr0<-datExpr
#save(datExpr0, datTraits, file="Samples_Traits_ALL.RData")

#Network construction and module detection
options(stringsAsFactors = FALSE)
#Choose a set of soft-thresholding powers
powers = c(seq(1,14,by=2), seq(15,30, by=0.5));
#Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, networkType="signed", verbose = 2) #want smallest value, closest to 0.9
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
#Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
    main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
    labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red") #This line corresponds to using an R^2 cut-off of h
#Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
    xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
    main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
softPower=14 #smallest value to plateau at ~0.85 
adjacency=adjacency(datExpr0, power=softPower,type="signed") 
#Translate the adjacency into topological overlap matrix and calculate the corresponding dissimilarity:
TOM= TOMsimilarity(adjacency,TOMType = "signed")
dissTOM= 1-TOM

#Each leaf corresponds to a gene, branches grouping together densely are interconnected, highly co-expressed genes
geneTree= flashClust(as.dist(dissTOM), method="average")
minModuleSize=50
dynamicMods= cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2, pamRespectsDendro= FALSE, minClusterSize= minModuleSize)
table(dynamicMods)
dynamicColors= labels2colors(dynamicMods)
sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang= 0.05, main= "Gene dendrogram and module colors")

#Merge modules whose expression profiles are very similar
MEList= moduleEigengenes(datExpr0, colors= dynamicColors,softPower = 14)
MEs= MEList$eigengenes
#Calculate dissimilarity of module eigenegenes
MEDiss= 1-cor(MEs)
#Cluster module eigengenes
METree= flashClust(as.dist(MEDiss), method= "average")
sizeGrWindow(7,6)
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")
MEDissThres= 0.3
abline(h=MEDissThres, col="red")
merge= mergeCloseModules(datExpr0, dynamicColors, cutHeight= MEDissThres, verbose =3)
mergedColors= merge$colors
mergedMEs= merge$newMEs
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)
moduleColors= mergedColors
colorOrder= c("grey", standardColors(50))
moduleLabels= match(moduleColors, colorOrder)-1
MEs=mergedMEs

#Relating modules to traits and finding important genes
options(stringsAsFactors = FALSE);
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)
table(moduleColors)

#Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

#Represent module trait correlations as a heatmap
sizeGrWindow(10,6)
#Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                        signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
             xLabels = names(datTraits),
             yLabels = names(MEs),
             ySymbols = names(MEs),
             colorLabels = FALSE,
             colors = blueWhiteRed(50),
             textMatrix = textMatrix,
             setStdMargins = FALSE,
             cex.text = 0.8,
             zlim = c(-1,1),
             main = paste("Module-trait relationships"))
             
#Gene relationship to trait and important modules
weight = as.data.frame(datTraits$Fe); #change to your trait name
names(weight) = "Fe"
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr0, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");


#Heatmap of module expression with bar plot of eigengene
sizeGrWindow(8,7);
which.module="yellow" #pick module of interest
ME=MEs[, paste("ME",which.module, sep="")]
genes=datExpr0[,moduleColors==which.module ] 
par(mfrow=c(2,1), mar=c(0.3, 5.5, 5, 2))
plotMat(t(scale(genes) ),nrgcols=30,rlabels=F, clabels=rownames(genes), rcols=which.module,)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="sample")

#Output ME by sample
meyellow<-ME #get ME for yellow module above
meout<-data.frame(cbind(rownames(datExpr),meyellow,memidnightblue))

write.csv(meout,"MEbySample_MidnightblueBlueMods.csv",quote=F,row.names=F)

#Connect Gene ID to functions, save file with p values. Outputs gene significance results, module colors and KOs
annot = read.csv(file = "kodef.csv"); #KO definition to KO
dim(annot)
names(annot)
probes = names(datExpr0)
probes2annot = match(probes, annot$KO)

geneInfo0 = data.frame(def = probes,
                       desc = annot$def[probes2annot],
                       class = annot$name[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
write.csv(geneInfo0, 'TPM_WGCNA_021319_allsamples.csv')

#Note: The blue color module was changed to the color black, and turquoise to white for clarity during visualizations

#KEGG enrichment
library(clusterProfiler)
a<-read.csv('blue.csv',header=F) #Read in KOs in module of interest
b <- a[,1] #Convert to vector
x <- enrichKEGG(b, organism='ko', keyType='kegg', universe = bg)
#write.csv(x, 'blue_enriched.csv')
bg<-read.csv('bg.csv',header=F) 
bg<-bg[,1]#all dino KOs in modules
aa<-read.csv('turquoise.csv',header=F) #Read in KOs in module of interest
bb <- aa[,1]
xx <- enrichKEGG(bb, organism='ko', keyType='kegg', universe = bg)
#write.csv(xx, 'turquoise_enriched.csv')
barplot(xx, colorBy = "p.adjust", showCategory = 50)

a<-read.csv('turquoise_enriched.csv')
library(ggplot2)
ggplot(data=a, aes(x= reorder (Description, GeneRatio), y=GeneRatio)) + xlab("KEGG Pathway Description") + ylab("Ratio in Module") + theme_minimal() +
  geom_bar(colour="black", fill="#0437e0", width=.8, stat="identity") + coord_flip()
  guides(fill=FALSE)
