even_indexes<-seq(2,42,2)
fpkm<-fpkms[,even_indexes]
rownames(fpkm)<- c(fpkms$gene_id)
head(fpkm[1:10,1:10])
colnames(fpkm)<-c((samples$ID_seq))
aa<-samples$ID_seq
aa[order(aa),]
samples<-samples[order(samples$ID_seq),]

write.table(fpkm, file="fpkm.txt", sep=" ",quote = F, col.names = T, row.names = T)   
############
library(VennDiagram) 
library(dplyr)
library(WGCNA)
library(DESeq2)
loga<-function(x){log2(x+1)}
fpkm <- read.csv("C:/Users/carlo/Desktop/cosasdoctorao/rstudio/tutorial wgcna/fpkm.txt", sep="")
dir<-"C:/Users/carlo/Desktop/cosasdoctorao/rstudio/countgenes/"
dir
samples<- read.table(file.path(dir, "datos2.txt"), header = TRUE)
samples<-samples %>% 
  filter(time == "post")
fpkm<- fpkm %>%
  select_if(colnames(fpkm) %in% samples$ID_seq)
fpkm<-loga(fpkm)
fpkm <-fpkm[rowSums(fpkm)>1,]
head(fpkm[1:10,1:10])
femData<-fpkm
options(stringsAsFactors = FALSE)
head(femData)
dim(femData);
names(femData)

datExpr0 = as.data.frame(t(femData[, ]));####si es el originas 1:8, para los mios 1
head(datExpr0[1:10,1:10])
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}
sampleTree = hclust(dist(datExpr0), method = "average")  
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)
abline(h = 160, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 260, minSize = 1)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
namesSamples <- row.names(datExpr)
namesSamples
trais <- read.delim2("C:/Users/carlo/Desktop/cosasdoctorao/rstudio/tutorial wgcna/trais.txt")
trat<-merge(x=trais, y=samples ,by.x="nombre", by.y="ID_CNAG", all.y = F)
trat
traitData<-cbind.data.frame(trat[,6],trat[,2:5])
dim(traitData)
names(traitData)
head(traitData)
trassa<-factorizeNonNumericColumns(traitData$clase)
trass<-as.numeric(trassa$data)
traitData$data<- trass
traitData = traitData[, -c(2)];
head(traitData)
allTraits<-traitData
head(allTraits[1:10,])
dim(allTraits)
names(allTraits)
# Form a data frame analogous to expression data that will hold the clinical traits.
#datExpr<-datExpr0
femaleSamples = rownames(datExpr);
traitRows = match(femaleSamples, allTraits$`trat[, 6]`);
datTraits = allTraits[traitRows, -1];
datTraits<-as.data.frame(datTraits)
rownames(datTraits) = allTraits[traitRows, 1];
head(datTraits)
datTraits <- datTraits[, -c(2)]
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
datTraits
traitColors = numbers2colors(datTraits, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
#################################################
##############construccion network###############
#################################################
powers = c(c(1:10), seq(from = 12, to=30, by=2))
# Call the network topology analysis function1
##guardo esto y lo hago en el servidor pork no tengo potencia
save(datExpr, datTraits, file = "wgcnafpkm.RData")
lnames = load(file = "wgcnafpkm.RData")
enableWGCNAThreads(8)
sft = pickSoftThreshold(datExpr,RsquaredCut=0.85, powerVector = powers, blockSize = 15000, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
sft
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
net = blockwiseModules(datExpr, power = 22,
                       TOMType = "unsigned", minModuleSize = 30,
                       maxBlockSize = 15000,
                       reassignThreshold = 0, mergeCutHeight = 0.1,
                       numericLabels = F, pamRespectsDendro = F,
                       saveTOMs = F,
                       saveTOMFileBase = "prueba18000",
                       verbose = 3)
lnames = load(file = "lanet.RData")
# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms;
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits);
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = F,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = F,
               textAdj = c(0.1, 0.1),
               cex.text = 0.6,
               cex.lab = 0.5,
               cex.lab.x = 0.8,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
# Define variable weight containing the weight column of datTrait
weight = as.data.frame(datTraits$porcgrasa);
names(weight) = "weight"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="")
module = "lightgreen"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for fat %",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
names(datExpr)
names(datExpr)[moduleColors=="lightgreen"]
annot = read.delim("C:/Users/carlo/Desktop/cosasdoctorao/rstudio/glist-ram1.0r101.txt", header=FALSE);
dim(annot)
names(annot)
probes = names(datExpr)
probes2annot = match(probes, annot$V4)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.
# Create the starting data frame
geneInfo0 = data.frame(substanceBXH = probes,
                       geneSymbol = annot$V5[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, weight, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.weight));
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = "geneInfopruebaconclase.csv")
head(geneInfo)

probes = names(datExpr)
head(probes)
probes2annot = match(probes, annot$V4)
# Get the corresponding Locuis Link IDs
allLLIDs = annot$V5[probes2annot];
head(allLLIDs)
# $ Choose interesting modules
intModules = c("lightgreen", "white", "salmon")
for (module in intModules)
{
  # Select module probes
  modGenes = (moduleColors==module)
  # Get their entrez ID codes
  modLLIDs = allLLIDs[modGenes];
  # Write them into a file
  fileName = paste("LocusLinkIDs-", module, ".txt", sep="");
  write.table(as.data.frame(modLLIDs), file = fileName,
              row.names = FALSE, col.names = FALSE)
}
# As background in the enrichment analysis, we will use all probes in the analysis.
fileName = paste("LocusLinkIDs-all.txt", sep="");
write.table(as.data.frame(allLLIDs), file = fileName,
            row.names = FALSE, col.names = FALSE)
head(geneModuleMembership)
getgenes<-function(resul){
  genspp<-as.data.frame(row.names(resul))
  colnames(genspp)<-c("V1")
  glist <- read.delim("C:/Users/carlo/Desktop/cosasdoctorao/rstudio/glist-ram1.0r101.txt", header=FALSE)
  biem<-merge(x = genspp, y = glist, by.x = "V1", by.y = "V4", all.y =F)
  return(as.data.frame(biem$V5)) }
genes<-getgenes(geneTraitSignificance)
junto<-cbind.data.frame(genes$`biem$V5`,moduleColors)
head(junto)
genesaroa <- read.table("C:/Users/carlo/Desktop/cosasdoctorao/rstudio/genesaroa.txt", quote="\"", comment.char="")
X<-split(junto, junto$moduleColors)
write.table(junto, file = "genycolor.txt",
            row.names = FALSE, col.names = FALSE,quote = F)
for (x in colores){
  modulo<-genycolor %>% filter(genycolor$V2 == x )
  
  venn.plot <- venn.diagram(
    x = list(
      modulo = modulo$V1,
      paperAroa = genesaroa$V1
    ),
    filename = paste0(x,"modul.jpg"),
    col = "black",
    fill = c("goldenrod1", "darkorange1"),
    alpha = 0.50,
    cat.col = c("goldenrod1", "darkorange1"),
    cat.cex = 1.5,
    cat.fontface = "bold",
    margin = 0.05
  )
  
}


colores<-unique(junto$moduleColors)
modulo<-genycolor %>% filter(genycolor$V2 == x )
coloresaroa<- merge.data.frame(genesaroa$V1, genycolor$V1)
coloresaroa<-merge(x=genesaroa, y=genycolor ,by.x= "V1", by.y="V1", all.y = F)
coloresaroa<-cbind.data.frame(coloresaroa$V1,coloresaroa$V2.y)                   
head(coloresaroa)
write.table(coloresaroa, file = "coloresaroa.txt",
            row.names = FALSE, col.names = FALSE,quote = F)

a <- table(coloresaroa$`coloresaroa$V2.y`)
as<-as.data.frame(a)
write.table(as, file = "tablecoloresaroa.txt",
            row.names = FALSE, col.names = FALSE,quote = F)
######
numeros<-as.data.frame(table(junto$moduleColors))
numeros
colnames(numeros)<-c("")
tablecoloresaroa <- read.table("C:/Users/carlo/Desktop/cosasdoctorao/rstudio/tutorial wgcna/tablecoloresaroa.txt", quote="\"", comment.char="")
tablecoloresaroa
totalvsaroa<-merge(numeros,tablecoloresaroa, by.x = "Var1",by.y = "V1", all.x = T)

totalvsaroa
totalvsaroa[is.na(totalvsaroa)] <- 0
genes<-colSums(totalvsaroa[,2:3])
genes
avc<-as.data.frame(t(as.matrix(totalvsaroa)))
colnames(avc)<-avc[1,]
avc<- avc[-1,]
avc
install.packages("varhandle")
library(varhandle)
chisq.test(avc)
chisq.test(as.matrix(totalvsaroa[,-1]))
col
chisq.test(avc[,-1])
cva<-as.data.frame(t(as.matrix(avc)))
cva
cva<-unfactor(cva)
cva$Freq<-(as.integer(cva$Freq)) 
cva$esperado<-cva$Freq/14382*154
cva$V2<-as.numeric(cva$V2)
cva$chi<-(cva$V2-cva$esperado)^2/cva$esperado
cva
str(cva)
brown<-cva[3,]
suma<-sum(cva$chi)
suma
pchisq(suma, df=29, ncp = 0, lower.tail = FALSE, log.p = FALSE)
pchisq(1.469583e+01, df=1, ncp = 0, lower.tail = TRUE, log.p = FALSE)
cvaa<-cva[-3,]
cvaa
cvaa2<-as.data.frame(t(colSums(cvaa)))
cvaa2
cvaa2$esperado<-cvaa2$Freq/14382*154
cvaa2$chi<-(cvaa2$V2-cvaa2$esperado)^2/cvaa2$esperado
junt<-rbind.data.frame(cvaa2,brown)
junt
suma2<-sum(junt$chi)
suma2
cva$pvalor<-pchisq(cva$chi, df=1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
cva$padju<-p.adjust(cva$pvalor, method = "fdr", n = length(cva$pvalor))
pchisq(suma2, df=1, ncp = 0, lower.tail = FALSE, log.p = FALSE)
chicuadrao<-function(resul){
  
  genspp<-as.data.frame(row.names(resul))
  colnames(genspp)<-c("V1")
  glist <- read.delim("C:/Users/carlo/Desktop/cosasdoctorao/rstudio/glist-ram1.0r101.txt", header=FALSE)
  biem<-merge(x = genspp, y = glist, by.x = "V1", by.y = "V4", all.y =F)
  return(as.data.frame(biem$V5)) }
cva
sum(cva$chi)
junto
marron<-junto %>% filter( moduleColors == "brown")
marron<-(unique(marron[,1]))
#Webgstalt


library(WebGestaltR)
BiocManager::install("WebGestaltR")
Genesmarron<-WebGestaltR(enrichMethod="ORA", organism="hsapiens", 
                                      enrichDatabase="geneontology_Biological_Process",enrichDatabaseFile=NULL, 
                                      enrichDatabaseType=NULL,enrichDatabaseDescriptionFile=NULL,interestGeneFile=NULL, 
                                      interestGene=marron,interestGeneType="genesymbol",collapseMethod="mean",referenceGeneFile=NULL,
                                      referenceGene=NULL,referenceGeneType="genesymbol",referenceSet="genome", minNum=10, maxNum=500,
                                      fdrMethod="BH",sigMethod="fdr",fdrThr=0.05,topThr=10,dNum=20,perNum=1000,
                                      lNum=20,is.output=TRUE,outputDirectory=getwd(),projectName="upregularedAssaf_BP",keepGSEAFolder=FALSE,
                                      methodType="R",dagColor="binary",hostName="http://www.webgestalt.org/")
category<-rep ("Biological_Process", times= nrow(Genesmarron))
Genesmarron_BP_2<-cbind(Genesmarron,category)
Genesmarron_BP_2

darkgreen<-junto %>% filter( moduleColors == "darkgreen")
darkgreen<-(unique(darkgreen[,1]))
genesverde<-WebGestaltR(enrichMethod="ORA", organism="hsapiens", 
                         enrichDatabase="geneontology_Biological_Process",enrichDatabaseFile=NULL, 
                         enrichDatabaseType=NULL,enrichDatabaseDescriptionFile=NULL,interestGeneFile=NULL, 
                         interestGene=darkgreen,interestGeneType="genesymbol",collapseMethod="mean",referenceGeneFile=NULL,
                         referenceGene=NULL,referenceGeneType="genesymbol",referenceSet="genome", minNum=10, maxNum=500,
                         fdrMethod="BH",sigMethod="fdr",fdrThr=0.05,topThr=10,dNum=20,perNum=1000,
                         lNum=20,is.output=TRUE,outputDirectory=getwd(),projectName="darkgreen",keepGSEAFolder=FALSE,
                         methodType="R",dagColor="binary",hostName="http://www.webgestalt.org/")
lightcyan<-junto %>% filter( moduleColors == "lightcyan")
lightcyan<-(unique(lightcyan[,1]))
genescyan<-WebGestaltR(enrichMethod="ORA", organism="hsapiens", 
                        enrichDatabase="geneontology_Biological_Process",enrichDatabaseFile=NULL, 
                        enrichDatabaseType=NULL,enrichDatabaseDescriptionFile=NULL,interestGeneFile=NULL, 
                        interestGene=lightcyan,interestGeneType="genesymbol",collapseMethod="mean",referenceGeneFile=NULL,
                        referenceGene=NULL,referenceGeneType="genesymbol",referenceSet="genome", minNum=10, maxNum=500,
                        fdrMethod="BH",sigMethod="fdr",fdrThr=0.05,topThr=10,dNum=20,perNum=1000,
                        lNum=20,is.output=TRUE,outputDirectory=getwd(),projectName="lightcyan",keepGSEAFolder=FALSE,
                        methodType="R",dagColor="binary",hostName="http://www.webgestalt.org/")
orange<-junto %>% filter( moduleColors == "orange")
orange<-(unique(orange[,1]))
orange
genesorange<-WebGestaltR(enrichMethod="ORA", organism="hsapiens", 
                       enrichDatabase="geneontology_Biological_Process",enrichDatabaseFile=NULL, 
                       enrichDatabaseType=NULL,enrichDatabaseDescriptionFile=NULL,interestGeneFile=NULL, 
                       interestGene=orange,interestGeneType="genesymbol",collapseMethod="mean",referenceGeneFile=NULL,
                       referenceGene=NULL,referenceGeneType="genesymbol",referenceSet="genome", minNum=5, maxNum=500,
                       fdrMethod="BH",sigMethod="fdr",fdrThr=0.9,topThr=10,dNum=20,perNum=1000,
                       lNum=20,is.output=TRUE,outputDirectory=getwd(),projectName="orange",keepGSEAFolder=FALSE,
                       methodType="R",dagColor="binary",hostName="http://www.webgestalt.org/")
genesorange

genemarron_MF<-WebGestaltR(enrichMethod="ORA", organism="hsapiens", 
                                       enrichDatabase="geneontology_Molecular_Function",enrichDatabaseFile=NULL, 
                                       enrichDatabaseType=NULL,enrichDatabaseDescriptionFile=NULL,interestGeneFile=NULL, 
                                       interestGene=marron,interestGeneType="genesymbol",collapseMethod="mean",referenceGeneFile=NULL,
                                       referenceGene=NULL,referenceGeneType="genesymbol",referenceSet="genome", minNum=10, maxNum=500,
                                       fdrMethod="BH",sigMethod="fdr",fdrThr=0.05,topThr=10,dNum=20,perNum=1000,
                                       lNum=20,is.output=TRUE,outputDirectory=getwd(),projectName="marron_MF",keepGSEAFolder=FALSE,
                                       methodType="R",dagColor="binary",hostName="http://www.webgestalt.org/")
category<-rep ("Molecular_F", times= nrow(genemarron_MF))
genemarron_MF_2<-cbind(genemarron_MF,category)
genesmarronbp1<-Genesmarron_BP_2 %>% select(1,2,4,9,11)
genesdarkgreen1<-genesverde %>% select(1,2,4,9,11)
genescyan1<-genescyan %>% select(1,2,4,9,11)
lightgreen<-junto %>% filter( moduleColors == "lightgreen")
lightgreen<-(unique(lightgreen[,1]))
lightgreen
geneslightgreen<-WebGestaltR(enrichMethod="ORA", organism="hsapiens", 
                         enrichDatabase="geneontology_Biological_Process",enrichDatabaseFile=NULL, 
                         enrichDatabaseType=NULL,enrichDatabaseDescriptionFile=NULL,interestGeneFile=NULL, 
                         interestGene=lightgreen,interestGeneType="genesymbol",collapseMethod="mean",referenceGeneFile=NULL,
                         referenceGene=NULL,referenceGeneType="genesymbol",referenceSet="genome", minNum=5, maxNum=500,
                         fdrMethod="BH",sigMethod="fdr",fdrThr=0.05,topThr=10,dNum=20,perNum=1000,
                         lNum=20,is.output=TRUE,outputDirectory=getwd(),projectName="lightgreen",keepGSEAFolder=FALSE,
                         methodType="R",dagColor="binary",hostName="http://www.webgestalt.org/")
geneslightgreen1<-geneslightgreen %>% select(1,2,4,9,11)


white<-junto %>% filter( moduleColors == "white")
white<-(unique(white[,1]))
white
geneswhite<-WebGestaltR(enrichMethod="ORA", organism="hsapiens", 
                             enrichDatabase="geneontology_Biological_Process",enrichDatabaseFile=NULL, 
                             enrichDatabaseType=NULL,enrichDatabaseDescriptionFile=NULL,interestGeneFile=NULL, 
                             interestGene=white,interestGeneType="genesymbol",collapseMethod="mean",referenceGeneFile=NULL,
                             referenceGene=NULL,referenceGeneType="genesymbol",referenceSet="genome", minNum=4, maxNum=500,
                             fdrMethod="BH",sigMethod="fdr",fdrThr=0.05,topThr=10,dNum=20,perNum=1000,
                             lNum=20,is.output=TRUE,outputDirectory=getwd(),projectName="white",keepGSEAFolder=FALSE,
                             methodType="R",dagColor="binary",hostName="http://www.webgestalt.org/")
geneswhite1<-geneswhite %>% select(1,2,4,9,11)

red<-junto %>% filter( moduleColors == "red")
red<-(unique(red[,1]))
red
genesred<-WebGestaltR(enrichMethod="ORA", organism="hsapiens", 
                        enrichDatabase="geneontology_Biological_Process",enrichDatabaseFile=NULL, 
                        enrichDatabaseType=NULL,enrichDatabaseDescriptionFile=NULL,interestGeneFile=NULL, 
                        interestGene=red,interestGeneType="genesymbol",collapseMethod="mean",referenceGeneFile=NULL,
                        referenceGene=NULL,referenceGeneType="genesymbol",referenceSet="genome", minNum=4, maxNum=500,
                        fdrMethod="BH",sigMethod="fdr",fdrThr=0.05,topThr=10,dNum=20,perNum=1000,
                        lNum=20,is.output=TRUE,outputDirectory=getwd(),projectName="red",keepGSEAFolder=FALSE,
                        methodType="R",dagColor="binary",hostName="http://www.webgestalt.org/")
genesred1<-genesred %>% select(1,2,4,9,11)

salmon<-junto %>% filter( moduleColors == "salmon")
salmon<-(unique(salmon[,1]))
salmon
genessalmon<-WebGestaltR(enrichMethod="ORA", organism="hsapiens", 
                        enrichDatabase="geneontology_Biological_Process",enrichDatabaseFile=NULL, 
                        enrichDatabaseType=NULL,enrichDatabaseDescriptionFile=NULL,interestGeneFile=NULL, 
                        interestGene=salmon,interestGeneType="genesymbol",collapseMethod="mean",referenceGeneFile=NULL,
                        referenceGene=NULL,referenceGeneType="genesymbol",referenceSet="genome", minNum=4, maxNum=500,
                        fdrMethod="BH",sigMethod="fdr",fdrThr=0.05,topThr=10,dNum=20,perNum=1000,
                        lNum=20,is.output=TRUE,outputDirectory=getwd(),projectName="salmon",keepGSEAFolder=FALSE,
                        methodType="R",dagColor="binary",hostName="http://www.webgestalt.org/")
genessalmon1<-genessalmon %>% select(1,2,4,9,11)
