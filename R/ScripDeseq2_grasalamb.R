
setwd("~/grasa_prueba/")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")
BiocManager::install ("tximportData")
library("tximportData")

BiocManager::install ("tximport")
BiocManager::install ("readr")
library(tximport)
library(readr)

setwd("~/grasa_prueba/")
dir<-"~/grasa_prueba/"
dir
samples <- read.table(file.path(dir, "samples.txt"), header = TRUE)
files <- file.path(dir, "rsem", paste0(samples$ID, ".counts.genes.results"))
names(files) <-paste0(samples$ID)
file.exists(files)
all(file.exists(files))
txi.rsem<-tximport(files, type="rsem")
head(txi.rsem$counts)
names(txi.rsem)

BiocManager::install("pheatmap")
library(pheatmap)
jpeg(filename="Heatmap.jpg",width=7,height=7, units="in",res=500)
pheatmap(cor(log10(txi.rsem$counts+1)))
dev.off()
mm<-cor(log10(txi.rsem$counts+1))
which(mm == min(mm), arr.ind = TRUE)

##Analisis de expresion diferencial con DESeq2

#install DESeq2

BiocManager::install ("DESeq2")
library (DESeq2)

#Tabla para DESeq2
colData2<-data.frame(breed=c(rep("assaf",8),rep("churra",6)), age=c(rep("less25",8), rep("more25",3),rep("less25",3)),row.names = colnames(txi.rsem$counts))

#generar un DESeqDataSet
txi.rsem$length[txi.rsem$length == 0] <- 1
dds <- DESeqDataSetFromTximport(txi.rsem, colData2, design = ~breed+age)
dds<- DESeq (dds)
resultsNames(dds)

#to obtain regularized log-transformed values
rlog.DESeq.sumExp<- rlog(dds, blind = T)
rlog.norm.counts<- assay (rlog.DESeq.sumExp)
meanSdPlot (rlog.norm.counts)
dds <- dds[ rowSums(counts(dds)) > 1, ]
pheatmap(cor(log10((counts(dds)+1))))
#-------##Correlation and PCA ##

distance.m_rlog <- as.dist (1- cor ((rlog.norm.counts), method = "spearman"))
jpeg(filename="2-Dendrogram & PCA.jpg",width=20,height=20, units="in",res=500)
par(mfrow=c(2,1))
plot (hclust(distance.m_rlog), labels= colnames (rlog.norm.counts), main="rlog transformed read counts\ndistance: Pearson correlation")
pc <- prcomp (t(rlog.norm.counts))
scores = as.data.frame(pc$x)
plot(scores$PC1, scores$PC2, main= "PCA of seq.depth normalized\n and rlog-transdormed read counts", col= colData(dds) [,1])
dev.off()

#hcluster Y PCA de CO vs LU
par(mfrow=c(1,1))
distance.m_rlog <-as.dist (1- cor ((rlog.norm.counts), method = "pearson"))
jpeg(filename=" distance_LUvsCO.jpg",width=10,height=10, units="in",res=500)
plot (hclust(distance.m_rlog), labels= colnames (rlog.norm.counts), main="rlog transformed read counts\ndistance: Pearson correlation")
dev.off()
pc <- prcomp (t(rlog.norm.counts))
scores = as.data.frame(pc$x)
summary(pc)

install.packages ("ggplot2")
library(ggplot2)
scores<- data.frame (scores, Condition= c("assaf","assaf","assaf","assaf", "assaf","assaf","assaf","assaf","churra", "churra", "churra", "churra","churra", "churra"))
x<-ggplot(data = scores, aes(x = PC1, y = PC2, label = rownames(scores))) +
  geom_hline(yintercept = 0, colour = "gray65") +
  geom_vline(xintercept = 0, colour = "gray65") +
  geom_text(aes(label = rownames(scores), colour=Condition, fontface="bold")) +
  scale_colour_manual(values=c("#000000", "#999999")) +
  ggtitle("PCA of seq.depth normalized\n and rlog-transformed read counts") +
  theme(plot.title = element_text(face="bold",hjust = 0.5)) +
  theme(panel.background = element_rect(fill = 'white'))
tiff(filename=" PCA2_byn.tiff",width=7,height=7, units="in",res=500)
x
dev.off()

cor(rlog.norm.counts,method="pearson")

plot(scores$PC1, scores$PC3, main= "PCA of seq.depth normalized\n and rlog-transdormed read counts", col= colData(dds) [,1])

pcFOvsCO <- prcomp (t(rlog.norm.counts))
scoresFOvsCO = as.data.frame(pcFOvsCO$x)
plot(scoresFOvsCO$PC1, scoresFOvsCO$PC2, main= "PCA of seq.depth normalized\n and rlog-transdormed read counts", col= colData(dds) [,1])

rlog.DESeq.sumExp<- rlog(dds)
jpeg(filename="3- PCA.jpg",width=7,height=7, units="in",res=500)
P <- plotPCA(rlog.DESeq.sumExp, intgroup= "breed")
P<- P+ theme_bw() +ggtitle ("Rlog transformed counts")
print(P)
dev.off()


##Differential expression analysis: differentDays
res<-results(dds, name="breed_churra_vs_assaf", alpha = 0.05)
head(res)
summary(res)
resOrdered<-res[order(res$padj),] ##ordering by p-value
na.omit(resOrdered) ##eliminating na
resOrdered<-subset(resOrdered, padj<0.05) ##extracting genes with padj<0.05
write.table(resOrdered, file="resOrdered.txt", sep='\t') ## writing the table
summary(resOrdered)

resOrdered2<-as.data.frame(resOrdered)
gene_id <- rownames(resOrdered2)
rownames(resOrdered2) <- NULL
resOrdered3 <- cbind(gene_id,resOrdered2)


gene_symbol_names<- read.table('~/trabajo2017/Grasa_cordero_paper/Grasa_r88_RSEM/quant_genes/Ovis_aries_v3.1-r88.corretalionEnsemble_Genesymbol_2', header=TRUE, sep= "\t")
resOrdered3_genesymbol<-merge(resOrdered3,gene_symbol_names,by=c("gene_id"))
write.table(resOrdered3_genesymbol, file="resOrdered_gene_symbol.txt", sep='\t') ## writing the table

resOrdered3_genesymbol_WebgStalt<-subset(resOrdered3_genesymbol,select=c("gene_symbol","log2FoldChange"))
colnames(resOrdered3_genesymbol_WebgStalt) <- c("genes", "scores")
resOrdered3_genesymbol_WebgStalt_ordered<-resOrdered3_genesymbol_WebgStalt[order(resOrdered3_genesymbol_WebgStalt$scores),]

resOrdered_menor0<- subset(resOrdered3_genesymbol_WebgStalt_ordered, scores<0)
resOrdered_menor0_gene_symbol<-resOrdered_menor0$genes
resOrdered_menor0_gene_symbol<-as.vector(resOrdered_menor0_gene_symbol)
write.table(resOrdered_menor0_gene_symbol, file="resOrdered_menor0_gene_symbol_WebGstalt.txt", sep='\t')

resOrdered_mayor0<- subset(resOrdered3_genesymbol_WebgStalt_ordered, scores>0)
resOrdered_mayor0_gene_symbol<-resOrdered_mayor0$genes
resOrdered_mayor0_gene_symbol<-as.vector(resOrdered_mayor0_gene_symbol)
write.table(resOrdered_mayor0_gene_symbol, file="resOrdered_mayor0_gene_symbol_WebGstalt.txt", sep='\t')

ADIPOQ <- plotCounts(dds, gene="ENSOARG00000020509", intgroup="breed", 
                returnData=TRUE)
library("ggplot2")
ggplot(ADIPOQ, aes(x=breed, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))

ADTRP <- plotCounts(dds, gene="ENSOARG00000014611", intgroup="breed", 
                     returnData=TRUE)
ggplot(ADTRP, aes(x=breed, y=count)) + 
  geom_point(position=position_jitter(w=0.1,h=0)) + 
  scale_y_log10(breaks=c(25,100,400))


#Webgstalt
install.packages("pkgmaker")
install.packages("rjson")
install.packages("data.table")
install.packages("PythonInR")
install.packages("parallel")
install.packages("doParallel")
install.packages("foreach")

biocLite ("WebGestaltR")
library("WebGestaltR")



outputDirectory<-getwd()

##Enrichment analysis downregulatedGenes


GenesupregulatedAssaf_BP<-WebGestaltR(enrichMethod="ORA", organism="hsapiens", 
            enrichDatabase="geneontology_Biological_Process",enrichDatabaseFile=NULL, 
            enrichDatabaseType=NULL,enrichDatabaseDescriptionFile=NULL,interestGeneFile=NULL, 
            interestGene=resOrdered_menor0_gene_symbol,interestGeneType="genesymbol",collapseMethod="mean",referenceGeneFile=NULL,
            referenceGene=NULL,referenceGeneType="genesymbol",referenceSet="genome", minNum=10, maxNum=500,
            fdrMethod="BH",sigMethod="fdr",fdrThr=0.05,topThr=10,dNum=20,perNum=1000,
            lNum=20,is.output=TRUE,outputDirectory=getwd(),projectName="upregularedAssaf_BP",keepGSEAFolder=FALSE,
            methodType="R",dagColor="binary",hostName="http://www.webgestalt.org/")
category<-rep ("Biological_Process", times= 18)
GenesupregulatedAssaf_BP_2<-cbind(GenesupregulatedAssaf_BP,category)

GenesupregulatedAssaf_MF<-WebGestaltR(enrichMethod="ORA", organism="hsapiens", 
                                      enrichDatabase="geneontology_Molecular_Function",enrichDatabaseFile=NULL, 
                                      enrichDatabaseType=NULL,enrichDatabaseDescriptionFile=NULL,interestGeneFile=NULL, 
                                      interestGene=resOrdered_menor0_gene_symbol,interestGeneType="genesymbol",collapseMethod="mean",referenceGeneFile=NULL,
                                      referenceGene=NULL,referenceGeneType="genesymbol",referenceSet="genome", minNum=10, maxNum=500,
                                      fdrMethod="BH",sigMethod="fdr",fdrThr=0.05,topThr=10,dNum=20,perNum=1000,
                                      lNum=20,is.output=TRUE,outputDirectory=getwd(),projectName="upregularedAssaf_MF",keepGSEAFolder=FALSE,
                                      methodType="R",dagColor="binary",hostName="http://www.webgestalt.org/")
category<-rep ("Molecular_Function", times= nrow(GenesupregulatedAssaf_MF))
GenesupregulatedAssaf_MF_2<-cbind(GenesupregulatedAssaf_MF,category)

GenesupregulatedAssaf_CC<-WebGestaltR(enrichMethod="ORA", organism="hsapiens", 
                                      enrichDatabase="geneontology_Cellular_Component",enrichDatabaseFile=NULL, 
                                      enrichDatabaseType=NULL,enrichDatabaseDescriptionFile=NULL,interestGeneFile=NULL, 
                                      interestGene=resOrdered_menor0_gene_symbol,interestGeneType="genesymbol",collapseMethod="mean",referenceGeneFile=NULL,
                                      referenceGene=NULL,referenceGeneType="genesymbol",referenceSet="genome", minNum=10, maxNum=500,
                                      fdrMethod="BH",sigMethod="fdr",fdrThr=0.05,topThr=10,dNum=20,perNum=1000,
                                      lNum=20,is.output=TRUE,outputDirectory=getwd(),projectName="upregularedAssaf_CC",keepGSEAFolder=FALSE,
                                      methodType="R",dagColor="binary",hostName="http://www.webgestalt.org/")
category<-rep ("Cellular_Component", times= nrow(GenesupregulatedAssaf_CC))
GenesupregulatedAssaf_CC_2<-cbind(GenesupregulatedAssaf_CC,category)

GO_GenesupregulatedAssaf_results<-rbind(GenesupregulatedAssaf_BP_2,GenesupregulatedAssaf_MF_2,GenesupregulatedAssaf_CC_2)
GO_GenesupregulatedAssaf_results<-subset(GO_GenesupregulatedAssaf_results, select=c("geneset","description","C","O","E","R","PValue","FDR","OverlapGene_UserID","category"))
#???write.table(GO_GenesupregulatedAssaf_results, "GO_GenesupregulatedAssaf_results.txt")

##Enrichment analysis upregulatedGenes

GenesupregulatedChurra_BP<-WebGestaltR(enrichMethod="ORA", organism="hsapiens", 
                                       enrichDatabase="geneontology_Biological_Process",enrichDatabaseFile=NULL, 
                                       enrichDatabaseType=NULL,enrichDatabaseDescriptionFile=NULL,interestGeneFile=NULL, 
                                       interestGene=resOrdered_mayor0_gene_symbol,interestGeneType="genesymbol",collapseMethod="mean",referenceGeneFile=NULL,
                                       referenceGene=NULL,referenceGeneType="genesymbol",referenceSet="genome", minNum=10, maxNum=500,
                                       fdrMethod="BH",sigMethod="fdr",fdrThr=0.05,topThr=10,dNum=20,perNum=1000,
                                       lNum=20,is.output=TRUE,outputDirectory=getwd(),projectName="upregularedChurra_BP",keepGSEAFolder=FALSE,
                                       methodType="R",dagColor="binary",hostName="http://www.webgestalt.org/")
category<-rep ("Biological_Process", times= nrow(GenesupregulatedChurra_BP))
GenesupregulatedChurra_BP_2<-cbind(GenesupregulatedChurra_BP,category)

GenesupregulatedChurra_MF<-WebGestaltR(enrichMethod="ORA", organism="hsapiens", 
                                       enrichDatabase="geneontology_Molecular_Function",enrichDatabaseFile=NULL, 
                                       enrichDatabaseType=NULL,enrichDatabaseDescriptionFile=NULL,interestGeneFile=NULL, 
                                       interestGene=resOrdered_mayor0_gene_symbol,interestGeneType="genesymbol",collapseMethod="mean",referenceGeneFile=NULL,
                                       referenceGene=NULL,referenceGeneType="genesymbol",referenceSet="genome", minNum=10, maxNum=500,
                                       fdrMethod="BH",sigMethod="fdr",fdrThr=0.05,topThr=10,dNum=20,perNum=1000,
                                       lNum=20,is.output=TRUE,outputDirectory=getwd(),projectName="upregularedChurra_MF",keepGSEAFolder=FALSE,
                                       methodType="R",dagColor="binary",hostName="http://www.webgestalt.org/")
#category<-rep ("Molecular_Function", times= nrow(GenesupregulatedChurra_MF))
#GenesupregulatedChurra_MF_2<-cbind(GenesupregulatedChurra_MF,category)

GenesupregulatedChurra_CC<-WebGestaltR(enrichMethod="ORA", organism="hsapiens", 
                                       enrichDatabase="geneontology_Cellular_Component",enrichDatabaseFile=NULL, 
                                       enrichDatabaseType=NULL,enrichDatabaseDescriptionFile=NULL,interestGeneFile=NULL, 
                                       interestGene=resOrdered_mayor0_gene_symbol,interestGeneType="genesymbol",collapseMethod="mean",referenceGeneFile=NULL,
                                       referenceGene=NULL,referenceGeneType="genesymbol",referenceSet="genome", minNum=10, maxNum=500,
                                       fdrMethod="BH",sigMethod="fdr",fdrThr=0.05,topThr=10,dNum=20,perNum=1000,
                                       lNum=20,is.output=TRUE,outputDirectory=getwd(),projectName="upregularedChurra_CC",keepGSEAFolder=FALSE,
                                       methodType="R",dagColor="binary",hostName="http://www.webgestalt.org/")
category<-rep ("Cellular_Component", times= nrow(GenesupregulatedChurra_CC))
GenesupregulatedChurra_CC_2<-cbind(GenesupregulatedChurra_CC,category)

GO_GenesupregulatedChurra_results<-rbind(GenesupregulatedChurra_BP_2,GenesupregulatedChurra_CC_2)
GO_GenesupregulatedChurra_results<-subset(GO_GenesupregulatedChurra_results, select=c("geneset","description","C","O","E","R","PValue","FDR","OverlapGene_UserID","category"))
#write.table(GO_GenesupregulatedChurra_results, "GO_GenesupregulatedChurra_results.txt")

breed<-rep("Assaf", times= nrow(GO_GenesupregulatedAssaf_results))
GO_GenesupregulatedAssaf_results<-cbind(GO_GenesupregulatedAssaf_results,breed)
breed<-rep("Churra", times= nrow(GO_GenesupregulatedChurra_results))
GO_GenesupregulatedChurra_results<-cbind(GO_GenesupregulatedChurra_results,breed)

GO_results<-rbind(GO_GenesupregulatedAssaf_results,GO_GenesupregulatedChurra_results)
GO_results<-subset(GO_results, select=c("breed","category","geneset","description","C","O","E","R","PValue","FDR","OverlapGene_UserID"))
write.table(GO_results, file="GO_results.txt",sep="\t")

library(foreign)
install.packages("xlsx")
library("xlsx")
write.xlsx(GO_results, file="GO_results.xlsx")

##Check the age of the animals


DESeq.age<-DESeq(dds, test="LRT", reduced=~breed)
#Function for extracting RNA-Seq tables with the results from the DE_analysis

table_results_gene_symbol<- function(x){
  x<-x[order(x$padj),] ##ordering by p-value
  na.omit(x) ##eliminating na
  x<-subset(x, padj<0.05) ##extracting genes with padj<0.05
  
  x<-as.data.frame(x)
  gene_id<- rownames(x)
  rownames(x) <- NULL
  x <- cbind(gene_id,x)
  
  
  gene_symbol_names<- read.table('Ovis_aries_v3.1-r88.corretalionEnsemble_Genesymbol_2',header=TRUE, sep= "\t")
  x_genesymbol<-merge(x,gene_symbol_names,by=c("gene_id"))
  return(x_genesymbol)
}

Resultados_age_reducedModel<-results(DESeq.age, alpha=0.05)
summary(Resultados_age_reducedModel)
genes_age<-table_results_gene_symbol(Resultados_age_reducedModel)
write.table(genes_age, file="resultados_genes_age_modeloltr_gene_symbol.txt", sep='\t') ## writing the table


#Comparison with the genes found DE by breed

DESeq.breed<-DESeq(dds, test="LRT", reduced=~age)
Resultados_breed_reducedModel<-results(DESeq.breed, alpha=0.05)
summary(Resultados_breed_reducedModel)
genes_breed<-table_results_gene_symbol(Resultados_breed_reducedModel)

#comparison of results:

genes_articulo<-read.table(file="resOrdered_gene_symbol.txt", sep='\t', header=T )

##Comparison

genes_articulo_gs<-as.vector(genes_articulo[,8])
genes_breed_lrt_gs<-as.vector(genes_breed[,8])
genes_age_lrt_gs<-as.vector(genes_age[,8])

biocLite ("limma")
library ("limma")
universe <- sort(unique(c(genes_articulo_gs, genes_breed_lrt_gs)))
Counts <- matrix(0, nrow=length(universe), ncol=2)
# Populate the said matrix
for (i in 1:length(universe)) {
  Counts[i,1] <- universe[i] %in% genes_articulo_gs
  Counts[i,2] <- universe[i] %in% genes_breed_lrt_gs
}
colnames(Counts) <- c("genes_articulo_gs", "genes_breed_lrt_gs")
cols<-c("Red", "Green")
jpeg(filename="VennDiagram_articlevslrt.jpg",width=10,height=10, units="in",res=500)
vennDiagram(vennCounts(Counts), circle.col=cols,cex= c(1.5, 1.5, 1))
dev.off()

universe1 <- sort(unique(c(genes_articulo_gs, genes_age_lrt_gs)))
Counts1 <- matrix(0, nrow=length(universe1), ncol=2)
# Populate the said matrix
for (i in 1:length(universe1)) {
  Counts1[i,1] <- universe1[i] %in% genes_articulo_gs
  Counts1[i,2] <- universe1[i] %in% genes_age_lrt_gs
}
colnames(Counts1) <- c("genes_articulo_gs", "genes_age_lrt_gs")
cols<-c("Red", "Green")
jpeg(filename="VennDiagram_articlevslrtage.jpg",width=10,height=10, units="in",res=500)
vennDiagram(vennCounts(Counts1), circle.col=cols,cex= c(1.5, 1.5, 1))
dev.off()

universe2 <- sort(unique(c(genes_breed_lrt_gs, genes_age_lrt_gs)))
Counts2 <- matrix(0, nrow=length(universe2), ncol=2)
# Populate the said matrix
for (i in 1:length(universe2)) {
  Counts2[i,1] <- universe2[i] %in% genes_breed_lrt_gs
  Counts2[i,2] <- universe2[i] %in% genes_age_lrt_gs
}
colnames(Counts2) <- c("genes_breed_lrt_gs", "genes_age_lrt_gs")
cols<-c("Red", "Green")
jpeg(filename="VennDiagram_lrtbreedvslrtage.jpg",width=10,height=10, units="in",res=500)
vennDiagram(vennCounts(Counts2), circle.col=cols,cex= c(1.5, 1.5, 1))
dev.off()

counts2<-Counts
row.names(counts2)<-universe
counts2<-as.data.frame(counts2)
counts3<-counts2[(counts2$genes_FO==1 & counts2$genes_COvsRes_modelCyB==1 & counts2$genes_COvsNores_modelCyB==1) ,]

#plot PCA
vsd <- vst(dds, blind=FALSE)
plotPCA(vsd, intgroup=c("breed", "age"), )



counts_articuloybreed<-Counts
row.names(counts_articuloybreed)<-universe
counts_articuloybreed<-as.data.frame(counts_articuloybreed)
counts_articuloybreedCommon<-counts_articuloybreed[(counts_articuloybreed$genes_articulo_gs==1 & counts_articuloybreed$genes_breed_lrt_gs==1) ,]
counts_articuloybreedCommon<-row.names(counts_articuloybreedCommon)
counts_articuloybreedCommon

counts_articuloyage<-Counts1
row.names(counts_articuloyage)<-universe1
counts_articuloyage<-as.data.frame(counts_articuloyage)
counts_articuloyageCommon<-counts_articuloyage[(counts_articuloyage$genes_articulo_gs==1 & counts_articuloyage$genes_age_lrt_gs==1) ,]
counts_articuloyageCommon<-row.names(counts_articuloyageCommon)
counts_articuloyageCommon


##lncRNA DE
genes_articulo_gs2<-as.vector(genes_articulo[,1])
lncRNA<- read.table(file = "~/trabajo2018/grasa_lamb/lncRNA2.txt",sep = "\t", header = T)
lncRNA<-as.vector(lncRNA[,1])
universe3 <- sort(unique(c(genes_articulo_gs2, lncRNA)))
Counts3 <- matrix(0, nrow=length(universe3), ncol=2)
# Populate the said matrix
for (i in 1:length(universe3)) {
  Counts3[i,1] <- universe3[i] %in% genes_articulo_gs2
  Counts3[i,2] <- universe3[i] %in% lncRNA
}
colnames(Counts3) <- c("genes_articulo_gs2", "lncRNA")
cols<-c("Red", "Green")
jpeg(filename="VennDiagram_articlevslncRNA.jpg",width=10,height=10, units="in",res=500)
vennDiagram(vennCounts(Counts3), circle.col=cols,cex= c(1.5, 1.5, 1))
dev.off()

counts_lncRNA<-Counts3
row.names(counts_lncRNA)<-universe3
counts_lncRNA<-as.data.frame(counts_lncRNA)
counts_lncRNACommon<-counts_lncRNA[(counts_lncRNA$genes_articulo_gs2==1 & counts_lncRNA$lncRNA==1) ,]
counts_lncRNACommon<-row.names(counts_lncRNACommon)
counts_lncRNACommon

##Comparison of deGenes with GENES in lncRNAregions
ncRNA_relatedGenes<- read.table(file = "~/trabajo2018/grasa_lamb/Grasa_r88_RSEM/Genes_withinlncRNA_regions2.txt",sep = "\t", header = T)
lncRNA_relatedGenes<-as.vector(ncRNA_relatedGenes[,1])
universe4 <- sort(unique(c(genes_articulo_gs2, lncRNA_relatedGenes)))
Counts4 <- matrix(0, nrow=length(universe4), ncol=2)
# Populate the said matrix
for (i in 1:length(universe4)) {
  Counts4[i,1] <- universe4[i] %in% genes_articulo_gs2
  Counts4[i,2] <- universe4[i] %in% lncRNA_relatedGenes
}
colnames(Counts4) <- c("genes_articulo_gs2", "lncRNA_relatedGenes")
cols<-c("Red", "Green")
jpeg(filename="VennDiagram_articlevslncRNA_relatedGenes.jpg",width=10,height=10, units="in",res=500)
vennDiagram(vennCounts(Counts4), circle.col=cols,cex= c(1.5, 1.5, 1))
dev.off()

counts_lncRNA_relatedgenes<-Counts4
row.names(counts_lncRNA_relatedgenes)<-universe4
counts_lncRNA_relatedgenes<-as.data.frame(counts_lncRNA_relatedgenes)
counts_lncRNA_relatedgenesCommon<-counts_lncRNA_relatedgenes[(counts_lncRNA_relatedgenes$genes_articulo_gs2==1 & counts_lncRNA_relatedgenes$lncRNA==1) ,]
counts_lncRNA_relatedgenesCommon<-row.names(counts_lncRNA_relatedgenesCommon)
counts_lncRNA_relatedgenesCommon<-as.vector(counts_lncRNA_relatedgenesCommon)

gene_symbol_names<- read.table('Ovis_aries_v3.1-r88.corretalionEnsemble_Genesymbol_2', header=TRUE, sep= "\t")
lncRNA_relatedGenes_gsymbol <- gene_symbol_names[gene_symbol_names$gene_id %in% counts_lncRNA_relatedgenesCommon, ]

#Test age independently
dds <- DESeqDataSetFromTximport(txi.rsem, colData2, design = ~age)

DESeq.age<-DESeq(dds)
Resultados_age<-results(DESeq.age, alpha=0.05)
resultsNames(DESeq.age)
summary(Resultados_age)
genes_age<-table_results_gene_symbol(Resultados_age)

genes_age_gs<-as.vector(genes_age[,8])

biocLite ("limma")
library ("limma")
universe <- sort(unique(c(genes_articulo_gs, genes_age_gs)))
Counts <- matrix(0, nrow=length(universe), ncol=2)
# Populate the said matrix
for (i in 1:length(universe)) {
  Counts[i,1] <- universe[i] %in% genes_articulo_gs
  Counts[i,2] <- universe[i] %in% genes_age_gs
}
colnames(Counts) <- c("genes_articulo_gs", "genes_age_gs")
cols<-c("Red", "Green")
jpeg(filename="VennDiagram_articlevsageindependent.jpg",width=10,height=10, units="in",res=500)
vennDiagram(vennCounts(Counts), circle.col=cols,cex= c(1.5, 1.5, 1))
dev.off()

counts_articuloyage<-Counts
row.names(counts_articuloyage)<-universe
counts_articuloyage<-as.data.frame(counts_articuloyage)
counts_articuloyageCommon<-counts_articuloyage[(counts_articuloyage$genes_articulo_gs==1 & counts_articuloyage$genes_age_gs==1) ,]
counts_articuloyageCommon<-row.names(counts_articuloyageCommon)
counts_articuloyageCommon

