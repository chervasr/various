#####funcion nombre de genes#######
getgenes<-function(resul){
  genspp<-as.data.frame(row.names(resul))
  colnames(genspp)<-c("V1")
  glist <- read.delim("C:/Users/carlo/Desktop/cosasdoctorao/rstudio/glist-ram1.0r101.txt", header=FALSE)
  biem<-merge(x = genspp, y = glist, by.x = "V1", by.y = "V4", all.y =F)
  return(as.data.frame(biem$V5)) }
###############getgenes2##########
getgenes2<-function(resul){
  genspp<-as.data.frame(row.names(resul))
  colnames(genspp)<-c("V1")
  glist <- read.delim("C:/Users/carlo/Desktop/cosasdoctorao/rstudio/glist-ram1.0r101.txt", header=FALSE)
  biem<-merge(x = genspp, y = glist, by.x = "V1", by.y = "V4", all.y =F)
  biem$V5 <- ifelse((biem$V5 == ""), biem$V1, biem$V5)
  return(as.data.frame(biem$V5)) }
############################
setwd("C:/Users/carlo/Desktop/cosasdoctorao/rstudio")
library(dplyr)
library("pheatmap")
#biocLite ("tximport")
#biocLite ("readr")
library(tximport)
library(readr)
##BiocManager::install("DESeq2")
dir<-"C:/Users/carlo/Desktop/cosasdoctorao/rstudio/retoinflamatorio"
dir
files
samples<- read.table(file.path(dir, "nombresreto.txt"), header = TRUE)
files<- file.path(dir, "counts", paste0(samples$File, "count.genes.results"))
#samples$x <- apply( samples[ , 3:4 ] , 1 , paste , collapse = "_" )
names(files) <-paste0(samples$File)
file.exists(files)
all(file.exists(files))
txi.rsem<-tximport(files, type="rsem")
head(txi.rsem$counts)
names(txi.rsem)
library(DESeq2)

trunk<-trunc(txi.rsem$counts)

dds <- DESeqDataSetFromMatrix(countData = trunk, colData = samples, design = ~Challenge)

dds <-dds[rowSums(counts(dds))>1,]
aaa<-collapseReplicates(dds, dds$ID_CNAG)

aaa1 <-DESeq(dds)
pheatmap(cor(log10(counts(aaa1)+1)))
rld<-rlog(aaa1)
plotPCA(rld, intgroup="Challenge",returnData= TRUE)

res<-results(aaa1, alpha=0.05)
resOrdered1<-res[order(res$padj),]
resOrdered1<-subset(resOrdered1, padj<0.05)
head(resOrdered1)
write.table(resOrdered, file="resOrderedprevspost.txt", sep="\t")
counts.normalized <- counts(aaa1, normalized= T)
log.norm.counts <- log2(counts.normalized +1)
DGEgenes<- rownames (resOrdered)
hm.mat_DGEgenes<-log.norm.counts[DGEgenes,]
#Heatmap con los genes DE
heatmap(hm.mat_DGEgenes)
####pre vs post###
dds <- DESeqDataSetFromMatrix(countData = trunk, colData = samples, design = ~ time)
plotCounts(dds, gene="ENSOARG00020012683", intgroup="time")
help("write.table")
library(dplyr)
####res_POS VS NORESPOS###
samplos<- samples %>% filter( time == "post")
files<- file.path(dir, "samples", paste0(samplos$ID_seq, "count.genes.results"))
names(files) <-paste0(samplos$ID_seq)
file.exists(files)
all(file.exists(files))
txi.rsem<-tximport(files, type="rsem")
head(txi.rsem$counts)
names(txi.rsem)
library(DESeq2)
trunk<-trunc(txi.rsem$counts)
dds <- DESeqDataSetFromMatrix(countData = trunk, colData = samplos, design = ~ Condition)

dds <-dds[rowSums(counts(dds))>1,]
aaa<-collapseReplicates(dds, dds$ID_CNAG)

aaa1 <-DESeq(aaa)
pheatmap(cor(log10(counts(aaa1)+1)))
rld<-rlog(aaa1)
plotPCA(rld, intgroup="x")

res<-results(aaa1, alpha=0.05)
resOrdered2<-res[order(res$padj),]
resOrdered2
resOrdered2<-subset(resOrdered, padj<0.20)
head(resOrdered2)
write.table(resOrdered, file="resOrderednoresposvsrespost.txt", sep="\t")
counts.normalized <- counts(aaa1, normalized= T)
log.norm.counts <- log2(counts.normalized +1)
DGEgenes<- rownames (resOrdered)
hm.mat_DGEgenes<-log.norm.counts[DGEgenes,]
#Heatmap con los genes DE
heatmap(hm.mat_DGEgenes)
##nores_pre vs res_post
library(dplyr)
samples$x <- apply( samples[ , 3:4 ] , 1 , paste , collapse = "_" )
samplos<- samples %>% filter( x == "nores_pre" | x == "res_post")
files<- file.path(dir, "samples", paste0(samplos$ID_seq, "count.genes.results"))
names(files) <-paste0(samplos$ID_seq)
file.exists(files)
all(file.exists(files))
txi.rsem<-tximport(files, type="rsem")
head(txi.rsem$counts)
names(txi.rsem)
library(DESeq2)
trunk<-trunc(txi.rsem$counts)
dds <- DESeqDataSetFromMatrix(countData = trunk, colData = samplos, design = ~ x)

dds <-dds[rowSums(counts(dds))>1,]
aaa<-collapseReplicates(dds, dds$ID_CNAG)

aaa1 <-DESeq(aaa)
pheatmap(cor(log10(counts(aaa1)+1)))
rld<-rlog(aaa1)
plotPCA(rld, intgroup="x")

res<-results(aaa1, alpha=0.05)
resOrdered3<-res[order(res$padj),]
resOrdered3
resOrdered3<-subset(resOrdered3, padj<0.05)
head(resOrdered3)
write.table(resOrdered, file="resOrderednonoresprevsrespost.txt", sep="\t")
counts.normalized <- counts(aaa1, normalized= T)
log.norm.counts <- log2(counts.normalized +1)
DGEgenes<- rownames (resOrdered)
hm.mat_DGEgenes<-log.norm.counts[DGEgenes,]
heatmap(hm.mat_DGEgenes)
plotCounts(dds, gene="ENSOARG00020000549", intgroup="x")
DGEgenes
install.packages("VennDiagram")
library(VennDiagram)
##res_pre vs nores_post
library(dplyr)

samples$x <- apply( samples[ , 3:4 ] , 1 , paste , collapse = "_" )
samplos<- samples %>% filter( x == "res_pre" | x == "nores_post")
files<- file.path(dir, "samples", paste0(samplos$ID_seq, "count.genes.results"))
names(files) <-paste0(samplos$ID_seq)
file.exists(files)
all(file.exists(files))
txi.rsem<-tximport(files, type="rsem")
head(txi.rsem$counts)
names(txi.rsem)
library(DESeq2)
trunk<-trunc(txi.rsem$counts)
dds <- DESeqDataSetFromMatrix(countData = trunk, colData = samplos, design = ~ x)

dds <-dds[rowSums(counts(dds))>1,]
aaa<-collapseReplicates(dds, dds$ID_CNAG)

aaa1 <-DESeq(aaa)
pheatmap(cor(log10(counts(aaa1)+1)))
rld<-rlog(aaa1)
plotPCA(rld, intgroup="x")

res<-results(aaa1, alpha=0.05)
resOrdered4<-res[order(res$padj),]
resOrdered4
resOrdered4<-subset(resOrdered4, padj<0.05)
head(resOrdered4)
write.table(resOrdered, file="resOrderedresprevsnorespost.txt", sep="\t")
counts.normalized <- counts(aaa1, normalized= T)
log.norm.counts <- log2(counts.normalized +1)
DGEgenes<- rownames (resOrdered)
hm.mat_DGEgenes<-log.norm.counts[DGEgenes,]
heatmap(hm.mat_DGEgenes)
plotCounts(dds, gene="ENSOARG00020000549", intgroup="x")

DGEgenes
install.packages("VennDiagram")
prevspos<-rownames(resOrdered1)
noresprevsrespos<-rownames(resOrdered3)
resprevsnorespos<-rownames(resOrdered4)
venn.plot <- venn.diagram(
  x = list(
    prevspos = prevspos,
    noresprevsrespos = noresprevsrespos,
    resprevsnorespos = resprevsnorespos
  ),
  filename = "venn.tiff",
  col = "black",
  fill = c("goldenrod1", "darkorange1", "seagreen3"),
  alpha = 0.50,
    cat.col = c("goldenrod1", "darkorange1", "seagreen3"),
  cat.cex = 1.5,
  cat.fontface = "bold",
  margin = 0.05
)
todos<-c(prevspos,noresprevsrespos)
tos<-todos[duplicated(todos)]
tos<-c(tos,resprevsnorespos)
tos<-tos[duplicated(tos)]
tos
install.packages("tidyverse")
library(tidyverse)
plotCounts(dds, gene="ENSOARG00020001966", intgroup="x")
####comparar prepost y norespreresp
prevspost <- read.delim("C:/Users/carlo/Desktop/cosasdoctorao/rstudio/resOrderedprevspost.txt")
noresprevsrespost <- read.delim("C:/Users/carlo/Desktop/cosasdoctorao/rstudio/resOrderednonoresprevsrespost.txt")
genspp<-row.names(prevspost)
gensdif<-row.names(noresprevsrespost)
todos<-c(genspp,gensdif)
tos<-todos[duplicated(todos)]
tos<-as.data.frame(tos)
####buscar nombre genes####
glist <- read.delim("C:/Users/carlo/Desktop/cosasdoctorao/rstudio/glist-ram1.0r101.txt", header=FALSE)
biem<-merge(x = tos, y = glist, by.x = "tos", by.y = "V4", all.y = F)
####compararlos solo los dos que se ven fuera####
dir<-"C:/Users/carlo/Desktop/cosasdoctorao/rstudio/countgenes"
dir
samples<- read.table(file.path(dir, "datos2.txt"), header = TRUE)
samples<-samples %>% 
  filter( samples$estos == "este")
samples
files<- file.path(dir, "samples", paste0(samples$ID_seq, "count.genes.results"))
samples$x <- apply( samples[ , 3:4 ] , 1 , paste , collapse = "_" )
names(files) <-paste0(samples$ID_seq)
file.exists(files)
all(file.exists(files))
txi.rsem<-tximport(files, type="rsem")
head(txi.rsem$counts)
names(txi.rsem)
trunk<-trunc(txi.rsem$counts)
dds <- DESeqDataSetFromMatrix(countData = trunk, colData = samples, design = ~ x)

dds <-dds[rowSums(counts(dds))>=1,]
aaa<-collapseReplicates(dds, dds$ID_CNAG)

aaa1 <-DESeq(aaa)
pheatmap(cor(log10(counts(aaa1)+1)))
rld<-rlog(aaa1)
plotPCA(rld, intgroup="x")

plotPCA(rld, intgroup="x",returnData= TRUE)

res<-results(aaa1, alpha=0.5)
resOrdered1<-res[order(res$padj),]
resOrdered1<-subset(resOrdered1, padj<0.05)
head(resOrdered1)
write.table(resOrdered1, file="2vs2t.txt", sep="\t")
counts.normalized <- counts(aaa1, normalized= T)
log.norm.counts <- log2(counts.normalized +1)
DGEgenes<- rownames (resOrdered)
hm.mat_DGEgenes<-log.norm.counts[DGEgenes,]
getgenes<-function(resul){
  genspp<-as.data.frame(row.names(resul))
  colnames(genspp)<-c("V1")
  glist <- read.delim("C:/Users/carlo/Desktop/cosasdoctorao/rstudio/glist-ram1.0r101.txt", header=FALSE)
  biem<-merge(x = genspp, y = glist, by.x = "V1", by.y = "V4", all.y =F)
  return(as.data.frame(biem$V5)) }
vs2t <- read.delim("C:/Users/carlo/Desktop/cosasdoctorao/rstudio/2vs2t.txt")
genes<-getgenes(vs2t)
write.table(genes,"nombregenes2vs2.txt", append = F, quote = F, row.names = F, col.names = F)
genesaroa <- read.table("C:/Users/carlo/Desktop/cosasdoctorao/rstudio/genesaroa.txt", quote="\"", comment.char="")
####comprobar con genes aroa####
juntos<-(c(genes$`biem$V5`,genesaroa$V1))
juntos<-juntos[duplicated(juntos)]
juntos
genesjunt<-unique(juntos)
write.table(genesjunt,"genesvsaora.txt", append = F, quote = F, row.names = F, col.names = F)
####diaframa de ven pa esto
venn.plot <- venn.diagram(
  x = list(
    g2vs2 = genes$`biem$V5`,
    genespaper = genesaroa$V1),
  filename = "g2vspaper.tiff",
  col = "black",
  fill = c("goldenrod1", "seagreen3"),
  alpha = 0.50,
  cat.col = c("goldenrod1", "seagreen3"),
  cat.cex = 1.5,
  cat.fontface = "bold",
  margin = 0.05
)

