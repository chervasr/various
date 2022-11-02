is_outlier<-function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}
install.packages("ggbeeswarm")
install.packages("ggrepel")
library(ggbeeswarm)
library(ggrepel)
install.packages("sctransform")
install.packages("reshape2")
library(reshape2)
library(sctransform)
install.packages("ggplot2")
library(dplyr)
library(ggplot2)

#####nuevo comienzo franceses#####
library(readxl)
Plasmafran <- read_excel("C:/Users/carlo/Desktop/cosasdoctorao/rstudio/franceses/UNILEONPlasmaresults02122020.xlsx")
UNILEONPlasmaresults02122020 <- read_excel("UNILEONPlasmaresults02122020.xlsx", 
                                           col_types = c("skip", "text", "text", 
                                                         "text", "text", "text", "numeric", 
                                                         "numeric", "numeric", "numeric", "numeric", "numeric", 
                                                         "numeric", "numeric", "numeric", "numeric", "numeric", 
                                                         "numeric", "numeric", "numeric"))
View(UNILEONPlasmaresults02122020)
grupos<-cbind.data.frame(UNILEON.Plasma.results.02122020_bea_forR$Animal_ID,UNILEON.Plasma.results.02122020_bea_forR$Tto_Recr,UNILEON.Plasma.results.02122020_bea_forR$Group_EBV)
colnames(grupos)<-c("Animal ID", "Challenge", "Response")
grupos<-unique(grupos)
plasmon<-merge(UNILEONPlasmaresults02122020,grupos, by= "Animal ID")
plasmon
#write.table(plasmon,"plasmaresult.txt", append = F, quote = F, row.names = F, col.names = T)
for (i in unique(plasmaresult$Time)){
agru<-plasmaresult %>% 
  filter(Time == paste0(i))
acv<-split(agru1, agru1$Challenge)
acv
agru1<-acv$C %>% 
    mutate(outlier=ifelse(is_outlier(IFng),Animal_ID,as.numeric(NA))) 
agru2<-acv$C %>% 
  mutate(outlier=ifelse(is_outlier(IFng),Animal_ID,as.numeric(NA)))
agru3<-rbind.data.frame(agru1,agru2)
agru3
print(agru1)
print(ggplot(agru1, aes(x=Challenge, y=IFng)) + 
  geom_boxplot(outlier.colour = NA)+
  ggtitle(paste0(i))+
    ggbeeswarm::geom_beeswarm(aes(color=IFng)) +
    ggrepel::geom_text_repel(data=. %>% filter(!is.na(outlier)), aes(label=Animal_ID))
)
}
bc<-cbind.data.frame(agru3$Time,agru3$outlier,aa)
bc<-bc[1,]
bc
inter<-names(agru[-c(1:5)])
plot_list = list()
plot_list2 = list()
for (aa in colnames(agru[c(6:19)])){
for (i in unique(plasmaresult$Time)){
  agru<-plasmaresult %>% 
    filter(Time == paste0(i))
  agru22<-cbind.data.frame(agru$Animal_ID,agru$Time,agru$Challenge,agru[[aa]])
  colnames(agru22)<-c("Animal_ID","Time","Challenge","marker")
  #print(agru22)
  acv<-split(agru22, agru22$Challenge)
  acv
  agru1<-acv$C %>% 
    mutate(outlier=ifelse(is_outlier(marker),Animal_ID,as.numeric(NA))) 
  agru2<-acv$NC %>% 
    mutate(outlier=ifelse(is_outlier(marker),Animal_ID,as.numeric(NA)))
  agru3<-rbind.data.frame(agru1,agru2)
  
  aaa<-cbind.data.frame(agru3$Time,agru3$outlier,aa)
  aaa<-na.omit(aaa)
  #write.table(aaa,file = "outliers.txt", append = T, row.names = F,col.names = F,quote = F)
  agru1<-agru22 %>% 
   group_by(Challenge) %>%
    mutate(outlier=ifelse(is_outlier(agru22[["marker"]]),Animal_ID,as.numeric(NA))) 
  #print(class(aa))
  #print(class(i))
  #print(agru1)
  p<-ggplot(agru3, aes(x=Challenge, y=marker)) + 
          geom_boxplot(outlier.colour = NA)+
          ggtitle(paste0(i),paste0(aa))+
          ggbeeswarm::geom_beeswarm(aes(color=marker)) +
          ggrepel::geom_text_repel(data=. %>% filter(!is.na(outlier)), aes(label=Animal_ID))
  plot_list[[i]] = p         
   
}
  plot_list2[[aa]] = plot_list  #print(bc)
}
pdf("plots.pdf")
for (i in 1:14) {
  print(plot_list2[[i]])
}
dev.off()

tabla<-table(outliers)
write.table(tabla,file = "tablaoutliers.txt", append = F,quote = F)
tablaanimal<-table(outliers$Animal)
write.table(tablaanimal,file = "tablaoutliersporanimal.txt", append = F,quote = F)
#####

data<- LPS_datos2
data<-data[,1:15]
data$cosa <- paste(data$TTO,data$VG)
data[2:4]<-NULL
data
desafioH<-data %>%
  filter(cosa == "DESAFÍO H")
desafioL<-data %>%
  filter(cosa == "DESAFÍO L")
controlH<-data %>%
  filter(cosa == "CONTROL H")
controlL<-data %>%
  filter(cosa == "CONTROL L")
data_transpose <- as.data.frame(t(as.matrix(desafioH)))
data_transpose$Names <- rownames(data_transpose)

data_transpose<-data_transpose[-13,]

colnames(data_transpose)[9]<-"temp"
colnames(data_transpose)<-data_transpose[1,]
data_transpose<-data_transpose[-1,]
df <- melt(data_transpose ,  id.vars = 'temp', variable.name = 'names')
data1
df$temp <- factor(df$temp, levels = unique(df$temp))
df$value<-as.numeric(df$value)
str(df)
ggplot(df, aes(x = temp, y = value, 
               group = names , color = names)) + geom_line() + ggtitle("Challenged High response")+labs(y="Temperature")+ facet_grid(names ~. ) +
  scale_y_continuous(limits=c(37.5, 40.5), breaks=c(38, 39, 40))
 ####♠desafio low####
data_transpose <- as.data.frame(t(as.matrix(desafioL)))
data_transpose$Names <- rownames(data_transpose)

data_transpose<-data_transpose[-13,]

colnames(data_transpose)[9]<-"temp"
colnames(data_transpose)<-data_transpose[1,]
data_transpose<-data_transpose[-1,]
df <- melt(data_transpose ,  id.vars = 'temp', variable.name = 'names')
data_transpose
df$temp <- factor(df$temp, levels = unique(df$temp))
df$value<-as.numeric(df$value)
str(df)
ggplot(df, aes(x = temp, y = value, 
               group = names , color = names)) + geom_line() + ggtitle("Challenged Low response")+labs(y="Temperature")+ facet_grid(names ~. ) +
  scale_y_continuous(limits=c(37.5, 40.5), breaks=c(38, 39, 40))
####♠control hidh####
data_transpose <- as.data.frame(t(as.matrix(controlH)))
data_transpose$Names <- rownames(data_transpose)

data_transpose<-data_transpose[-13,]

colnames(data_transpose)[8]<-"temp"
colnames(data_transpose)<-data_transpose[1,]
data_transpose<-data_transpose[-1,]
df <- melt(data_transpose ,  id.vars = 'temp', variable.name = 'names')
data_transpose
df$temp <- factor(df$temp, levels = unique(df$temp))
df$value<-as.numeric(df$value)
str(df)
ggplot(df, aes(x = temp, y = value, 
               group = names , color = names)) + geom_line() + ggtitle("Control High response")+labs(y="Temperature")+ facet_grid(names ~. ) +
  scale_y_continuous(limits=c(37.5, 40.5), breaks=c(38, 39, 40))
####♠control hidh####
data_transpose <- as.data.frame(t(as.matrix(controlL)))
data_transpose$Names <- rownames(data_transpose)

data_transpose<-data_transpose[-13,]

colnames(data_transpose)[8]<-"temp"
colnames(data_transpose)<-data_transpose[1,]
data_transpose<-data_transpose[-1,]
df <- melt(data_transpose ,  id.vars = 'temp', variable.name = 'names')
data_transpose
df$temp <- factor(df$temp, levels = unique(df$temp))
df$value<-as.numeric(df$value)
str(df)
ggplot(df, aes(x = temp, y = value, 
               group = names , color = names)) + geom_line() + ggtitle("Control Low response")+labs(y="Temperature")+ facet_grid(names ~. ) +
  scale_y_continuous(limits=c(37.5, 40.5), breaks=c(38, 39, 40))

############celulas#########
data<- LPS_datos2celulas
data<-data[,1:15]
data$cosa <- paste(data$TTO,data$VG)
data[2:3]<-NULL
data
desafioH<-data %>%
  filter(cosa == "DESAF?O H")
desafioL<-data %>%
  filter(cosa == "DESAF?O L")
controlH<-data %>%
  filter(cosa == "CONTROL H")
controlL<-data %>%
  filter(cosa == "CONTROL L")
data_transpose <- as.data.frame(t(as.matrix(desafioH)))
data_transpose$Names <- rownames(data_transpose)

data_transpose<-data_transpose[-9,]

colnames(data_transpose)[9]<-"cells"
colnames(data_transpose)<-data_transpose[1,]
data_transpose<-data_transpose[-1,]
df <- melt(data_transpose ,  id.vars = 'cells', variable.name = 'names')
df$cells <- factor(df$cells, levels = unique(df$cells))
df$value<-as.numeric(df$value)
str(df)
ggplot(df, aes(x = temp, y = value, 
               group = names , color = names)) + geom_line() + ggtitle("Challenged High response")+labs(y="Cell count")+ facet_grid(names ~. )+
  scale_y_continuous(breaks=c(0, 15000, 30000))
####♠desafio low####
data_transpose <- as.data.frame(t(as.matrix(desafioL)))
data_transpose$Names <- rownames(data_transpose)

data_transpose<-data_transpose[-9,]

colnames(data_transpose)[9]<-"cells"
colnames(data_transpose)<-data_transpose[1,]
data_transpose<-data_transpose[-1,]
df <- melt(data_transpose ,  id.vars = 'cells', variable.name = 'names')
data_transpose
df$cells <- factor(df$cells, levels = unique(df$cells))
df$value<-as.numeric(df$value)
str(df)
ggplot(df, aes(x = cells, y = value, 
               group = names , color = names)) + geom_line() + ggtitle("Challenged Low response")+labs(y="Cell count")+ facet_grid(names ~. ) +
  scale_y_continuous(breaks=c(0, 15000, 30000))
####♠control hidh####
data_transpose <- as.data.frame(t(as.matrix(controlH)))
data_transpose$Names <- rownames(data_transpose)

data_transpose<-data_transpose[-9,]

colnames(data_transpose)[8]<-"cells"
colnames(data_transpose)<-data_transpose[1,]
data_transpose<-data_transpose[-1,]
df <- melt(data_transpose ,  id.vars = 'cells', variable.name = 'names')
data_transpose
df$cells <- factor(df$cells, levels = unique(df$cells))
df$value<-as.numeric(df$value)
str(df)
ggplot(df, aes(x = cells, y = value, 
               group = names , color = names)) + geom_line() + ggtitle("Control High response")+labs(y="Cell count")+ facet_grid(names ~. ) +
  scale_y_continuous(breaks=c(0, 15000, 30000))
####♠control low####
data_transpose <- as.data.frame(t(as.matrix(controlL)))
data_transpose$Names <- rownames(data_transpose)

data_transpose<-data_transpose[-9,]

colnames(data_transpose)[8]<-"cells"
colnames(data_transpose)<-data_transpose[1,]
data_transpose<-data_transpose[-1,]
df <- melt(data_transpose ,  id.vars = 'cells', variable.name = 'names')
data_transpose
df$cells <- factor(df$cells, levels = unique(df$cells))
df$value<-as.numeric(df$value)
str(df)
ggplot(df, aes(x = cells, y = value, 
               group = names , color = names)) + geom_line() + ggtitle("Control Low response")+labs(y="Cell count")+ facet_grid(names ~. ) +
  scale_y_continuous(breaks=c(0, 15000, 30000))
outliers <- read.csv("C:/Users/carlo/Desktop/cosasdoctorao/rstudio/franceses/outliers.txt", sep="")
animalyhora<-table(outliers$Time, outliers$Animal)
write.table(animalyhora,file = "tablaoutliersanimalyhora.txt", append = F,quote = F, sep = "\t")
###############################
plasmaresultpp <- read.delim("C:/Users/carlo/Desktop/cosasdoctorao/rstudio/franceses/plasmaresultpp.txt")
data<-plasmaresultpp
data$cosa <- paste(data$Challenge,data$Response)
data[17:18]<-NULL
data[3:5]<-NULL
data$Time <- factor(data$Time, levels = unique(data$Time))
data
desafioH<-data %>%
  filter(cosa == "C H")
desafioL<-data %>%
  filter(cosa == "C L")
controlH<-data %>%
  filter(cosa == "NC H")
controlL<-data %>%
  filter(cosa == "NC L")
####
library(varhandle)
aaa<-unfactor(desafioH)
aaa[aaa$Time=="-24h"]<- "T1(-24h)"
str(aaa)
dx<-desafioH[,1:3]
df <- melt(dx ,  id.vars = 'Animal_ID', variable.name = 'names')
df <- melt(data_transpose ,  id.vars = 'cells', variable.name = 'names')
ggplot(data = dx, aes(x = Time, y = IFng, 
                      group = Animal_ID, color = "retoh")) +
  geom_point(data = desafioH , aes(x = Time, y = IFng, 
                                   group = Animal_ID, color = "retoL"))+
  geom_point(size = 2) +  # add points at the ends, size = 2
  geom_line()
scale_y_continuous(limits=c(0, 150))+
  geom_point(size = 2) +  # add points at the ends, size = 2
  geom_line() 

dx$Animal_ID<-as.character(dx$Animal_ID)

ggplot(dx, aes(x = Time, y = IFng, 
               group = Animal_ID , color = Animal_ID )) + geom_line() + ggtitle("Challenged High response")+labs(y="Temperature")+ facet_grid(Animal_ID ~. )# +
#scale_y_continuous(limits=c(37.5, 40.5), breaks=c(38, 39, 40))
plot_list = list()
plot_list2 = list()
#for (zz in c("desafioH","desafioL","controlH","controlL")){
for (aa in colnames(desafioH[c(3:16)])){
  #az<-get(paste0(zz))
  zz<-"controlL"
  az<-controlL
  dx<-cbind.data.frame(az[,1:2],az[[aa]])
  colnames(dx)<-c("Animal","Time","marker")
  dx$Animal<-as.character(dx$Animal)
  p<-ggplot(dx, aes(x = Time, y = marker, 
                    group = Animal , color = Animal )) + geom_line() + ggtitle(paste0(zz))+labs(y=paste0(aa))+ facet_grid(Animal ~. )# +
  #scale_y_continuous(limits=c(37.5, 40.5), breaks=c(38, 39, 40))
  print(p)
  plot_list[[aa]] = p         
  
}

pdf("controlL3.pdf")
for (i in 1:14) {
  print(plot_list[[i]])
}
dev.off()
data_transpose <- as.data.frame(t(as.matrix(desafioH)))
data_transpose$Names <- rownames(data_transpose)

data_transpose<-data_transpose[-13,]

for (aa in colnames(desafioH[c(3:16)])){
  #az<-get(paste0(zz))
  zz<-"controlL"
  az<-controlL
  dx<-cbind.data.frame(az[,1:2],az[[aa]])
  colnames(dx)<-c("Animal","Time","marker")
  dx$Animal<-as.character(dx$Animal)
  p<-ggplot(dx, aes(x = Time, y = marker, 
                    group = Animal , color = Animal )) + geom_line() + ggtitle(paste0(zz))+labs(y=paste0(aa))+ facet_grid(Animal ~. )# +
  #scale_y_continuous(limits=c(37.5, 40.5), breaks=c(38, 39, 40))
  print(p)
  plot_list[[aa]] = p         
  
}
###################
plasmaresultpp <- read.delim("C:/Users/carlo/Desktop/cosasdoctorao/rstudio/franceses/plasmaresultpp.txt")
data<-plasmaresultpp
data$cosa <- paste(data$Challenge,data$Response)
data[17:18]<-NULL
data$Time<-as.factor(data$Time)
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}
plotcondesv<-list()
for (aa in colnames(data[c(3:16)])){
  df2 <- data_summary(data, varname=paste0(aa), 
                      groupnames=c("cosa","Time"))
  head(df2)
  library(ggplot2)
  p<- ggplot(df2, aes(x=Time, y=get(paste0(aa)), group=cosa, color=cosa))+ facet_grid(cosa ~. ) + 
    labs(y = paste0(aa))+
    labs(colour = "Groups")+
    geom_line() +
    geom_point()+
    geom_errorbar(aes(ymin=get(paste0(aa))-sd, ymax=get(paste0(aa))+sd), width=.2,
                  position=position_dodge(0.05))
  plotcondesv[[aa]] = p
  print(p)
}
pdf("plotsdesvi45.pdf")
for (i in 1:14) {
  print(plotcondesv[[i]])
}
dev.off()
print(plotcondesv[[1]])
#######
library(dplyr)
data<-plasmaresultpp
data$Response<-NULL
p2<-ggplot(df2, aes(x=Time, y=IFng,group=Challenge,color=Challenge))+
  geom_point()+geom_line()+scale_color_manual(values=c('#0000FF','#FF0000'))+
  #geom_errorbar(aes(ymin=SCC-sd_SCC, ymax=SCC+sd_SCC), width=.1)+
  labs(title="Evolution of SCC accross the sampling time points",x="Time point", y = "SCC")#+
#theme_classic()
p2
df2 <- data_summary(data, "IFng", 
                    groupnames=c("Challenge","Time"))
head(df2)
plotcnc<-list()
for (aa in colnames(data[c(3:16)])){
  df2 <- data_summary(data, varname=paste0(aa), 
                      groupnames=c("Challenge","Time"))
  head(df2)
  library(ggplot2)
  p<- ggplot(df2, aes(x=Time, y=get(paste0(aa)), group=Challenge, color=Challenge))+ #facet_grid(cosa ~. ) + 
    labs(y = paste0(aa))+
    labs(colour = "Groups")+
    geom_line() +
    geom_point()+
    geom_errorbar(aes(ymin=get(paste0(aa))-sd, ymax=get(paste0(aa))+sd), width=.2,
                  position=position_dodge(0.05))
  plotcnc[[aa]] = p
  print(p)
}
plotcnc[["IFng"]]
pdf("plotscnc.pdf")
for (aa in colnames(data[c(3:16)])) {
  print(plotcnc[[aa]])
}
dev.off()
plotcnc<-list()
for (aa in colnames(data[c(3:16)])){
  df2 <- data_summary(data, varname=paste0(aa), 
                      groupnames="Time")
  head(df2)
  library(ggplot2)
  p<- ggplot(df2, aes(x=Time, y=get(paste0(aa)),group = 1))+ #facet_grid(cosa ~. ) + 
    labs(y = paste0(aa))+
    geom_line() +
    geom_point(color='darkblue')+
    geom_point()+
    geom_errorbar(aes(ymin=get(paste0(aa))-sd, ymax=get(paste0(aa))+sd), width=.2,
                  position=position_dodge(0.05))
  plotcnc[[aa]] = p
  print(p)
}

pdf("lineassolo1.pdf")
for (aa in colnames(data[c(3:16)])) {
  print(plotcnc[[aa]])
}
dev.off()
ifng<-cbind.data.frame(data[,1:3])
res.aov2 <- aov(len ~ supp + dose, data = ifng)
ifng2<-as.data.frame((t(ifng)))
head(ifng2)
colnames(ifng2)<-ifng2[2,]
ifng2<- ifng2[-c(1, 2), ]
abcd<-cbind.data.frame(as.data.frame(t(ifng2[,1:24])),as.data.frame(t(ifng2[,25:48])),as.data.frame(t(ifng2[,49:72])),as.data.frame(t(ifng2[,73:96])),as.data.frame(t(ifng2[,97:120])),as.data.frame(t(ifng2[,121:144])),as.data.frame(t(ifng2[,145:168])),as.data.frame(t(ifng2[,169:192])),as.data.frame(t(ifng2[,193:216])),as.data.frame(t(ifng2[,217:240])),as.data.frame(t(ifng2[,241:264])))
head(abcd[,1:11])

colnames(abcd)<-unique(data$Time)
colnames(abcd)<-c("T01","T02","T03","T04","T05","T06","T07","T08","T09","T10","T11")


anovaIFgn<-aov(lm( IFng~ Time+Challenge, data=data))
summary(anovaIFgn)
TukeyHSD(anovaIFgn)
listanova<-list()
listatukey<-list()
for (aa in colnames(data[c(3:16)])){
  anova<-aov(lm( get(paste0(aa))~ Time*Challenge, data=data))
  listanova[[aa]]<-anova
  print(anova)
  tukey<-TukeyHSD(anova)
  print(tukey)
  listatukey[[aa]]<-tukey
  
}
plot(listatukey[["IFng"]])
lapply(listatukey, write, "listatukeyconchallenge.txt", append=TRUE, ncolumns=10000)
lapply(listatukey, function(x) write.table( data.frame(x), 'listatukey.csv'  , append= T, sep=',' ))
for (aa in colnames(data[c(3:16)])){
  capture.output(aa, file = "listatukeyconchallenge.txt", append = TRUE)
  capture.output(listatukey[[aa]], file = "listatukeyconchallenge.txt", append = TRUE)
}
#####quitar outlyers###
library(dplyr)
conout<-data %>%
  filter(Animal_ID != "61506" | Animal_ID != "61521" | Animal_ID != "61526" )
listanova1<-list()
listatukey1<-list()
for (aa in colnames(conout[c(3:16)])){
  anova<-aov(lm( get(paste0(aa))~ Time*Challenge, data=conout))
  listanova1[[aa]]<-anova
  print(anova)
  tukey1<-TukeyHSD(anova)
  print(tukey1)
  listatukey1[[aa]]<-tukey1
}
for (aa in colnames(data[c(3:16)])){
  capture.output(aa, file = "listatukeysinoutliers.txt", append = TRUE)
  capture.output(listatukey[[aa]], file = "listatukeysinoutliers.txt", append = TRUE)
}

plotcnc<-list()
for (aa in colnames(conout[c(3:16)])){
  df2 <- data_summary(conout, varname=paste0(aa), 
                      groupnames=c("Time","Challenge"))
  head(df2)
  library(ggplot2)
  p<- ggplot(df2, aes(x=Time, y=get(paste0(aa)),group = Challenge, colour= Challenge))+ #facet_grid(cosa ~. ) + 
    labs(y = paste0(aa))+
    geom_line() +
    geom_point(color='darkblue')+
    geom_point()+
    geom_errorbar(aes(ymin=get(paste0(aa))-sd, ymax=get(paste0(aa))+sd), width=.2,
                  position=position_dodge(0.05))
  plotcnc[[aa]] = p
  print(p)
}

pdf("sinoutliers.pdf")
for (aa in colnames(data[c(3:16)])) {
  print(plotcnc[[aa]])
}
dev.off()
###### high vs low####
data$Challenge<-NULL
listanova<-list()
listatukey<-list()
for (aa in colnames(data[c(3:16)])){
  anova<-aov(lm( get(paste0(aa))~ Time*Response, data=data))
  listanova[[aa]]<-anova
  print(anova)
  tukey<-TukeyHSD(anova)
  print(tukey)
  listatukey[[aa]]<-tukey
  
}
for (aa in colnames(data[c(3:16)])){
  capture.output(aa, file = "listatukeyresponse.txt", append = TRUE)
  capture.output(listatukey[[aa]], file = "listatukeyresponse.txt", append = TRUE)
}
plotcnc<-list()
colnames(data$Response)<-"Breading value"
for (aa in colnames(data[c(3:16)])){
  df2 <- data_summary(data, varname=paste0(aa), 
                      groupnames=c("Time","Breeding.Value"))
  head(df2)
  library(ggplot2)
  p<- ggplot(df2, aes(x=Time, y=get(paste0(aa)),group = Breeding.Value, colour= Breeding.Value))+ #facet_grid(cosa ~. ) + 
    labs(y = paste0(aa))+
    geom_line() +
    geom_point(color='darkblue')+
    geom_point()+
    labs(fill = "Breeding value")+
    theme(axis.text.x = element_text(size=5))+
    geom_errorbar(aes(ymin=get(paste0(aa))-sd, ymax=get(paste0(aa))+sd), width=.2,
                  position=position_dodge(0.05))
  plotcnc[[aa]] = p
  print(p)
}
pdf("breedingvalue1.pdf")
for (aa in colnames(data[c(3:16)])) {
  print(plotcnc[[aa]])
}
dev.off()
