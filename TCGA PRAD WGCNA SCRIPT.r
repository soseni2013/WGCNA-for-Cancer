#rm(list=ls(all=TRUE))
#setwd("C:/Users/sosen/Documents")
#getwd()
#outputfigs = "C:/Users/sosen/Documents"
#traits<-read.csv("Pan_traits.csv", header = TRUE, row.names = 1)
#cleanDat<-read.csv("Pancancer_TGCA.csv", header = TRUE, row.names = 1, check.names = FALSE)
#head(cleanDat)
#head(traits)
#dim(cleanDat)
#dim(traits)
#rownames(traits)==colnames(cleanDat)
## Declare as outliers those samples which are more than sdout sd above the mean connectivity based on the chosen measure
#sdout=3

rm(list=ls(all=TRUE))
setwd("C:/Users/sosen/Documents")
getwd()
outputfigs = "C:/Users/sosen/Documents"
load("finaloutput_for_Saheed_06_06_2020.RData")
numericMeta2<-read.csv("NumericMeta2.csv", header = TRUE, row.names = 1)
head(numericMeta2)
dim(numericMeta2)
head(cleanDat)
dim(cleanDat)
rownames(numericMeta2)==colnames(cleanDat)
## Declare as outliers those samples which are more than sdout sd above the mean connectivity based on the chosen measure
sdout=3

library(WGCNA)
library(NMF)
library(igraph)
library(ggplot2)
library(RColorBrewer)
library(Cairo)
library("doParallel")
library("biomaRt")
library("NMF")
library("plotly")
library("stringr")
library(boot)
install.packages("callr")
packageVersion("callr")
library("callr")


if(!exists("numericMeta")) numericMeta<-NumericMeta2
head(numericMeta)
dim (numericMeta)
for (repeated in 1:5) {
  normadj <- (0.5+0.5*bicor(cleanDat,use="pairwise.complete.obs")^2)
  ## Calculate connectivity
  netsummary <- fundamentalNetworkConcepts(normadj)
  ku <- netsummary$Connectivity
  z.ku <- ku-(mean(ku))/sqrt(var(ku))
  ## Declare as outliers those samples which are more than sdout sd above the mean connectivity based on the chosen measure
  outliers <- (z.ku > mean(z.ku)+sdout*sd(z.ku))|(z.ku < mean(z.ku)-sdout*sd(z.ku))
  print(paste0("There are ",sum(outliers)," outlier samples based on a bicor distance sample network connectivity standard deviation above ",sdout,".  [Round ",repeated,"]"))
  targets.All=numericMeta
  cleanDat <- cleanDat[,!outliers] 
  numericMeta <- targets <- targets.All[!outliers,]
  outliers.All<-c(outliers.All,outliers)
} #repeat 5 times
print(paste0("There are ",sum(outliers.All)," total outlier samples removed in ",repeated," iterations:"))
names(which(outliers.All))
outliersRemoved<-names(which(outliers.All))
dim(outliersRemoved)
dim(outliers.All)
LThalfSamples<-length(colnames(cleanDat))/2
LThalfSamples<-LThalfSamples - if ((length(colnames(cleanDat)) %% 2)==1) { 0.5 } else { 1.0 }
IndexHighMissing<-rowsRemoved<-zeroVarRows<-vector()
temp2<-as.data.frame(cleanDat[which(rowSums(as.matrix(is.na(cleanDat)))>LThalfSamples),])
if (ncol(temp2)==1) {
  temp2<-t(temp2)
  rownames(temp2)=rownames(cleanDat)[which(rowSums(as.matrix(is.na(cleanDat)))>LThalfSamples)]
}
if (nrow(temp2)>0) { IndexHighMissing=which(rowSums(as.matrix(is.na(cleanDat)))>LThalfSamples); rowsRemoved<-rownames(cleanDat)[IndexHighMissing]; cleanDat<-cleanDat[-IndexHighMissing,]; }
dim(cleanDat)
cleanDat.unreg<-cleanDat
dim(numericMeta)
library("doParallel")
parallelThreads=15
clusterLocal <- makeCluster(c(rep("localhost",parallelThreads)),type="SOCK")
registerDoParallel(clusterLocal)
allowWGCNAThreads() 
powers <- seq(4,12,by=1)
sft <- pickSoftThreshold(t(cleanDat),blockSize=nrow(cleanDat)+1000,   #always calculate power within a single block (blockSize > # of rows in cleanDat)
                         powerVector=powers,
                         corFnc="bicor",networkType="signed")
jpeg(file="savingsah_plottablesft.jpeg")
tableSFT<-sft[[2]]
plot(tableSFT[,1],tableSFT[,2],xlab="Power (Beta)",ylab="SFT R?")
dev.off()
tableSFT<-sft[[2]]
#plot(tableSFT[,1],tableSFT[,2],xlab="Power (Beta)",ylab="SFT R?") 
#dev.off()
plot(tableSFT[,1],tableSFT[,2],xlab="Power (Beta)",ylab="SFT R?") 
head(sft)
head(tableSFT)
write.csv(sft, file = "sft.csv")
write.csv(tableSFT, file = "tableSFT.csv")
powers <- seq(9,20,by=0.5)
sft <- pickSoftThreshold(t(cleanDat),blockSize=nrow(cleanDat)+1000,   #always calculate power within a single block (blockSize > # of rows in cleanDat)
                         powerVector=powers,
                         corFnc="bicor",networkType="signed")
tableSFT<-sft[[2]]
dim(sft)
head(sft)
head(tableSFT)
write.csv(sft, file = "sft1.csv")
write.csv(tableSFT, file = "tableSFT1.csv")
#power=8
power=10
net <- blockwiseModules(t(cleanDat),power=power,deepSplit=2,minModuleSize=100,
                        mergeCutHeight=0.15,TOMdenom="mean", #detectCutHeight=0.9999,                        #TOMdenom="mean" may get more small modules here.
                        corType="bicor",networkType="signed",pamStage=TRUE,pamRespectsDendro=TRUE,
                        verbose=3,saveTOMs=FALSE,maxBlockSize=nrow(cleanDat)+1000,reassignThresh=0.05)
nModules<-length(table(net$colors))-1
modules<-cbind(colnames(as.matrix(table(net$colors))),table(net$colors))
orderedModules<-cbind(Mnum=paste("M",seq(1:nModules),sep=""),Color=labels2colors(c(1:nModules)))
modules<-modules[match(as.character(orderedModules[,2]),rownames(modules)),]
as.data.frame(cbind(orderedModules,Size=modules))
net.ds2<-net
FileBaseName=paste0("MyNetworkDescription_power_",power,"_PAMstageTRUE")
CairoPDF(file="3Sa.GlobalNetworkPlots-FileBaseName.pdf",width=12,height=10.5)
MEs<-tmpMEs<-data.frame()
MEList = moduleEigengenes(t(cleanDat), colors = net$colors)
MEs = orderMEs(MEList$eigengenes)
colnames(MEs)<-gsub("ME","",colnames(MEs)) #let's be consistent in case prefix was added, remove it.
rownames(MEs)<-rownames(numericMeta)
numericIndices<-unique(c( which(!is.na(apply(numericMeta,2,function(x) sum(as.numeric(x))))), which(!(apply(numericMeta,2,function(x) sum(as.numeric(x),na.rm=T)))==0) ))
geneSignificance <- cor(sapply(numericMeta[,numericIndices],as.numeric),t(cleanDat),use="pairwise.complete.obs")
rownames(geneSignificance) <- colnames(numericMeta)[numericIndices]
geneSigColors <- t(numbers2colors(t(geneSignificance),,signed=TRUE,lim=c(-1,1),naColor="black"))
rownames(geneSigColors) <- colnames(numericMeta)[numericIndices]
plotDendroAndColors(dendro=net$dendrograms[[1]],
                    colors=t(rbind(net$colors,geneSigColors)),
                    cex.dendroLabels=1.2,addGuide=TRUE,
                    dendroLabels=FALSE,
                    groupLabels=c("Module Colors",colnames(numericMeta)[numericIndices]))
head(MEs)
tmpMEs <- MEs #net$MEs
colnames(tmpMEs) <- paste("ME",colnames(MEs),sep="")
MEs[,"grey"] <- NULL
tmpMEs[,"MEgrey"] <- NULL
plotEigengeneNetworks(tmpMEs, "Eigengene Network", marHeatmap = c(3,4,2,2), marDendro = c(0,4,2,0),plotDendrograms = TRUE, xLabelsAngle = 90,heatmapColors=blueWhiteRed(50))
head(numericMeta)
Grouping<-numericMeta$PFS
Grouping[numericMeta$PFS==1]<-"Progression"  #only necessary if Group was numerically encoded; does nothing if the Grouping vector has no numeric values
Grouping[numericMeta$PFS==0]<-"Non_progression"
#Grouping[numericMeta$Group==3]<-"Age.Dx"  #there is one of these lines for each numeric value set in the Group column of the traits.csv file
#Grouping[numericMeta$Group==4]<-"Stg4"
head(numericMeta)
#regvars <- data.frame(as.factor( numericMeta$PROGRESSION_FREE_STATUS ), as.numeric(numericMeta$TUMOR_TYPE))
regvars <- data.frame(as.factor( numericMeta$PFS), as.numeric(numericMeta$CANCER_RELAPSE))
colnames(regvars) <- c("PFS","CANCER_RELAPSE")
#colnames(regvars) <- c("PROGRESSION_FREE_STATUS","TUMOR_TYPE")
lm1 <- lm(data.matrix(MEs)~PFS,data=regvars)
pvec <- rep(NA,ncol(MEs))
for (i in 1:ncol(MEs)) {
  f <- summary(lm1)[[i]]$fstatistic ## Get F statistics
  pvec[i] <- pf(f[1],f[2],f[3],lower.tail=F) ## Get the p-value corresponding to the whole model
}
names(pvec) <- colnames(MEs)
kMEdat <- signedKME(t(cleanDat), tmpMEs, corFnc="bicor")
library(RColorBrewer)
MEcors <- bicorAndPvalue(MEs,numericMeta[,numericIndices])
moduleTraitCor <- MEcors$bicor
moduleTraitPvalue <- MEcors$p
textMatrix = apply(moduleTraitCor,2,function(x) signif(x, 2))
par(mfrow=c(1,1))
par(mar = c(6, 8.5, 3, 3));
cexy <- if(nModules>75) { 0.8 } else { 1 }
colvec <- rep("white",1500)
colvec[1:500] <- colorRampPalette(rev(brewer.pal(8,"BuPu")[2:8]))(500)
colvec[501:1000]<-colorRampPalette(c("white",brewer.pal(8,"BuPu")[2]))(3)[2] #interpolated color for 0.05-0.1 p
labeledHeatmap(Matrix = apply(moduleTraitPvalue,2,as.numeric),
               xLabels = colnames(numericMeta)[numericIndices],
               yLabels = paste0("ME",names(MEs)),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = colvec,
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               cex.lab.y= cexy,
               zlim = c(0,0.15),
               main = paste("Module-trait relationships\n bicor r-value shown as text\nHeatmap scale: Student correlation p value"),
               cex.main=0.8)
numericMetaCustom<-numericMeta[,numericIndices]
MEcors <- bicorAndPvalue(MEs,numericMetaCustom)
moduleTraitCor <- MEcors$bicor
moduleTraitPvalue <- MEcors$p
moduleTraitPvalue<-signif(moduleTraitPvalue, 1)
moduleTraitPvalue[moduleTraitPvalue > as.numeric(0.05)]<-as.character("")
textMatrix = moduleTraitPvalue; #paste(signif(moduleTraitCor, 2), " / (", moduleTraitPvalue, ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
#textMatrix = gsub("()", "", textMatrix,fixed=TRUE)
labelMat<-matrix(nrow=(length(names(MEs))), ncol=2,data=c(rep(1:(length(names(MEs)))),labels2colors(1:(length(names(MEs))))))
labelMat<-labelMat[match(names(MEs),labelMat[,2]),]
for (i in 1:(length(names(MEs)))) { labelMat[i,1]<-paste("M",labelMat[i,1],sep="") }
for (i in 1:length(names(MEs))) { labelMat[i,2]<-paste("ME",labelMat[i,2],sep="") }
xlabAngle <- if(nModules>75) { 90 } else { 45 }
par(mar=c(16, 12, 3, 3) )
par(mfrow=c(1,1))
bw<-colorRampPalette(c("#0058CC", "white"))
wr<-colorRampPalette(c("white", "#CC3300"))
colvec<-c(bw(50),wr(50))
labeledHeatmap(Matrix = t(moduleTraitCor)[,],
               yLabels = colnames(numericMetaCustom),
               xLabels = labelMat[,2],
               xSymbols = labelMat[,1],
               xColorLabels=TRUE,
               colors = colvec,
               textMatrix = t(textMatrix)[,],
               setStdMargins = FALSE,
               cex.text = 0.5,
               cex.lab.x = cexy,
               xLabelsAngle = xlabAngle,
               verticalSeparator.x=c(rep(c(1:length(colnames(MEs))),as.numeric(ncol(MEs)))),
               verticalSeparator.col = 1,
               verticalSeparator.lty = 1,
               verticalSeparator.lwd = 1,
               verticalSeparator.ext = 0,
               horizontalSeparator.y=c(rep(c(1:ncol(numericMetaCustom)),ncol(numericMetaCustom))),
               horizontalSeparator.col = 1,
               horizontalSeparator.lty = 1,
               horizontalSeparator.lwd = 1,
               horizontalSeparator.ext = 0,
               zlim = c(-1,1),
               main = "Module-trait Relationships\n Heatmap scale: signed bicor r-value", # \n (Signif. p-values shown as text)"),
               cex.main=0.8)
toplot <- MEs
colnames(toplot) <- colnames(MEs)
rownames(toplot) <- rownames(MEs)
toplot <- t(toplot)
pvec <- pvec[match(names(pvec),rownames(toplot))]
#rownames(toplot) <- paste(rownames(toplot),"\np = ",signif(pvec,2),sep="")
rownames(toplot) <- paste(orderedModules[match(colnames(MEs),orderedModules[,2]),1]," ",rownames(toplot),"  |  K-W p=",signif(pvec,2),sep="")
# add any traits of interest you want to be in the legend
#TUMOR_TYPE=as.numeric(numericMeta$tumor_type)
#TUMOR_TYPE[tumor_type==2]<-"F"
#TUMOR_TYPE[tumor_type==1]<-"M"
#metdat=data.frame(PROGRESSION_FREE_STATUS=Grouping, TUMOR_TYPE=tumor_type)
CANCER_RELAPSE=as.numeric(numericMeta$cancer_relapse)
CANCER_RELAPSE[cancer_relapse==1]<-"F"
CANCER_RELAPSE[cancer_relapse==0]<-"M"
metdat=data.frame(PFS=Grouping, CANCER_RELAPSE=cancer_relapse)
# set colors for the traits in the legend
heatmapLegendColors=list('Group'=c("midnightblue","dodgerblue","goldenrod","red"), #,"seagreen3","hotpink","purple"),
                         
                         #'CANCER_RELAPSE'=c("pink","dodgerblue"), #F, M
                         'CANCER_RELAPSE'=c("pink","dodgerblue"), #F, M
                         'Modules'=sort(colnames(MEs)))
library(NMF)
par(mfrow=c(1,1))
aheatmap(x=toplot, ## Numeric Matrix
         main="Plot of Eigengene-Trait Relationships - SAMPLES IN ORIGINAL, e.g. BATCH OR REGION ORDER",
         annCol=metdat,
         annRow=data.frame(Modules=colnames(MEs)),
         annColors=heatmapLegendColors,
         border=list(matrix = TRUE),
         scale="row",
         distfun="correlation",hclustfun="average", ## Clustering options
         cexRow=0.8, ## Character sizes
         cexCol=0.8,
         col=blueWhiteRed(100), ## Color map scheme
         treeheight=80,
         Rowv=TRUE, Colv=NA) ## Do not cluster columns - keep given order
aheatmap(x=toplot, ## Numeric Matrix
         main="Plot of Eigengene-Trait Relationships - SAMPLES CLUSTERED",
         annCol=metdat,
         annRow=data.frame(Modules=colnames(MEs)),
         annColors=heatmapLegendColors,
         border=list(matrix = TRUE),
         scale="row",
         distfun="correlation",hclustfun="average", ## Clustering options
         cexRow=0.8, ## Character sizes
         cexCol=0.8,
         col=blueWhiteRed(100), ## Color map scheme
         treeheight=80,
         Rowv=TRUE,Colv=TRUE) ## Cluster columns
colnames(numericMeta)[numericIndices]
table(Grouping)
par(mfrow=c(4,5)) #rows,columns on each page (try to choose to keep all plots for each module on one page or row(s), without squeezing too much in)
par(mar=c(5,6,4,2))
for (i in 1:nrow(toplot)) {
  boxplot(toplot[i,]~factor(Grouping,c("Progression","Non_progression")),col=colnames(MEs)[i],ylab="Eigengene Value",main=rownames(toplot)[i],xlab=NULL,las=2)  #you choose the order of groups for boxplots
  boxplot(toplot[i,]~factor(cancer_relapse,c("Relapse","No_relapse")),col=colnames(MEs)[i],ylab="Eigengene Value",main=rownames(toplot)[i],xlab=NULL,las=2)  #you choose the order of groups for boxplots
  verboseScatterplot(x=numericMeta[,"T_STAGE"],y=toplot[i,],xlab="Pathology T-Stage (2-7)",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=colnames(MEs)[i],pch=19)
  verboseScatterplot(x=numericMeta[,"N_STAGE"],y=toplot[i,],xlab="Pathology N-Stage (0=no; 1=yes)",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=colnames(MEs)[i],pch=19)
  verboseScatterplot(x=numericMeta[,"RADIOTHERAPY"],y=toplot[i,],xlab="Prior Radiotherapy (0=no; 1=yes)",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=colnames(MEs)[i],pch=19)
  verboseScatterplot(x=numericMeta[,"GLEASON"],y=toplot[i,],xlab="Gleason Score (2-3)",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=colnames(MEs)[i],pch=19)
  verboseScatterplot(x=numericMeta[,"RACE"],y=toplot[i,],xlab="Race (1=European-American; 2=African-American, 3=Asian-American)",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=colnames(MEs)[i],pch=19)
  verboseScatterplot(x=numericMeta[,"CANCER_RELAPSE"],y=toplot[i,],xlab="Cancer Relapse (0=no; 1=yes)",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=colnames(MEs)[i],pch=19)
  #verboseScatterplot(x=numericMeta[,"Dead"],y=toplot[i,],xlab="Deceased (0=no; 1=yes)",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=colnames(MEs)[i],pch=19)
  #verboseScatterplot(x=numericMeta[,"days.to.death"],y=toplot[i,],xlab="Days to Death",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=colnames(MEs)[i],pch=19)
  #verboseScatterplot(x=numericMeta[,"gender"],y=toplot[i,],xlab="Gender (1=male; 2=female)",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=colnames(MEs)[i],pch=19)
  #verboseScatterplot(x=numericMeta[,"TUMOR_TYPE"],y=toplot[i,],xlab="TUMOR_TYPE at Sampling (1-2)",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=colnames(MEs)[i],pch=19)
  #verboseScatterplot(x=numericMeta[,"in_treatment"],y=toplot[i,],xlab="Number of Days in Treatment",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=colnames(MEs)[i],pch=19)
  #verboseScatterplot(x=numericMeta[,"Melphalan"],y=toplot[i,],xlab="Treatment: Melphalan",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=colnames(MEs)[i],pch=19)
  #verboseScatterplot(x=numericMeta[,"Melphalan"],y=toplot[i,],xlab="Treatment: Melphalan",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=colnames(MEs)[i],pch=19)
  #verboseScatterplot(x=numericMeta[,"Melphalan"],y=toplot[i,],xlab="Treatment: Melphalan",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=colnames(MEs)[i],pch=19)
  #verboseScatterplot(x=numericMeta[,"Melphalan"],y=toplot[i,],xlab="Treatment: Melphalan",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=colnames(MEs)[i],pch=19)
  #verboseScatterplot(x=numericMeta[,"Melphalan"],y=toplot[i,],xlab="Treatment: Melphalan",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=colnames(MEs)[i],pch=19)
  #verboseScatterplot(x=numericMeta[,"NoTreatment"],y=toplot[i,],xlab="Treatment: None (--)",ylab="Eigengene",abline=TRUE,cex.axis=1,cex.lab=1,cex=1,col=colnames(MEs)[i],pch=19)
  for (i in 1:7) frame()
}
dev.off()

orderedModulesWithGrey=rbind(c("M0","grey"),orderedModules)
kMEtableSortVector<-apply( as.data.frame(cbind(net$colors,kMEdat)),1,function(x) if(!x[1]=="grey") { paste0(paste(orderedModulesWithGrey[match(x[1],orderedModulesWithGrey[,2]),],collapse=" "),"|",round(as.numeric(x[which(colnames(kMEdat)==paste0("kME",x[1]))+1]),4)) } else { paste0("grey|AllKmeAvg:",round(mean(as.numeric(x[-1],na.rm=TRUE)),4)) } ) 
kMEtable=cbind(c(1:nrow(cleanDat)),rownames(cleanDat),net$colors,kMEdat,kMEtableSortVector)[order(kMEtableSortVector,decreasing=TRUE),]
