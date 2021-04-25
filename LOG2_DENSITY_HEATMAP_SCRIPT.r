if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("limma")

library(limma)
DEnShino<-matplot(SU2C_2019_RNASeq)

NormShino<- normalizeCyclicLoess(SU2C_2019_RNASeq)

# to remove or hide first gene Id column 
SU2C_2019_RNASeq$Hugo_Symbol=NULL

# to log2 transform data
Logshino<-log2(SU2C_2019_RNASeq)

#to normalize data
Norm2shino<- normalizeQuantiles(Logshino)

is.na(Logshino)<- sapply(Logshino, is. infinite) # replace inf with NA, problem solved
#it replaces both negative and positive infinity

#to gerid of NAS
na.omit(Logshino)

Den1shino<-matplot(Logshino)
DEn2Shino<-matplot(Norm2shino)


SU2C_2019_RNASeq$Hugo_Symbol= NULL
NormCGA<- normalizeBetweenArrays(Pancancer_TCGA_RNASeq_2018)
LogTCGA<-

Pancancer_TCGA_RNASeq_2018$Hugo_Symbol=NULL
LogShinoTCGA_2018<-log2(Pancancer_TCGA_RNASeq_2018)
NormShinoTCGA<-normalizeQuantiles(LogShinoTCGA_2018)
rlang::last_error()


write.csv(Logshino,"shinoSU2Clogtransformed.rdata")
write.csv(Norm2shino,"shinoSUC2CNormalized.rdata")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq")

SURVIVAL ANALYSIS

##################################Density Plot for IRAK1 in Datasets####################################################################################
#To make Density plots for IRAK1 in TCGA (INFLAMMAHEATMAP_TCGA has been log2 transformed and normalized already)
#Kernel Density Plot
getwd()
setwd()
library(ggplot2)
library(dplyr)
library(hrbrthemes)
library(viridis)
library(tidyr)

p1 <- ggplot(data=diamonds, aes(x=IRAKS_DENSITY_PLOT_TCGA, group=cut, fill=cut)) +
  geom_density(adjust=1.5) +
  theme_ipsum()
#Stack Density plot
p <- ggplot(data=IRAKS_DENSITY_PLOT_TCGA, aes(x= "IRAK1", group=cut, fill=cut)) +
  geom_density(adjust=1.5, position="fill") + polygon(d, col="red", border="black")
  theme_ipsum()
########Redo Density Plot########worked really well############
  head(IRAKS_DENSITY_PLOT_TCGA)
  density(IRAKS_DENSITY_PLOT_TCGA$IRAK1)
  plot(density(IRAKS_DENSITY_PLOT_TCGA$IRAK1))
  library(ggplot2)
  ggplot(IRAKS_DENSITY_PLOT_TCGA, aes(IRAK1)) + geom_density()
  ggplot(IRAKS_DENSITY_PLOT_TCGA, aes(IRAK1)) + geom_density() + labs(x= "IRAK Log2(FPKM+1)", y= "Distribution of Gene Expression")
  ggplot(IRAKS_DENSITY_PLOT_TCGA, aes(IRAK1)) + geom_density(fill="pink") + labs(x= "IRAK1 Log2(FPKM+1)", y= "Distribution of Gene Expression")
  ggplot(IRAKS_DENSITY_PLOT_TCGA, aes(IRAK1)) + geom_density(fill="pink", alpha = 0.5) + labs(x= "IRAK1 Log2(FPKM+1)", y= "Distribution of Gene Expression")
  ggplot(IRAKS_DENSITY_PLOT_TCGA, aes(IRAK1)) + geom_density(fill="pink", alpha = 0.5) + labs(x= "IRAK1 Log2(FPKM+1)", y= "Distribution of Gene Expression") + geom_vline(xintercept = 12.18, size=1.5, color="red", linetype = "dashed")
  
  ##############This works for IRAK1 MAGENTA GENE NEIGBHORS HEATMAP################
  library(limma)
  library(ggplot2)
  library(dplyr)
  library(plotly)
  library(heatmap3)
  
  M<-data.matrix(IRAK1_MAGENTA_NEIGBHORS)
  heatmap3(M, scale = "column", balanceColor=TRUE, main = "IRAK Neigbhoring Genes Co-expression Profiles (Magenta Module) in Localized PCa Patients", Rowv = TRUE, Colv = TRUE, cexRow = 1.0, margins = c(8,8)) 
  
  heatmap3(M, scale = "column", balanceColor=TRUE, main = "Co-expression Profiles of IRAK1 Neigbhors within the Magenta Module", Rowv = TRUE, Colv = TRUE, cexRow = 1.0, margins = c(8,8)) 
  
  heatmap3(M, scale = "column", balanceColor=TRUE, main = "", Rowv = TRUE, Colv = TRUE, cexRow = 1.0, margins = c(8,8)) 
  
  
  ##############This works for IRAK1-4 HEATMAP################
  library(limma)
  library(ggplot2)
  library(dplyr)
  library(plotly)
  library(heatmap3)
  
  X<-data.matrix(IRAKS_DENSITY_PLOT_TCGA)
  heatmap3(X, scale = "row", balanceColor=TRUE, main = "Expression Profiles of IRAK Genes in Localized PCa Patients", Rowv = TRUE, Colv = TRUE, cexRow = 1.0, margins = c(8,8)) 
  
X<-data.matrix(IRAKS_DENSITY_PLOT_TCGA)
heatmap3(X, scale = "row", balanceColor=TRUE, main = "", Rowv = TRUE, Colv = TRUE, cexRow = 1.0, margins = c(8,8)) 
                    
#NEPC IRAK1-4 only Heatmap

Z<-data.matrix(IRAKS_DENSITY_PLOT_NEPC)
heatmap3(Z, scale = "row", balanceColor=TRUE, main = "Expression Profile of IRAK Genes in NEPC Patients", Rowv = TRUE, Colv = TRUE, cexRow = 1.0, margins = c(8,8)) 

Z<-data.matrix(IRAKS_DENSITY_PLOT_NEPC)
heatmap3(Z, scale = "row", balanceColor=TRUE, main = "", Rowv = TRUE, Colv = TRUE, cexRow = 1.0, margins = c(8,8)) 

#CRPC IRAK1-4 only Heatmap
y<-data.matrix(IRAK1.4_SU2C_CRPC)
heatmap3(y, scale = "row", balanceColor=TRUE, main = "Expression Profile of IRAK Genes in Metastatic/CRPC Patients", Rowv = TRUE, Colv = TRUE, cexRow = 1.0, margins = c(8,8)) 

y<-data.matrix(IRAK1.4_SU2C_CRPC)
heatmap3(y, scale = "row", balanceColor=TRUE, main = "", Rowv = TRUE, Colv = TRUE, cexRow = 1.0, margins = c(8,8)) 

########################basic Density Plot#########worked####
 library(ggplot2)
  head(IRAKS_DENSITY_PLOT_TCGA) 
 density(IRAKS_DENSITY_PLOT_TCGA$IRAK1)
 ggplot(IRAKS_DENSITY_PLOT_TCGA$IRAK1) + geom_density(fill= "red", alpha = 0.5) + labs (x="", y= "Distribution of IRAK1 Gene Expression")    
#####################################################################################      
        density(IRAKS_DENSITY_PLOT_TCGA$IRAK2))
 library(ggplot2)
 ggplot(IRAKS_DENSITY_PLOT_TCGA, aes(IRAKS_DENSITY_PLOT_TCGA)) + geom_density()
 
 
 
 ggplot(x=newirak, aes(fill=genenames, color=genenames)) + 
   geom_density(alpha = 0.5)+
     scale_fill_viridis(discrete = T)+
     scale_color_viridis(discrete = T)+
     xlab("IRAKs Log2 FPKM+1")+
     ylab("Distribution of IRAK Gene Expression")
     
   
   aes(IRAKS_DENSITY_PLOT_TCGA$IRAK1,IRAKS_DENSITY_PLOT_TCGA$IRAK2))+ 
                
 
 ggplot(IRAKS_DENSITY_PLOT_TCGA, aes(IRAK1, fill= as.factor((IRAKS_DENSITY_PLOT_TCGA$IRAK2) + geom_density(alpha = 0.5) + labs (x="IRAK1 Log2 (FPKM+1)", y= "Distribution of IRAK1 Gene Expression"), geom_vline(xintercept = 12.18, size=2, color="red", linetype="dashed")
 plot(density(IRAKS_DENSITY_PLOT_TCGA$IRAK3))
ggplot(IRAKS_DENSITY_PLOT_TCGA, aes(IRAK3)) + geom_density(fill= "red", borders = "black", alpha = 0.5) + labs (x="IRAK3 Log2 (FPKM+1)", y= "Distribution of IRAK1 Gene Expression")
density.default(x = IRAKS_DENSITY_PLOT_TCGA$IRAK3)
density.default(x = IRAKS_DENSITY_PLOT_TCGA$IRAK1)
density.default(x = IRAKS_DENSITY_PLOT_TCGA$IRAK4)
IRAK1_Dplot <- density(TCGA_IRAK1_Density)
plot(IRAK1_Dplot, main="RNA log2 (FPKM+1)")
polygon(d, col="red", border="black")
plot(IRAK1_Dplot, cex = 0.1)

TCGA_IRAK1_Density <- as.numeric(INFLAMMAHEATMAP_TCGA$IRAK1)
IRAK1_Dplot <- geom_density(TCGA_IRAK1_Density)
plot(IRAK1_Dplot, cex = 0.1)

# Make the histogram
data <- as.numeric(IRAKS_DENSITY_PLOT_TCGA)
data %>%
  filter(IRAK1)
  ggplot(aes(x=IRAK1)) +
  geom_density(fill="#69b3a2", color="#e9ecef", alpha=0.8) +
  ggtitle("Density Plot of IRAK1 Expression (Log2 (FPKM+1)") +
  theme_ipsum()



##################################################################################################################################
getwd()
ShinoStart<-Pancancer_TCGA_RNA_Seq_v2_expression_2018_
ShinoStart$Entrez_Gene_Id=NULL
ShinoStart$Hugo_Symbol=NULL
library(limma)
NormShino<-normalizeBetweenArrays(ShinoStart)
View(Log2Shino)<-log2(ShinoStart+1) 
View(LogShino)<- as.data.frame(Log2Shino, na.rm = TRUE)
View()

Log2TCGA_2018 <- log2(Pancancer_TCGA_RNASeq_2018+1)
write.csv(TCGA_log2TCGA_2018, "TCGADATA_logged_Data")

ShinoStart2<-ShinoStart+1
log2Shino<-log2(ShinoStart2)
write.csv(log2Shino,"shino.rdata")
library(quantro)
?cbind()
View(log2Shino)
library(ggplot2)
library()
Denshino<-matplot(shino)




write.csv(logShinoTCGA_2018,"TCGA_2018_Logged2")

library(ggplot2)
HeaTCGAShino <-data.matrix(LogTCGAShino) 
heatmap(HeaTCGAShino, main = "Genes in Local PCa Patients", Rowv = TRUE, Colv = NA, cexRow = 0.5, margins = c(10, 12))
#Make heatmap using ggplot2
?heatmap
heatmap.2(HeaTCGAShino, main = "Sample", trace = "none", margins =  c(10,12))
x<-data.matrix(INFLAMMAHEATMAP)          
?heatmap
heatmap(x, main = "Sample",trace = "none", margins = c(10,12), cexRow = 0.5, Colv=NULL, Rowv=NULL, keep.dendro = FALSE, key = TRUE)
?colorRampPalette
x<-data.matrix(INFLAMMAHEATMAP)
yb <- col = colorRampPalette(c("darkgreen","black","red"))
heatmap(x, main = "Sample", Rowv = NULL, showColDendro = FALSE, showRowDendro = TRUE, cexRow = 0.5, margins = c(10, 12))


X<-data.matrix(INFLAMMAHEATMAP)
heatmap(X, main = "Inflammation Genes in Local PCa Patient", Rowv = TRUE, Colv = NA, cexRow = 0.5, margins = c(10, 12))
my_color<-colorRampPalette(c("red","yellow", "green")) (n =299)
heatmap(X, col = my_color (299), main = "Inflammation Genes in Local PCa Patients", Rowv = TRUE, Colv = TRUE, cexRow = 1, margins = c(10, 12), colorkey = TRUE)
legend(x = "topleft", legend = c("min", "ave", "max"), fill = colorRampPalette(c("yellow", "red"))(3)

library(reshape2)
data_melt<-melt(X)
library(ggplot2)
ggp <- ggplot(data_melt, aes(Var1, Var2) + geom_tile(aes(fill = value))
ggp + scale_fill_gradient(low = "green", high = "black"),space = "Lab" guide = "colourbar", aesthetics = "fill") 


library(plotly)

plot_ly(z = , type = "heatmap")
plot_ly(z = X, colorscale = "rgb", type = "heatmap")

#To creat Heatmap in R for Log2 Transformed TCGA PRAD_2018 dataset
#I Called Inflammation-associated genes (30) by creating new csv file with name: INFLAMMAHEATMAP
#Convert to matrix and run heatmap (see code below). #obviously the dendograms could be removed (Rowv or Colv = FALSE)
X<-data.matrix(INFLAMMAHEATMAP)
heatmap(X, main = "Inflammation Genes in Local PCa Patients", Rowv = TRUE, Colv = TRUE, cexRow = 1, margins = c(10, 12))

#To change the default color of the heatmap to my_colors. I did the following
X<-data.matrix(INFLAMMAHEATMAP)
heatmap(X, main = "Inflammation Genes in Local PCa Patients", notecol="black", density.info="none", trace = "none", Rowv = TRUE, Colv = TRUE, cexRow = 1, margins = c(10, 10))

X<-data.matrix(INFLAMMAHEATMAP)
my_palette <- colorRampPalette(c("red", "yellow", "green")) (n=299)
col_breaks = c(seq(-1, 0, length=100), seq(0, 5, length=100), seq(5-20, length=100)) 
heatmap(X, main = "Inflammation Genes in Local PCa Patients", notecol="black", density.info="none", trace="none", Rowv = TRUE, Colv = TRUE, cexRow = 1, margins = c(10, 12))

#To use the heatmap3 for visualizing gene expression
library(heatmap3)
X<-data.matrix(INFLAMMAHEATMAP)
my_palette <- colorRampPalette(c("navy", "White", "firebrick"))(1024)
heatmap3(X, main = "Inflammation Genes in Local PCa Patients", Rowv = TRUE, Colv = TRUE, cexRow = 1, margins = c(10, 12), col=my_palette) 


#SU2C_METASTATIC CRPC DATA (LOG TRANSFORMATION AND RECALL)
getwd()
library(limma)
SU2C_MET<- SU2C_2019_RNASeq+1
logSU2C_MET<-log2(SU2C_MET)
View(logSU2C_MET)
write.csv(logSU2C_MET, "OSENI_R_STUDIO_DATA.rdata")

##########################################HEATMAPS FOR DATATSETS#############################################################################################################
#THIS REALLY WORKED FOR THE HEATMAP_TCGA LOCAL PCA...SO HAPPY
X<-data.matrix(INFLAMMAHEATMAP_TCGA)
heatmap3(X, scale = "column", balanceColor=TRUE, main = "Inflammation-associated Genes in Localized PCa Patients", Rowv = TRUE, Colv = TRUE, cexRow = 1, margins = c(8, 8)) 

#THIS REALLY WORKED FOR THE HEATMAP_SU2C...SO HAPPY
Y<-data.matrix(INFLAMMAHEATMAP_SU2C_MET)
heatmap3(Y, scale = "column", balanceColor=TRUE, main = "Inflammation-associated Genes in CRPC Patients", Rowv = TRUE, Colv = TRUE, cexRow = 1, margins = c(8, 8)) 

#NEPC DATASET (LOG TRANSFORMATION AND RECALL)
library(limma)
NEPC<- NEPC_RNASEQ+5
logNEPC<-log2(NEPC)
View(logNEPC)
write.csv(logNEPC,"NEPC_Logged2.rdata")


#THIS REALLY WORKED FOR THE HEATMAP_NEPC...SO HAPPY
Z<-data.matrix(INFLAMMAHEATMAP_NEPC)
heatmap3(Z, scale = "column", balanceColor=TRUE, main = "Inflammation-associated Genes in NEPC Patients", Rowv = TRUE, Colv = TRUE, cexRow = 1, margins = c(8, 8)) 
#########################################################################################################################################################################################

# For Differential Gene Expression with Deseq2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
browseVignettes("DESeq2")
library(DESeq2)


View(WGCNA_MODULES_TCGA)
library(limma)
install.packages("BiocManager")
BiocManager::install("WGCNA")
library(WGCNA)

library(readr)
WGCNA_ALL_MODULES_TCGA <- read_csv("WGCNA_ALL_MODULES_TCGA.csv")
View(WGCNA_ALL_MODULES_TCGA)

WGCNA_YY<- data.matrix(WGCNA_ALL_MODULES_TCGA)



(topGens, gene_modules, file_prefix = NULL, plotdims= c(9,9), KME_cut=0.75, fc_cut=log2(1.5),x_lim=NULL, y_lim=NULL, gene_labs=FALSE, x_cut=0, y_cut=0, point_order="random")

install.packages("dplyr")
install.packages("ggplot2")
install.packages("tidyr")
install.packages("readr")
install.packages("stringr")
install.packages("knitr")
install.packages("rmarkdown")


library(dplyr)
library(ggplot2)
library(tidyr)
library(readr)
library(stringr)
library(knitr)
library(rmarkdown)

#I want to use WGCNA on the SU2C Dataset_just trying my luck
CRPC2_MET<-REFINED.CRPC_SU2C_RNASEQ_DATASET.+1
logCRPC2_MET<-log2(CRPC2_MET)
View(logCRPC2_MET)
write.csv(logCRPC2_MET,"CRPC_MET_Logged2")
#Next is the Start of the WGCNA anaysis
library(WGCNA)
options(stringsAsFactors = FALSE)

datExpr0 = as.data.frame(t(logCRPC2_MET[, -c(1:8)]));
#Check for genes and samples with too many missing values
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
#Remove offending genes
datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
sampleTree = hclust(dist(datExpr0), method = "average");
sizeGrWindow(12,9)
sampleTree = hclust(dist(datExpr0), method = "average");
par(cex = 0.5);
par(mar = c(0,5,3,0))
plot(sampleTree, main = "Sample clustering to detect outliers in metastatic/CRPC patients dataset", sub="", xlab="", cex.lab = 2,
     cex.axis = 2, cex.main = 2)

CRPC_WGCNA = rownames(datExpr0);
traitRows = match(CRPC_WGCNA, REFINED.CRPC_SU2C_CLINICAL_TRAITS_DATASET.);
datTraits = REFINED.CRPC_SU2C_CLINICAL_TRAITS_DATASET.[traitRows, -1];
rownames(datTraits) = REFINED.CRPC_SU2C_CLINICAL_TRAITS_DATASET.[traitRows, 1];
collectGarbage();

# For NEPC Dataset
#I want to use WGCNA on the NEPC Dataset_just trying my luck

#Next is the Start of the WGCNA anaysis


NEPC_datExpr0 = as.data.frame(t(logNEPC[, -c(1:8)]));
#Check for genes and samples with too many missing values
gsg = goodSamplesGenes(NEPC_datExpr0, verbose = 3);
gsg$allOK
#Remove offending genes
NEPC_datExpr0 = NEPC_datExpr0[gsg$goodSamples, gsg$goodGenes]
sampleTree = hclust(dist(NEPC_datExpr0), method = "average");
sizeGrWindow(12,9)
sampleTree = hclust(dist(NEPC_datExpr0), method = "average");
par(cex = 0.5);
par(mar = c(0,5,3,0))
plot(sampleTree, main = "Sample clustering to detect outliers in NEPC patients dataset", sub="", xlab="", cex.lab = 2.0,
     cex.axis = 2.0, cex.main = 2)

#TCGA_WGCNA_SAMPLE CLUSTERING
library(limma)
view(logNEPC)
library(WGCNA)
options(stringsAsFactors = FALSE)
TCGA_datExpr0 = as.data.frame(t(LogTCGA[, -c(1:8)]));
#Check for genes and samples with too many missing values
gsg = goodSamplesGenes(TCGA_datExpr0, verbose = 3);
gsg$allOK
#Remove offending genes
TCGA_datExpr0 = TCGA_datExpr0[gsg$goodSamples, gsg$goodGenes]
sampleTree = hclust(dist(TCGA_datExpr0), method = "average");
sizeGrWindow(12,9)
sampleTree = hclust(dist(NEPC_datExpr0), method = "average");
par(cex = 0.5);
par(mar = c(0,5,3,0))
plot(sampleTree, main = "Sample clustering to detect outliers in NEPC patients dataset", sub="", xlab="", cex.lab = 2.0,
     cex.axis = 2.0, cex.main = 2)

#I want to use WGCNA on the SU2C INFLAMMATION Dataset_just trying my luck
CRPC2_MET<-REFINED.CRPC_SU2C_RNASEQ_DATASET.+1
logCRPC2_MET<-log2(CRPC2_MET)
View(logCRPC2_MET)
write.csv(logCRPC2_MET,"CRPC_MET_Logged2")
#Next is the Start of the WGCNA anaysis
library(WGCNA)
options(stringsAsFactors = FALSE)
library(limma)
datExpr0 = as.data.frame(t(INFLAMMAHEATMAP_SU2C_MET[, -c(1:8)]));
#Check for genes and samples with too many missing values
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
#Remove offending genes
datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
sampleTree = hclust(dist(datExpr0), method = "average");
sizeGrWindow(12,9)
sampleTree = hclust(dist(datExpr0), method = "average");
par(cex = 0.5);
par(mar = c(0,5,3,0))
plot(sampleTree, main = "Sample clustering to detect outliers in CRPC patients inflammation dataset", sub="", xlab="", cex.lab = 2,
     cex.axis = 2, cex.main = 2)
CRPC_Traits <- REFINED.CRPC_SU2C_CLINICAL_TRAITS_DATASET.
CRPC_WGCNA = rownames(datExpr0)
traitRows = match(CRPC_WGCNA,CRPC_Traits$PSA..ng.ml.)
datTraits = CRPC_Traits[traitRows, -1]
rownames(datTraits) = make.names(traitRows, 1, unique=TRUE)
rownames(datTraits) = CRPC_Traits[traitRows, 1]
collectGarbage()
sampleTree2 = hclust(dist(datExpr0), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE)
plotDendroAndColors(sampleTree2, traitColors,groupLabels = names(datTraits), main = "Sample dendrogram and trait heatmap")


#REVISED CODES FOR CRPC_WGCNA
#I want to use WGCNA on the SU2C Dataset_just trying my luck

library(limma)
CRPC2_MET<-REFINED.CRPC_SU2C_RNASEQ_DATASET.+1
logCRPCnew_MET<-log2(CRPC2_MET)
View(logCRPCnew_MET)
write.csv(logCRPCnew_MET,"CRPC_MET2_Logged2")
#Next is the Start of the WGCNA anaysis
library(WGCNA)
options(stringsAsFactors = FALSE)

datExpr0 = as.data.frame(t(logCRPCnew_MET[, -c(1:8)]));
#Check for genes and samples with too many missing values
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
#Remove offending genes
datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
sampleTree = hclust(dist(datExpr0), method = "average");
sizeGrWindow(12,9)
sampleTree = hclust(dist(datExpr0), method = "average");
par(cex = 0.5);
par(mar = c(0,5,3,0))
plot(sampleTree, main = "Sample clustering to detect outliers in metastatic/CRPC patients dataset", sub="", xlab="", cex.lab = 2,
     cex.axis = 2, cex.main = 2)

CRPC_WGCNA = rownames(datExpr0);
traitRows = match(CRPC_WGCNA,CRPC_Traits$SAMPLE_ID)
datTraits = CRPC_Traits[traitRows, -1]
rownames(datTraits) = CRPC_Traits[traitRows, 1]
#ANOTHER RUN
CRPC_WGCNA = rownames(datExpr0);
traitRows = match(CRPC_WGCNA,CRPC_Traits$SAMPLE_ID)
library(dplyr)
datTraits2 <- data.frame (CRPC_Traits)
distinct_datTraits <- dplyr:: distinct(datTraits2)
rownames(distinct_datTraits) = datTraits2[traitRows, 1]

  distinct_(datTraits2, keep_all = TRUE)
# getting the indices of duplicate elements in CRPC_Traits
dups <- which(duplicated(CRPC_Traits))
View(dups)

length(dups)
#prints the elements in CRPC-Traits that are occuring more than once
dupp <-unique(dups)
View(dupp)
rownames (datTraits) = CRPC_Traits(unique(traitRows,(incomparables=FALSE))
rownames(datTraits) = make.names(CRPC_Traits)
rownames(datTraits2) = datTraits2[traitRows, 1]
collectGarbage()
sampleTree2 = hclust(dist(datExpr0), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE)
plotDendroAndColors(sampleTree2, traitColors,groupLabels = names(datTraits), main = "Sample dendrogram and trait heatmap")


#Re WGCNA of CRPC Dataset After removing 4 duplicates (i.e. from 270 samples to 266 samples) 
library(limma)
CRPC2_MET<-DUBREMOVED_CRPC_SU2C_RNASEQ.+1
logCRPCnew_MET<-log2(CRPC2_MET)
View(logCRPCnew_MET)
write.csv(logCRPCnew_MET,"CRPC_MET_DUBREMOVED_Logged2")
#Next is the Start of the WGCNA anaysis
library(WGCNA)
options(stringsAsFactors = FALSE)

datExpr0 = as.data.frame(t(logCRPCnew_MET[, -c(1:8)]));
#Check for genes and samples with too many missing values
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK
#Remove offending genes
datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
sampleTree = hclust(dist(datExpr0), method = "average");
sizeGrWindow(12,9)
sampleTree = hclust(dist(datExpr0), method = "average");
par(cex = 0.5);
par(mar = c(0,5,3,0))
plot(sampleTree, main = "Sample clustering to detect outliers in metastatic/CRPC patients dataset", sub="", xlab="", cex.lab = 2,
     cex.axis = 2, cex.main = 2)
CRPC_WGCNA = rownames(datExpr0);
CRPC_Traits <-data.frame(DUBREMOVED_CRPC_SU2C_CLINICAL_TRAITS.)

traitRows = match(CRPC_WGCNA,CRPC_Traits$SAMPLE_ID)
datTraits = CRPC_Traits[traitRows, (na.rm = TRUE)]
rownames(datTraits) = CRPC_Traits[traitRows, 1]
collectGarbage()
sampleTree2 = hclust(dist(datExpr0), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE)
plotDendroAndColors(sampleTree2, traitColors,groupLabels = names(datTraits), main = "Sample dendrogram and trait heatmap")
                                   
##############################################################################################################################################################################
#Sample WGCNA and Cluster for TCGA PCa Data
library(readr)
Pancancer_TCGA_RNASeq_2018 <- read_csv("OSENI_R_TUTORIAL/Pancancer_TCGA_RNASeq_2018.csv", 
#To skip a column (Entrez_Gene_Id) in the dataset
col_types = cols(Entrez_Gene_Id = col_skip()))
View(Pancancer_TCGA_RNASeq_2018)
TCGA1_RNASEQ<-na.omit(Pancancer_TCGA_RNASeq_2018)
library(limma)
TCGA_RNASEQ<- as.numeric(TCGA1_RNASEQ +1)
TCGA_RNASEQ <- as.data.frame(TCGA1_RNASEQ +1)

TCGA_RNASEQ <- as.data.frame(Pancancer_TCGA_RNASeq_2018 +1, na.rm = TRUE)
log2TCGA_RNASEQ <- log2(TCGA_RNASEQ)
view(log2TCGA_RNASEQ)
write.csv(log2TCGA_RNASEQ,"log2TCGA_RNASEQ.RDATA")
#Nextis the WGCNA anaysis
library(WGCNA)
options(stringsAsFactors = FALSE)

datExpr0 = as.data.frame(t(log2TCGA_RNASEQ[, -c(1:8)]));

# I want to Check for genes and samples with too many missing values
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

# Next, I want to remove offending genes
datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]

#Next I want to cluster my samples to check for outliers
sampleTree = hclust(dist(datExpr0), method = "average");
sizeGrWindow(12,9)
sampleTree = hclust(dist(datExpr0), method = "average");
par(cex = 0.5);
par(mar = c(0,5,3,0))
plot(sampleTree, main = "Sample clustering to detect outliers in localized PCa patients dataset", sub="", xlab="", cex.lab = 2,
     cex.axis = 2, cex.main = 2)

#############################################################################################################################################################################
#6/17/2020

#Restarting the WGCNA for CRPC_SU2C

#I am going to import both RNASEQ and Clinical Traits to my R environment
#My new data was generated after removing 4 duplicates as reported in my previous trials

#Re WGCNA of CRPC Dataset After removing 4 duplicates (i.e. from 270 samples to 266 samples) 
library(limma)
CRPC2_MET<-DUBREMOVED_CRPC_SU2C_RNASEQ.+1
logCRPCnew_MET<-log2(CRPC2_MET)
View(logCRPCnew_MET)

#Nextis the WGCNA anaysis
library(WGCNA)
options(stringsAsFactors = FALSE)

datExpr0 = as.data.frame(t(logCRPCnew_MET[, -c(1:8)]));

# I want to Check for genes and samples with too many missing values
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

# Next, I want to remove offending genes
datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]

#Next I want to cluster my samples to check for outliers
sampleTree = hclust(dist(datExpr0), method = "average");
sizeGrWindow(12,9)
sampleTree = hclust(dist(datExpr0), method = "average");
par(cex = 0.5);
par(mar = c(0,5,3,0))
plot(sampleTree, main = "Sample clustering to detect outliers in metastatic/CRPC patients dataset", sub="", xlab="", cex.lab = 2,
     cex.axis = 2, cex.main = 2)

datExpr = datExpr0
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#I am recalling the RNASeQ and loading the clinical trait dataset

traitData<- DUBREMOVED_CRPC_SU2C_CLINICAL_TRAITS.
dim(traitData)
names(traitData)

CRPC_Samples = rownames(datExpr)
traitRows = match(CRPC_Samples, traitData)
datTraits = traitData[traitRows, -1]

row.names(datTraits) = traitData[traitRows, 1]
collectGarbage();

# To check for duplicates
dups <- which(duplicated(CRPC_Samples))
View(dups)
length(dups)
dupp <-unique(dups)
View(dupp)
#prints the elements in CRPC-Traits that are occuring more than once




# To re-cluster samples
sampleTree2 = hclust(dist(datExpr0), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE)
plotDendroAndColors(sampleTree2, traitColors,groupLabels = names(datTraits), main = "Sample dendrogram and trait heatmap")

# To save data outcome
save(datExpr, datTraits, file = "CRPC_WGCNA_OUTPUT1.RData")

###############Creating network consessus and modules###########################################
getwd()
library(WGCNA)
options(stringsAsFactors = FALSE);
lnames = load(file = "CRPCC_WCGNA_OUTPUT1.RData");
lnames
nSets = checkSets(multiExpr)$nSets

# Choose a set of soft-thresholding powers
powers = c(seq(4,10,by=1), seq(12,20, by=2));
# Initialize a list to hold the results of scale-free analysis
powerTables = vector(mode = "list", length = nSets);
# Call the network topology analysis function for each set in turn
for (set in 1:nSets)
  powerTables[[set]] = list(data = pickSoftThreshold(multiExpr[[set]]$data, powerVector=powers,
                                                     verbose = 2)[[2]]);
collectGarbage();
# Plot the results:
colors = c("black", "red")
# Will plot these columns of the returned scale free analysis tables
plotCols = c(2,5,6,7)
colNames = c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
             "Max connectivity");
# Get the minima and maxima of the plotted points
ylim = matrix(NA, nrow = 2, ncol = 4);
for (set in 1:nSets)
{
  for (col in 1:length(plotCols))
    2
  {
    ylim[1, col] = min(ylim[1, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
    ylim[2, col] = max(ylim[2, col], powerTables[[set]]$data[, plotCols[col]], na.rm = TRUE);
  }
}
# Plot the quantities in the chosen columns vs. the soft thresholding power
sizeGrWindow(8, 6)
#pdf(file = "Plots/scaleFreeAnalysis.pdf", wi = 8, he = 6)
par(mfcol = c(2,2));
par(mar = c(4.2, 4.2 , 2.2, 0.5))
cex1 = 0.7;
for (col in 1:length(plotCols)) for (set in 1:nSets)
{
  if (set==1)
  {
    plot(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         xlab="Soft Threshold (power)",ylab=colNames[col],type="n", ylim = ylim[, col],
         main = colNames[col]);
    addGrid();
  }
  if (col==1)
  {
    text(powerTables[[set]]$data[,1], -sign(powerTables[[set]]$data[,3])*powerTables[[set]]$data[,2],
         labels=powers,cex=cex1,col=colors[set]);
  } else
    text(powerTables[[set]]$data[,1], powerTables[[set]]$data[,plotCols[col]],
         labels=powers,cex=cex1,col=colors[set]);
  if (col==1)
  {
    legend("bottomright", legend = setLabels, col = colors, pch = 20) ;
  } else
    legend("topright", legend = setLabels, col = colors, pch = 20) ;
}
dev.off();

#Network construction and Module detection for CRPCs
net = blockwiseConsensusModules(
  multiExpr, power = 6, minModuleSize = 30, deepSplit = 2,
  pamRespectsDendro = FALSE,
  mergeCutHeight = 0.25, numericLabels = TRUE,
  minKMEtoStay = 0,
  saveTOMs = TRUE, verbose = 5)

#Generating Eingengenes Construct
consMEs = net$multiMEs;
moduleLabels = net$colors;
# Convert the numeric labels to color labels
moduleColors = labels2colors(moduleLabels)
consTree = net$dendrograms[[1]];
#Plotting Output
sizeGrWindow(8,6);
#pdf(file = "Plots/ConsensusDendrogram-auto.pdf", wi = 8, he = 6)
plotDendroAndColors(consTree, moduleColors,
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Consensus gene dendrogram and module colors")
dev.off()
save(consMEs, moduleLabels, moduleColors, consTree, file = "CRPC_Consensus-NetworkConstruction-auto.RData")
