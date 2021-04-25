###################################################################
#Anova Plot and Volcano Plots
###################################################################
rm(list=ls(all=TRUE))
setwd("C:/Users/sosen/Documents")
getwd()
#numericMeta2<-read.csv("NumericMeta2.csv", header = TRUE, row.names = 1)
#head(numericMeta2)
#dim(numericMeta2)

load("finaloutput_for_Saheed_06_06_2020.RData")
dim(cleanDat)
dim(numericMeta)
#write.csv(cleanDat,"C:\\Users\\sosen\\Documents\\OSENI_R_STUDIO_DATA\\NCleaned_TCGAloggedData.csv", row.names = TRUE)
#write.csv(numericMeta,"C:\\Users\\sosen\\Documents\\OSENI_R_STUDIO_DATA\\Cleaned_TCGA_Clinical_TriatsData.csv", row.names = FALSE)


#Grouping<-numericMeta$PFS
#Grouping[numericMeta$PFS==1]<-"Progression" #only necessary if Group was numerically encoded; does nothing if the Grouping vector has no numeric values
#Grouping[numericMeta$PFS==0]<-"Non-progression"


Grouping<-numericMeta$CANCER_RELAPSE
Grouping[numericMeta$CANCER_RELAPSE==0]<-"No_PCa_Relapse"  #only necessary if Group was numerically encoded; does nothing if the Grouping vector has no numeric values
Grouping[numericMeta$CANCER_RELAPSE==1]<-"PCa_Relapse"

#Grouping<-numericMeta$RADIOTHERAPY
#Grouping[numericMeta$RADIOTHERAPY==0]<-"Non-exposed_Radiotherapy"  #only necessary if Group was numerically encoded; does nothing if the Grouping vector has no numeric values
#Grouping[numericMeta$RADIOTHERAPY==1]<-"Exposed_Radiotherapy"

length(Grouping)==length(colnames(cleanDat))
table(net$colors)
outputfigs=getwd()
FileBaseName = "9Modules"
groupsToTest=c("ALL") #unique(numericMeta$Group)

ANOVAoutList<-list()
#for (caseSubset in groupsToTest) {
#  cleanDat.sub<-cleanDat[,match(rownames(numericMeta)[numericMeta$Group==caseSubset],colnames(cleanDat))]
#  Grouping.sub <- numericMeta$Region.Treatment[numericMeta$Group==caseSubset] #character strings
caseSubset="ALL"
cleanDat.sub<-cleanDat
Grouping.sub<-Grouping

outfile = paste0("New-ANOVA_diffEx-",caseSubset,"-",FileBaseName)
data = as.data.frame(cbind(colnames(cleanDat.sub), Grouping.sub, t(cleanDat.sub)))
colnames(data)[1:2]<-c("CODE","SampleType") # just changes the 1st 2 rows to CODE & SampleType
#test run gets column headers for output
i=3
aov<-aov(data[,i]~SampleType, data=data) #SampleType is a factor 
#data2 <- as.numeric(as.character(colnames(data[3:31393])))
#dffinalSep$PC1 <- as.numeric(as.character(dffinalSep$PC1))
#aov<-aov(data2[,i]~SampleType, data2=data2)


anovaresult<-anova(aov)
tuk <- TukeyHSD(aov)
tukresult1<-data.frame(tuk$SampleType) #ASSUMES NO PAIRWISE COMPARISONS ARE MISSING FOR FIRST cleanDat.sub protein (ROWS OF THIS data frame)--this is the template for comparison columns in ANOVAoutList[[caseSubset]]
j=length(rownames(tukresult1))
comparisonList<-rownames(tukresult1)

line = c(paste("Protein", "F-Value", "Pr(>F)", sep=","))
for (a in 1:length(comparisonList)) {
  line=c(paste(line,comparisonList[a],sep=","))
}
for (a in 1:length(comparisonList)) {
  line=c(paste(line,paste0("diff ",comparisonList[a]),sep=","))
}

## Fast code with apply
dataclipped <- data[, 3:ncol(data)] # as.matrix(as.numeric(data[,3:ncol(data)]))
SampleType <- data$SampleType
ANOVAoutList[[caseSubset]] <- apply(dataclipped, 2, function(x) {
  if(length(unique(SampleType[which(!is.na(x))]))<2) { #<2, handles "only one contrast level..." 
    x=data.frame(x=rep(1,length(x)), SampleType = SampleType)  #handles too many missing values
  } else {
    x <- data.frame(x = as.double(x), SampleType = SampleType) #as.double(x) instead of x corrects:   Error in lm.fit(x, y,... NA/NaN/Inf in 'y'
  }
  aov <- aov(x ~ SampleType, data = x)
  anovaresult <- anova(aov)
  tuk <- TukeyHSD(aov)
  tukresult <- data.frame(tuk$SampleType)
  if (length(rownames(tukresult)) == length(rownames(tukresult1))) {
    c(anovaresult$F[1], anovaresult$Pr[1], as.vector(tukresult[, "p.adj"]), as.vector(tukresult[, "diff"]))
  } else {
    tukresult <- tukresult[match(rownames(tukresult1), rownames(tukresult)), ]
    c(anovaresult$F[1], anovaresult$Pr[1], as.vector(tukresult[, "p.adj"]), as.vector(tukresult[, "diff"]))
  }
})
ANOVAoutList[[caseSubset]] <- t(ANOVAoutList[[caseSubset]])
ANOVAcols <- as.vector(data.frame(do.call("rbind", strsplit(as.character(line), "[,]"))))
ANOVAcols <- ANOVAcols[2:length(ANOVAcols)]
if (length(unique(SampleType))==2) { #for single pairwise comparison, essentially Pr(>F) from ANOVA is equivalent to T-test result (except 1 vs 2 non-NA measurements are allowed)
  ANOVAoutList[[caseSubset]][,3] <- p.adjust(ANOVAoutList[[caseSubset]][,2],method="BH") #get BH FDR for ANOVA/T-test p values
  ANOVAoutList[[caseSubset]][,2:3]<-ANOVAoutList[[caseSubset]][,c(3,2)]
  ANOVAcols[2] <- "FDR (BH)"
}
colnames(ANOVAoutList[[caseSubset]]) <- ANOVAcols
ANOVAoutList[[caseSubset]] <- as.data.frame(ANOVAoutList[[caseSubset]])


## Add network module colors
ANOVAoutList[[caseSubset]]$NETcolors <- net$colors

#*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+
#write.csv(ANOVAoutList[[caseSubset]],file=paste0(outputtabs,"/",outfile,".csv"))
#*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+*+

#  cat(paste0("Finished ANOVA for subset ",caseSubset,". [",length(groupsToTest)," total groups to test]\n"))
#}
ANOVAout<-ANOVAoutList[["ALL"]]
#write.csv(ANOVAout, file = "Module9_ANOVAout.csv")
colnames(ANOVAout)

## Generate Volcano plots, both PDF and interactive (mouseover) HTML for selected pairwise comparisons in ANOVA-Tukey/T-test stats. 
###################################################################################################################################
library(ggplot2)

#** this loop would allow for multiple different ANOVA groupings corresponding to multiple DFs in ANOVAoutList, instead of the single DF ANOVAout
#** for (caseSubset in groupsToTest) { #*+*+*+*+*
#** baseNameVolcanoes <- paste0(caseSubset,".",FileBaseName)
#** ANOVAout<-ANOVAoutList[[caseSubset]]

baseNameVolcanoes <- paste0(FileBaseName) #** paste0(FileBaseName,".",caseSubset)


# We want to note the positions of the important genes in each volcano
colnames(ANOVAout) # choose a column or columns by column number(s) (counting from left) below (column numbers point to testIndex below)

#Configuration parameters for this code block  (stop here, edit below parameters based on console output showing columns)
#*****************************





#Configuration parameters for this code block  (stop here, edit below parameters based on console output showing columns)
#*****************************
testIndexMasterList <- c(3)  #What columns (column numbers) of ANOVAout (T-test or Tukey p values) do you want volcanoes for? 
flip <- as.vector(c())      #Edit this vector, to flip sign of this/these comparison(s); e.g. put WT as denominator for log2(ratio) if you see it first in the colnames(ANOVAout) you just selected
#flip <- as.vector(c(3))      #Edit this vector, to flip sign of this/these comparison(s); e.g. put WT as denominator for log2(ratio) if you see it first in the colnames(ANOVAout) you just selected
cutoff <- log2(1)            # log2(1) means NO change minimum to be counted in the volcano bookends; log2(1.25) for 25% FC min.
sigCutoff <- 0.05            # p value cutoff for Volcano significant hit counting; dashed line at -log10(sigCutoff)

useNETcolors<-TRUE  #Use WGCNA colors if TRUE, otherwise red=up, green=down, and blue=no change beyond cutoffs
#splitColors<-FALSE  #Do not Make a separate volcano for each color of spot?
splitColors<-TRUE  #Make a separate volcano for each color of spot?

BIGspots <- c("IRAK1","IRAK2","IRAK3","IRAK4") #c('TARDBP','C9ORF72','SCA36','SCA3','SQSTM1','SFPQ','MSN','PLEC','HEPACAM') #which gene symbols should be checked for and highlighted as bigger spots, if any?
#*****************************
cutoff  #shows what your cutoff for log2(FC) calculates as


n <- nrow(ANOVAout)
dexComps <- list()
iter <- length(testIndexMasterList) + 1
comparisonIDs <- data.frame(dfVariable = rep(NA, length(testIndexMasterList)), Comparison = rep(NA, length(testIndexMasterList)))
numberOfNonComparisonColumns=length(colnames(ANOVAout)) - length(which(grepl("diff ",colnames(ANOVAout))))*2
numComp <- (length(colnames(ANOVAout)) - numberOfNonComparisonColumns) / 2 # of columns separating comparisons from matched column of log2(diffs), i.e. # of comparisons
for (i in testIndexMasterList) {
  iter <- iter - 1
  # dexRows<-which(ANOVAout[,i]<sigCutoff) #choose rows where the DEX p<sigCutoff
  comparisonIDs[iter, ] <- as.vector(c(paste0("dexTargets.", gsub("-", ".", colnames(ANOVAout)[i])), paste0(as.character(gsub("-", " vs ", colnames(ANOVAout)[i])))))
  dexComps[[comparisonIDs[iter, 1]]] <- ANOVAout
  if (!is.na(match(i, flip))) {
    dexComps[[comparisonIDs[iter, 1]]][, i + numComp] <- -1 * as.numeric(dexComps[[comparisonIDs[iter, 1]]][, i + numComp])
    comparisonIDs[iter, 2] <- gsub("(*.*) vs (*.*)", "\\2 vs \\1", comparisonIDs[iter, 2]) # flip label "vs" in comParisonIDs$Comparison[iter]
  }
}
comparisonIDs # list element names and Logical comparisons for those retrievable Dex measurements in the list elements
ls(dexComps) # list elements are dataframes with the DEX entries for that comparison


## volcano plots with (module or other) colors, SPLITTABLE on COLOR

pointColorsVectorListForPlots <- list()
volcListModColorsWeb <- list()
volcListModColors <- list()
dfListModColors <- list()

iter <- length(testIndexMasterList) + 1
for (testIndex in testIndexMasterList) {
  iter <- iter - 1
  df <- eval(parse(text = "dexComps[[comparisonIDs$dfVariable[iter]]]"))
  cat(paste0("Processing ANOVA column ", testIndex, " (", comparisonIDs$Comparison[iter], ") for volcano...\n"))
  # correct 0 Tukey pValues to ANOVA p (in column 2); it's better than taking -log10 of 0 in the next step
  df[which(df[, testIndex] == 0), testIndex] <- as.numeric(df[which(df[, testIndex] == 0), 2])
  
  ## Check if ANOVA pVal is Significant and above FC cutoff defined above. Thresholds are used to set volcano point colors
  df[, testIndex][is.na(df[, testIndex])] <- 1 # p=0.9999 instead of NA
  df[, testIndex + numComp][is.na(df[, testIndex + numComp])] <- 0 # log2(difference)=0 instead of NA
  
  df$negLogP <- -log10(as.numeric(df[, testIndex]))
  
  df$threshold1 <- as.numeric(rep(0, n))
  ## Any COMPARISON SIGNIFICANT (uses ANOVA p in column 2 of df instead of Tukey p): # for (i in 1:n) { if (abs(as.numeric(df[i,testIndex+numComp]))<cutoff | df[i,2]>sigCutoff ) {df$threshold1[i]=3} else { if (df[i,testIndex+numComp]<cutoff) {df$threshold1[i]=2} else {df$threshold1[i]=1}} }
  for (i in 1:n) {
    if (abs(as.numeric(df[i, testIndex + numComp])) < cutoff | as.numeric(df[i, testIndex]) > sigCutoff) {
      df$threshold1[i] <- 3
    } else {
      if (as.numeric(df[i, testIndex + numComp]) < cutoff) {
        df$threshold1[i] <- 2
      } else {
        df$threshold1[i] <- 1
      }
    }
  }
  df$threshold1 <- as.factor(df$threshold1)
  
  df$Symbol <-   do.call("rbind", strsplit(as.character(rownames(df)), "[|]"))[, 1]
  
  ## Color Interesting Gene Product Spots DIFFERENTLY as 4th color if doing blue/red/green (no module colors) -- (4=gold1 below)
  if(useNETcolors) { 
    df$color1<-df$NETcolors
    df$threshold2 <- as.numeric(df$threshold1)
    df$threshold2[match(intersect(df$Symbol, BIGspots), df$Symbol)] <- 4 
  } else {
    df$color1 <- as.numeric(df$threshold1)
    df$color1[match(intersect(df$Symbol, BIGspots), df$Symbol)] <- 4
    df$threshold2 <- as.numeric(df$threshold1)
    df$threshold2[match(intersect(df$Symbol, BIGspots), df$Symbol)] <- 4 
  }
  df$color1 <- as.factor(df$color1)
  
  if(useNETcolors) {
    df$size1 <- as.numeric(df$threshold2)
  } else {
    df$size1 <- as.numeric(df$color1)
  }
  df$size1[df$size1 < 4] <- 2.5
  df$size1[df$size1 == 4] <- 6
  
  #df$color2: actual spot color as.character() ; df$color1 is a factorable dummy
  if(useNETcolors) {
    df$color2<-as.character(df$color1)
    df$color2[df$threshold2 == 4] <- "gold1" #for BIGspots
  } else {
    df$color2 <- as.numeric(df$color1)
    df$color2[df$color2 == 1] <- "darkred"
    df$color2[df$color2 == 2] <- "darkgreen"
    df$color2[df$color2 == 3] <- "dodgerblue"
    df$color2[df$color2 == 4] <- "gold1" #for BIGspots
  }                              #gold1 is also the 435th, last WGCNA unique color
  
  #df$color3 is outline color, where outlined pch symbols (21) used
  df$color3 <- df$color2
  df$color3[df$color3 == "gold1"] <- "black" #for BIGspots outline
  df$pch <- as.numeric(df$threshold2)
  df$pch[df$pch < 4] <- 16 # unfilled circles (use color2)
  df$pch[df$pch == 4] <- 21 # filled, outlined circles (border uses color3)
  
  #put gold1 back to module color for fill of BIGspots if (useNETcolors)
  if (useNETcolors) { df$color2[df$color2=="gold1"] <- as.character(df$color1[df$color2=="gold1"]) }
  
  df <- df[order(df$size1, decreasing = FALSE), ] # puts larger dots on top (at bottom of df)
  
  
  #splitColors TRUE/FALSE: make one volcano with all colors (FALSE), or make volcanoes for each color (TRUE)
  # SPLIT DATA FRAME FOR VOLCANO PLOT BY COLORS (if multiple eachColorSplit items)
  df.AllColors <- df
  eachColorSplit <- if (splitColors) {
    unique(df.AllColors$NETcolors)
  } else {
    c("allcolors")
  }
  for (eachColor in eachColorSplit) {
    if (splitColors) {
      df.oneColor <- df.AllColors[which(df.AllColors$NETcolors == eachColor), ]
      df.oneColor$color1 <- factor(df.oneColor$NETcolors)
    }
    else {
      df.oneColor <- df.AllColors
    } # end if (splitColors)
    
    names(df.oneColor)[testIndex + numComp] <- "xdata" # x=as.numeric(df.oneColor[,testIndex+numComp])
    list_element <- paste(comparisonIDs$dfVariable[iter], eachColor, sep = ".") # colnames(df)[testIndex]
    pointColorsVectorListForPlots[[list_element]] <- data.frame(color1 = factor(as.integer(df.oneColor$color1)), color2 = as.character(df.oneColor$color2), color3 = as.character(df.oneColor$color3), size = as.numeric(df.oneColor$size), pch = as.numeric(df.oneColor$pch)) #*** df.oneColor$NETcolors
    
    
    volcano1 <- ggplot(data = df.oneColor, aes(x = xdata, y = negLogP, color = color1, text = Symbol)) +
      # scale_colour_manual(values = unique(data.frame(col1=df.oneColor$color1,col2=df.oneColor$color2))[order(unique(data.frame(col1=df.oneColor$color1,col2=df.oneColor$color2))[,1]),2] )+ #THIS COLOR(S) IS LOOKED UP ACTIVELY DURING GRAPHICS OUTPUT IN THE VARIABLE, SO WE'VE USED A LIST ELEMENT THAT IS NEVER CHANGED
      geom_point(aes(fill = pointColorsVectorListForPlots[[list_element]][, "color3"]), alpha = 0.66, size = pointColorsVectorListForPlots[[list_element]]$size, pch = pointColorsVectorListForPlots[[list_element]][, "pch"], color = pointColorsVectorListForPlots[[list_element]][, "color3"]) +
      theme(legend.position = "none") +
      xlim(c(min(as.numeric(df.oneColor[, testIndex + numComp])), max(as.numeric(df.oneColor[, testIndex + numComp])))) + ylim(c(0, max(df.oneColor$negLogP))) +
      xlab(as.expression(bquote("Difference, log"[2] ~ .(comparisonIDs$Comparison[iter])))) + # colnames(df.oneColor)[testIndex]
      ylab(as.expression(bquote("-log"[10] ~ "p value"))) +
      theme(axis.title.x = element_text(size = rel(1.8), angle = 00)) +
      theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
      
      geom_hline(yintercept = 1.30103, linetype = "dashed", color = "black", size = 1.2) +
      #geom_text(aes(0,1.30103,label = 1.30103, vjust = -1))+
      geom_vline(xintercept = cutoff, linetype = "dashed", color = "black", size = 1.2) +
      geom_vline(xintercept = -cutoff, linetype = "dashed", color = "black", size = 1.2) +
      annotate("text", x = min(as.numeric(df.oneColor[, testIndex + numComp])) / 2, y = max(df.oneColor$negLogP) * .95, size = 5, label = paste0("Downregulated: ", bquote(.(length(which(as.numeric(df.oneColor$threshold1) == 2)))))) +
      annotate("text", x = max(as.numeric(df.oneColor[, testIndex + numComp])) / 2, y = max(df.oneColor$negLogP) * .95, size = 5, label = paste0("Upregulated: ", bquote(.(length(which(as.numeric(df.oneColor$threshold1) == 1)))))) +
      
      theme(
        # axis.text = element_text(size = 14),
        # legend.key = element_rect(fill = "navy"),
        # legend.background = element_rect(fill = "white"),
        # legend.position = c(0.14, 0.80),
        panel.grid.major = element_line(color = "darkgrey", linetype = "dashed"),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white")
      )
    
    
    # web version doesn't use as.expression! plotly fails with those, so we rebuild the volcano for the web.
    volcanoweb <- ggplot(data = df.oneColor, aes(x = xdata, y = negLogP, color = color1, text = Symbol)) +
      scale_colour_manual(values = unique(data.frame(col1 = df.oneColor$color1, col2 = df.oneColor$color2))[order(unique(data.frame(col1 = df.oneColor$color1, col2 = df.oneColor$color2))[, 1]), 2]) + # THIS COLOR(S) IS LOOKED UP ACTIVELY BY PLOTLY IN THE VARIABLE, SO WE'VE USED A LIST ELEMENT THAT IS NEVER CHANGED
      # scale_y_continuous(breaks = seq(0, 8, by = 1))+
      # geom_point(aes(fill=pointColorsVectorListForPlots[[list_element]][,"color2"]), alpha=0.66, size=pointColorsVectorListForPlots[[list_element]]$size, pch=16, color=pointColorsVectorListForPlots[[list_element]][,"color3"]) + #pch=pointColorsVectorListForPlots[[list_element]][,"pch"] not using variable symbol types (21 for outlined circle only)
      geom_point(alpha = 0.66, size = pointColorsVectorListForPlots[[list_element]]$size, pch = 16) + # pch=pointColorsVectorListForPlots[[list_element]][,"pch"] just uses the higher pch code in the web render.
      theme(legend.position = "none") +
      xlim(c(min(as.numeric(df.oneColor[, testIndex + numComp])), max(as.numeric(df.oneColor[, testIndex + numComp])))) + ylim(c(0, max(df.oneColor$negLogP))) +
      xlab(paste0("Difference, log2 ", comparisonIDs$Comparison[iter])) +
      ylab(paste0("-log10 p value")) +
      theme(axis.title.x = element_text(size = rel(1.8), angle = 00)) +
      theme(axis.title.y = element_text(size = rel(1.8), angle = 90)) +
      
      geom_hline(yintercept = 1.30103, linetype = "dashed", color = "black", size = 1.2) +
      #geom_text(aes(0,1.30103,label = 1.30103, vjust = -1))+
      #geom_hline(yintercept = 1.0, linetype = "dashed", color = "black", size = 1.2) +
      # geom_text(aes(0,1.30103,label = 1.30103, vjust = -1))+
      
      geom_vline(xintercept = cutoff, linetype = "dashed", color = "black", size = 1.2) +
      geom_vline(xintercept = -cutoff, linetype = "dashed", color = "black", size = 1.2) +
      annotate("text", x = min(as.numeric(df.oneColor[, testIndex + numComp])) / 2, y = max(df.oneColor$negLogP) * .95, size = 5, label = paste0("Downregulated: ", bquote(.(length(which(as.numeric(df.oneColor$threshold1) == 2)))))) +
      annotate("text", x = max(as.numeric(df.oneColor[, testIndex + numComp])) / 2, y = max(df.oneColor$negLogP) * .95, size = 5, label = paste0("Upregulated: ", bquote(.(length(which(as.numeric(df.oneColor$threshold1) == 1)))))) +
      
      theme(
        # axis.text = element_text(size = 14),
        # legend.key = element_rect(fill = "navy"),
        # legend.background = element_rect(fill = "white"),
        # legend.position = c(0.14, 0.80),
        panel.grid.major = element_line(color = "darkgrey", linetype = "dashed"),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white")
      )
    
    
    volcListModColors[[list_element]] <- volcano1
    volcListModColorsWeb[[list_element]] <- volcanoweb
    
    print(volcano1) # prints to active output (separate page)
    rm(volcano1)
    rm(volcanoweb)
    dfListModColors[[list_element]] <- df.oneColor
  } # closes for(eachColor...
} # closes for(testIndex...



if(splitColors) { dir.create(file.path(outputfigs, "/SplitVolcano/")) }

# Print to PDFs, one per color (per comparison, if multiple)
iter <- length(testIndexMasterList) + 1
for (testIndex in testIndexMasterList) {
  iter <- iter - 1
  for (eachColor in eachColorSplit) {
    list_element <- paste(comparisonIDs$dfVariable[iter], eachColor, sep = ".")
    df.oneColor <- dfListModColors[[list_element]]
    if(splitColors) {
      file <- paste0(outputfigs,"/SplitVolcano/","Volcano_", eachColor,"-",gsub(" ", "_", comparisonIDs$Comparison[iter]),"-",baseNameVolcanoes,".pdf")
    } else {
      file <- paste0(outputfigs,"/","Volcano_", eachColor,"-",gsub(" ", "_", comparisonIDs$Comparison[iter]),"-",baseNameVolcanoes,".pdf")
    }
    pdf(file = file, colormodel="srgb", height = 8, width = 8)
    par(mfrow = c(1, 1))
    par(mar = c(6, 8.5, 3, 3))
    print(volcListModColorsWeb[[list_element]])
    #     volcListModColors[[list_element]])    #bigspots colors not handled correctly for PDF output of modColored spots:
    dev.off()
  }
}

## Print to html volcano plots (one per module colors and per comparison, if applicable)
library(plotly)

iter <- length(testIndexMasterList) + 1
for (testIndex in testIndexMasterList) {
  iter <- iter - 1
  for (eachColor in eachColorSplit) {
    list_element <- paste(comparisonIDs$dfVariable[iter], eachColor, sep = ".")
    df.oneColor <- dfListModColors[[list_element]]
    webPlot <- ggplotly(volcListModColorsWeb[[list_element]])
    if(splitColors) {
      tempfilename=paste0(outputfigs,"/SplitVolcano/","HTMLvolcano_",eachColor,"-",gsub(" ", "_", comparisonIDs$Comparison[iter]),"-",baseNameVolcanoes,".html")
    } else {
      tempfilename=paste0(outputfigs,"/","HTMLvolcano_",eachColor,"-",gsub(" ", "_", comparisonIDs$Comparison[iter]),"-",baseNameVolcanoes,".html")
    }
    htmlwidgets::saveWidget(webPlot, tempfilename, selfcontained = TRUE, libdir = "delete.me")
  }
}


