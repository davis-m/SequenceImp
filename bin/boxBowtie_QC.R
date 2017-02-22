
###########################################################################

#This script is part of the Kraken framework which aims to facilitate RNA 
#sequence analysis in a streamlined and efficient manner.
#Copyright (C) 2011 2012 2013 EMBL - European Bioinformatics Institute 

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program (Please see the COPYING file for details).
#If not, see <http://www.gnu.org/licenses/>.

#Please send bug reports to:
#kraken@ebi.ac.uk

###########################################################################

library(gplots)

# NOTE: gplots reporting requirements for additional libraries to run read.xls()

##################

maxPlotExp <- 1
debugme <- FALSE  # Print debugging information

args <- R.utils::commandArgs(asValue=TRUE)

if (!is.null( args$input)){
   inputDir <- args$input
}else{
   stop ("Require an input directory containing Bowtie QC tables : --input=<DIRECTORY>")
}

if (!is.null(args$basicID)){
   basicID <- args$basicID
}else{
   stop ("Require a base name to uniquely identify the results pdf: --basicID=<STRING>")
}

if (!is.null(args$grangesClass)){
   grangesObj <- args$grangesClass
}else{
   stop ("Require the GRanges annotation object specifiying the annotation classes of the relevant Ensembl version : --grangesClass=<FILE.RData>")
}

if (!is.null(args$output)){
   outDir <- args$output
}else{
   stop ("Require an output directory into which the QC plots will be deposited: --output=<DIRECTORY>")
}

if(!is.null(args$maxHits)){
   maxPlotExp <- args$maxHits
}
#stopifnot(FALSE)

if(!is.null(args$chrTypes)){
   chromosomes <- args$chrTypes
}else{
   stop ("Require a file containing information concerning toplevel designation of chromosomes within the Ensembl API (eg. chromosome/supercontig) : --chrTypes=<FILE>")
}

if(!is.null(args$procSummary)){
   savedReads <- args$procSummary
}else{
   stop ("Require a summary of the total number of processed reads (--procSummary=.../total.saved.trimmit.tab).")
}

if (!is.null(args$debug)){
#   debugme <- args$debug
   debugme <- TRUE
}

if (debugme){write("\n### DEBUGGING ###\n\n",stderr())}

#stopifnot(FALSE)



#########################################################################################
### Judges the number of samples based on the number of ".bowtie.output.bwt.gz" files ###
#########################################################################################

sampleNo <- length(dir(inputDir, pattern = ".bowtie.conv.GR.RData"))

##################
### Set up PDF ###
##################

outputFile <- paste(outDir,"/",basicID,"_Bowtie_qc.pdf", sep="")

pdf(file = outputFile ,height=(4.6+8*sampleNo), width=9.6)
#dev.new(height=(4.6+8*sampleNo), width=9.6)

layoutStuff <- c(1,2,3)
curCount <- 3
for (k in seq_len(sampleNo)){
   layoutStuff <- c(layoutStuff, seq(curCount+1, curCount+3))   
   curCount <- curCount+3
}

for (i in seq_len(sampleNo)){
   layoutStuff <- c(layoutStuff, rep.int(curCount+1, 3))
   curCount <- curCount+1
}

layoutStuff <- c(layoutStuff, rep.int((curCount+1),3))

layout(matrix(layoutStuff, nrow = (2+2*(sampleNo)),
           ncol= 3, byrow = TRUE))

par(omi=c(0.3,0.3,0.3,0.3))

##################################################
### Mappable Reads vs. Uniquely mappable reads ###
##################################################
write("Plotting mappable reads against those that map uniquely",stderr())

#stopifnot(FALSE)

ug <- par(mar=c(5,5,4,2))
mapSummary <- dir(inputDir, pattern = "bowtie.total.mapping.tab", full.names=TRUE)

if(length (mapSummary) != 1){
   stop("Detected the wrong number of 'bowtie.total.mapping.tab' files in input directory")
}

mappableTable <- as.matrix(read.table(mapSummary[1], sep="\t", as.is=TRUE, header=TRUE))
mapCols <- c("blue", "red")
names(mapCols) <- rownames(mappableTable)

barpoints <- barplot(mappableTable, col= mapCols, xlab= "Samples", beside=TRUE, main = "The number of mappable reads for each sample",cex.axis = 0.75, cex.main=0.75, las=2, xaxt='n')
title(ylab= "Read Number", line=4, cex = 0.75)
axis(1, at=apply(barpoints,2,mean), las=2, cex.axis=0.8, lty=0,  labels = colnames(mappableTable))

### Read in processed read number

savedTable <- as.matrix(read.table(savedReads, sep="\t", as.is=TRUE, row.names=1))

rownames(barpoints) <- rownames(mappableTable)
colnames(barpoints) <- colnames(mappableTable)

for(thisone in colnames(barpoints)){
   for(thatone in rownames(barpoints)){
      text(barpoints[thatone, thisone], as.numeric(mappableTable[thatone, thisone])/2,paste(as.character(mappableTable[thatone, thisone])," (",signif(as.numeric(mappableTable[thatone, thisone])*100/as.numeric(savedTable[thisone,1]),digits=3),'%)', sep=""), srt=90, cex=0.5, col='white')
   }
}



par(ug)
plot(NA, xlim=c(0,1), ylim=c(0,1),axes=F, xlab="", ylab="")
legend("left",pt.cex = 0.75, cex = 0.75, fill=mapCols, legend=names(mapCols))

######################################################################
### Proportion of reads mapping to a specified number of locations ###
######################################################################

write("Plotting reads against the number of locations to which they map",stderr())

multiMapTabs <- list()
multiMapFiles <- dir(inputDir, pattern = "bowtie.hitFreq.tab", full.names=TRUE)

for (multiMapName in multiMapFiles){
   multiMapTabs[[sub(".+\\.([[:alpha:]]+)\\.bowtie.+","\\1",multiMapName)]] <-  as.matrix(read.table(multiMapName, sep="\t", as.is=TRUE, header=TRUE))
}

maximumLocs <- max(unlist(lapply(multiMapTabs, function(z){
   max(as.numeric(rownames(z)))
}), use.names=FALSE))

multiDist <- matrix(data=0, nrow=maximumLocs, ncol=length(names(multiMapTabs)),dimnames=list("rows"=seq_len(maximumLocs), "columns"=names(multiMapTabs)))

for(tableChoice in names(multiMapTabs)){
   for (rowDef in rownames(multiMapTabs[[tableChoice]])){
      multiDist[rowDef,tableChoice] <- multiMapTabs[[tableChoice]][rowDef,"Frequency"]
   } 
}

multiCols <- rainbow(ncol(multiDist))
names(multiCols) <- colnames(multiDist)
#stopifnot(FALSE)
plot(multiDist[,1],type = "o", ylim=c(0,100),col= multiCols[colnames(multiDist)[1]], axes=FALSE, 
      main=paste("Percentage of processed reads mapping to a specified\nnumber of loci given that reads reporting more than", 
      maxPlotExp,"\nloci are discarded"), xlab="Loci to which a read maps", ylab="Percentage of reads", cex.main=0.75)

if(ncol(multiDist) > 1){
   write("Loci information being plotted for multiple samples",stderr())
   for (zzz in 2:ncol(multiDist)){
      lines(multiDist[,zzz],type = "o", col = multiCols[colnames(multiDist)[zzz]])
   }
}else{
   write("Loci information being plotted for single sample", stderr())
}

legend("topright", pt.cex=0.75, cex=0.75, fill=multiCols, legend=names(multiCols))

axis(1,at=1:nrow(multiDist),labels=rownames(multiDist))
axis(2,at=seq.int(0,100,by=10), labels = seq.int(0,100,by=10), cex.axis = 0.75, las =2)

########################################
### Distribution of reads to classes ###
########################################

################################################
### Function for plotting barpots of classes ###
################################################

bowtieSummary <- function (distribution, cols, title, ylimit){

  # shuffledCols <- c()
  # for (i in c(1:length(cols))) {
  #    k <- i * 5
  #    j <- (k %% length(cols))+1 ### R strings indexed beginning at 1
  #    #print(class(j))
  #    shuffledCols[i] <-  cols[j]
  #    #print (paste(i,j,k,shuffledGeneCols[i], sep =" "))
  # }
   #OR:
   #
   #bRb <- unlist(lapply(seq_along(aRb), function (x) { aRb[1 + ((x * 3) %% length(aRb))] }))
 
   oldpar <- par(mar=c(7.5,4.1,4.1,2.1))
   barplot(distribution,col = cols, main = title, beside = FALSE, ylab="Proportion of reads", ylim=c(0,ylimit),las=2, cex.main=0.75, cex.axis = 0.75, cex.names=0.5)
   legend("right",pt.cex= 0.75, cex = 0.75 , fill=cols, legend=rownames(distribution))
   par(oldpar)
}



##############################
### Class plotting process ###
##############################

write("Plotting the number of reads overlapping with specified annotation classes",stderr())

load(file = grangesObj)

readDistFile <- dir(inputDir, pattern = "bowtie.class.summary.RData", full.names=TRUE)

if(length (readDistFile) != 1){
   stop("Detected the wrong number of 'bowtie.class.summary.RData' files in input directory")
}

#readDistTable <- as.matrix(read.table(readDistFile[1], sep="\t", as.is=TRUE, header=TRUE))
load(readDistFile)

#stopifnot(FALSE)

Other_annotation <- list()
for (proportionalise in names(sumValues)){
   tableTotal <- sum(sumValues[[proportionalise]]['total',])
   sumValues[[proportionalise]] <- as.matrix(sumValues[[proportionalise]][-(which(rownames(sumValues[[proportionalise]]) == 'total')),])
   sumValues[[proportionalise]] <- sumValues[[proportionalise]]/tableTotal
   whichothers <- (apply(sumValues[[proportionalise]],1,sum)<0.01) # Added 170714 - Number of classes no longer known. Limit plot.
   Other_annotation[[proportionalise]] <- sumValues[[proportionalise]][whichothers,]
   sumValues[[proportionalise]]  <- sumValues[[proportionalise]][! whichothers,, drop=FALSE] 
}

#cycle <- 0

repeatCols <- c("orange","green")
geneCols <- c("purple","yellow")
yaxisLimit <- max(unlist(lapply(sumValues,function(x){max(rowSums(x))})))
repeatOrder <- c()
geneOrder <- c()

for (sampleChoice in names(sumValues)){
   
   # Altered 170714 - Total annotation classes no longer controlled. Need to limit.

   readDistTable <- sumValues[[sampleChoice]]
   otherDistTable <- Other_annotation[[sampleChoice]]
   
   if (! is.null(classRecords[["Repeats"]])){
      allClasses <- c(classRecords[["Repeats"]], classRecords[["GeneClasses"]])
   }else{
      write("No repeat annotation found - comparing to genes only\n",stderr())
      allClasses <- classRecords[["GeneClasses"]]
   }

#   if (any(!rownames(readDistTable)%in% allClasses) | any(! allClasses %in% rownames(readDistTable) )){stop("Bowtie distribution table and expected classes do not align")}
   if (any(!rownames(readDistTable)%in% allClasses)){stop("Bowtie distribution table and expected classes do not align")}
   if (any(!rownames(otherDistTable)%in% allClasses)){stop("Bowtie distribution table and expected classes do not align")}

   ### GENE PLOTS:
   #readDistRepeats <- matrix(data=0, nrow = length(classRecords[["Repeats"]]), ncol = ncol(readDistTable), dimnames = list("row" = classRecords[["Repeats"]],"col" = colnames(readDistTable))) # Altered 170714 - Total annotation classes no longer controlled. Need to limit.
   genesDistTable <- readDistTable[rownames(readDistTable)%in%classRecords[["GeneClasses"]],,drop=FALSE]
   remaining_gene_classes <- apply(otherDistTable[rownames(otherDistTable) %in% classRecords[["GeneClasses"]],],2,sum)
   readDistGenes <- rbind(genesDistTable,t(as.data.frame(remaining_gene_classes))[,colnames(genesDistTable),drop=FALSE])

   #readDistGenes <- matrix(data=0, nrow = length(classRecords[["GeneClasses"]]), ncol = ncol(readDistTable), dimnames = list("row" = classRecords[["GeneClasses"]],"col" = colnames(readDistTable)))
  
#   if (cycle == 0){ # Removed 170714 - different class numbers for each set
   if(nrow(readDistGenes)>1){
      readDistGenes <- readDistGenes[order(rowSums(readDistGenes), decreasing=TRUE),]
   }
   #repeatOrder <- rownames(readDistRepeats) 
   #geneOrder <- rownames(readDistGenes)
#   }
   geneTitle <- paste("Reads mapping to genic features as a proportion\nof the reads that map to an annotation class for", sampleChoice)
   bowtieSummary(t(readDistGenes), geneCols, geneTitle, yaxisLimit)
   
   ### REPEAT PLOTS:
   repeatTitle <- paste("Reads mapping to repeat features as a proportion\nof the reads that map to an annotation class for",sampleChoice)
   if (! is.null(classRecords[["Repeats"]])){
      repeatDistTable  <- readDistTable[rownames(readDistTable)%in%classRecords[["Repeats"]],,drop=FALSE]
      remaining_repeat_classes <- apply(otherDistTable[rownames(otherDistTable) %in% classRecords[["Repeats"]],],2,sum)
      readDistRepeats <- rbind(repeatDistTable,t(as.data.frame(remaining_repeat_classes))[,colnames(repeatDistTable),drop=FALSE])
      
      if(nrow(readDistRepeats)>1){
         readDistRepeats <- readDistRepeats[order(rowSums(readDistRepeats), decreasing=TRUE),] 
      }
      
      bowtieSummary(t(readDistRepeats), repeatCols, repeatTitle, yaxisLimit)
   }else{
      plot(NA, main=repeatTitle, cex = 0.5, xlim=c(0,1), ylim=c(0,1), axes=F, xlab="", ylab="",cex.main=0.75)
      text(0.5,0.5,"Repeat annotation not found for\ngenome and version specified")
   }  

   plot(NA, xlim=c(0,1), ylim=c(0,1), axes=F, xlab="", ylab="")

#   cycle <- cycle +1

}
#stopifnot(FALSE)



#repeatCols <- colorpanel(nrow(readDistRepeats), low="Red", high="DarkSlateGray3")


#geneCols <- colorpanel(nrow(readDistGenes), low="yellow", high="MediumPurple3")



###############################
### Chromosome Distribution ###
###############################

write("Plotting the distribution of reads accross chromosomes",stderr())

chrmFiles <- dir(inputDir, pattern = ".bowtie.chromLoc.tab", full.names=TRUE)
chrInfo <- read.table(chromosomes, header=TRUE, sep="\t", row.names = 1, colClasses = c("character","character"))
toplevels <- split(rownames(chrInfo), chrInfo$Class)
#stopifnot(FALSE)
for (chrmSel in chrmFiles){
   chrmTable <- read.table(chrmSel, sep="\t", as.is=TRUE, check.names = FALSE)
  
   ### Organising chromsomes for plotting
   chrPlot <- chrmTable[, colnames(chrmTable) %in% toplevels[["chromosome"]]]
   chrPlotNames <- colnames(chrPlot)   ### Sorting chromosome numbers numerically
   rr <- regexpr("^\\d+$", chrPlotNames) > 0
   chrPlotNames[rr] <- as.character(sort(as.numeric(chrPlotNames[rr])))  
   chrPlotNames[!rr] <- sort(chrPlotNames[!rr])  
   chrPlot<- chrPlot[,chrPlotNames]

   otherChr <- paste(names(toplevels)[!names(toplevels) == "chromosome"], collapse="/")
   
   if (nchar(otherChr) > 0){
      chrPlot[,otherChr] <- c(0,0,0,0)
   
      for (rowstuff in rownames(chrPlot)){
         chrPlot[rowstuff,otherChr ] <- sum(chrmTable[rowstuff,!colnames(chrmTable) %in% toplevels[["chromosome"]]])
      }
   }

   chrPlot[c("Mappable_AntiSense","UniquelyMappable_AntiSense"),] <- -chrPlot[c("Mappable_AntiSense","UniquelyMappable_AntiSense"),]
   
   ### Sorting plot colours - NOTE: order of names for colours and column names must correspond!
   plot_colours <- c("red","blue", "orange", "purple")
   names(plot_colours) <- c("Mappable, Positive Strand", "Uniquely Mappable, Positive Strand", "Mappable, Negative Strand", "Uniquely Mappable, Negative Strand")
   chrPlot <- chrPlot[ c("Mappable_Sense","UniquelyMappable_Sense","Mappable_AntiSense","UniquelyMappable_AntiSense"),]

   ### Sorting yaxis
   if (max(chrPlot) > 0 & min(chrPlot) < 0 ){
      chrPoints <- unique(c(round(seq((ceiling(min(chrPlot))-1),0, 1)),round(seq(0,ceiling(max(chrPlot)),1))))
   }else if(max(chrPlot) > 0 & min(chrPlot) >= 0){
      chrPoints <- seq(0,ceiling(max(chrPlot)),1)
   }else if(max(chrPlot) <= 0 & min(chrPlot) < 0){
      chrPoints <- seq((ceiling(min(chrPlot))-1),0, 1)
   }

   ### Making the plot
   barTitle <- paste(sub(".+\\.([[:alpha:]]+)\\.bowtie.+","\\1",chrmSel),": The percentage of processed reads mapping to each chromosome on either the positive or negative strands")
   barplot(as.matrix(chrPlot) , beside=TRUE, col=plot_colours,ylim =c(min(chrPoints),max(chrPoints)), ylab="Percentage of processed reads" , xlab="Ensembl chromosomes" ,main=barTitle, axes = FALSE, cex.axis = 0.75, cex.main=0.75, las = 2)

   axis(side = 2, at = chrPoints, labels=abs(chrPoints), cex.axis=0.75, las=2)

}
plot(NA, xlim=c(0,1), ylim=c(0,1), axes=F, xlab="", ylab="")
legend("top", horiz = TRUE, legend = names(plot_colours), fill = plot_colours, ,pt.cex=0.75, cex=0.75 )

mtext(basicID, outer=TRUE,cex=1)

write("Bowtie QC plots complete",stderr())

devname <- dev.off()
