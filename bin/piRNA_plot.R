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
library(IRanges)

### TO DO:
# Present all parametres in plot - Note max mapping numbers for reads mapping to canonical sequence but not genome.

debugme <- FALSE 

args <- R.utils::commandArgs(asValue=TRUE)

if (!is.null( args$input)){
   inputDir <- args$input
}else{
   stop ("Require an input directory containing Bowtie QC tables : --input=<DIRECTORY>")
}

if (!is.null(args$output)){
   outDir <- args$output
}else{
   stop ("Require an output directory into which the QC plots will be deposited: --output=<DIRECTORY>")
}

if (!is.null(args$basicID)){
   basicID <- args$basicID
}else{
   stop ("Require a base name to uniquely identify the results pdf: --basicID=<STRING>")
}

if (!is.null(args$repName)){
   repName <- args$repName
}else{
   stop ("Require a repeat name for which to organise the results pdf: --repName=<STRING>")
}

if (!is.null(args$genMap)){
   genomeSummary <- args$genMap
}else{
   stop ("Require a file summarising all reads which map to the genome with repeat analysis criteria: --genMap=<FILE>")
}

if(!is.null(args$alignRec)){
   alignRecFile <- args$alignRec
}else{
   stop ("Require an alignment record: --alignRec=<FILE>")
}

if (!is.null(args$debug)){
   debugme <- TRUE
}

if (debugme){write("\n### DEBUGGING ###\n\n",stderr())}



### Output file
outFile <- paste(outDir,"/",basicID,".",repName,".mapping_QC.pdf", sep ="")


#######################
### If no alignment ###
#######################
allAlignments <- read.table(file=alignRecFile, sep="\t", header=TRUE, row.names=1, as.is=TRUE)
if (sum(allAlignments$Alignments)==0){
   pdf( file = outFile )
   par(omi=c(0.3,0.3,0.3,0.3))
   plot(NA,xlim=c(0,1),ylim=c(0,1),xaxt="n",yaxt="n",xlab="",ylab="",bty="n")
   mtext(paste(basicID,repName), outer=TRUE, cex=1) 
   text(x=0.5,y=0.5,labels="No alignments found",font=2)
   dev.off()
   quit(save="no",status=0)
}

########################
### Arranging Inputs ###
########################

repSummary <- dir(inputDir, pattern = "bowtie.total.mapping.tab", full.names=TRUE)
if(length (repSummary) != 1){
   stop("Detected the wrong number of 'bowtie.total.mapping.tab' files in input directory")
}

#genomeSummary <- dir(inputDir, pattern = "genome_mapped.txt", full.names=TRUE)
#if(length (genomeSummary) != 1){
#   stop("Detected the wrong number of 'genome_mapped.txt' files in input directory")
#}

lengthSummary <- dir(inputDir, pattern = ".length_distrib.tab", full.names=TRUE)
if(length (lengthSummary) != 1){
   stop("Detecting the wrong number of 'length_distrib.tab files in input directory'")
}

coveragePattern <- paste(repName,".coverage.RData", sep = "")
coverageSummary <- dir(inputDir, pattern = coveragePattern, full.names=TRUE)
if(length (coverageSummary) != 1){
   stop(paste("Detected the wrong number of",coveragePattern,"files in input directory"))
}

overlapPattern <- paste(repName,".pingpong.tab", sep="") 
overlapSummary <- dir(inputDir, pattern = overlapPattern, full.names=TRUE)
if(length (overlapSummary) != 1){
   stop(paste("Detected the wrong number of",overlapPattern,"files in input directory"))
}

# stopifnot(FALSE)

#####################################################
### Load coverage data to determine sample number ###
#####################################################

load(coverageSummary[1])
sampleNo <- length(names(coverage))

##################
### Set up PDF ###
##################

write ("Creating QC PDF",stderr())


pdf( file = outFile ,height=(9.6+(sampleNo*3)), width=7.6)

layoutStuff <- c(1,2,3,4,5,6)
for (i in seq_len(sampleNo)){
   layoutStuff <- c(layoutStuff,rep(i+6, 2))
}

layout(matrix(layoutStuff, nrow = 3+sampleNo, ncol= 2, byrow=TRUE))
par(omi=c(0.3,0.3,0.3,0.3))

############################
### Repeat mapping reads ###
############################

write ("Recording number of repeat mapping reads",stderr())

repTable <- read.table(repSummary[1], sep="\t", as.is=TRUE, header=TRUE)
genomeTable <- read.table(genomeSummary[1], sep="\t",row.names=1, header=TRUE)

# Add back samples with no read alignments
missing_samples <- rownames(allAlignments[which(allAlignments$Alignments == 0),,drop=FALSE])
if(length(missing_samples)>0){
   if (any(missing_samples %in% colnames(repTable))){stop("Repeat information present for sample lacking alignments")}
   missing_sample_data<- as.data.frame(matrix(data=0,nrow=2,ncol=length(missing_samples),dimnames=list("rows"=rownames(repTable),"cols"=missing_samples)))
   repTable <- cbind(repTable,missing_sample_data)
}

barplot(main = "Canonical repeat mapping reads",height=as.matrix(repTable), beside=TRUE,col=c('green','red'), ylab='Read Number',xlab='Sample')
plot(NA, xlim=c(0,1), ylim=c(0,1),axes=F, xlab="", ylab="")
legend('topleft', rownames(repTable), fill=c('green', 'red'))

################################
### Read length distribution ###
################################

write ("Plotting repeat-mapping read length distributions",stderr())

lenPreTable <- read.table(lengthSummary[1], sep="\t", as.is=TRUE, header=TRUE, check.names = FALSE)
lenTable <- as.data.frame(matrix(data=0, ncol=length(min(as.numeric(colnames(lenPreTable))):max(as.numeric(colnames(lenPreTable)))), nrow=nrow(lenPreTable), 
                   dimnames=list("row"=rownames(lenPreTable),"col"=as.character(min(as.numeric(colnames(lenPreTable))):max(as.numeric(colnames(lenPreTable)))))))
lenTable[,colnames(lenPreTable)] <- lenPreTable
plotcols <- rainbow(nrow(lenTable))
names(plotcols) = rownames(lenTable)
plot(as.numeric(lenTable[1,])~as.numeric(colnames(lenTable)), main="Length distribution of reads mapping\nto the canonical sequence", type = "o",xlab = 'Read length', 
      ylab='Read number', col = plotcols[rownames(lenTable[1,,drop=FALSE])], ylim= c(0,max(as.matrix(lenTable))), bty = 'n',pch=20)
if(nrow(lenTable) > 1){
   write("Length distribution plot being drawn for multiple samples",stderr())
   for (selected in 2:nrow(lenTable)){
      lines(as.numeric(lenTable[selected,])~as.numeric(colnames(lenTable)),type = "o", col = plotcols[rownames(lenTable[selected,,drop=FALSE])],pch=20)
   }
}else{
   write("Single sample: Only single line for the length distribution data",stderr())
}
plot(NA, xlim=c(0,1), ylim=c(0,1),axes=F, xlab="", ylab="")
legend('topleft', legend=names(plotcols), fill = plotcols)

######################
### Ping Pong Plot ###
######################

write ("Plotting sense-antisense read overlap lengths",stderr())

overlapPreTable <- read.table(overlapSummary[1], sep="\t",row.names=1, header=TRUE,check.names=FALSE)
if (! all(is.na(overlapPreTable))){
   overlapTable <- as.data.frame(matrix(data=0, ncol=length(min(as.numeric(colnames(overlapPreTable))):max(as.numeric(colnames(overlapPreTable)))), nrow=nrow(overlapPreTable),
                      dimnames=list("row"=rownames(overlapPreTable),"col"=as.character(min(as.numeric(colnames(overlapPreTable))):max(as.numeric(colnames(overlapPreTable)))))))
   overlapTable[,colnames(overlapPreTable)] <- overlapPreTable
   
   plot(as.numeric(overlapTable[1,])~as.numeric(colnames(overlapTable)),main="Sense vs. antisense nucleotide overlap of reads\nmapping to canonical repeat", type = "l",
      xlab = 'Sense-antisense nucleotide overlaps', ylab='Proportion', col = plotcols[rownames(overlapTable[1,])], ylim= c(0,max(as.matrix(overlapTable))), bty = 'n')
   
   if(nrow(overlapTable) > 1){
      write("Ping-pong plot being drawn for multiple samples",stderr())
      for (zzz in 2:nrow(overlapTable)){
         lines(as.numeric(overlapTable[zzz,])~as.numeric(colnames(overlapTable)),type = "l", col = plotcols[rownames(overlapTable[zzz,])])
      }
   }else{
      write("Single sample: Only single line for the ping-pong data",stderr())
   }
   plot(NA, xlim=c(0,1), ylim=c(0,1),axes=F, xlab="", ylab="")
   legend('topleft', legend=names(plotcols), fill = plotcols)
}else{
   plot(NA, xlim=c(0,1), ylim=c(0,1),axes=F, xlab="", ylab="", main="Sense vs. antisense nucleotide overlap of reads\nmapping to canonical repeat")
   text(x=0.5, y=0.5, labels="No overlaps found for samples in this sequencing lane")
   plot(NA, xlim=c(0,1), ylim=c(0,1),axes=F, xlab="", ylab="")
}
################################
### Mapping accross consensi ###
################################

write ("Plotting read coverage across canonical repeat",stderr())

replen <- length(coverage[[1]][[1]][[1]])
cols <- c("sense"='red', 'antisense'='blue') 

axislim <- signif(1.1 * (max(unlist(lapply(coverage,function(x){lapply(x,function(y){max(y)})})))),2)
steps <- signif(axislim/10,1)

yaxisPoints <- unique(c(seq(axislim,0,-1 * steps),seq(-axislim,0,steps)))
yaxisLabs <- abs(yaxisPoints)


for (bar in names(coverage)){
   plot(main=paste("Coverage of canonical repeat for the",bar,"sample\ncorrected to the number of genome mapping reads"), NA, ylim =c(-axislim,axislim), xlim= c(0,replen) , ylab="Coverage per million" , xlab="Repeat position", axes = FALSE, cex.axis = 0.75, las = 2)
   axis(side = 2, labels= yaxisLabs, at=yaxisPoints, cex.axis=0.75, las=2)
   axis(side = 1)
   if (length(coverage[[bar]][['Sense']][[1]]) != replen){stop(paste("The sense strand coverage vectors for",bar,"are of a different length to other coverage vectors"))}
   lines(as.numeric(coverage[[bar]][['Sense']][[1]]),type="l", col = cols['sense'])
   if (length(coverage[[bar]][['Antisense']][[1]]) != replen){stop(paste("The antisense strand coverage vectors for",bar,"are of a different length to other coverage vectors"))}
   lines(-(as.numeric(coverage[[bar]][['Antisense']][[1]])),type="l", col = cols['antisense'])
   legend('topright', legend=names(cols), fill=cols)
   legend('bottomright',paste(bar, genomeTable[bar,'Reads'], sep=": "),title= "Total reads that map to the genome at least once", bty='n')
   abline(h=0)
}
mtext(paste(basicID,repName), outer=TRUE, cex=1)
devname <- dev.off()
