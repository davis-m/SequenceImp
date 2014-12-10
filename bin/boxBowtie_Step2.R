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

library(GenomicRanges)
debugme <- FALSE  # Print debugging information


########################
### Argument Sifting ###
########################

args <- R.utils::commandArgs(asValue=TRUE)

if(!is.null(args$annot)){
   annotationFile <- args$annot
}else{
   stop("Require a GRangeList of annotation: --annot=<GRANGESLIST.RData>")
}

if(!is.null(args$input)){
   inputDir <- args$input
}else{
   stop("Require a specified input directory containing bowtie GRanges objects: --input=<DIRECTORY>")
}

if(!is.null(args$output)){
   outputDir <- args$output
}else{
   stop("Require a specified output directory for depositing class summary objects: --output=<DIRECTORY>")
}

if (!is.null(args$debug)){
   debugme <- TRUE
}
if (debugme){write("\n### DEBUGGING ###\n\n",stderr())}


##############
### Inputs ###
##############

inputGRs <- dir(inputDir, pattern= ".bowtie.conv.GR.RData", full.names=TRUE)

if(debugme){write(paste("Recognised input files:", inputGRs),stderr())}
write(paste("Performing comparison to the annotation file:", annotationFile,"\n"),stderr())

####################################
### Loading the Annotation Files ###
####################################

load(annotationFile)

sumValues <- list()

for (selectedGR in inputGRs){
   
   write(paste("Finding overlaps with annotation for",selectedGR),stderr())
#   stopifnot(FALSE)

   #########################################
   ### Loading the Bowtie GRanges object ###
   #########################################
   
   load(selectedGR)
   
   ############################################################################
   ### Only comparing overlaps with the first base of the mapped shortreads ###
   ############################################################################
   
   sampleGR <- flank(GR,width=-1,start=T)
   
   ###############################################################
   ### Find the overlaps and total the appropriate read counts ###
   ###############################################################
   
   barthing <- sub(".*\\.(([[:alnum:]])+).bowtie\\.conv\\.GR\\.RData$",'\\1',selectedGR)
   sumValues[[barthing]] <- data.frame(row.names = names(EnsemblAnnotReduced[['sense']]))

   for(strand in names(EnsemblAnnotReduced)){
#      stopifnot(FALSE)
      tempZero <- rep(0,length(names(EnsemblAnnotReduced[[strand]])))
      names(tempZero) <- names(EnsemblAnnotReduced[[strand]]) 
      overlapMat <- findOverlaps(EnsemblAnnotReduced[[strand]],sampleGR) # Subject is the Bowtie range matrix, Query is the GRangeLst of all of the different RNA classes
      sumPerCategory <- sapply(split(as.numeric(values(sampleGR)[,"count"])[subjectHits(overlapMat)], queryHits(overlapMat)),sum)
      tempZero[as.numeric(names(sumPerCategory))] <- sumPerCategory
      columnName <- base::paste(barthing,strand, sep="_")
      sumValues[[barthing]][names(tempZero),columnName] <- tempZero
   }
   sampleSums <- apply(sumValues[[barthing]],2,sum)
   for(choice in colnames(sumValues[[barthing]])){
      sumValues[[barthing]]['total',choice]<-sampleSums[choice]
   }
   
}

###############################
### WRITING RESULTS TO FILE ###
###############################

write("\nWriting results to file",stderr())

outFile <- paste(outputDir,"bowtie.class.summary.RData", sep ="")

save(sumValues,file = outFile)

