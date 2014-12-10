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


### Transforming the raw sample bowtie results into GRanges class for subseqeunt alignments.

library("ShortRead")

SampleGRs <- list()

maxbowtie <- 0 # Maximum matches reported by Bowtie
samFormat <- FALSE
step <- 'genomic' # Want this script to work for both alignments to the genome and specified features will default to 'genomic' to ensure it is compatible with earlier versions.
debugme <- FALSE # Print debugging information

#print (class(samFormat))

########################
### ARGUMENT SIFTING ###
########################


args <- R.utils::commandArgs(asValue=TRUE)

if(!is.null(args$step)){
   step <- args$step
   write(paste("\nConfigured Bowtie conversion to deal with",step,"alignments"),stderr())
}else{
   write(paste("\nConfigured Bowtie conversion to deal with",step,"alignments"),stderr())
}

if(!is.null(args$input)){
   inputDir <- args$input 
}else{
   stop("Require a specified input directory: --input=<DIRECTORY>")
}

if(!is.null(args$output)){
   outputDir <- args$output
}else{
   stop("Require a specified output directory: --output=<DIRECTORY>")
}

if(!is.null(args$chrLen)){
   chrLenFile <- args$chrLen
}else{
   stop("Require a specified chromosome length file --chrLen=<FILE>")
}

if(!is.null(args$maxHits)){
   maxbowtie <- as.integer(args$maxHits)
}else{
   write("\nNo specification of the maximum allowed Bowtie Hits (--maxHits=<INTEGER>). Can not check if a limit has been reached",stderr())
}

if(!is.null(args$procSummary)){
   procSummary <- args$procSummary
}else{
   stop("Require a PROCESSED summary file (total.saved.trimmit.tab) --procSummary=<FILE>")
}

if (!is.null(args$sam)){
   samFormat <- TRUE
   write("\nProcessing Bowtie results in SAM format\n",stderr())
}else{
   write("\nProcessing Bowtie results in standard Bowtie format\n",stderr())
}

if (!is.null(args$debug)){
   debugme <- TRUE
}

if (debugme){write("\n### DEBUGGING ###\n\n",stderr())}

#stopifnot(FALSE)

###############################
### Identifying input files ###
###############################

if(samFormat){
   library(Rsamtools)    ### Only require sam tools if SAM format to be used
   inputFileSet <- dir(inputDir,pattern = "bowtie.unique.output.sort.bam$" ,full.names=TRUE)
}else{
   inputFileSet <- dir(inputDir,pattern = "bowtie.output.bwt.gz$" ,full.names=TRUE)
}
if(debugme){write(inputFileSet, stderr())}

####################################
### Reading in Chromosome Length ###
####################################

chromLengthTab <- read.table(chrLenFile, sep="\t", as.is=TRUE, header = TRUE,colClasses = c("character", "integer"))
chromLengthTab <- chromLengthTab[order(chromLengthTab[,1]),]
chrLens <- chromLengthTab[,2]
names(chrLens) <- chromLengthTab[,1]

#print(chrLens)

#stopifnot(FALSE)

############################################
### Reading in the processed read totals ###
############################################

procTotals <- read.table(procSummary, sep = "\t",row.names=1, as.is = TRUE, header = FALSE,col.names = c("barcode","total"), colClasses = c("character", "integer"))

#######################################################
### Creating tables to collate all QC data together ###   ######################## UNIVERSAL
#######################################################

if(step == "genomic"){
   mappingComparison <- as.data.frame(matrix(data=0, nrow=2, ncol = length(inputFileSet), dimnames = list("rows" = c("Mappable", "UniquelyMappable"), "col" = sub(".+\\.([[:alpha:]]+)\\.bowtie.+","\\1",inputFileSet) )))
}else if (step == "repeat"){
   mappingComparison <- as.data.frame(matrix(data=0, nrow=2, ncol = length(inputFileSet), dimnames = list("rows" = c("Mappable", "UniquelyMappable"), "col" = sub(".+\\.([[:alpha:]]+)\\.[[:alnum:]_]+\\.bowtie.+","\\1",inputFileSet) )))
}else{
   stop("Unrecognised 'step' specified")
}


mapFile <- paste(outputDir,"/bowtie.total.mapping.tab",sep="") 

####################################################
### Converting each file to GenomicRanges object ###
####################################################

for(inputFile in inputFileSet){
#   stopifnot(FALSE)
   write(paste("Working with file:", inputFile),stderr())

   ######################
   ### Arrange output ###
   ######################
   inputBasePattern <- paste(inputDir,"/(.+)\\.bowtie.+", sep="")
   fileBase <- as.character(sub(inputBasePattern,"\\1", inputFile)) 

   outFile <- paste(outputDir,"/",fileBase,".bowtie.conv.GR.RData", sep = "")
   hitFREQFile <- paste(outputDir,"/",fileBase,".bowtie.hitFreq.tab", sep = "")
   
   if(step == 'genomic'){
      chromASSIGNFile <- paste(outputDir,"/",fileBase,".bowtie.chromLoc.tab", sep = "")   
   }

   if(debugme){write(paste("Bowtie GRanges Output:", outFile, sep=""),stderr())}
   if(debugme){write(paste("Bowtie Hit Frequency Output:", hitFREQFile, sep=""),stderr())}

   ###################################
   ### Reading in SAM format files ###
   ###################################
   #stopifnot(FALSE)
   if(samFormat){
      what <- c("qname","rname","strand","pos","qwidth")
      flag <- scanBamFlag(isUnmappedQuery=FALSE)
      param <- ScanBamParam(what=what, flag=flag)
      bamHITS <- scanBam(inputFile,param=param)[[1]]
      input <- data.frame(bamHITS,stringsAsFactors=FALSE) 
      colnames(input) <- c("Read","Seq","Strand","Start","Width")
      #stopifnot(FALSE)
      if (! all (levels(input$Seq) %in% names(chrLens))){stop("There is a discrepency between the chromosome length file and the chromosomes represented in the SAM file")}
      #if (! all (names(chrLens) %in% levels(input$Seq))){stop("There is a discrepency between the chromosome length file and the chromosomes represented in the SAM file")}
  }else{ 

   # indels
   # cigar format - appears to be incorrect in Bowtie SAM files
   # use of qwidth for alignment width - if gapped alignments are unsupported then width will always equal the length of the alignment.
   # soft-clipping - Bowtie does not support gapped alignments - is scoft clipping considered a gap?

   #################################
   ### Reading in bowtie results ###
   #################################

      input <- read.table(gzfile(inputFile),header=FALSE, sep = "\t", 
           col.names = c("Read", "Strand", "Seq", "Start", "AltLocs", "Mismatch"), 
           colClasses = c("character","character","character","numeric","numeric","character")) ### Note offset = 0 for coordinates.
      input$Width <- as.numeric(sub(".+_w(\\d+)_x.+","\\1",input$Read)) ### Check R regular expressions - what if there is no counts included? Chance of confusion due to names?
   }
   if(debugme){write(paste("The bowtie input has:",nrow(input),"rows and:",ncol(input),"columns.", sep = " "),stderr())} 
   ############################################################################################
   ### Rearranging bowtie table columns to include addditional data found in sequence names ###
   ############################################################################################

   input$Counts <- as.numeric(sub(".+_x(\\d+)$","\\1",input$Read))

   # Loci number based on sequence name frequency
   locNumber <- base::table(input$Read)
   input$Locations <- as.integer(locNumber[input$Read])
   # input$Locations <- input$AltLocs +1 #### Depricated 230511 to replace with more standard method which considers the number of each sequence name instead. 
  
   ######################################################################################################  
   ### Check whether read mapping may have exceeded maximum reporting instances - Loss of information ###  ### Altered (180511) to become a simple check that Bowtie is not exceeding its -m boundaries or that the AltLocs column isn't being weird.
   ######################################################################################################
   
   if(debugme){write("Checking the number of hits reported by Bowtie to ensure Bowtie discarding multiple matches correctly",stderr())}

   #system.time(
   if(! maxbowtie == 0){
      if(debugme){write(paste("Checking to see if the maximum hit number (",maxbowtie,") has been exceeded in error"),stderr())}
      if(max(input$Locations) > maxbowtie){
         stop("WARNING - Bowtie may have exceeded maximum reportable matches in some cases!")
      }else{
         if(debugme){write("Bowtie does not appear to have exceeded its limits of maximum reportable matches",stderr())}
      }
   }
   #)
   #input$Read <- sapply(input$Read,function(x){strsplit(x,"\\.")[[1]][1]})
   
   ###############################################################
   ### Saving aspects of the Bowtie result for plotting and QC ###
   ###############################################################
   #stopifnot(FALSE)
   
   write("Recording Bowtie mapping loci counts for QC",stderr())
   
   ### Recognising the current barcode
#   stopifnot(FALSE)
   if(step == "genomic"){
      barnicle <- as.character(sub(".+\\.([[:alpha:]]+)\\.bowtie.+","\\1",inputFile))
   }else if(step == "repeat"){
      barnicle <- as.character(sub(".+\\.([[:alpha:]]+)\\.[[:alnum:]_]+\\.bowtie.+","\\1",inputFile)) 
   }else{
      stop("Unrecognised 'step' specified")
   }

   nonDups <- !duplicated(input$Read)

   ReadLocs   <- Rle(input$Locations)
   ReadCounts <- Rle(input$Counts)

   nonDupLocs <- ReadLocs[nonDups]
   nonDupCounts <- ReadCounts[nonDups]

   readLocsSummary  <- as.data.frame(matrix(data=NA, nrow=max(nonDupLocs), ncol =1, dimnames = list("rows" = seq_len(max(nonDupLocs)), "col" = c("Frequency") )))

   #system.time(
   for (x in rownames(readLocsSummary)){
      readLocsSummary[x,"Frequency"] <- sum(nonDupCounts[nonDupLocs==x])
   }
   #)
   readLocsSummary <- (readLocsSummary/procTotals[barnicle,"total"])*100

   write.table(readLocsSummary, file = hitFREQFile, row.names = TRUE, quote =FALSE, sep = "\t")

   ### Summary of unique vs. multiple mappings

   mappingComparison["Mappable",barnicle] <- sum(nonDupCounts)
   mappingComparison["UniquelyMappable", barnicle] <- sum(nonDupCounts[nonDupLocs==1])
   
   ###########################################################################################
   ### Original dataframe manipulation - Removed - When dataframe gets too big - UNWIELDY ####
   ###########################################################################################


   ### Distribution of reads over locations - Converted to a percentage of processed reads
#   readLocs <- input[!duplicated(input$Read),c("Read","Locations","Counts")] 
#   readLocsSummary  <- as.data.frame(matrix(data=NA, nrow=max(readLocs[,"Locations"]), ncol =1, dimnames = list("rows" = seq_len(max(readLocs[,"Locations"])), "col" = c("Frequency") )))   
#   
#   system.time(
#   for (x in rownames(readLocsSummary)){
#      readLocsSummary[x,"Frequency"] <- sum(readLocs[readLocs$Locations == x,"Counts"])
#   }
#   )
#   
#   readLocsSummary <- (readLocsSummary/procTotals[barnicle,"total"])*100
#   write.table(readLocsSummary, file = hitFREQFile, row.names = TRUE, quote =FALSE, sep = "\t")
#   
#   ### Summary of unique vs. multiple mappings
#
#   mappingComparison["Mappable",barnicle] <- sum(readLocs$Counts) 
#   mappingComparison["UniquelyMappable", barnicle] <- sum(readLocs[readLocs$Locations == 1,"Counts"])
#


   
   
   ######################################
   ### Converting to a GRanges Output ###
   ######################################
   #stopifnot(FALSE)
   write("Converting mapping tables to GRanges results",stderr())

   if(samFormat){
      # Already 1-based in the case of sam format - no alteration to start required 
      GR  <-   GRanges(
         seqnames   = Rle(factor(input$Seq,levels=names(chrLens))),
         ranges     = IRanges(start=input$Start, width=input$Width),
         strand     = Rle(factor(input$Strand)),
         seqlengths = chrLens,
         locations  = input$Locations,
         count      = input$Counts/(input$Locations),
         name       = input$Read
      )    
   }else{
      # Note that the Start index must be corrected as Bowtie is 0-start based! Add one to first base
      GR  <-   GRanges(
         seqnames   = Rle(factor(input$Seq,levels=names(chrLens))),
         ranges     = IRanges(start=input$Start+1, width=input$Width),
         strand     = Rle(factor(input$Strand)),
         seqlengths = chrLens,
         locations  = input$Locations,
         count      = input$Counts/(input$Locations),
         name       = input$Read
      )
   }


   ##############################################
   ### Collecting Chromosome mapping of reads ###
   ##############################################
   ### NOTE: Faster conducted on GRanges objects than dataframes
   
   if (step == "genomic"){
      write("Assigning reads to chromosomes",stderr())
   
      chromosomeSummary <- as.data.frame(matrix(data=0, nrow=4, ncol = length(chrLens), dimnames = list("rows" = c("Mappable_Sense", "UniquelyMappable_Sense", "Mappable_AntiSense", "UniquelyMappable_AntiSense"), "col" =  names(chrLens))))
      
      chrNames <- colnames(chromosomeSummary)
   
      #system.time(c(
      chromosomeSummary["Mappable_Sense",] <- sapply(chrNames, function(x){
                  MSsum <- sum(elementMetadata(GR[seqnames(GR)==x & strand(GR)=="+",])[,"count"])
                  return(MSsum) 
      })
      chromosomeSummary["UniquelyMappable_Sense",] <- sapply(chrNames, function(x){
                  UMSsum <- sum(elementMetadata(GR[seqnames(GR)==x & strand(GR)=="+" & elementMetadata(GR)[,"locations"] == 1,])[,"count"])
                  return(UMSsum)
      })
      chromosomeSummary["Mappable_AntiSense",] <- sapply(chrNames, function(x){
                  MASsum <- sum(elementMetadata(GR[seqnames(GR)==x & strand(GR)=="-",])[,"count"])
                  return(MASsum) 
      })
      chromosomeSummary["UniquelyMappable_AntiSense",] <- sapply(chrNames, function(x){
                  UMASsum <- sum(elementMetadata(GR[seqnames(GR)==x & strand(GR)=="-" & elementMetadata(GR)[,"locations"] == 1,])[,"count"])
                  return(UMASsum)
      })
   
      chromosomeSummary <- (chromosomeSummary/procTotals[barnicle,"total"])*100   
      #))   
   #      chromosomeSummary["UniquelyMappable_Sense", chrName] <- sum(as.numeric(apply(input,1, function(x){
   #         if (x["Strand"] == "+" && x["Locations"]==1 && x["Seq"]==chrName){
   #            return(as.numeric(x["Counts"]))
   #         }else{
   #            return(0)
   #         }
   #      })))
   #
   #      chromosomeSummary["UniquelyMappable_AntiSense", chrName] <- sum(as.numeric(apply(input,1, function(x){
   #         if (x["Strand"] == "-" && x["Locations"]==1 && x["Seq"]==chrName){
   #            return(as.numeric(x["Counts"]))
   #         }else{
   #            return(0)
   #         }
   #      })))
   #
   #      chromosomeSummary["Mappable_Sense", chrName] <- sum(as.numeric(apply(input,1, function(x){
   #         if (x["Strand"] == "+" && x["Seq"]==chrName){
   #            return(as.numeric(x["Counts"])/as.numeric(x["Locations"]))
   #         }else{
   #            return(0)
   #         }
   #      })))
   #
   #
   #      chromosomeSummary["Mappable_AntiSense", chrName] <- sum(as.numeric(apply(input,1, function(x){
   #         if (x["Strand"] == "-" && x["Seq"]==chrName){
   #            return(as.numeric(x["Counts"])/as.numeric(x["Locations"]))
   #         }else{
   #            return(0)
   #         }
   #      }))) 
   #   }
      write("Writing results to file",stderr())
   #system.time(
      write.table(chromosomeSummary, file=chromASSIGNFile, row.names = TRUE, quote = FALSE, sep= "\t")
      #)
   }else{write("Writing results to file",stderr())}
   ### Saving the new Granges construct ###
   save(GR, file = outFile)
   
}

write.table(mappingComparison, file = mapFile, row.names = TRUE, quote = FALSE, sep = "\t")

rm(input)
gc()

