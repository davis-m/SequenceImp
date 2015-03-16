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


### A script to take a small RNA sequencing GRanges result and convert it into a miRBase sequence count table
### This will produce a table including the depth of unique reads and duplicate reads which map to specified features
### The table will also include feature lengths and mapped library sizees to allow ease of compatability with BaySeq.

library(GenomicRanges)

debugme <- FALSE
separate <- FALSE # 020413
proportional <- FALSE # 030413
collapse_type <- ""

args <- R.utils::commandArgs(asValue=TRUE)

if(!is.null(args$annot)){
   annotationFile <- args$annot
}else{
   stop("Require a list of GRanges objects corresponding to miRNA annotation: --annot=<LIST.RData>")  # A list containing the miRNA mature and precursor annotation
}

if(!is.null(args$inDir)){
   inDir <- args$inDir
}else{
   stop("Require a specified input bowtie directory containing mapped read GRanges objects: --inDir=<DIRECTORY/>")
}

if(!is.null(args$outDir)){
   outDir <- args$outDir
}else{
   stop("Require a directory for the output: --outDir=<DIRECTORY/>")
}

if(!is.null(args$overlap)){
   minOverlap <- as.integer(args$overlap)
}else{
   stop("Require a minimum overlap value for reads and annotation: --overlap=<INTEGER>")
}

if (!is.null(args$debug)){
   debugme <- TRUE
}

if (!is.null(args$nomerge)){
   separate <- TRUE
}

if (!is.null(args$proportional)){
   proportional <- TRUE
}

# 170414 - Allow the collapse of miRNA counts based on sequence rather than ID.   
if (!is.null(args$collapsetype)){     
   collapse_type <- args$collapsetype
}else{
   if(! separate){
      stop("Require a specified criteria for collapsing/merging miRNA counts (--collapse_type=<mature_id/sequence>)")
   }
}

if (!is.null(args$treatOverlap)){
   treatOverlap <- args$treatOverlap
}else{
   stop("Require a specified behaviour if miRNA loci overlapping '--treatOverlap=<VALUE>'")
}

if(is.na(minOverlap)){
   stop("--overlap was not an integer")
}

if (treatOverlap != "ignore" && treatOverlap != "remove" && treatOverlap != "merge"){
   stop("Acceptable actions if miRNA loci overlap are 'ignore', 'remove' and 'merge'")
}

if (debugme){write("\n### DEBUGGING ###\n",stderr())}

#######################################
### Summarise parameters being used ###
#######################################
write(paste(),stderr())

write(paste("Hunting for overlaps of",minOverlap,"nucleotides between reads and annotation"), stderr())
write(paste("Treating overlapping miRNA loci with the following protocol:",treatOverlap),stderr())

if(proportional){
   write(paste("Ignoring non-miRNA alignments and dividing depth of multi-mapping reads between miRNAs to reflect the depth of uniquely mapped reads"),stderr())
}else{
   write(paste("Dividing reads between multiple alignments equally"),stderr())
}

if(separate){
   write("All loci considered independently and attributed a unique ID",stderr())
   if(nchar(collapse_type) > 0){
      write(paste("WARNING: Ignoring specified collapse type:",collapse_type),stderr())
   }
}else{
   if(collapse_type == "mature_id"){
      write("miRNAs sharing a miRBase identifier merged to single count",stderr())
   }else if(collapse_type == "sequence"){
      write("miRNAs sharing an identical sequence in miRBase merged to single count",stderr())
   }else{
      stop("Require a recognised method with which to merge miRNA counts.")
   }
}

##########################
### Loading annotation ###
##########################

write(paste("\nLoading miRBase annotation:",annotationFile),stderr())
load(annotationFile) # miRBaseGRs

miRs <- miRBaseGRs[["mature"]]         ### Still contains multiple entries for mature miRNAs from several precursors.
#pres <- miRBaseGRs[["precursor"]]     ### Removed 150714 Removed for now - room to return later.

#################################################################
### Clean up overlapping loci in the miRBase miRNA annotation ### # 040413
#################################################################

### 020413: Added to form unique IDs #! NOTE: In miRBase precursor names are not always unique!
### Added IDs to miRs and DF to allow them to be merged and compared.

write("\nIdentifying overlapping miRBase annotation",stderr())

values(miRs)$uniqueID <- paste(elementMetadata(miRs)[,"precursor"],".",elementMetadata(miRs)[,"mature"],":",as.character(seqnames(miRs)),":",start(flank(miRs, -1)), sep="")  # Unique ID consists of: precursor ID, mature ID, chromosome name and start coordinate (strand specific).

if(any(duplicated(values(miRs)$uniqueID))){
   stop("miRBase identifiers indistinguishable when creating a unique ID for mature locus")
}

knownOverlaps <- c()

miROverlaps <- findOverlaps(miRs, ignoreSelf=TRUE)

overlappingMiRs <- miRs[unique(queryHits(miROverlaps))]

if (length(miROverlaps) > 0){
   if (treatOverlap == "ignore"){
      write(paste("Identified overlapping miRNAs:", paste(values(overlappingMiRs)$uniqueID, collapse = ", "), sep = " "), stderr())   
      write("Proceeding with no alterations to miRNA annotation as instructed",stderr())
   }else if (treatOverlap == "remove"){
      write(paste("Identified overlapping miRNAs:", paste(values(overlappingMiRs)$uniqueID, collapse = ", "), sep = " "), stderr())   
      write("Removing overlapping mature miRNA coordinates from the miRNA annotation as instructed",stderr())
      
      miRs <- miRs[-unique(queryHits(miROverlaps))]
   
   } else if (treatOverlap == "merge"){
      write(paste("Identified overlapping miRNAs:", paste(values(overlappingMiRs)$uniqueID, collapse = ", "), sep = " "), stderr())   
      write("Merging overlapping mature miRNA coordinates from the miRNA annotation as instructed",stderr())
      
      reducedRegions <- reduce(overlappingMiRs)
      namingRefs <- findOverlaps(reducedRegions, overlappingMiRs)
      namingRefs <- split(subjectHits(namingRefs),queryHits(namingRefs))
      
      values(reducedRegions)$mature <- as.character(NA)
      values(reducedRegions)$precursor <- as.character(NA)
      values(reducedRegions)$uniqueID <- as.character(NA)
   
      for(thisOverlap in names(namingRefs)){
         values(reducedRegions[as.numeric(thisOverlap)])$mature <- paste(values(overlappingMiRs[as.numeric(namingRefs[[thisOverlap]])])$mature ,collapse="//")
         values(reducedRegions[as.numeric(thisOverlap)])$precursor <- paste(values(overlappingMiRs[as.numeric(namingRefs[[thisOverlap]])])$precursor ,collapse="//")
         values(reducedRegions[as.numeric(thisOverlap)])$uniqueID <- paste(values(overlappingMiRs[as.numeric(namingRefs[[thisOverlap]])])$uniqueID ,collapse="//")
      }
   
      miRs <- miRs[-unique(queryHits(miROverlaps))]
      miRs <- c(miRs,reducedRegions)

   }else{
      write(paste("Identified overlapping miRNAs:", paste(values(overlappingMiRs)$uniqueID, collapse = ", "), sep = " "), stderr())   
      stop(paste("Unrecognised behaviour requested when overlaps identified:", treatOverlap))
   }
}

### YUP!!

#stopifnot(FALSE)
############################################
### Set up dataframe for evaluating data ###
############################################

write("\nConstructing dataframe for collating mature miRNA depth",stderr())

miRCountFrame <- as.data.frame(matrix(data=NA, nrow= length(elementMetadata(miRs)[,"mature"]), ncol = 3, dimnames= list("rows"= NULL, "cols"= c("Mature", "Precursor", "Width"))))     # The data frame comntains a space for each mature miRNA locus

miRCountFrame[,"Mature"] <- elementMetadata(miRs)[,"mature"]
miRCountFrame[,"Precursor"] <- elementMetadata(miRs)[,"precursor"]
miRCountFrame[,"Width"] <- width(miRs)


if (separate){    ### Added 020413
   write("Maintaining each mature miRNA locus as a separate entity in results table",stderr())
   

   miRCountFrame$Unique_ID <- values(miRs)$uniqueID
   if(any(duplicated(miRCountFrame$Unique_ID))){stop("Duplicate 'unique' IDs found for miRNAs while constructing dataframe")}
   rownames(miRCountFrame) <- miRCountFrame[,"Unique_ID"]

   miRCountFrame$Sequence <- values(miRs)$sequence

   ### Add a row to the table to record the total read depth for the samples
   
   libraryDepth<- as.data.frame(matrix(data=c("Mapped_library_depth",NA,NA,NA,NA), nrow=1, ncol=5, dimnames=list("rows"=c("Mapped_library_depth"), "cols"=c("Mature","Precursor","Width","Unique_ID","Sequence"))))
   miRCountFrame <- rbind(miRCountFrame,libraryDepth)


   #Example:
   #                                                 Mature    Precursor Width                                  Unique_ID                Sequence
   #mmu-mir-6341.mmu-miR-6341:1:12426016       mmu-miR-6341 mmu-mir-6341    23       mmu-mir-6341.mmu-miR-6341:1:12426016 CAGUGCAAUGAUAUUGUCACUAU
   #mmu-mir-206.mmu-miR-206-5p:1:20679017    mmu-miR-206-5p  mmu-mir-206    23      mmu-mir-206.mmu-miR-206-5p:1:20679017 ACAUGCUUCUUUAUAUCCUCAUA
   #mmu-mir-206.mmu-miR-206-3p:1:20679055    mmu-miR-206-3p  mmu-mir-206    22      mmu-mir-206.mmu-miR-206-3p:1:20679055  UGGAAUGUAAGGAAGUGUGUGG
   #mmu-mir-133b.mmu-miR-133b-5p:1:20682797 mmu-miR-133b-5p mmu-mir-133b    22    mmu-mir-133b.mmu-miR-133b-5p:1:20682797  GCUGGUCAAACGGAACCAAGUC
   #mmu-mir-133b.mmu-miR-133b-3p:1:20682834 mmu-miR-133b-3p mmu-mir-133b    22    mmu-mir-133b.mmu-miR-133b-3p:1:20682834  UUUGGUCCCCUUCAACCAGCUA
   #mmu-mir-30a.mmu-miR-30a-5p:1:23272274    mmu-miR-30a-5p  mmu-mir-30a    22      mmu-mir-30a.mmu-miR-30a-5p:1:23272274  UGUAAACAUCCUCGACUGGAAG

}else{
   if(collapse_type == "mature_id"){   
      write("Merging mature miRNA loci that share the miRBase ID into a single entity in results table",stderr())
   
      #######################################################################################
      ### Collect information concerning mature miRNA expression from multiple precursors ###
      #######################################################################################
      
      summaryPre <- list() # Collate precursor information for all mature names  
      summaryLen <- list() # Also collect miRNA length information per mature name
      
      for(i in seq_len(nrow(miRCountFrame))){
         if(is.null(summaryPre[[miRCountFrame[i,"Mature"]]])){
            summaryPre[[miRCountFrame[i,"Mature"]]] <- miRCountFrame[i,"Precursor"] 
            summaryLen[[miRCountFrame[i,"Mature"]]] <- miRCountFrame[i,"Width"]
         }else{
            summaryPre[[miRCountFrame[i,"Mature"]]] <- c(summaryPre[[miRCountFrame[i,"Mature"]]],miRCountFrame[i,"Precursor"]) # paste together all hairpin names for identicle miRNA mature names
            if(summaryLen[[miRCountFrame[i,"Mature"]]] != miRCountFrame[i,"Width"]){
               stop("miRBase mature miRNAs have two lengths associated with the same name!") # For each mature sequence ensure all widths are identical.
            }
         }
      }
   
      summaryPre <- sapply(summaryPre,function(x){paste(unique(x),collapse=";")})      

      ### Remove duplicate mature miRNA names (Merge counts from multiple loci).
      
      miRCountFrame <- miRCountFrame[!duplicated(miRCountFrame[,"Mature"]),]
      rownames(miRCountFrame) <- miRCountFrame[,"Mature"]
      miRCountFrame[,"Precursor"] <- unlist(summaryPre[rownames(miRCountFrame)],use.names=FALSE)
      ### Add a row to the table to record the total read depth for the samples
      
      libraryDepth<- as.data.frame(matrix(data=c("Mapped_library_depth",NA,NA), nrow=1, ncol=3, dimnames=list("rows"=c("Mapped_library_depth"), "cols"=c("Mature","Precursor","Width"))))
      miRCountFrame <- rbind(miRCountFrame,libraryDepth)
   
      #Example:
      #                                                       Mature                     Precursor Width
      #mmu-miR-206-5p                                 mmu-miR-206-5p                   mmu-mir-206    23
      #mmu-miR-206-3p                                 mmu-miR-206-3p                   mmu-mir-206    22
      #mmu-miR-133b-5p                               mmu-miR-133b-5p                  mmu-mir-133b    22
      #mmu-miR-133b-3p                               mmu-miR-133b-3p                  mmu-mir-133b    22
      #mmu-miR-30a-5p                                 mmu-miR-30a-5p                   mmu-mir-30a    22
      #mmu-miR-30a-3p                                 mmu-miR-30a-3p                   mmu-mir-30a    22
      #mmu-miR-669d-5p                               mmu-miR-669d-5p   mmu-mir-669d;mmu-mir-669d-2    22
      #Mapped_library_depth                     Mapped_library_depth                          <NA>  <NA>
   }else if(collapse_type == "sequence"){
      write("Merging mature miRNA loci that share the miRBase mature sequence into a single entity in results table",stderr())
      
      miRCountFrame$Sequence <- values(miRs)$sequence
      
      summaryPre <- list() # Collate precursor information for all sequences 
      summaryMat <- list() # Collate mature information for all sequences
      summaryLen <- list() # Also collect miRNA length information per mature name
 
      for(i in seq_len(nrow(miRCountFrame))){
         if(is.null(summaryPre[[miRCountFrame[i,"Sequence"]]])){
            summaryPre[[miRCountFrame[i,"Sequence"]]] <- miRCountFrame[i,"Precursor"]
            summaryMat[[miRCountFrame[i,"Sequence"]]] <- miRCountFrame[i,"Mature"]
            summaryLen[[miRCountFrame[i,"Sequence"]]] <- miRCountFrame[i,"Width"]
         }else{
            summaryPre[[miRCountFrame[i,"Sequence"]]] <- c(summaryPre[[miRCountFrame[i,"Sequence"]]],miRCountFrame[i,"Precursor"]) # collect together all hairpin names for identicle miRNA mature sequences
            summaryMat[[miRCountFrame[i,"Sequence"]]] <- c(summaryMat[[miRCountFrame[i,"Sequence"]]],miRCountFrame[i,"Mature"]) # collect together all hairpin names for identicle miRNA mature sequences
            if(summaryLen[[miRCountFrame[i,"Sequence"]]] != miRCountFrame[i,"Width"]){
               stop("miRBase mature miRNAs have two lengths associated with the same sequence!") # For each mature sequence ensure all widths are identical.
            }
         }
      }

      summaryPre <- sapply(summaryPre,function(x){paste(unique(x),collapse=";")})
      summaryMat <- sapply(summaryMat,function(x){paste(unique(x),collapse=";")})

      ### Remove duplicate mature miRNA names (Merge counts from multiple loci).
 
      miRCountFrame <- miRCountFrame[!duplicated(miRCountFrame[,"Sequence"]),]
      rownames(miRCountFrame) <- miRCountFrame[,"Sequence"]
      miRCountFrame[,"Mature"] <- summaryMat[rownames(miRCountFrame)]
      miRCountFrame[,"Precursor"] <- summaryPre[rownames(miRCountFrame)]
      ### Add a row to the table to record the total read depth for the samples
 
      libraryDepth<- as.data.frame(matrix(data=c("Mapped_library_depth",NA,NA,NA), nrow=1, ncol=4, dimnames=list("rows"=c("Mapped_library_depth"), "cols"=c("Mature","Precursor","Width","Sequence"))))
      miRCountFrame <- rbind(miRCountFrame,libraryDepth)

      #Example
      #                                 Mature    Precursor Width Sequence
      #CAGUGCAAUGAUAUUGUCACUAU    mmu-miR-6341 mmu-mir-6341    23 CAGUGCAAUGAUAUUGUCACUAU
      #ACAUGCUUCUUUAUAUCCUCAUA  mmu-miR-206-5p  mmu-mir-206    23 ACAUGCUUCUUUAUAUCCUCAUA
      #UGGAAUGUAAGGAAGUGUGUGG   mmu-miR-206-3p  mmu-mir-206    22 UGGAAUGUAAGGAAGUGUGUGG
      #GCUGGUCAAACGGAACCAAGUC  mmu-miR-133b-5p mmu-mir-133b    22 GCUGGUCAAACGGAACCAAGUC
      #UUUGGUCCCCUUCAACCAGCUA  mmu-miR-133b-3p mmu-mir-133b    22 UUUGGUCCCCUUCAACCAGCUA
      #UGUAAACAUCCUCGACUGGAAG   mmu-miR-30a-5p  mmu-mir-30a    22 UGUAAACAUCCUCGACUGGAAG

   }else{
      stop("Unrecognised collapse type when producing table")
   }
}

# YUP!

###################################
### Identify all files in inDir ###
###################################

inputFiles <- list.files(inDir, pattern = "bowtie\\.conv\\.GR\\.RData$", full.names = TRUE)
if (length(inputFiles)==0){stop(paste("No Bowtie GRanges (RData) files found in ",inDir,". Is this the correct directory?",sep= " "))}
withoutDir <- as.character(sub(".+\\/","",inputFiles[1]))
filePrefix <- as.character(sub("([^\\.]+\\.).+$","\\1", withoutDir)) ### WORK ON THIS

for (input in inputFiles ){
   ##########################
   ### Reading input data ###
   ##########################
   
   write(paste("\nLoading input file",input,sep=" "),stderr())

   load(input)    # GR
   
   ########################
   ### Arranging output ###
   ########################
   
   inputBasePattern <- as.character(sub("(.+)\\.bowtie.+","\\1", input))
   fileBase <- as.character(sub(".+\\/","",inputBasePattern))
   barcode <- as.character(sub(".+\\.([[:alpha:]]+)$","\\1", fileBase))

   ###################################################
   ### Find overlaps between Bowtie and annotation ###
   ###################################################
   
   write("Checking the overlap between reads and miRNA annotation",stderr())

   miRoverlapMat <- findOverlaps(miRs,GR,minoverlap=minOverlap)   # Apply minimum overlap. Found between each loci and reads.
      
   
   if (any(duplicated(subjectHits(miRoverlapMat)))){              ### Added 040413 - Identify alignments associated with multiple miRNA loci and react accordingly.
      if(treatOverlap == 'merge' || treatOverlap == 'remove'){
      
         conflictedAlignments <- subjectHits(miRoverlapMat[which(duplicated(subjectHits(miRoverlapMat)))])
         write("Alignments spanning multiple loci using specified minimum overlap:",stderr())
         write(paste( unique(values(miRs[unique(queryHits(miRoverlapMat[which(subjectHits(miRoverlapMat) %in% conflictedAlignments)]))])$mature), collapse = "\n"),stderr())
         stop("Increase the 'overlap' parameter in the Configuration file to ressolve ambiguity")
      
      }else if(treatOverlap == "ignore"){
         
         conflictedAlignments <- subjectHits(miRoverlapMat[which(duplicated(subjectHits(miRoverlapMat)))])
         uniqueConflicts <- unique(values(miRs[unique(queryHits(miRoverlapMat[which(subjectHits(miRoverlapMat) %in% conflictedAlignments)]))])$uniqueID)   
         
         if (all(uniqueConflicts %in% values(overlappingMiRs)$uniqueID)){
            ### NOTE: No guarantee that the alignments are overlapping pairs from this list that overlap each other.
            write("WARNING: Where an alignment overlaps multiple loci, these loci themselves are amongst those known to overlap another miRNA. Consider increasing the 'overlap' parameter in the Configuration file and merging or removing overlapping miRNA loci as this may ressolve conflicts.",stderr())
            write("WARNING: Beware of potential double counting.", stderr())
         }else{
            write("Alignments spanning multiple loci using specified minimum overlap. These loci are not all in the list of miRNAs known to overlap other loci:",stderr())
            write(paste( unique(values(miRs[unique(queryHits(miRoverlapMat[which(subjectHits(miRoverlapMat) %in% conflictedAlignments)]))])$mature), collapse = "\n"),stderr())
            stop("Increase the 'overlap' parameter in the Configuration file to ressolve ambiguity")
         }
            
      }else{
         stop("Unrecognised treatOverlap parameter in miR_table.R script")
      }
   }else{
      write("No alignments found to overlap multiple miRNA loci",stderr())
   }

   # YUP!


   ################################################################################
   ### Create columns in table within which to store depths of mature sequences ###
   ################################################################################
   
   mapCol <- paste(barcode,"_Read_Depth" , sep ="")
   uniqCol <- paste(barcode,"_Non_Redundant_Reads" , sep ="")

   miRCountFrame[,mapCol] <- rep.int(0, nrow(miRCountFrame))
   miRCountFrame[,uniqCol] <- rep.int(0, nrow(miRCountFrame))
   
   ######################
   ### Fill the frame ###
   ######################
   
   write("\nCalculating read depth and number of non-redundant reads mapping to each locus",stderr())
   #stopifnot(FALSE)
   
   if (proportional){

      # 030413 - Add a large section here to separate counts of multi-mappers by unique mappers. 
      # Work in an entirely miRNA centric manner to allow unique depth for all potential multimapping loci to be estimated.
      write("Dividing multi-mapping reads by a ratio reflecting an incremented unique read depth",stderr())      
      
      values(GR)$raw_count <- as.integer(sub("\\w+_x(\\d+)$","\\1",values(GR)$name,perl=TRUE)) # Distribute reads between miRNA loci alignments only - Max depth of reads to redivide between miRNA loci

      write("Identifying reads associated with a single miRNA",stderr())
      duplReadNames <- unique(values(GR[subjectHits(miRoverlapMat)])[duplicated(values(GR[subjectHits(miRoverlapMat)])$name),"name"])
      uniqReadNames <- unique(values(GR[subjectHits(miRoverlapMat)])[(! values(GR[subjectHits(miRoverlapMat)])$name %in% duplReadNames),"name"])
      
      write("Assigning depth of uniquely associated reads to miRNA",stderr())
      uniqOverlaps <- miRoverlapMat[which(values(GR[subjectHits(miRoverlapMat)])$name %in% uniqReadNames)]
      uniqueDepth <- sapply(split(as.numeric(values(GR)[,"raw_count"])[subjectHits(uniqOverlaps)], queryHits(uniqOverlaps)),sum)
      names(uniqueDepth) <- elementMetadata(miRs)[as.numeric(names(uniqueDepth)),"uniqueID"]

      uniqueNumber <- sapply(split(subjectHits(uniqOverlaps), queryHits(uniqOverlaps)), length)
      names(uniqueNumber) <- elementMetadata(miRs)[as.numeric(names(uniqueNumber)),"uniqueID"]

      write("Assigning depth of reads associated with multiple miRNAs", stderr())
      duplOverlaps <- miRoverlapMat[which(values(GR[subjectHits(miRoverlapMat)])$name %in% duplReadNames)]
      DupMiRNAsNames <- unique(values(miRs[queryHits(duplOverlaps)])$uniqueID)
      
      ### Create a matrix of all of the uniquely mapping reads that map to the loci occupied by multimapping reads.
      uniqDepthMatrix <- rep(0,length(DupMiRNAsNames))
      names(uniqDepthMatrix) <- DupMiRNAsNames
      uniqNumbMatrix <- uniqDepthMatrix # Create a matrix to contain the divided, uniquely mapped, non-redundant read numbers.
      
      DupliDistribDepth <- uniqDepthMatrix # Create a vector to record the divided multimapping read counts.
      DupliDistribNumber <- uniqDepthMatrix # Create a vector to record the divided multimapping non-redundant read counts.

      uniqDepthMatrix[names(uniqDepthMatrix)%in%names(uniqueDepth)] <- uniqueDepth[names(uniqDepthMatrix)[names(uniqDepthMatrix)%in%names(uniqueDepth)]]
      uniqNumbMatrix[names(uniqNumbMatrix)%in%names(uniqueNumber)] <- uniqueNumber[names(uniqNumbMatrix)[names(uniqNumbMatrix)%in%names(uniqueNumber)]]
      # Add 1 to the matrix in case low counts accross multiple sites leads to 100% allocation of multimappers to a single locus.
      uniqDepthMatrix <- uniqDepthMatrix +1
      uniqNumbMatrix <- uniqNumbMatrix +1

      ### Calculate the multiple loci for each read
      duplLoci <- split(queryHits(duplOverlaps), values(GR[subjectHits(duplOverlaps)])$name)

      # YUP!

      ### Cycle through the reads and divide multimappers according to the uniquely mapping depths
      for(depthSplit in names(duplLoci)){
         targetmiRuniqDepth <- uniqDepthMatrix[values(miRs[duplLoci[[depthSplit]]])$uniqueID]
         totalTargetDepth <- sum(targetmiRuniqDepth)
         depthToSplit <- sub("\\w+_x(\\d+)$","\\1",depthSplit,perl=TRUE)
         dividedDepths <- sapply(names(targetmiRuniqDepth), function(x){
            return((targetmiRuniqDepth[[x]]/totalTargetDepth)*as.integer(depthToSplit))
         },USE.NAMES = FALSE)
         names(dividedDepths) <- names(targetmiRuniqDepth) 
         DupliDistribDepth[names(dividedDepths)] <- DupliDistribDepth[names(dividedDepths)]+ dividedDepths
      
         targetmiRuniqNumb <- uniqNumbMatrix[values(miRs[duplLoci[[depthSplit]]])$uniqueID]
         totalTargetNumb <- sum(targetmiRuniqNumb)
         dividedNumb <- sapply(names(targetmiRuniqNumb), function(x){
            return(targetmiRuniqNumb[[x]]/totalTargetNumb)
         },USE.NAMES = FALSE)
         names(dividedNumb) <- names(targetmiRuniqNumb)
         DupliDistribNumber[names(dividedNumb)] <- DupliDistribNumber[names(dividedNumb)]+ dividedNumb
      }

      # YUP!

      if(!separate){
         names(miRs) <- values(miRs)$uniqueID
         if(collapse_type == "mature_id"){
            write("Merging miRNA loci according to the miRBase mature miRNA ID",stderr())
            names(uniqueDepth) <- values(miRs[names(uniqueDepth)])$mature
            names(DupliDistribDepth) <- values(miRs[names(DupliDistribDepth)])$mature
            names(uniqueNumber) <- values(miRs[names(uniqueNumber)])$mature
            names(DupliDistribNumber) <- values(miRs[names(DupliDistribNumber)])$mature
         }else if(collapse_type == "sequence"){
            write("Merging miRNA loci according to the mature sequence",stderr())
            names(uniqueDepth) <- values(miRs[names(uniqueDepth)])$sequence
            names(DupliDistribDepth) <- values(miRs[names(DupliDistribDepth)])$sequence
            names(uniqueNumber) <- values(miRs[names(uniqueNumber)])$sequence
            names(DupliDistribNumber) <- values(miRs[names(DupliDistribNumber)])$sequence
         }else{
            stop("Unrecognised collapse type")
         }
         uniqueDepth <- tapply(uniqueDepth,names(uniqueDepth),sum)
         DupliDistribDepth <- tapply(DupliDistribDepth,names(DupliDistribDepth),sum)
         uniqueNumber <- tapply(uniqueNumber,names(uniqueNumber),sum)
         DupliDistribNumber <- tapply(DupliDistribNumber,names(DupliDistribNumber),sum)
      }
      
      miRCountFrame[names(uniqueDepth),mapCol] <- uniqueDepth
      miRCountFrame[names(DupliDistribDepth),mapCol] <-  miRCountFrame[names(DupliDistribDepth),mapCol] + DupliDistribDepth
      miRCountFrame[names(uniqueNumber),uniqCol] <- uniqueNumber
      miRCountFrame[names(DupliDistribNumber),uniqCol] <-  miRCountFrame[names(DupliDistribNumber),uniqCol] + DupliDistribNumber
      miRCountFrame[,mapCol] <-  round(miRCountFrame[,mapCol],1)
      miRCountFrame[,uniqCol] <-  round(miRCountFrame[,uniqCol],1)
      miRCountFrame <- miRCountFrame[base::order(miRCountFrame[,mapCol], decreasing=TRUE),]
      
      write("Calculating Mapped_library_depth associated with miRNAs", stderr())
      miRCountFrame["Mapped_library_depth",mapCol] <- sum(miRCountFrame[,mapCol])  # Note: total row is 0 at this point
      miRCountFrame["Mapped_library_depth",uniqCol] <- sum(miRCountFrame[,uniqCol])  # Note: total row is 0 at this point


   }else{

      countDepth <- sapply(split(as.numeric(values(GR)[,"count"])[subjectHits(miRoverlapMat)], queryHits(miRoverlapMat)),sum)  # Calculate the rounded redundant read depth for all loci
      
      readAssociation  <- split(subjectHits(miRoverlapMat), queryHits(miRoverlapMat)) # Split reads between loci
      countNumber <- lapply(readAssociation,function(x){1/elementMetadata(GR)[x,"locations"]}) # Assume each read has a depth of 1 and divide this between mapped loci to find non-redudnant read depth
      countNumber <- sapply(countNumber,sum)
      
      
      if(separate){   ### 020413 - Added ability to populate table by either merging based on miRBase mature ID or not (using unique ID of mature and precursor IDs with Chromosome and Start positions)
         names(countDepth) <- elementMetadata(miRs)[as.numeric(names(countDepth)),"uniqueID"]
         names(countNumber) <- elementMetadata(miRs)[as.numeric(names(countNumber)),"uniqueID"]
      }else{
         if(collapse_type == "mature_id"){
            write("Merging miRNA loci according to the miRBase mature miRNA ID",stderr())
            names(countDepth) <- elementMetadata(miRs)[as.numeric(names(countDepth)),"mature"]
            countDepth <- tapply(countDepth,names(countDepth),sum)  # Summed the miRNA counts for each unique locus and labelled them with non-unique miRNA mature sequence names -> Finally summed
            names(countNumber) <- elementMetadata(miRs)[as.numeric(names(countNumber)),"mature"]
            countNumber <- tapply(countNumber,names(countNumber),sum) # Non-redundant read count for each locus summed depending upon mature miRNA name
         }else if(collapse_type == "sequence"){
            write("Merging miRNA loci according to the mature sequence",stderr())
            names(countDepth) <- elementMetadata(miRs)[as.numeric(names(countDepth)),"sequence"]
            countDepth <- tapply(countDepth,names(countDepth),sum)  # Summed the miRNA counts for each unique locus and labelled them with non-unique miRNA mature sequence names -> Finally summed
            names(countNumber) <- elementMetadata(miRs)[as.numeric(names(countNumber)),"sequence"]
            countNumber <- tapply(countNumber,names(countNumber),sum) # Non-redundant read count for each locus summed depending upon mature miRNA name
         }else{
            stop("Unrecognised collapse type")
         }
      }

      miRCountFrame[names(countNumber),uniqCol] <- round(countNumber,1) # 020413 - Added rounding to Count Number - rounded at the end
      miRCountFrame[names(countDepth),mapCol] <- round(countDepth,1)
 
      miRCountFrame <- miRCountFrame[base::order(miRCountFrame[,mapCol], decreasing=TRUE),]
      
      ####################
      ### Library Size ###       ### Calculate and append the library sizes (redundant and non-redundant)
      ####################
      
      write("Calculating Mapped_library_depth for reads mapping to genome",stderr())
      
      miRCountFrame["Mapped_library_depth",mapCol] <- sum(elementMetadata(GR)$count)
      miRCountFrame["Mapped_library_depth",uniqCol] <- sum(1/(elementMetadata(GR)$locations))
   }
}   

   
#####################
### Write to file ###
#####################

write("\nWriting read depth summaries to file",stderr())
write.table(miRCountFrame, file=paste(outDir,"/",filePrefix,"mature.counts.txt",sep=""), quote=FALSE, sep="\t", na = "-", row.names=FALSE )

