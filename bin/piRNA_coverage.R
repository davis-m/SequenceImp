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


library('GenomicRanges')
library('IRanges')

coverage <- list()
overlap <- list()
debugme <- FALSE

args <- R.utils::commandArgs(asValue=TRUE)

if (is.null(args$inputDir)){stop('Need an input directory containing alignment GRanges data')}else{inDir <- args$inputDir}
if (is.null(args$outputDir)){stop('Need an output directory into which data will be deposited')}else{outDir <- args$outputDir}
if (is.null(args$mapDepth)){stop('Need a table containing the mapped read depth for samples against the genome to calculate read density')}else{mappedTable <- args$mapDepth}
if (is.null(args$repeatName)){stop('Require a repeat to be specified')}else{repName <- args$repeatName} 
if (!is.null(args$debug)){debugme <- TRUE}

if (debugme){write("\n### DEBUGGING ###\n\n",stderr())}

maps <- read.table(mappedTable, sep = "\t",row.names=1, as.is = TRUE, header = TRUE, col.names = c("barcode","total"), colClasses = c("character", "integer"))
nameScan <- paste(repName,".bowtie.conv.GR.RData", sep = "")
inputGRs <- dir(inDir, pattern= nameScan, full.names=TRUE)
lengthDepth <- list()

covOut <- paste(outDir,"/",repName,".coverage.RData", sep = "")
pingOut <- paste(outDir,"/",repName,".pingpong.tab", sep = "")
lenOut <- paste(outDir,"/",repName,".length_distrib.tab", sep = "")

for (thisFile in inputGRs){
   barcode <- as.character(sub(".+\\.([[:alpha:]]+)\\.[[:alnum:]_]+\\.bowtie.+","\\1",thisFile))
   load(thisFile)

   write(paste("Determining coverage for", barcode),stderr())

   ########################################
   ### Calculate the coverage of repeat ###
   ########################################

   #stopifnot(FALSE)
   GRSense <- GR[strand(GR)=='+']
   GRAntiSense <- GR[strand(GR)=='-']
   SenseCoverage <- coverage(GRSense, weight=round(as.numeric(values(GRSense)$count/(maps[barcode,'total']/1000000)),15))
   AntiSenseCoverage <- coverage(GRAntiSense, weight=round(as.numeric(values(GRAntiSense)$count/(maps[barcode,'total']/1000000)),15))

   coverage[[barcode]] <- list()
   coverage[[barcode]][['Sense']] <- SenseCoverage
   coverage[[barcode]][['Antisense']] <- AntiSenseCoverage
   #stopifnot(FALSE)

   ######################################################
   ### Length distribution of reads mapping to repeat ###
   ######################################################
   
   write("Calculating repeat associated read length distribution",stderr())

   readLength <- as.data.frame(list('length'=width(GR),'depth'= elementMetadata(GR)$count))
   splitLength <- split(readLength, readLength$length)
   lengthDepth[[barcode]] <- sapply(splitLength,function(x){sum(x$depth)})

   ##################################
   ### Calculate pingpong overlap ###
   ##################################
   
   write("Calculating sense-antisense overlap for reads",stderr())

   #stopifnot(FALSE)
   if(length(GRAntiSense)>0){
      # Convert antisense to the same strand as sense to find overlaps.
      strand(GRAntiSense) <- '+'
      # Reduce the antisense reads to a single base (Last base equivalent to the 5' end).
      GRAntiSenseShort <- flank(GRAntiSense,start=F,-1)
      # Find overlaps: GRSense = query, GRAntiSenseShort = subject.
      mat <- findOverlaps(GRSense,GRAntiSenseShort)
      # Split AntiSense sequences according to the sense transcripts that overlap with their last base.
      tL=split(GRAntiSense[subjectHits(mat)],queryHits(mat))
      # Select only the sense sequences with at least one overlap.
      tA=GRSense[unique(queryHits(mat))]

      # May not have any overlaps: Will depend on parameters
      if(length(tA) > 0){
         # Overlap the selected unique sense transcripts with the list of antisense transcripts they overlap with. 
         intersectGR=pintersect(tA,tL)
         # Split the depths of the antisense reads in the identicle manner to the overlaps with the sense reads
         tC=split(as.numeric(values(GRAntiSense[subjectHits(mat)])$count),queryHits(mat)) # Must be a list of the same length as tA
         # Multiply the antisense depths by the corresponding sense depths.
         tC=mapply('*',tC,as.numeric(values(tA)$count),SIMPLIFY=FALSE) # Overlap of a single size is affecting object class of TC and leading to crash?
                                                        # Single size, 2 overlaps = matrix, single size, 1 overlap = numberic 

         #Calculate the length of the overlaps for each read.
         widths=unlist(width(intersectGR))
         #Pair the overlaps with the count depths. 
         control = data.frame(widths=widths,counts=unlist(tC))
         #Split the dataframe into widths
         t=split(control,control$widths)
         #For each width sum up the counts
         counts=sapply(t,function(x){sum(x$counts)})
         #Recalculate as a proportion of the total overlap counts
         counts=counts/sum(counts)
         overlap[[barcode]] <- counts
      }
   }
}

#################################################
### Fill in the gaps in the overlap dataframe ###
#################################################

write("\nOrganising the writing of data to file",stderr())

# All samples (multiplexed) must have an entry even if there is no overlap
expectedSamples <- sapply(inputGRs,function(x){as.character(sub(".+\\.([[:alpha:]]+)\\.[[:alnum:]_]+\\.bowtie.+","\\1",x))})

if(length(names(overlap))>0){
   maxOverlap <- max(unlist(lapply(overlap, function(x){max(as.numeric(names(x)))})))
   overlapFrame <- as.data.frame(matrix(data=0, nrow=length(expectedSamples), ncol=maxOverlap, dimnames=list('row' = expectedSamples,"col"=seq_len(maxOverlap))))
   
   for(selSample in names(overlap)){
      overlapFrame[selSample,names(overlap[[selSample]])] <- overlap[[selSample]]
   }
}else{
   # A special case if there are no overlaps for any sample
   overlapFrame <- as.data.frame(matrix(data=NA, nrow=length(inputGRs), ncol=1, dimnames=list('row' = expectedSamples,"col"="overlaps")))
}

##################################################
### Summarise length distribution to dataframe ###
##################################################

largest <- max(unlist(lapply(lengthDepth, function(x){as.numeric(names(x))})))
smallest <- min(unlist(lapply(lengthDepth, function(x){as.numeric(names(x))})))
lengthRange <- seq(smallest,largest)
lengthTable <- as.data.frame(matrix(data=0, ncol = length(lengthRange), nrow = length(names(lengthDepth)), dimnames = list("rows"=names(lengthDepth), "cols"=as.character(lengthRange))))    ### Fill in the gaps with 0s.
for(barCall in names(lengthDepth)){
   lengthTable[barCall,names(lengthDepth[[barCall]])] <- lengthDepth[[barCall]]
}

save(coverage, file=covOut)
write.table(overlapFrame, file=pingOut ,quote=FALSE, sep="\t")
write.table(lengthTable, file=lenOut, quote=FALSE, sep="\t")

