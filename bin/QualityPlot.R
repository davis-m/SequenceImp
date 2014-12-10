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


### Script intended to take the reports from tuck and transalate them into pretty plots.
library(RColorBrewer)

### Parameters and input files

minimumSize <- -1 # DEFAULT - Minimum size for length plot cutoff           
maximumSize <- -1 # DEFAULT                                                
zero <- FALSE     # Consider 0 count when assessing length plot ymax cutoff                              # Add a flag for intial QC vs trimmit QC plots
debugme <- FALSE  # Print debugging information
amipaired <- FALSE # Reformat plotting to suit paired file system (tally hasn't been run at this point...)

### Generic parametres

args <- R.utils::commandArgs(asValues=TRUE)

if (is.null(args$dataDir)) {
   stop("I require a path to the data: --dataDir=directory/")
} else { 
   inPath <- args$dataDir                                                                                # Input - Reaper output
}

if (is.null(args$basicID)){
   stop("Require a basic lane identifier to distinguish pdfs --basicID=string")                          # basic ID -  a unique prefix for pdf file nomenclature
}else{
   basicID <- args$basicID
}

if (is.null(args$outDir)) {
   stop("I require a path to a suitable destination for my qc output: --outDir=directory/")
} else {
   outPath <- args$outDir                                                                                # Output - QC pdf
}

if (!is.null(args$noZero)){
   zero <- TRUE                                                                                 # Plot zero parameter passed from seq_imp 
}

if (!is.null(args$debug)){
   debugme <- TRUE
}

if (!is.null(args$paired)){
   amipaired <- TRUE
}

if (!is.null(args$maxSize)){
   maximumSize <- as.numeric(args$maxSize)
}
if (!is.null(args$minSize)){
   minimumSize <- as.numeric(args$minSize)
}

if (is.null(args$plotfunc)){
   stop("I require the script containing the plotting functions: --plotfunc=script")
}else{
   plotSubroutines <- args$plotfunc
}

### Stage specific arguments

if (is.null(args$stage)){
   stop("I require a --stage=No. flag to specify the round of qc plotting to be conducted")
}else{
   stage <-as.numeric(args$stage)
}

if (stage == 1){
   if (is.null(args$multi)){
      stop("I require a --multi=(0/1) to specify whether the samples are multiplexed with barcodes ")
   }else{
      multiplexed <- as.numeric(args$multi)
   }
}

if (stage == 2){
   if (is.null(args$saved)){
      stop("I require --saved=<FILE> to determine the number of reads that passed processing ")
   }else{
      allsaved <- args$saved
   }
}

if (debugme){write("\n### DEBUGGING ###\n\n",stderr())}

#####################################################################

#### Dodgy script for estimating platform from minimumQC score and correcting quality scores to Phred
# DEPRICATED 060111 - Not a fool proof way to determine the scoring method used.

#platform <- function(qcTemplate){
#   minex <- min(qcTemplate) # Guesses platform based upon the minimum qc score in matrix
#   if(minex < 59){
#      qcTemplate <- qcTemplate - 33
#   }else{
#      if(minex < 64){
#         qcTemplate <- round(10 * log10((10^((qcTemplate-64)/10))+1)) 
#      }else{
#         qcTemplate <- qcTemplate - 64
#      }
#   }
#   return(qcTemplate)
#}

#####################################################################

##########################################
### FUNCTIONS FOR PLOTTING ###############
##########################################

source(plotSubroutines)

### Plots:

if (stage == 1){

   ##################################################################
   #### Identifying all of the files for qc plotting - STAGE ONE ####
   ##################################################################
   
   repPattern <- "report"                                                                                   # Distinct pattern for reaper outputs
   sstatPattern <- "sumstat"                                                                                  # Reaper summary file distinct pattern
   
   reportFiles <- list.files(inPath, pattern= repPattern)
   sepVector <- sub("^[[:alnum:]\\_]+\\.([[:alpha:]]+)\\..+","\\1",reportFiles)                                # Isolating barcodes
   groupFiles <- split(reportFiles, as.factor(sepVector))
   
   sstatFile <- list.files(inPath, pattern = sstatPattern)
   
   #############################
   ### MAKING THE PLOTS ########
   #############################
   
   ##########################
   ### Stage one QC plots ###
   ##########################
   
   totalReads <- list()
   outFile <- paste(outPath,"/",basicID,"_Reaper_qc.pdf", sep = "")
   
   pdf(outFile, width = 15, height = (length(names(groupFiles))+1)*3)  
   
   plotNumb <- length(names(groupFiles))* 5
   
   if(debugme){write(names(groupFiles),stderr())}
   
   layout(matrix(c(seq_len(plotNumb),rep.int(plotNumb+1, 2), rep.int(plotNumb+2,3)), nrow = length(names(groupFiles))+1, 
          ncol= 5, byrow = TRUE), widths= c(2,2,2,2,2))
   
   par(omi=c(0.3,0.3,0.3,0.3))
  

   write("Beginning plots",stderr())
   #stopifnot(FALSE)   
   for(barcode in names(groupFiles)){
   
      write(paste("Beginning plots for",barcode),stderr())

   ####################################
   ### BASE COMPOSITION INPUT PLOTS ###
   ####################################
   
      # Input composition plot
   
      basePlotFdirty <- grep("\\.input\\.nt", groupFiles[[barcode]], value = TRUE)                             # Select file
      basePlotFdirtyPath <- paste(inPath, basePlotFdirty, sep ="/")                                             # Crate path
      if(debugme){write(paste("Path to the base plot directory:",basePlotFdirtyPath,sep= " "),stderr())}
      baseCompTitleDirty <- paste("Nucleotide frequency per base for\nthe",barcode,"sample input", sep = " ")  # Create title
  
      #print(paste("Plotting base composition before cleaning:", basePlotFdirtyPath))

   #   stopifnot(FALSE)
      dummy <- baseCompPlot(basePlotFdirtyPath,baseCompTitleDirty,debugme)                                                      # Plot
      
   #   print("Plot 1 complete")
   
   
   ##############################
   ### QUALITY INPUT PLOT #######
   ##############################
   
   
      qualFile <- grep("\\.input\\.q",groupFiles[[barcode]], value = TRUE)
      qualFilePath <- paste(inPath, qualFile, sep = "/")
      qualTitle <- paste("Quality per cycle for the",barcode,"\nsample - displayed as raw ASCII values", sep =" ")
      #stopifnot(FALSE)
      dummy <- qualityPlot(qualFilePath, qualTitle,debugme)
  
   #   print("Plot 2 complete")
   
      
   ###############################################
   ### NUCLEOTIDE COMPOSITION FOLLOWING REAPER ###
   ###############################################
      
      if(barcode == "total"){                                                                                # 'total' does not have a .clean file
         basePlotFclean <- grep("\\.miss\\.nt", groupFiles[[barcode]], value = TRUE)
         basePlotFcleanPath <- paste(inPath, basePlotFclean, sep ="/")
         baseCompTitleClean <- paste("Nucleotide frequency per base for the reads with\nno correct barcodes following Reaper")
         dummy <- baseCompPlot(basePlotFcleanPath,baseCompTitleClean,debugme)
      }else{
         basePlotFclean <- grep("\\.clean\\.nt", groupFiles[[barcode]], value = TRUE)
         basePlotFcleanPath <- paste(inPath, basePlotFclean, sep ="/")
         baseCompTitleClean <- paste("Nucleotide frequency per base for\nthe",barcode,"sample following Reaper", 
                            sep = " ")
         dummy <- baseCompPlot(basePlotFcleanPath,baseCompTitleClean,debugme)
      }
   
   #   print("Plot 3 complete")
   
   ##############################
   ### COMPLEXITY CUTOFF PLOT ###
   ##############################
      #stopifnot(FALSE)   
      if(barcode == "total"){
          plot.new()                                                                                           # 'total' does not have a complexity plot
      }else{
         complFile <- grep("\\.clean\\.trinucl",groupFiles[[barcode]], value = TRUE)
         complFfn <- paste(inPath, complFile, sep = "/")
         complTitle <- paste("Cleaned read complexity - The proportion\nof reads removed by a specified\ntrinucleotide complexity threshold")
         dummy <- complexityPlot(complFfn,complTitle,debugme)
      }
   ##############################                                                                            
   ### LENGTH PLOT ##############
   ##############################
   
      if(barcode == "total"){
         plot.new()                                                                                            # 'total' does not have a length distribution plot
      }else if(amipaired){
         simpleLenPlotF <- grep("\\.clean\\.len", groupFiles[[barcode]], value = TRUE)
         simpleIggle <- paste(inPath, simpleLenPlotF, sep ="/")
         simpleLenTitle <- paste("Read length plot for",barcode,"sample\nfollowing Reaper", sep = " ")
         totalReads[[barcode]] <- simplelengthPlot(simpleIggle, simpleLenTitle, zero ,debugme) ## Adding up reads for piechart of read counts per barcode in lane
      }else{
         lenPlotF <- grep("\\.clean\\.annotlen", groupFiles[[barcode]], value = TRUE)
         iggle <- paste(inPath, lenPlotF, sep ="/")
         lenTitle <- paste("Read length plot for",barcode,"sample\nfollowing Reaper", sep = " ")
         totalReads[[barcode]] <- lengthPlot(iggle, lenTitle, maximumSize, minimumSize, zero ,debugme) ## Adding up reads for piechart of read counts per barcode in lane
      }
   }
   
   
   ##############################
   #### LANE SUMMARY PLOTS ######
   ##############################
   
   
   write("Beginning general plots",stderr())

   sstatPath <-  paste(inPath, sstatFile, sep="/")
   readUse <- read.table(sstatPath , header = FALSE, sep = "=", row.names = 1) # Loading the file
   
   
   # Sanity checks! - DO NOT USE 'doubles' IN IF STATEMENTS - INTEGERS - BIT LIMIT
   
   withBarcodes <- sum(as.integer(unlist(totalReads)))                                             # Total read count from all of the barcoded reads - from .annotlen file - 
   if (as.integer(readUse["total_accepted",]) != withBarcodes){                                    # Output of Reaper recorded in .sumstat file (total reads accepted) - Thes
      stop ("The total reads in the uniquify .annot file is not the same as the Reaper output")
   }
   
   # All sets of discarded reads (with reasoning) should add up to total
   
   if (sum(as.integer(readUse[grep("^discarded", rownames(readUse)),1])) != as.integer(readUse["total_discarded",1])){                 # Must be a better way to write this w
      stop ("The total discarded reads in the Reaper .sumstat file is not the same as the total of the individual discarded classes")
   }
   
   # The breakdown of reads should add up to the total
   
   if (as.integer(readUse["total_input",]) != sum(as.integer(readUse[c(grep("^discarded", rownames(readUse), value = TRUE),"total_error", "total_accepted" ),]))){
      stop ("The total_input in the .sumstat file is not the sum of the other output classes")
   }
   

   sumstat_plot(readUse)

   
   
   ###############################################################################################
   #### Summarisng the split of total reads in the uniquified files to each available barcode ####
   ###############################################################################################



   if(multiplexed != 0){
      totalReadsPlot <- as.matrix(lapply(totalReads,function(x){x/1000000}))
      #pin = c(1,3 )
      old <- par(mar=c(5.1,20,4.1,15)) # Querky R rule - will store old settings before setting new.
      #print(old)
      colBarSplit <- terrain.colors(nrow(totalReadsPlot))
      #print(colBarSplit)
      barplot(cbind(totalReadsPlot,NA), col = colBarSplit, ylab = "Number of reads (millions)",
         main= "Number of non-discarded reads which\nbelong to each sample", ylim=c(0,withBarcodes/1000000*1.2), cex.main = 1)
   
      legend("bottomright",horiz = FALSE, bty="o",cex = 0.5, pt.cex = 0.5, bg="white",fill=colBarSplit,legend = paste(rownames(totalReadsPlot),"- Million reads:", totalReadsPlot[,1]))
      par(old)                                                                                                      
   #pie(unlist(totalReads), main = "Pie chart of the number of reads for whom\nthe barcode was recognised")         
   }else{                                                                                                           
      write("The lane does not contain multiplexed samples so the barcode split plot is skipped",stderr())        
   }
   
   mtext(basicID, outer=TRUE, cex=1)
   
   write("Plots complete", stderr())
   
   devname2 <- dev.off()
   
   #################################
   ### END of Stage one QC plots ###
   #################################


} else {

   if (stage == 2){
   
      write("Plotting length summary",stderr())
      
      ##########################################################
      ### Plotting QC plots for processed and filtered reads ###
      ##########################################################

      totalReads <- list()

      repPattern <- "report"                                                                                   # Distinct pattern for reaper outputs
     
      savedTotals <- read.table(allsaved, sep = "\t",row.names=1, as.is = TRUE, header = FALSE,col.names = c("barcode","total"), colClasses = c("character", "integer"))

      reportFiles <- list.files(inPath, pattern= repPattern)
      sepVector <- sub("^[[:alnum:]\\_]+\\.([[:alpha:]]+)\\..+","\\1",reportFiles)                                # Isolating barcodes
      groupFiles <- split(reportFiles, as.factor(sepVector))
      
      #############################
      ### MAKING THE PLOTS ########
      #############################
      #stopifnot(FALSE)      
      outFile <- paste(outPath,"/",basicID,"_Processed_reads_qc.pdf", sep = "")
      
      pdf(outFile, width = 5, height = length(names(groupFiles))*5)  
      
      plotNumb <- length(names(groupFiles))
      
      if(debugme){write(names(groupFiles),stderr())}
      
      layout(matrix(seq_len(plotNumb), nrow = length(names(groupFiles)), 
             ncol= 1, byrow = TRUE), widths= lcm(10.5))
      
      par(omi=c(0.3,0.3,0.3,0.3))

      for(barcode in names(groupFiles)){
         
         write(paste("Sample:",barcode),stderr())

         ################################
         ### LENGTH DISTRIBUTION PLOT ###
         ################################
   
         lenPlotF <- grep("\\.clean\\.processed\\.annotlen", groupFiles[[barcode]], value = TRUE)
         iggle <- paste(inPath, lenPlotF, sep ="/")
         lenTitle <- paste("Read length plot for",barcode,"\nfollowing filtering\nReads saved:",savedTotals[barcode,"total"], sep = " ")
         totalReads[[barcode]] <- lengthPlot(iggle, lenTitle, maximumSize, minimumSize, zero ,debugme) ## Adding up reads for piechart of read counts per barcode in lane
      }
      
      write("Summary plots complete",stderr())
      
      mtext(basicID, outer=TRUE, cex=1)

      devname <- dev.off()

   }
}



