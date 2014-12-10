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

library(RColorBrewer)

### This script contains the functions for Reaper QC plotting

######################################
#### Plotting base composition #######
######################################

baseCompPlot <- function(inBaseFile, title, debuggleBase){
   if (debuggleBase){write(paste("Plotting base composition for",inBaseFile),stderr())}
   baseTable <- read.delim(inBaseFile, quote = "", row.names = 1,
                           colClasses = c("character", "numeric", "numeric", "numeric" ,"numeric", "numeric"))

   if (nrow(baseTable) == 0){
      plot(NA, main=title, cex = 0.5, xlim=c(0,1), ylim=c(0,1), axes=F, xlab="", ylab="")
      text(0.5,0.5,"No reads found to match these criteria")
      return()
   }

   #print("Here")
   baseTable <- baseTable[,c("A","T","G","C","N")]
   rowSum <- apply(baseTable,1,sum)
   baseTableFreqs <- as.data.frame(matrix(data = NA, nrow = nrow(baseTable), ncol = ncol(baseTable),
                            dimnames = list(row = rownames(baseTable), col = colnames(baseTable))))
   #print("No here")
   for(i in rownames(baseTable)){
      #print(paste("Rowname 1:",i,"Rowname 2:", rownames(baseTable[i,])))
      baseTableFreqs[i,] <- baseTable[i,]/rowSum[i]  # Convert base frequencies per cycle to base proportions
   }
   baseTableFreqs <- t(baseTableFreqs)
   baseTableFreqs <- baseTableFreqs[c("A","T","G","C","N"),]                                                # Assume these are the only bases we need to deal with - CHECK
   cols <- c("palegreen3","palevioletred3","palegoldenrod","paleturquoise3","grey")
   names(cols) <- rownames(baseTableFreqs)
   #print(paste("Ahahaha here:", baseTableFreqs[,1]))
   barplot(as.matrix(baseTableFreqs), ylim=c(0,1.1), col=cols, beside=FALSE, border="black", space=0,
           xlab="Base frequency by cycle", ylab="Frequency", las=3, cex.axis=0.7, cex.names=0.6,
           main=title, cex.main = 1, xaxt='n')
   #print("Hello")
   abline(h=seq(0,1,0.2), lty=3, lwd=1, col="darkgrey")
   legend("top",legend=rownames(baseTableFreqs), fill=cols, bty="o", bg="white", horiz=TRUE,cex = 0.5, pt.cex = 0.5)
   axis(side=1,at=seq.int(0,ncol(baseTableFreqs),5),cex.axis= 0.7)
}



##################################################################
#### Quality score plotting - Plot 10%, 50% and 90% quantiles ####
##################################################################

qualityPlot <- function(inQualFile, titleQ, debuggleQual){
      if (debuggleQual){write(paste("Plotting sequencing quality for",inQualFile),stderr())}
      quants2plot <- c("q10","q50","q90")                                                                   # Selected quantiles plotted by function
      qualTable <- read.delim(inQualFile, quote = "", row.names = 1,
                             colClasses = c("character","numeric","numeric","numeric",
                             "numeric","numeric","numeric","numeric"))


      if (nrow(qualTable) == 0){
         plot(NA, main=titleQ, xlim=c(0,1), cex = 0.5, ylim=c(0,1), axes=F, xlab="", ylab="")
         text(0.5,0.5,"No reads found to match these criteria")
         return()
      }
                      #qualTable <- platform(qualTable)  ### Probably needs removal? - Calls platform function above to translate quality scores - DEPRICATED due to the potential for error
      qualTransl <- c(
                        q0 = "0% Quantile",
                        q10 = "10% Quantile",
                        q25 = "25% Quantile",
                        q50 = "50% Quantile",
                        q75 = "75% Quantile",
                        q90 = "90% Quantile",
                        q100 = "100% Quantile"
      )

      qualCols <- brewer.pal(ncol(qualTable[,quants2plot]),"RdPu")                                          # Assign a palette
      plot(NULL, xlim=c(1,nrow(qualTable)), ylim=c(min(qualTable)*0.9,max(qualTable)*1.1),
           xlab="Quality by cycle", ylab="Quality",
           main = titleQ, cex.main = 1, cex.axis=0.7)
      matlines(as.matrix(qualTable[,quants2plot]), type="l", col = qualCols, lty=1,lwd =2,bty = "n")
      legend("bottomright", legend=qualTransl[colnames(qualTable[,quants2plot])],lty = 1, col=qualCols,
           bty="n", lwd = 2, cex = 0.5, pt.cex = 0.5)

}

#################################################################
####### Plotting read complexity thresholds against length ######
#################################################################

complexityPlot <- function(complexityFile, titleC, debuggleComp){

   if (debuggleComp){write(paste("Plotting read complexity for",complexityFile),stderr())}

   complexity <- read.table(file=complexityFile, sep="\t", colClasses= "integer", header=TRUE)

   if (nrow(complexity) == 0){
      plot(NA, main=titleC, xlim=c(0,1), cex = 0.5, ylim=c(0,1), axes=F, xlab="", ylab="")
      text(0.5,0.5,"No reads found to match these criteria")
      return()
   }

   colnames(complexity) <- c("0-9","10-19","20-29","30-39","40-49","50-59","60-69","70-79","80-89","90-99","100+")

   complexityMatrix <- as.data.frame(matrix(data = NA, nrow = nrow(complexity), ncol = ncol(complexity),
                                     dimnames = list(row = rownames(complexity), col = colnames(complexity))))

   rowSums <- apply(complexity,1,sum)
   for(j in rownames(complexity)){
      complexity[j,] <- complexity[j,]/rowSums[j]
      elements2sum <- c()
      for (k in rev(seq_len(length(complexity[j,])))){
         #print (paste("Element:",argh[j,k]))
         elements2sum <- c(elements2sum, k)
         #print (elements2sum)
         complexityMatrix[j,k] <- sum(complexity[j, elements2sum])
      }
   }

   cols2plot <- c("10-19","20-29","30-39","40-49","50-59","60-69","70-79","80-89","90-99")
   compPlotLabels <- c("10", "20", "30", "40", "50", "60", "70","80", "90")                  # Each column is up to and including the sequences within the specified bine. Therefore a threshold of 20 will include all reads up to and including those with a score of 19...

   triCols <- brewer.pal(ncol(complexityMatrix[,cols2plot]),"RdYlGn")

   plot(NULL, xlim=c(0,nrow(complexityMatrix)),ylim=c(0,1),
       xlab="Sequence length (nt)", ylab="Proportion removed by threshold",
       cex.main = 1, cex.axis=0.7,xaxt="n", main = titleC)

   matlines(as.matrix(complexityMatrix[,cols2plot]), type="l",col=triCols, lty=1,lwd =2,bty = "n")

   legend("topleft", legend=compPlotLabels,lty = 1, col=triCols,
                  bty="n", lwd = 2, cex = 0.5, pt.cex = 0.5)

   axis(side=1,at=seq.int(1,nrow(complexityMatrix),5), label = rownames(complexityMatrix[seq.int(1,nrow(complexityMatrix), 5),]),cex.axis= 0.7)
}


####################################################################################
#### Potting read length distributions - More complex plot from uniquify output ####
####################################################################################


lengthPlot <- function(lenFile, titleL, upperL, lowerL, zeroPass, debuggleLen){

   if (debuggleLen){write(paste("Plotting length distribution plot for",lenFile),stderr())}

   lenTable <- as.matrix( read.table(lenFile, quote = "", row.names = 1, header = TRUE))                   # Pass an .annotlen table
   if (nrow(lenTable) == 0){
      plot(NA, main=titleL, xlim=c(0,1), cex = 0.5, ylim=c(0,1), axes=F, xlab="", ylab="")
      text(0.5,0.5,"No reads found to match these criteria")
      return(0)
   }

   if (upperL >= 0 && upperL > lowerL ){
   ### Checking parametres - Maximum and minimum cutoffs must be within data size distribution
      if (lowerL < 0){
         lowerL <- 0
         write(paste("Minimum size cutoff set to",lowerL), stderr())
      }
      if (upperL > max(as.numeric(rownames(lenTable)))){
         upperL <- max(as.numeric(rownames(lenTable)))
         write(paste("Maximum size cutoff too large, so corrected to",upperL), stderr())
      }
   }else if(upperL < 0 && lowerL >=0 ){
         upperL <- max(as.numeric(rownames(lenTable)))
         write(paste("No maximum boundary specified"), stderr())
   }else{
      write ("Plotting the entire length distribution plot", stderr())
      lowerL <- 0
      upperL <- max(as.numeric(rownames(lenTable)))
   }

   lenTableMills <- lenTable/1000000
   selectedSizes <- lowerL:upperL+1                                                                         # Adjustment- Row "0" is the first row of the matrix etc (Base 1).

#   print(paste("Min:",lowerL))

   complMatrix <- matrix(data = 0, nrow=max(as.numeric(rownames(lenTableMills)))+1, ncol=ncol(lenTableMills), dimnames=list("row" = as.character(0:max(as.numeric(rownames(lenTableMills)))), "columns"=colnames(lenTableMills)))
   selectedMatrix <- complMatrix
   nonSelectedMatrix <- complMatrix

#   print(class(rownames(lenTableMills)))

   complMatrix[as.character(rownames(lenTableMills)),] <- lenTableMills                                     # Fills in any gaps left in the length table
   selectedMatrix[selectedSizes,] <- complMatrix[selectedSizes,]                                            # Defining the desired sizes and splitting atrix for plotting
   nonSelectedMatrix[-selectedSizes,] <- complMatrix[-selectedSizes,]

#   print(paste("Base numbers",rownames(complMatrix)))
#   print(paste("The maximum number of bases is:",max(as.numeric(rownames(lenTable)))))
#   print(paste("Plot parametres:",lowerL,"and",upperL,sep = " "))
#   print(selectedSizes)
#   print(complMatrix[0,])

   readSets <- c("Unique_Low_Complexity","Unique_High_Complexity","Repeated_Low_Complexity","Repeated_High_Complexity")

   selectedMatrix <- t(selectedMatrix[,readSets])
   nonSelectedMatrix <- t(nonSelectedMatrix[,readSets])

   if (zeroPass) {                                                                                          # Defining ymax on the '0' length reads depending upon 'zero' parameter - NOTE global parameter used in subroutine
      ymaxLen = max(complMatrix[-1,"Total"])
   }else{
      ymaxLen = max(complMatrix[,"Total"])
   }

   lenCols <- c("red","darkblue","orange","lightblue")
   names(lenCols) <- readSets

   br <- barplot(selectedMatrix[readSets,], ylim = c(0,1.3*ymaxLen),las=3, col = lenCols[readSets],
                  xaxt = 'n', main = titleL,
                  ylab="Number of reads (millions)", border=NA,cex.axis=0.7,  cex.names=0.6, cex.main = 1)
   axis(cex.axis= 0.7,side=1,at=br[seq.int(1,ncol(selectedMatrix),5)], labels = colnames(selectedMatrix)[seq.int(1,ncol(selectedMatrix),5)], line=1)
   title(xlab = "Read length following Reaper", mgp = c(3.5, 1, 0))

   # Defining the legend depending upon whether there were any low complexity reads found or not -> At step=reaper they are not expected
   legendSums <- apply(complMatrix[,readSets],2,sum)!=0
   if (sum(legendSums[c("Unique_Low_Complexity","Repeated_Low_Complexity")]) == 0){
      legendLabels <- c("Unique_High_Complexity" = "Unique reads", "Repeated_High_Complexity" = "Repeated reads")
   }else{
      legendLabels <- c("Unique_Low_Complexity"="Unique, low-complexity reads","Unique_High_Complexity"="Unique, high-complexity reads","Repeated_Low_Complexity"="Repeated, low-complexity reads","Repeated_High_Complexity"="Repeated, high-complexity reads")
   }

   legend("top", legend = legendLabels, fill = lenCols[names(legendLabels)], bty="o", bg="white",cex = 0.5, pt.cex = 0.5, ncol = 2 )
   barplot(nonSelectedMatrix,xaxt = 'n', xlab="", ylab="", border=NA, col=lenCols, density=40, axes=FALSE, add=TRUE, axisnames=FALSE)

   return(sum(lenTable[,"Total"]))
}


################################################
#### Summarising the .sumstat read use file ####
################################################

sumstat_plot <- function(readUse_pass){
#   stopifnot(FALSE) 
   readUseMills <- readUse_pass/1000000

   destinyPlot <- readUseMills[c(grep("^discarded", rownames(readUseMills), value = TRUE), "total_error", "total_accepted"),,drop=FALSE]        # Subset .sumstat to desired 

   colReadUse <- rainbow(nrow(destinyPlot[grep("^discarded", rownames(readUseMills)),,drop=FALSE]))
   colReadUse <- c(colReadUse, "white", "black")
   names(colReadUse) <- rownames(destinyPlot)

   labelStuff <- rep.int(NA, length(rownames(destinyPlot)))
   names(labelStuff) <- rownames(destinyPlot)
   labelStuff["total_accepted"] <- paste("Total Accepted Reads (Reads: ",round(readUseMills["total_accepted",1],2)," mill, ",round(readUse_pass["total_accepted",1]*100/readUse_pass["total_input",1],digits=1),"%)", sep = "")

   pie(destinyPlot[,1],col = colReadUse[rownames(destinyPlot)], cex.main = 1, labels = labelStuff[rownames(destinyPlot)],
          main = paste("Pie chart displaying those reads with\nrecognised barcodes and the basis upon which\nothers were discarded\nTotal Reads:",readUse_pass['total_input',],  sep = " "))

   destinyPlotLeg <- destinyPlot[destinyPlot[,1] >0,,drop=FALSE]
   colReadUseLeg <- colReadUse[rownames(destinyPlotLeg)]

   legend("bottomright",legend=paste(names(colReadUseLeg),". Million Reads: ",format(destinyPlot[names(colReadUseLeg),1]),sep=""), fill=colReadUseLeg, bty="o", bg="white", horiz=FALSE,cex = 0.6, pt.cex = 0.5)


}


#############################################################################
#### Plotting read length distributions - Simple plot from Tuck output ######
#############################################################################

simplelengthPlot <- function(simpleLenFile, simpleTitleL,  simpleZeroPass, simpleDebuggleLen){
   #print (paste("Read length file: ", simpleLenFile))   

   simpleLenTable <- read.delim(simpleLenFile, quote = "", row.names = 1, colClasses = c("character","numeric"))
   simpleLenTable$millions <- simpleLenTable$count/1000000
   
   if (simpleZeroPass) {
      #print ("Ignoring reads of length 0")                                                                                          # Defining ymax on the '0' length reads depending upon 'zero' parameter - NOTE global parameter used in
      simpleYmaxLen = max(simpleLenTable[-1,"millions"])
   }else{
      simpleYmaxLen = max(simpleLenTable[,"millions"])
   }
   
   simpleLenTable <- t(simpleLenTable)
  

   #print (paste("The y axis maximum is: ",simpleYmaxLen))

   pointies <- barplot(simpleLenTable["millions",], ylim = c(0,1.2*simpleYmaxLen),las=3, 
            main = simpleTitleL, xaxt='n',
            ylab="Number of reads (millions)",col = "orange", border=NA,cex.axis=0.7,  cex.names=0.7, cex.main = 1)
   axis(cex.axis= 0.7,side=1,at=pointies[seq.int(1,ncol(simpleLenTable),5)], labels = colnames(simpleLenTable)[seq.int(1,ncol(simpleLenTable),5)], line=1)
   title(xlab = "Read length following Reaper", mgp = c(3.5, 1, 0))

   #print(paste("Total accepted reads: ", sum(simpleLenTable["count",])))
   return(sum(simpleLenTable["count",])) ### Return value to allow comparison of total reads with each barcode
}
