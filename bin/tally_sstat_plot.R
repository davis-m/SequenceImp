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

translation <- list()                                             # This will store the conversion data for interpreting 'left' and 'right' files

args <- R.utils::commandArgs(asValues=TRUE)

if (is.null(args$output)){
   stop("Require  an output file: --output=Ffn")                     # This is the PDF to which the QC data will be printed
}else{
   outfile <- args$out
}

if (is.null(args$left)){
   stop("Require a 'left' sample name: --left=string")            # This is the sample name for the 'left' sequence file as interpreted by Tally
}else{
   translation[["left"]] <- args$left
}

if (is.null(args$right)){
   stop("Require a 'right' sample name: --right=string")          # This is the sample name for the 'right' sequence file as interpreted by Tally 
}else{
   translation[['right']] <- args$right
}

if (is.null(args$input)){
   stop("Require an input file: --input=Ffn")                        # This is the input Tally paired end .sstat for interpretation 
}else{
   inFile <- args$input
}


### Read in sstat table

ug <- read.table(inFile, header = FALSE, sep = "=", row.names = 1)

discarded = list()

### Pull out the reasons for discarding a read

discarded[['left']] <- ug[grep("^discarded_left_", rownames(ug), value = TRUE),,drop=TRUE]
names(discarded[['left']]) <- sub("left_","", grep("^discarded_left_", rownames(ug), value=TRUE))
discarded[['right']] <- ug[grep("^discarded_right_", rownames(ug), value = TRUE),,drop=TRUE]
names(discarded[['right']]) <- sub("right_","", grep("^discarded_right_", rownames(ug), value=TRUE))

### Convert report into a dataframe

### Rownames:
reasons <- unique(c(names(discarded[['left']]), names(discarded[['right']])))

if (! "passed_unique" %in% rownames(ug)){stop("'passed_unique' missing from the Tally .sumstat files")}
if (! "passed_total" %in% rownames(ug)){stop("'passed_total' missing from the Tally .sumstat files")}
saved <- grep("passed",rownames(ug), value=TRUE)

if (length(saved) != 2){
   stop("There are too many 'passed' categories for the paired Tally QC plots. Check the Tally .sumstat file for 2 'passed' categories")
}

summary <- as.data.frame(matrix(data=NA, nrow = (length(reasons)+2), ncol =2, dimnames=list(rows=c(reasons,saved), cols=c(unlist(translation)))))

### Populate table with read info

for (partner in names(translation)){
#   print(paste("Current pair:", partner,"Column:", translation[[partner]]))
   for(currentSel in names(discarded[[partner]])){
#      print(paste("Current reasoning:", currentSel))
      summary[currentSel,translation[[partner]]]  <- discarded[[partner]][currentSel]
   }
   summary[saved, translation[[partner]]] <- ug[saved,]
}

### Convert "left" and "right" into a case specific nomencalture

rownames(summary) <- sub("left", paste("_",translation[["left"]],sep="") , rownames(summary))
rownames(summary) <- sub("right", paste("_",translation[["right"]],sep="") , rownames(summary))

### Work out unique vs redundant reads

summary <- apply(summary, 2, function(x){
   x[is.na(x)] <- 0
   x["passed_total"] <- x["passed_total"] - x["passed_unique"]
   return(x)
})

### Reorganise table to plot in correct orientation

rownames(summary)[which(rownames(summary) == "passed_total")] <- "passed_duplicates"
multicopy <- which(rownames(summary) == "passed_duplicates")
allalone <- which(rownames(summary) == "passed_unique")
summaryplot <- summary[-c(multicopy, allalone),]
summaryplot <- rbind( summary[c(allalone,multicopy),],summaryplot)

### Colours

colours <- rainbow(nrow(summaryplot[grep("^discarded", rownames(summaryplot)),,drop=FALSE]))
colours <- c( "black", "white", colours)
names(colours) <- rownames(summaryplot)


### Plot to PDF

pdf(outfile, width = 9, height = 4.6)
layout(matrix(c(1,2), nrow = 1,
          ncol= 2, byrow = TRUE), widths= c(3,3))

par(omi=c(0.3,1,0.3,1)) ### Change shape of graph
barplot(summaryplot,cex.axis=0.6,las=1,  cex.names=0.6, cex.main = 1, col=colours, ylab = "Reads", xlab = "Fastq File", main = "Filtering of reads and the removal\nof redundancy")
plot.new()

par(omi=c(0.3,0.3,0.3,0.3))
legend("top",legend=paste(names(colours),". Reads: ",translation[["left"]],": ",summaryplot[names(colours),translation[["left"]]],". ",translation[["right"]],": ", summaryplot[names(colours),translation[["right"]]],sep=""), fill=colours, bty="o", bg="white", horiz=FALSE,cex = 0.4, pt.cex = 0.5)

mtext(paste(translation[["left"]],"and",translation[["right"]],sep=" "),outer=TRUE,cex=1)
dev.off()
