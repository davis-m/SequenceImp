#!/usr/bin/env bash

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

set -e

references=false
# Required version of samtools
samver="1.3"

function clean_up {
   if [[ $? == 0 ]]; then
      echo -e "Systems check completed successfully\n" 1>&2
   else
      echo -e "System checks failed\n" 1>&2
   fi
}


while getopts :s:rh opt
do
    case "$opt" in
    s)
      step=$OPTARG
      ;;
    r)
      references=true
      ;;
    h)
      cat <<EOU
-s pipeline_step
EOU
      exit 0
      ;;
    :) echo "Flag $opt needs argument" 1>&2
        exit 1;;
    ?) echo "Flag $opt unknown" 1>&2
        exit 1;;
   esac
done

if [[  $step == "" ]]; then
   echo "Need step! (see -s)" 1>&2
   false
fi

if [[ $step == "complete" || $step == "reaper" || $step == "filter" || $step == "align" || $step == "features" ]]; then
   if [[ $step == "complete" || $step == "reaper" ]]; then
      if ! reaper -h > /dev/null 2>&1; then
         echo -e "\nCould not find reaper. To install reaper and for associated documentation please visit http://www.ebi.ac.uk/~stijn/reaper/. You will need to ensure the programme is installed in the file path" 1>&2
         false
      else
         echo -e "\nreaper found.\nVersion: " $(reaper --version) 1>&2
      fi
   fi
   
   if [[ $step == "complete" || $step == "reaper" || $step == "filter" ]]; then
      if ! tally -h > /dev/null 2>&1; then
         echo -e "\nCould not find tally. To install tally and for associated documentation please visit http://www.ebi.ac.uk/~stijn/reaper/." 1>&2
         false
      else
         echo -e "\ntally found.\nVersion: " $(tally --version) 1>&2
      fi
   fi
   
   if [[ $step == "complete" || $step == "align" || $step == "features" ]]; then
      if ! bowtie --version > /dev/null 2>&1; then
         echo -e "\nCould not find bowtie. Please download and install bowtie (http://bowtie-bio.sourceforge.net/index.shtml)." 1>&2
         false
      else
         echo -e "\nbowtie found.\nVersion: " $(bowtie --version) 1>&2
         if $references; then
            echo -e "Reference:\nUltrafast and memory-efficient alignment of short DNA sequences to the human genome.\nLangmead B., Trapnell C., Pop M. and Salzberg S.\nGenome Biology (2009), 10:R25." 1>&2
         fi
      fi
   fi
   
   if [[ $step == "complete" || $step == "align" || $step == "features" ]]; then
      if ! which samtools > /dev/null 2>&1; then
         echo -e "\nCould not find samtools. Please download and install samtools (v$samver+) (http://htslib.org//)." 1>&2
         false
      else
         # Check samtools version installed here.
 
         if  ! samtools --version > /dev/null 2>&1; then
            echo -e "\nCould not identify samtools version (samtools --version). Please use samtools v$samver+" 1>&2
            false
         fi

         currentver=$(samtools --version | head -n1 | cut -d " " -f2)

         if [[ ! $currentver =~ ^[[:digit:].]+$ ]]; then
            echo -e "\nSamtools version parsing has failed"
            false
         fi

         # Sorts in ascending order using version sorting
      
         if [[ "$(printf "$samver\n$currentver" | sort -V | head -n1)" == $samver ]]; then
            echo -e "\nsamtools found.\nVersion: $currentver"
         else
            echo -e "\nNeed a more recent version of samtools (v$samver+).\nCurrent version: $currentver"
            false
         fi
 

         if $references; then
            echo -e "Reference:\nThe Sequence alignment/map (SAM) format and SAMtools.\nLi H.*, Handsaker B.*, Wysoker A., Fennell T., Ruan J., Homer N., Marth G., Abecasis G., Durbin R. and 1000 Genome Project Data Processing Subgroup.\nBioinformatics  (2009), 25:2078-9." 1>&2
         fi
      fi
   fi
   
   if [[ $step == "complete" || $step == "reaper" || $step == "filter" || $step == "align" || $step == "features" ]]; then
      if ! R --version > /dev/null 2>&1; then 
         echo -e "\nCould not find R. Please download and install R (http://www.r-project.org/)." 1>&2
         false
      else
         echo -e "\nR found.\nVersion " $(R --version) 1>&2
      fi
   fi
   
   if [[ $step == "complete" || $step == "reaper" || $step == "filter" || $step == "align" || $step == "features" ]]; then
      if ! R --vanilla > /dev/null 2>&1 <<EOR
      library(GenomicRanges)
      library(gplots)
      library(IRanges)
      library(RColorBrewer)
      library(ShortRead)
      library(R.utils)
EOR
      then
         echo -e "\nCould not find all of the required R libraries. Require GenomicRanges, gplots, IRanges, RColorBrewer, ShortRead, R.utils and their dependencies." 1>&2
         false
      else
         echo -e "\nAll the required R libraries are present." 1>&2
         if $references; then
            echo -e "References:" 1>&2
            echo -e "R.utils: Various programming utilities. Bengtsson H." 1>&2
            echo -e "RColorBrewer: ColorBrewer palettes. Neuwirth E." 1>&2
            echo -e "GenomicRanges: Representation and manipulation of genomic intervals. Aboyoun P., Pages H. and Lawrence M." 1>&2
            echo -e "gplots: Various R programming tools for plotting data. Warnes G.R." 1>&2
            echo -e "IRanges: Infrastructure for manipulating intervals on sequences. Pages H., Aboyoun P. and Lawrence M." 1>&2
            echo -e "ShortRead: a Bioconductor package for input, quality assessment and exploration of high-throughput sequence data.\nMorgan M., Anders S., Lawrence M., Aboyoun P., Pages H. and Gentleman R. Bioinformatics (2009), 25:2607-2608.\n" 1>&2
         fi
      fi
   fi
else
   echo "Unrecognised step specified to system_check.sh" 1>&2
   exit 1
fi

