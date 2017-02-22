#!/usr/local/bin/perl

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

# This script will take an index to add (1 or 2), a format (fasta of fastq) and read file.
# It will subsequently parse the read file and append the index to each read ID (/1 or /2).

use strict;
use warnings;
use Getopt::Long;

my $readFile = "";
my $format = "";
my $index = 0;
my $outFile = "";

if(!GetOptions(
   "reads=s" => \$readFile,
   "format=s" => \$format,
   "index=i" => \$index,
   "out=s" => \$outFile
)){
   die "Could not parse the command line options\n";
}

die "Need a file of read sequences (--reads)\n" if length($readFile)==0;
die "Need a format for the read file (--format)\n" if length($format)==0;
die "Need an index integer to append to read IDs (--index)\n" if ($index==0);
die "Need an output file for formatted read sequences (--out)\n" if length($outFile)==0;

die "The read file does not seem to exist: $readFile\n" if (! -s $readFile);
die "Format must be either fasta OR fastq\n" if (($format ne 'fasta') && ($format ne 'fastq'));

my $lineCount = 0;

open(IN, "gzip -cdf $readFile |") || die "Could not open the read file: $readFile\n";
open(OUT, "| gzip -c1 > $outFile") || die "Could not opent the output file: $outFile\n";

while(<IN>){
   chomp;
   my $fastLine = $_;
   $lineCount++;

   if($format eq "fastq"){
      if(($lineCount%4)==1){
         die "Expecting an @ sign: line $lineCount\n" if (!($fastLine=~/^@/));
         $fastLine="$fastLine/$index";
      }else{
         die "@ sign on the wrong line: line $lineCount\n" if ((($lineCount%4)!=0) && ($fastLine=~/^@/)); # NOTE: Quality line can start with an @
         if(($lineCount%4)==2){
            die "Sequence not as expected (ATGCN): line $lineCount\n" if (!($fastLine=~/^[AGCTN]*$/));
         } elsif(($lineCount%4)==3){
            die "Quality ID not as expected: line $lineCount\n" if (!($fastLine=~/^\+$/));
         } elsif(($lineCount%4)==0){
            die "Quality line contains spaces: line $lineCount\n" if ($fastLine=~/\s/);
         }
      }
   }elsif($format eq "fasta"){
      if(($lineCount%2)==1){
         die "Expecting a > sign: line $lineCount\n" if (!($fastLine=~/^>/));
         $fastLine="$fastLine/$index";
      }else{
         die "> sign on the wrong line: line $lineCount\n" if ($fastLine=~/^>/);
         die "Sequence not as expected (ATGCN): line $lineCount\n" if (!($fastLine=~/^[AGCTN]*$/));
      }
   }else{
      die "Format is unacceptable\n";
   }

   print OUT "$fastLine\n";
}

close(IN) || die "Could not close the read file: $readFile\n";
close(OUT) || die "Could not close the output file: $outFile\n";



