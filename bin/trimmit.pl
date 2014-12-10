#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long;

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

my $inputDir= "";
my $outputDir="";
my $trim5p   = 0;
my $trim3p   = 0;
my $triCut  =  -1;
my $top     = -1;
my $bot     = -1;
my @inFiles     ;
my $outputTOTAL ;
my $do_tally = 0;
my $tallyFfn = "";
my $debug = 0;

#######################################
### Organising the input parimetres ###
#######################################

GetOptions( "inDir=s" => \$inputDir,     # Specify input directory containing all input files
            "outDir=s" => \$outputDir,   # Specify output directory
            "fiveP=i"  => \$trim5p,      # Specify 5' nucleotides to be trimmed
            "threeP=i" => \$trim3p,      # Specify 3' nucleotides to be trimmed
            "tri=i"    => \$triCut,      # Specify a trinucleotide complexity cutoff
            "upper=i"  => \$top,         # Specify an upper bound for read length
            "lower=i"  => \$bot,         # Specify a lower bound for read length
            "tally=s"  => \$tallyFfn,    # Specify tally full file name
            "debug"    => \$debug        # Print extra debugging info
         );

if (length $inputDir == 0){
   die "Need a an input directory: --inDir==<DIRECTORY>\n";
}

if (length $outputDir == 0){
   die "Need a an output directory: --outDir=<DIRECTORY>\n";
}

if (length $tallyFfn ==0){
   die "Need a full file name for tally: --tally=<FFN>\n";
}

die "The input directory does not exist!" unless -e $inputDir;
die "The input directory is not a directory!" unless -d $inputDir;

die "The output directory does not exist!" unless -e $outputDir;
die "The output directory is not a directory!" unless -d $outputDir;

@inFiles = <$inputDir/*clean.uniquified.fa.gz>;
print STDERR "@inFiles\n" if $debug;

if ($triCut >= 0){
   print STDERR "Removing reads equal to or above a trinucleotide threshold of: $triCut\n";
}else{
   print STDERR "No trinucleotide cut off specified\n";
}

if ($top >= 0){
   print STDERR "Removing reads more than $top nucleotides in length\n";
}else{
   print STDERR "No maximum read length specified\n";
}

if ($bot >= 0) {
   print STDERR "Removing reads less than $bot nucleotides in length\n";
}else{
   print STDERR "No minimum read length specified\n";
}

print STDERR "Trimming $trim5p nucleotides from the 5' end of reads\n";
print STDERR "Trimming $trim3p nucleotides from the 3' end of reads\n";

die "Cannot have negative trim parameter" unless $trim5p >= 0 && $trim3p >= 0;
$do_tally = 1 if $trim5p > 0 || $trim3p > 0;

################################################
### Opening file to recorded total read data ###
################################################

$outputTOTAL = "$outputDir/total.saved.trimmit.tab";
open(OUTTOTAL, "> $outputTOTAL") || die "Not able to create $outputTOTAL : $!\n";

##########################################
### Parse the fasta files into trimmit ###
##########################################

print STDERR "\nBeginning to cycle through processed read files\n";

foreach my $inputFfn (@inFiles){


   print STDERR "\nCleansing file: $inputFfn\n";

   my $readNum     = 0;
   my $readSave    = 0;
   my $wasteTri    = 0;
   my $wasteLong   = 0;
   my $wasteShort  = 0;
   my $base           ;
   my $outputFASTA    ;
   my $outputLEN      ;
   my %unique    = ();
   my %repeated  = ();
   my %LCDiscardedUnique = ();
   my %LCDiscardedRepeat = ();
   my %total     = ();
   my $barcode        ; 

   if ($inputFfn =~ /^$inputDir\/(\S+)\.clean\.uniquified\.fa\.gz$/){
      $base = $1;                                  # Storing basic read file ID
   } else {
      die "The input file is not of the expected format!\n";
   }

   ######### HARVEST THE BARCODE SEQUENCE

   $outputFASTA = "$outputDir/$base\.clean\.processed\.fa\.gz";
   $outputLEN   = "$outputDir/$base\.report\.clean\.processed\.annotlen";

   open(IN, "gzip -cdf $inputFfn| cat|") || die "No $inputFfn : $!\n"; 

   if ($do_tally) {
      print STDERR "Will reuniquifying sequences with tally following clipping as bases will have been removed from the processed reads\n";
      open(OUTFASTA, "| $tallyFfn -i - -record-format '%R%t%X%n' -format '>tr%I_t%T_w%L_x%C%n%R%n' -o $outputFASTA") || die "Not able to create $outputFASTA : $!\n";
      # "$ID"."_t$triN"."_w$trimWidth"."_x$frequency\n$trimSeq\n";
   }
   else {
      open(OUTFASTA, "| gzip -c1 > $outputFASTA") || die "Not able to create $outputFASTA : $!\n";
   }
   open(OUTLEN, "> $outputLEN") || die "Not able to create $outputLEN : $!\n";

   while(<IN>){
      chomp;
      
      my $trimSeq   = "";
      my @fastaID       ;
      my $ID        = "";
      my $triN      = "";
      my $width     = "";
      my $frequency = "";
      my $sequence  = "";
      my $trimWidth     ;

      if(!/^>/){
         die "Faulty Fasta format - missing '>' in Uniquify output\n";
      }

      #######################################################################
      ### Splitting the read IDs into components and loading the sequence ###
      #######################################################################

   #   print"$_\n";
      @fastaID = split /_/, $_, -1 ;

      #print "@fastaID\n";

      $ID    = shift @fastaID;
      $triN  = shift @fastaID;
      $triN  =~ s/[a-zA-Z]//;
      $width = shift @fastaID;
      $width =~ s/[a-zA-Z]//;
      $frequency = shift @fastaID;
      $frequency =~ s/[a-zA-Z]//;

      $readNum += $frequency;

      chomp($sequence = <IN>);

      ######################################
      ### Apply cutoff criteria to reads ###
      ######################################

      #print "ID: $ID\tTriNuclScore: $triN\tWidth: $width\tFrequency: $frequency\nSequence: $sequence\n";
      
      ######################################
      ### Check trinucleotide complexity ###
      ######################################

      if ($triCut >= 0){
         if ($triN >= $triCut){
            #      print "Read $ID removed based on trinucleotide score\n";
            $wasteTri += $frequency;
            $LCDiscardedUnique{$width}++;
            $LCDiscardedRepeat{$width} += ($frequency - 1);
            $total{$width} += $frequency;
            next;
         }
      }

      ###################################
      ### Check if reads are too long ###
      ###################################

      if ($top >= 0){
         if ($width > $top){
            #      print "Read $ID removed based on maximum read length cutoff\n";
            $wasteLong += $frequency; 
            $total{$width} += $frequency;
            $unique{$width}++;
            $repeated{$width} += ($frequency -1) ;
            next;
         }
      }
      
      ####################################
      ### Check if reads are too short ###
      ####################################

      if ($bot >= 0){
         if ($width < $bot){
            #      print "Read $ID removed based on minimum read length cutoff\n";
            $wasteShort += $frequency;
            $total{$width} += $frequency;
            $unique{$width}++;
            $repeated{$width} += ($frequency -1) ;
            next;
         }
      }

     
      ############################################
      ### Trim nucleotodes from 5' and 3' ends ###
      ############################################

      $trimSeq = substr($sequence, $trim5p);   
      $trimSeq = substr($trimSeq, 0, -$trim3p) unless $trim3p == 0;
      

      ###############################################
      ### Write processed sequences to FASTA file ###
      ###############################################

      #print "Trimmed sequence: $trimSeq\n";
      $readSave += $frequency;
      $total{$width} += $frequency;
      $unique{$width}++;
      $repeated{$width} += ($frequency -1) ;

      $trimWidth = $width - $trim5p - $trim3p;
      if ($do_tally) {
         print OUTFASTA "$trimSeq\t$frequency\n";
      }
      else {
         print OUTFASTA "$ID"."_t$triN"."_w$trimWidth"."_x$frequency\n$trimSeq\n";
      }

   }

   #######################################################################
   ### Record length of reads along with those discarded by complexity ###
   #######################################################################

   # Note: Length and uniqueness are determined based in input sequences rather than trimmed sequences

   print OUTLEN "Read_Length\tTotal\tUnique_Low_Complexity\tUnique_High_Complexity\tRepeated_Low_Complexity\tRepeated_High_Complexity\n"; 
   foreach my $theLengths (sort { $b <=> $a } keys %total ){
       print OUTLEN "$theLengths\t$total{$theLengths}\t";

       if (defined $LCDiscardedUnique{$theLengths}){
          print OUTLEN "$LCDiscardedUnique{$theLengths}\t";
       }else{
          print OUTLEN "0\t";
       }

       if (defined $unique{$theLengths}){
          print OUTLEN "$unique{$theLengths}\t";
       }else{
          print OUTLEN "0\t";
       }

       if (defined $LCDiscardedRepeat{$theLengths}){
          print OUTLEN "$LCDiscardedRepeat{$theLengths}\t";
       }else{
          print OUTLEN "0\t";
       }

       if (defined $repeated{$theLengths}){
          print OUTLEN "$repeated{$theLengths}\n";
       }else{
          print OUTLEN "0\n";
       }
   }

   ################
   ### Clean up ###
   ################

   print STDERR "Read Number: $readNum\nReads Saved: $readSave\nReads discarded, too short: $wasteShort\nReads discarded, too long: $wasteLong\nReads discarded, complexity: $wasteTri\n";

   ###################################################
   ### WRITE TOTAL SAVED READS TO FILE PER BARCODE ###
   ###################################################
   
   die "Cannot find barcodes in file name: $base" unless $base =~ /\.(\w+)$/;
   $barcode = $1;
   
   print OUTTOTAL "$barcode\t$readSave\n";

   ################
   ### CLEAN UP ###
   ################

   close IN ;
   close OUTFASTA;
   close OUTLEN;
}

print STDERR "\nCompleted read clipping and filtering\n";

close OUTTOTAL ;


__DATA__
>1_t0.6_w22_x149810
TGAGGTAGTAGGTTGTATGGTT
>2_t0.3_w22_x103703
TGAGGTAGTAGGTTGTATAGTT
>3_t0.9_w21_x51504
CACCCGTAGAACCGACCTTGC
>4_t0.5_w22_x45737
TGAGGTAGTAGATTGTATAGTT
