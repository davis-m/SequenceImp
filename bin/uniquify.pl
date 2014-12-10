#!/usr/bin/env perl

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

## Collapsing duplicates from a single file 

use warnings;
use strict;
use Getopt::Long;

print STDERR "Beginning to uniquify\n";

##########################
### Programme settings ###
##########################

#my $fn = pop;
my $base = "";
my $FASTAoutext = ".uniquified.fa.gz";             # Extension of unique read output name
my $LenTABoutext = ".annotlen";                    # Extension of length output name
my $TriTABoutext = ".trinucl";
my $logext = ".log";                               # Log file extension
my $dataDir = "test";                              # Default input and output directory
my $cutoffLC = -1;                                # Default low complexity cutoff for length plot calls - Note - by default no effective cutoff applied
my $oclock = localtime();
my $maximum = -1;                                  # Default maximum length for unique read output - Default out of bounds - no cutoff applied - Not max must be more than min or ignored
my $minimum = -1;                                  # Default minimum length for unique read output - Default out of bounds - no cutoff applied
my $binSize = 10;
my $tallyffn = undef;
my $debug = 0;
my $fastq = "";
my $inFormat ="";

GetOptions( "low=i" => \$cutoffLC,                 # Not used in the current software. LC cutoff dealt with by the trimmit.pl program.
            "input_directory=s" => \$dataDir,
            "maxSize=i" => \$maximum,              # Not used in the current version of the software. Size selection dealt with by trimmit.pl
            "tallyffn=s" => \$tallyffn,
            "minSize=i" => \$minimum,
            "inFormat=s" => \$inFormat,            # If tallyFFn supplied, this will gove the accompanying input format
            "fastq"     => \$fastq,                # Input format flag. If tallyFFn not supplied this will be used to specify perl parsing.
            "debug"     => \$debug,   
         );


###################################
### Gathering all Data Files ######
###################################

my @files = <$dataDir/*>;                    
my @readFiles = grep /\.clean\.gz$/, @files;       # Cleaned read output files from reaper
#print "@files\n";
print STDERR "Input files: @readFiles\n" if $debug;

##########################################################################################
### Reading and checking the input file and preparing output files for each input file ###
##########################################################################################
foreach my $fn (@readFiles){

   my %sequences;
   my %uniqueLowCompl;
   my %repeatedLowCompl;
   my %uniqueHighCompl;
   my %repeatedHighCompl;
   my %seqLenSummary;
   my @lenComplSum;

   if ($fn =~ /^$dataDir\/(\S+)\.clean\.gz$/){
      $base = $1;                                  # Storing basic read file ID
   } else {
      die "The input file is not of the expected format!";
   }

   # Opening irrelevant files

   my $FASTAoutput = "$dataDir/$base.clean$FASTAoutext"; 
   my $lengthTABoutput = "$dataDir/$base.report.clean$LenTABoutext";
   my $trinuclTABoutput = "$dataDir/$base.report.clean$TriTABoutext";
   my $logFile = "$dataDir/$base.uniquify$logext";

   open(LOG, "> $logFile") || die "Not able to create $logFile : $!";

#####
   if (defined($tallyffn)) {
      print STDERR "Using the Tally to collapse duplicate sequences\n";
      open(IN, "$tallyffn -i $fn -o - -record-format $inFormat -format %R%t%C%t%T%n | gunzip |") || die "no tally";   # NOTE:There is a maximum read limit to sequence sorting after which the sample will be split in half and each half sorted seperately
      # tally gives us <READ>tab<COUNT>tab<TRINT>
   }
   else {
      print STDERR "Duplicate reads will be collapsed using a perl based system\n";
      open(IN, "gzip -cdf $fn| cat|") || die "No $fn : $!";
   }
   open(OUTFASTA, "| gzip -c1 > $FASTAoutput") || die "Not able to create $FASTAoutput : $!";
   open(OUTTAB, "> $lengthTABoutput") || die "Not able to create $lengthTABoutput : $!";
   open(OUTTRI, "> $trinuclTABoutput") || die "Not able to create $trinuclTABoutput : $!";

   print STDERR "Initiating uniquify for $fn. See: $logFile for details\n";
   print LOG "$oclock\n\n";
   print LOG "Input file: $fn\n";
   print LOG "FastA output file: $FASTAoutput\n";
   print LOG "Tab Length output file: $lengthTABoutput\n\n";
   
   if($maximum < 0 || $maximum < $minimum){
      if($cutoffLC < 0){
         print LOG "Settings:\nNo low complexity threshold applied\n";
      }else{
         print LOG "Settings:Low complexity cutoff: $cutoffLC\n";
      }
      print LOG "No read length cutoff applied\n";
      print LOG "Maximum or minimium read lengths conflict or are out of bounds" if ($maximum != -1 && $minimum != -1);
   }else{
      if($cutoffLC < 0){
         print LOG "Settings:\nNo low complexity threshold applied\n";
      }else{
         print LOG "Settings:Low complexity cutoff: $cutoffLC\n";
      }
      print LOG "Maximum read length: $maximum\nMinimum read length: $minimum\n";
   }

   ##################################################
   #### Read in file and sum identical sequences #### 
   ##################################################

   if (defined($tallyffn)) {
      my $which = 1;
      while (<IN>) {                               # FASTA input

         chomp;
         my ($baseSeq, $count, $trint) = split("\t", $_);

         die "There is a format error in the input file" unless defined($trint);

         my $seqLen = length($baseSeq);
         
         #################################################################
         ### Collating trinucleotide scores relative to length buckets ###
         #################################################################

         my $bin = ($trint - ($trint % $binSize))/$binSize;
         $bin = 10 if $bin >10;         
         $lenComplSum[$seqLen][$bin] += $count;

         #################################################################################
         ### Idenifying replicates and collating length distribution and annotation ######
         #################################################################################
               
         # print "$baseSeq has a length of $seqLen\n";
         $seqLenSummary{$seqLen} += $count;

         if ( $cutoffLC >= 0){            # NOTE: Default cutoff is -1. Therefore no cutoff for complexity applied until positive value given.
            if ($trint >= $cutoffLC){     # Record sequence length in a hash table - Broken up into unique and repeated reads and low and high complexity
               $uniqueLowCompl{$seqLen}++;          
            }else{
               $uniqueHighCompl{$seqLen}++;
            }
            if ($trint >= $cutoffLC){
               $repeatedLowCompl{$seqLen} += ($count -1);
            }else{
               $repeatedHighCompl{$seqLen} += ($count -1);
            }
         }else{
            $uniqueHighCompl{$seqLen}++;
            $repeatedHighCompl{$seqLen} += ($count -1);
         }

         if ($maximum < 0 || $maximum < $minimum){
            print OUTFASTA ">$which"."_t$trint"."_w$seqLen"."_x$count\n$baseSeq\n";
            $which++;
         }elsif($seqLen <= $maximum && $seqLen >= $minimum){
            print OUTFASTA ">$which"."_t$trint"."_w$seqLen"."_x$count\n$baseSeq\n";
            $which++;
         }
      }
   }
   else {
      my $total = 0;
      while (<IN>) {                               # FASTA input

         chomp;
         if (!$fastq){
            if(!/^>/){
               die "Faulty Fasta format - missing '>' in Reaper output";
            }
         }else{
            if(!/^@/){
               die "Faulty Fastq format - missing @ in Reaper output";
            }
         }

         my @lineSplit = split "\t", $_, -1;       # split without stripping trailing white space - important if pop is to be used.
         my $read_tri = pop @lineSplit;            # Extract tri nucleotide scores from read identifier
         my %annotation = split "=", $read_tri;
         my $uniqID = shift @lineSplit;            # Store ID
         
         chomp(my $baseSeq = <IN>);                # Second line is the base sequence
         die "There is a sequence missing within the input file" unless ( defined $baseSeq);
         
         my $seqLen = length($baseSeq);
         
         #################################################################
         ### Collating trinucleotide scores relative to length buckets ###
         #################################################################

         my $bin = ($annotation{tri} - ($annotation{tri} % $binSize))/$binSize;
         $bin = 10 if $bin >10;         
         $lenComplSum[$seqLen][$bin]++;

         #################################################################################
         ### Idenifying replicates and collating length distribution and annotation ######
         #################################################################################
               
         # print "$baseSeq has a length of $seqLen\n";
         $seqLenSummary{$seqLen}++;

         my $IDreadinfo = "$baseSeq"."__$annotation{tri}";      
         
         if(!exists($sequences{$IDreadinfo})){
            if ( $cutoffLC >= 0){
               if ($annotation{tri} >= $cutoffLC){     # Record sequence length in a hash table - Broken up into unique and repeated reads and low and high complexity
                  $uniqueLowCompl{$seqLen}++;          # NOTE: by default cutoffLC is 101: therefor cutoff isn't applied as maximum triscore is 100- effective cutoff applied by trimmit.pl
               }else{
                  $uniqueHighCompl{$seqLen}++;
               }
            }else{
               $uniqueHighCompl{$seqLen}++;
            }
         }else{
            if ( $cutoffLC >= 0){
               if ($annotation{tri} >= $cutoffLC){
                  $repeatedLowCompl{$seqLen}++;
               }else{
                  $repeatedHighCompl{$seqLen}++;
               }
            }else{
               $repeatedHighCompl{$seqLen}++;
            }
         }

         $sequences{$IDreadinfo}++;                   # Record sequence frequency in a hash table - Sequence is the key
         
         if ($fastq){
            chomp(my $QualID = <IN>);
            die "Faulty Fastq format - missing '+' in Reaper output\n" if ($QualID !~ /^\+/);
            chomp(my $QualScore = <IN>);
         }
      }

      #######################
      #### Sanity  check ####
      #######################
      foreach my $numb (values %sequences){
         $total += $numb;                          # Record of total reads
      }
      print "The total reads is $total\n" if $debug;
   }


   ###################################################################################
   ### Producing a table summarising read length data and trinucleotide complexity ###
   ###################################################################################
   
   print STDERR "Recording trinucleotide complexity in relation to read length\n";

   print OUTTRI "0-9\t10-19\t20-29\t30-39\t40-49\t50-59\t60-69\t70-79\t80-89\t90-99\t100+\n";
   
   for my $row (0..$#lenComplSum){
      print OUTTRI "$row\t";
      if(!defined $lenComplSum[$row]){
         print OUTTRI "0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";
         next;
      }
      #print "\n\n$row\n";
      #print "@{$lenComplSum[$row]}\n";
      for my $column (0..$#{$lenComplSum[$row]}){
         if (!defined $lenComplSum[$row][$column]){
            $lenComplSum[$row][$column] = 0;
         }
      }
      #print "@{$lenComplSum[$row]}\n";
      my $data = join "\t", @{$lenComplSum[$row]};
      print OUTTRI "$data";
      my $tabs = 11-scalar(@{$lenComplSum[$row]});
      if ($tabs > 0){
         print OUTTRI "\t";
         if ($tabs > 1){
            for my $padding (1..($tabs -1)){
               print OUTTRI "0\t";
            }
         }
         print OUTTRI "0";
      }
      print OUTTRI "\n";
   } 

   ####################################################################
   #### Creating a tab delineated file of Read lengths for plotting ###
   ####################################################################

   print STDERR "Recording read length distribution data\n";

   print OUTTAB "Read_Length\tTotal\tUnique_Low_Complexity\tUnique_High_Complexity\tRepeated_Low_Complexity\tRepeated_High_Complexity\n";

   foreach my $theLengths (sort { $b <=> $a } keys %seqLenSummary ){
       print OUTTAB "$theLengths\t$seqLenSummary{$theLengths}\t";

       if (defined $uniqueLowCompl{$theLengths}){
          print OUTTAB "$uniqueLowCompl{$theLengths}\t";
       }else{
          print OUTTAB "0\t";
       }

       if (defined $uniqueHighCompl{$theLengths}){
          print OUTTAB "$uniqueHighCompl{$theLengths}\t";
       }else{
          print OUTTAB "0\t";
       }

       if (defined $repeatedLowCompl{$theLengths}){
          print OUTTAB "$repeatedLowCompl{$theLengths}\t";
       }else{
          print OUTTAB "0\t";
       }
    
       if (defined $repeatedHighCompl{$theLengths}){
          print OUTTAB "$repeatedHighCompl{$theLengths}\n";
       }else{
          print OUTTAB "0\n";
       }
   }

###############################################
#### Creating a FASTA File of Unique Reads ####
###############################################

   if (!defined($tallyffn)) {
      my $which = 1;
      foreach my $element (sort { $sequences{$b} <=> $sequences{$a} } keys %sequences){
         # print "Width cutoffs: Max-$maximum Min-$minimum\n";
         
         my @seqStuff = split /__/ , $element ;
         my $relevantSeq = shift @seqStuff;
         my $relevantTri = shift @seqStuff;

         my $width = length($relevantSeq);
         #     print "Element width: $width\n";                                              # Alter - insertion of an if statement to only screen for size if the maximum and minimum sizes are in bounds
         if ($maximum < 0 || $maximum < $minimum){
            print OUTFASTA ">$which"."_t$relevantTri"."_w$width"."_x$sequences{$element}\n$relevantSeq\n";
            $which++;
         }elsif($width <= $maximum && $width >= $minimum){
            print OUTFASTA ">$which"."_t$relevantTri"."_w$width"."_x$sequences{$element}\n$relevantSeq\n";
            $which++;
         }else{
            next;
         }
      }
   }

####################
### Closing Down ###
####################

   close(OUTTRI);
   close(OUTFASTA);
   close(OUTTAB);
   close(IN);
   close(LOG);
}

print STDERR "Collapse duplicate reads to unique sequences with associated read depths for trimmed reads in $dataDir\n";

__DATA__
GTATGCCGTCTT   tri=100
GCCCTAAGGTGAATTTTTTGGG  tri=85
GAATTGATCAGGACATAGA  tri=100
GGGGGTATAGCTCAGTGGCAGAGCATTTGACTCGT tri=81
GGTAAAATGGCTGAGTAAGCATTAGACTGTATCGT tri=84
GTAAACTGCTTAGATTTT   tri=93
GGGGGTGTAGCTCAGTGGTAGAGAGCGTGCTTG   tri=64
GGAGGGGGAAAAATCGTA   tri=68
GGGCGCTGTAGGCTGT  tri=71
GTAGGCGGGGCCGAGTTGGTAAGGCTATTGTTTTA tri=75
GGGCGCTGTAGGCTGT  tri=71
