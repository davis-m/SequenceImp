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



my $logFile     = "test/bowtie.log";
my $inDir      = "";
my $outDir     = "";
#my $genome      = "";
my $maxHits     = 10 ;
my $maxHitsOpt  = "";
my $chunk       = "";
my $chunkOpt    = "";
my $mismatch    = 0 ;
my $mismatchOpt = "";
my $stratum     = 0 ;      # strata always passed -strata from sequence_imp.pl
my $stratumOpt  = "";
my $genomeFile  = "";
my $inFile      = "";
#my $maxReports  = 0 ;     # Removed maxHits now defines -m option rather than -k option - 180511. -k option is determined by maxHits +1
my $maxHitsCorrection = 0;
my $maxRepsOpt  = "";
my $samFormat   = "";
my $samOpt      = "";
my $suppressOpt = "";
my $depthField  = "";
my @inDirList ;
my $debug = 0;


##############################
### COLLECTING OPTIONS #######
##############################

print STDERR "\nBeginning bowtie organisation\n";

### Insert --sam option / flag

GetOptions( "in=s"            => \$inDir, 
            "out=s"           => \$outDir, 
            "genome=s"        => \$genomeFile,        # genome Bowtie index basename (including directory)
            "log=s"           => \$logFile, 
            "chunkmbs=i"      => \$chunk,             # integer value to be passed to -chunkmbs option
            "number=i"        => \$maxHits,           # -m options of Bowtie - integers - also specifies -k 
            "strata"          => \$stratum,           # --strata option for Bowtie - ALWAYS PASSED --strata REMOVE BOWTIE CONFUSION
            "mismatches=i"    => \$mismatch,          # options to be passed to -v option of Bowtie
            "sam"             => \$samFormat,         # should Bowtie be run with SAM format ? 
            "opts=s"          => \$depthField,         # Full file name of script to add optional fields to SAM  - required if --sam specified
            "debug"           => \$debug,             # Prints extra output
         );
            
die "--number must be 1 or more. This variable sets -m filtering in Bowtie\n" if $maxHits < 1;           ####### ALTER: MaxHits will now be >= 1. Will report hits at 1 but won't know relative quantities of reads that hit twice.
die "--mismatches must be >= 0\n" if $mismatch < 0;

if(length($inDir) == 0){
   die "Need a path to a directory containing cleaned sequence files: --in=Directory/\n";
}

if (! -d $inDir){
   die "Aha! $inDir does not seem to be a directory please specify an appropriately constructed input directory\n";
}   

if(length($outDir) == 0){
   die "Need a path to a directory for bowtie output: --out=Directory/\n";
}

if (! -d $outDir){
   die "Aha! $outDir does not seem to be a directory please specify an appropriately constructed output directory\n";
}
if ($samFormat){
   die "Require the full file name of  depthField.pl if SAM format specified: --opts:<FULL FILE NAME>\n" if (! $depthField);
   die "depthField.pl does not seem to exist\n" if (! -s $depthField);
}

########################################
### Check appropriate index basename ###                                      #### Alter - now passed genome file itself rather than genome -regulated by sequence_imp.pl
########################################

print STDERR "Checking index passed: $genomeFile\n" if $debug;

die "Need a specified index basename to map the sequence reads to: --genome='bowtie index basename'" if length $genomeFile ==0;

#if ($genome eq "mouse"){
#   $genomeFile = "/nfs/research2/enright/Bowtie/Mus_musculus.NCBIM37.60.dna";
#}else{
#   die "Do not have access to the $genome genome"
#}

###################################
### Identifying all input files ###
###################################

print STDERR "Globbing input files:\n" if $debug;

@inDirList = <$inDir/*\.clean\.processed\.fa\.gz>;
print STDERR "@inDirList\n" if $debug;

########################################
### Report the parametres to be used ###
########################################

print STDERR "Bowtie parameters include:\n";
print STDERR "log file: $logFile\n";
print STDERR "input directory: $inDir\n";
print STDERR "output directory: $outDir\n";
print STDERR "genome index basename: $genomeFile\n";
$samFormat ? print STDERR "format: SAM\n" : print "format: standard Bowtie\n"; 

################################
### Organising a log file ######
################################

print STDERR "\nCreating log file\n";

if (!open LOG, "> $logFile"){
   die "Cannot open the logfile: $!\n";
}

############## Alter: maxHits will become the -m option (Will report reads with -m or fewer reportable hits), -k option will be maxHits +1. This will limit the Bowtie query and hopefully speed the process. 

$maxHitsCorrection = $maxHits + 1;
$maxHitsOpt  = "-k $maxHitsCorrection" ;

print STDERR "\n### WARNING: Discarding reads that are aligning to the genome more than $maxHits times to remove BOWTIE reporting bias\n";

#$maxReports = $maxHits-1  ;         ### Note: always include a -m option as Bowtie reports alignments sequentially and can introduce a strand bias if limited by -k
$maxRepsOpt = "-m $maxHits" ;
#$stratumOpt  = "--strata";

$chunkOpt = length $chunk != 0 ? "--chunkmbs $chunk" : "";
$mismatchOpt = "-v $mismatch" ;                            ### Sensible default set to prevent BOWTIE reporting biases.
$stratumOpt  = $stratum == 1 ? "--strata" : "";                ### Optional Strata option hardcoded by sequence_imp to --strata to avoid confusion
$samOpt = $samFormat ? "-S" : "";
$suppressOpt = $samFormat ? "" : "--suppress 5,6";

foreach my $inFile (@inDirList) {
   print STDERR "\nSending Bowtie request for $inFile\n";

   my $fileBase = "";
   my $outFile = "";

   if ($inFile =~ /^$inDir\/(\S+)\.clean\.processed\.fa\.gz$/){
      $fileBase = $1 ;
      if ($samFormat){
         $outFile = "$outDir/$fileBase.bowtie.temp.output.sam.gz";
      }else{   
         $outFile = "$outDir/$fileBase.bowtie.output.bwt.gz";                                  # Storing basic read file ID
      }
   }

   ####################################
   ### Constructing the system call ###
   ####################################
   
   ### Insert option to invoke --sam flag and retrieve sam format output

   my $samfilter = $samFormat ? " | samtools view -Sh -F 4 - " : ""; # If SAM format unmapped reads must be filtered out.

   my $syscall = "gzip -dcf $inFile | bowtie --time $mismatchOpt --best $maxHitsOpt $stratumOpt $chunkOpt $maxRepsOpt $samOpt -f $suppressOpt $genomeFile - $samfilter | gzip -c1 > $outFile";
   
   print STDERR "###########\n$syscall\n###########\n";
   
   my $sysreport = system($syscall);
   # print "Return value:$sysreport\n";

   die "syscall hasn't worked: $syscall\n" unless $sysreport==0;

   print STDERR "Bowtie is complete\n";

   #################################################
   ### SAM format file conversion to BAM format  ###
   #################################################

   if ($samFormat){
      my $samExpanded = "$outDir/$fileBase.bowtie.unique.output.sam.gz";
      my $bamFile = "$outDir/$fileBase.bowtie.unique.output.bam";
      my $bamSort = "$outDir/$fileBase.bowtie.unique.output.sort.bam";
      #my $indexIn = "$outDir/$fileBase.bowtie.unique.output.sort.bam";

      my $depthFieldCall = "$depthField $outFile $samExpanded";
      my $viewCall = "gzip -cdf $samExpanded | samtools view -bS -o $bamFile -";
      my $sortCall = "samtools sort $bamFile -T $outDir/$fileBase -o $bamSort";
      my $indexCall = "samtools index -b $bamSort";
      
      print STDERR "Adding optional fields to SAM format\n"; 
      print STDERR "$depthFieldCall\n" if $debug;
      die "Failed to add additional fields to sam $outFile\n" unless system($depthFieldCall)==0;
      
      print STDERR "Converting SAM to BAM format\n";
      print STDERR "$viewCall\n" if $debug;
      die "Failed to convert $samExpanded to bam format\n" unless system($viewCall)==0;
      
      print STDERR "Sorting BAM file\n";
      print STDERR "$sortCall\n" if $debug;
      die "Failed to sort $bamFile\n" unless system($sortCall)==0;
      
      print STDERR "Indexing BAM file\n";
      print STDERR "$indexCall\n" if $debug;
      die "Failed to index $bamSort\n" unless system($indexCall)==0;
   }

   my $oclock = localtime();
   print LOG "$oclock \nBowtie run on $inFile \nThis should produce $outFile\n The relevant system call is as follows:\n$syscall\n";
   
}

close LOG;


#bsub -M 10000 -R "rusage[mem=10000]" "gzip -dcf out.uniquified.fa.gz | bowtie --time -v 2 --best -k 101 --strata --chunkmbs 512 -
#f --suppress 5,6 --phred64-quals /nfs/research2/enright/Bowtie/Mus_musculus.NCBIM37.58.dna - | gzip -c > bowtie.output.bwt.gz"
