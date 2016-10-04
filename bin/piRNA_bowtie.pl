#!/usr/bin/env perl

use warnings;
use strict;
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


### TO DO:

my $debug = 0;
my $input = "";
my $repeat = "";
my $output = "";
my $maxHits = 0;
my $chunk = -1;
my $mismatches = 0;
#my $genome = "";
my $record = "";

my $depthFfn = "";

my @inDirList;
my $repeatName;
my $maxHitsCorrection;
my $maxHitsOpt;
my $maxRepsOpt;
my $chunkOpt;
my $mismatchOpt;

#######################
### Collect Options ###
#######################

GetOptions(
      "in=s"         => \$input,
      "out=s"        => \$output,
      "debug"        => \$debug,                # Print extra output to aid debugging
      "maxHits=i"    => \$maxHits,
      "mismatch=i"   => \$mismatches,
      "chunkmbs=i"   => \$chunk,
      "depthField=s" => \$depthFfn,
      "repeat=s"     => \$repeat,               # Repeat index basename for current repeat feature. Will pass one at a time.
#      "genomeIndex=s"=> \$genome,              # Genome mapping now calculated in piRNA_analysis.pl
      "record=s"     => \$record                # The file for recording alignment number for each Bowtie call
   );

print STDERR "\n### DEBUGGING ###\n\n" if $debug;

die "--maxHits must be 1 or more. This variable sets -m filtering in Bowtie\n" if $maxHits < 1;
die "--mismatches must be >= 0\n" if $mismatches < 0;
die "Need a repeat index basename for Bowtie to repeats: --repeat=basename\n" if(length($repeat) == 0);
die "Need a file into which the alignment number for each processed data file can be stored\n" if (length($record)==0);
#die "Need a genome index basename for Bowtie and calculation of normalisation factor: --genomeIndex\n" if(length($genome) == 0);

die "--depthField must be specified\n" if ! $depthFfn;
die "depthField does not seem to exist\n" if ! -s $depthFfn;

die "Need a path to a directory containing cleaned sequence files: --in=Directory/\n" if(length($input) == 0);
die "Aha! $input does not seem to be a directory please specify an appropriately constructed input directory\n" if (! -d $input);

die "Need a path to a directory for piRNA bowtie output: --out=Directory/\n" if(length($output) == 0);
die "Aha! $output does not seem to be a directory please specify an appropriately constructed output directory\n" if (! -d $output);

###################################
### Identifying all input files ###
###################################

print STDERR "Identifying processed read files\n";
print STDERR "Globbing input files:\n" if $debug;

@inDirList = <$input/*\.clean\.processed\.fa\.gz>;
print STDERR "@inDirList\n" if $debug;

#####################################
### Organising output directories ###
#####################################
#
#my $mappedOut = "$output/genome_mapped.txt";
#open(MAPPED, "> $mappedOut") || die "Could not open the file for recording genome mapping data\n";
#print MAPPED "Sample\tReads\n";

open (ALIGNED, "> $record") || die "Could not open the file for recording alignments: $record\n";
print ALIGNED "Sample\tAlignments\n";
close (ALIGNED) || die "Could not close the alignment record: $record\n";

################################
### Arranging Bowtie Options ###
################################
if ($repeat =~ /^\S*\/\w+\.([^\/]+)$/){
   $repeatName = $1;                # Subselect repeat name
}

print STDERR "\n### WARNING: Discarding reads that are aligning to the genome more than $maxHits times to remove BOWTIE reporting bias\n";

$maxHitsCorrection = $maxHits + 1;
$maxHitsOpt  = "-k $maxHitsCorrection" ;
$maxRepsOpt  = "-m $maxHits";

$chunkOpt    =  $chunk <= 0 ? "" : "--chunkmbs $chunk" ;
$mismatchOpt = "-v $mismatches" ;

### Set up file name (repeatName.alignments...)

my %barcode_record;

foreach my $inFile (@inDirList) {     # Arrange sequential input files (clean.processed.fa.gz) 
   print STDERR "\n### Collecting alignment info for: $inFile ###\n\n";
   
   my $sampleName = "";
   my $barcode = "";
   my $fileBase = "";
   my $outFile = "";
#   my $rpkmOut = "";
#   my $total = 0;

   #############################
   ### Organise output names ###
   #############################

   if ($inFile =~ /^$input\/(\S+)\.clean\.processed\.fa\.gz$/){
      $sampleName = $1;
#      $sampleName =~ /\.(\w+)$/ && ($barcode = $1);
      $fileBase = "$sampleName.$repeatName" ;
      $outFile = "$output/$fileBase.bowtie.temp.output.sam.gz" ;
#      $rpkmOut = "$output/$fileBase.rpkm.bowtie.temp.ebwt.gz";
   
      # Also work out the barcode for the sample for the alignment record
      if($sampleName=~/^\w+\.(\w+)$/){
         $barcode = $1;
         die ("Barcode recognised twice in processed read file by repeat Bowtie submission script") if exists($barcode_record{$barcode});
         $barcode_record{$barcode}++;
      }else{
         die "Can't work out the barcode for: $sampleName\n";
      }
   }else{
      die "Could not calculate the file base name from the file name: $inFile\n";
   }

   my $samExpanded = "$output/$fileBase.bowtie.unique.output.sam.gz";
   my $bamFile = "$output/$fileBase.bowtie.unique.output.bam";
   my $bamSort = "$output/$fileBase.bowtie.unique.output.sort.bam";
   #my $indexIn = "$output/$fileBase.bowtie.unique.output.sort.bam";
   
   ################################################################
   ### Arrange all systems calls for sequential file processing ###
   ################################################################

   my $syscall = "gzip -dcf $inFile | bowtie --time $mismatchOpt --best $maxHitsOpt --strata $chunkOpt $maxRepsOpt -S -f $repeat - | samtools view -Sh -F 4 - | gzip -c1 > $outFile";
   my $depthFieldCall = "$depthFfn $outFile $samExpanded";
   my $viewCall = "gzip -cdf $samExpanded | samtools view -bS -o $bamFile -";
   my $sortCall = "samtools sort $bamFile -T $output/$fileBase -o $bamSort";
   my $indexCall = "samtools index $bamSort";

   print STDERR "Systems calls:\n\n" if $debug;

   print STDERR "Bowtie call:\n$syscall\n\n" if $debug;
   print STDERR "Read depth call:\n$depthFieldCall\n\n" if $debug;
   print STDERR "Samtools 'views' call:\n$viewCall\n\n" if $debug;
   print STDERR "Samtools 'sort' call:\n$sortCall\n\n" if $debug;
   print STDERR "Samtools 'index' call:\n$indexCall\n\n" if $debug;

   print STDERR "\nSending Bowtie request\n";
   print STDERR "###########\n$syscall\n###########\n";
   die "Bowtie hasn't worked: $syscall\n" unless system($syscall)==0;
   print STDERR "Bowtie completed\n";
  
   open (BOWRES, "gzip -cd $outFile |") || die "Could not open the Bowtie SAM file: $outFile\n";
   
   my $any_alignments = 0;
   while (<BOWRES>){
      chomp;
      if (/^[^@]/){
         $any_alignments++;
      }else{
         next;
      }
   }
   close (BOWRES) || die "Could not close the connection to the SAM file: $outFile\n";

   open (RECORD, ">> $record") || die "Could not open the file to record alignments for $barcode: $record\n";
   print RECORD "$barcode\t$any_alignments\n";
   close (RECORD) || die "Could not close the record file after recording data for $barcode: $record\n";

   # If there are no alignments in the SAM file skip the remaining steps to avoid a crash
   if($any_alignments > 0){

      print STDERR "\nAdding optional fields to SAM format\n";
      die "Depthfield conversion hasn't worked\n" unless system($depthFieldCall)==0;   
      print STDERR "Depth fields added to files\n";
      
      print STDERR "\nConverting SAM to BAM format\n";
      die "Samtools 'view' (Bam conversion) hasn't worked\n" unless system($viewCall)==0;
      print STDERR "Files converted to BAM format\n";

      print STDERR "\nSorting BAM file\n";
      die "Samtools 'sort' hasn't worked\n" unless system($sortCall)==0;
      print STDERR "BAM files sorted\n";
      
      print STDERR "\nIndexing BAM file\n";
      die "Samtools 'index' hasn't worked\n" unless system($indexCall)==0;
      print STDERR "BAM files indexed\n";
   
   }else{
      print STDERR "No alignments found for $barcode. BAM conversion skipped.\n"
   }
#   my $rpkmCall = "gzip -dcf $inFile | bowtie --time $mismatchOpt $chunkOpt --suppress 2,3,4,5,6,7,8 -f $genome - | gzip -c1 > $rpkmOut";
#   die "Samtools 'index' hasn't worked\n" unless system($rpkmCall)==0; 
#   open(GENBOW,"gzip -cd $rpkmOut | cat |") || die "Could not open the genome aligned read file\n";
#   while(<GENBOW>){
      #   print "$_";
#      chomp;
#      /_x(\d+)/ && ($total += $1);
#   }
#   close GENBOW;
#   print MAPPED "$barcode\t$total\n";
   
   print STDERR "\n### Finished alignment data collection for $inFile ###\n"
}

#close MAPPED;

######## Write to file a file containing all of the mapped read number information for each sample



