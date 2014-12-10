#!/usr/bin/env perl

##########################################################################################################
### A short lightweight programme for summarising trinucleotide complexity scores for Paired end reads ###
##########################################################################################################

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



my $dataDir = "";
my $binSize = 10;
my $fastq = "";

#print STDERR "Summarising trinucleotide complexity for Paired end files\n";

if
(!GetOptions(

         "input_directory=s" => \$dataDir,
         "fastq"             => \$fastq

)){print STDERR "Failed to parse commandline options\n";exit(1);};

die "Need to define --input_directory containing the trimmed files to be summarised\n" if (length($dataDir) == 0);

my @files = <$dataDir/*>;
my @seq_files = grep /\.clean\.gz$/,@files;

#print "@seq_files\n";

foreach my $thisfile (@seq_files){
   
   my $base;
   my @complexity;
   #my %score_tally;

   if ($thisfile =~ /^($dataDir\/\S+)\.clean\.gz$/){
      $base = $1;                                  # Storing basic read file ID
      #     print "$base\n";
   } else {
      die "The input file is not of the expected format!";
   }

   open (IN, "gzip -cd $thisfile | cat |") || die "Unable to open the input directory\n";

   while (<IN>){
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

      #print "$_\t";
     
      my @namesplit = split("\t", $_); 
      die "FASTA name format incompatible with trinucl_complexity_light.pl\n" if scalar(@namesplit)!=3 ;
      my  $trint = $namesplit[2];

      die "Please check tri score correctly postitioned in read names.\n" if ($trint !~ /tri=/);

      my ($tag,$score) = split("=",$trint);
      #$score_tally{$score}++;

      chomp(my $sequence = <IN>);
      my $seq_length = length($sequence);
     
     
      #print "$score\t$sequence\t$seq_length\n";
      
      my $bin = ($score - ($score % $binSize))/$binSize;
      $bin = 10 if $bin >10;
      $complexity[$seq_length][$bin]++;

      if ($fastq){
         chomp(my $QualID = <IN>);
         die "Faulty Fastq format - missing '+' in Reaper output\n" if ($QualID !~ /^\+/);
         chomp(my $QualScore = <IN>);
      }
   }
   close(IN);

   #while ( my ($k,$v) = each %score_tally ) {
   #   print "$k => $v\n";
   #}

   my $trinucl_file = $base.".report.clean.trinucl";         #### Where is extension defined? Consistent with uniquify script. Maintenance.
   open(OUTTRI, "> $trinucl_file") || die "Not able to create $trinucl_file : $!";

   print STDERR "Recording trinucleotide complexity in relation to read length\n";
   print OUTTRI "0-9\t10-19\t20-29\t30-39\t40-49\t50-59\t60-69\t70-79\t80-89\t90-99\t100+\n";

   for my $row (0..$#complexity){
      print OUTTRI "$row\t";
      if(!defined $complexity[$row]){
         print OUTTRI "0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n";
         next;
      }
      #print "\n\n$row\n";
      #print "@{$complexity[$row]}\n";
      for my $column (0..$#{$complexity[$row]}){
         if (!defined $complexity[$row][$column]){
            $complexity[$row][$column] = 0;
         }
      }
      #print "@{$complexity[$row]}\n";
      my $data = join "\t", @{$complexity[$row]};
      print OUTTRI "$data";
      my $tabs = 11-scalar(@{$complexity[$row]});
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
   close(OUTTRI);
}
