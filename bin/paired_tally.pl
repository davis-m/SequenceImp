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

### NOTE: Five prime bases will be removed before uniquifying

use warnings;
use strict;
use Getopt::Long;

my $tallyFfn = "";
my $inONE = "";
my $inTWO = "";
my $outONE = "";
my $outTWO = "";
my $min = -1;
my $max = -1;
my $tri = -1;
my $five = -1;
my $minOPT = "";
my $maxOPT = "";
my $triOPT = "";
my $fiveOPT = "";
my $stat = "";
my $inFormat = "";
my $outFormat = "";
my $maxQual = "";
my $fastq = 0;

if(!GetOptions(
   "tally=s"      => \$tallyFfn,
   "inONE=s"      => \$inONE,
   "inTWO=s"      => \$inTWO,
   "outONE=s"     => \$outONE,
   "outTWO=s"     => \$outTWO,
   "lower=s"      => \$min,
   "upper=s"      => \$max,
   "tri=s"        => \$tri,
   "fiveP=s"      => \$five,
   "stat=s"       => \$stat,
   "inFormat=s"   => \$inFormat,
   "outFormat=s"  => \$outFormat,
   "fastq"        => \$fastq
)){
   print STDERR "Failed to parse options from the command line.\n";
   exit(1);
}

$minOPT = "-l $min" if($min   >= 0);
$maxOPT = "-u $max" if($max   >= 0);
$triOPT = "-tri $tri" if($tri   >= 0);
$fiveOPT = "-si $five" if($five >= 0);  ### Note this will remove bases before uniquifying
$maxQual = $fastq ? "--with-quality": "";

my $penguin = "$tallyFfn -i $inONE -j $inTWO -o $outONE -p $outTWO -sumstat $stat -record-format \'$inFormat\' -format \'$outFormat\' $maxQual $minOPT $maxOPT $triOPT $fiveOPT";
print STDERR"\n-----------------------\nThe following command is being passed to Tally:\n\n$penguin\n-----------------------\n\n";

die "Paired tally has failed ($penguin): $!\n" if system($penguin);
