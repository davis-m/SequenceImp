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


##########################
## Pipeline parametres ###
##########################

my $progname = "sequence_imp.pl";
my $help = 0;
my $debug = 0;

my $stage = "";
my $logStore = "sequence_imp.log";
my $rawdata = "";
my $metadata = "";
my $base = "reap";                                                                                             # Passed to reaper as a base name for the reaper output files
my $baseDir = "";
my $reaperRoute = "reaper";                                    # Default route for Reaper script
my $reaperExtra = "";
my $uniquifyRoute = "$ENV{SEQIMP_ROOT}/bin/uniquify.pl";       # Default route for uniquify script
my $qcRoute = "$ENV{SEQIMP_ROOT}/bin/QualityPlot.R";           # Default route for QualityPlot script
my $trimmitFfn = "$ENV{SEQIMP_ROOT}/bin/trimmit.pl";
my $paired_tally_Ffn = "$ENV{SEQIMP_ROOT}/bin/paired_tally.pl";       # Default route for paired end read uniquify script
my $tally_sstat_Ffn = "$ENV{SEQIMP_ROOT}/bin/tally_sstat_plot.R";     # Default route for paired tally sumstat plot script
my $addIndexFfn = "$ENV{SEQIMP_ROOT}/bin/add_read_index.pl";          # Default route for script to add read name indexes
my $tallyFfn = "tally";                                               # Default route to the tally program - note if passed to uniquify will trigger its use for uniquifying sequences in a more memory efficient fashion. If not Perl method will be used.
my $LCcutoff = "";
my $smallest = "";
my $biggest = "";
my $trimmitMax = "" ;
my $trimmitMin = "" ;
my $qcMax = "" ;
my $qcMin = "" ;
my $trimmitLC = "";
my $plotZeros = "";
#my $qcZeros = "";  090913
my $R_qcZeros = "";
my $geometry = "";
my $configDir = "$ENV{SEQIMP_ROOT}/config/reaper_config/";   # Default Reaper config directory
my $formatconfig = "$ENV{SEQIMP_ROOT}/config/seqimp_config/fasta_fastq_config.txt"; # Default format file : Currently this can't be altered.
my $fastq = "";
my $fiveTrim = "";
my $threeTrim = "";
my $reapData = "";
my $trimmitFive = "";
my $trimmitThree = "";
my $trimQCMax = "" ;
my $trimQCMin = "" ;
#my $trimqcZeros = "" ; 090913
my $R_trimqcZeros = "" ;
my $shipBowtieFfn = "$ENV{SEQIMP_ROOT}/bin/shipBowtie.pl";
my $depthFieldFfn = "$ENV{SEQIMP_ROOT}/bin/depthField.pl";
my $boxBowtieONEFfn = "$ENV{SEQIMP_ROOT}/bin/boxBowtie_Step1.R";
my $boxBowtieTWOFfn = "$ENV{SEQIMP_ROOT}/bin/boxBowtie_Step2.R";
my $boxBowtieQCFFn = "$ENV{SEQIMP_ROOT}/bin/boxBowtie_QC.R";
my $miRTableFfn = "$ENV{SEQIMP_ROOT}/bin/miR_table.R";
my $kraken_plot = "$ENV{SEQIMP_ROOT}/bin/kraken_plots.R";
my $trinucl_light = "$ENV{SEQIMP_ROOT}/bin/trinucl_complexity_light.pl";
my $chunkmbs    = 0 ;
my $genome      = ""; 
my $misMatches  = 0 ;
#my $stratum     = 0 ; # Option removed to control for bowtie biases - Now always set to "ON" - Check below.
my $hitNo       = 10;
my $ensversion  = -1 ;
my $mirversion  = -1 ;
my $bowtieLog   = "";
my $annotDir    = "";
my $versionDir ; 
my $perlUnq     = "";
my $fType       = "";
my $overlap     = 0 ;
my $paired      = 0 ;

my $bowtie_sam  = "";

my $configRoute;
my @reaperConfig;
my $tallyHo;
my @files;
my $uniqueCall;
my $qcplotCall;
my $considerBarcode;

my $geoConfig = "$ENV{SEQIMP_ROOT}/config/seqimp_config/geometry_config.txt";       # The default geometry config file

my $repMisMatches = 0;
my $repHitNo = 10;
my $repChunkmbs = 0;
my $repversion = -1;

my $no_merge = 0;
my $proportional = 0;
my $annot_conflict = "";
my $collapse_method = "default";  # added 180714 to allow miRNAs to be collapsed on sequence as well as ID.

my $sampleTag = "";

my $piBowtieFfn = "$ENV{SEQIMP_ROOT}/bin/piRNA_bowtie.pl";
my $piCovFfn = "$ENV{SEQIMP_ROOT}/bin/piRNA_coverage.R";
my $piPlotFfn = "$ENV{SEQIMP_ROOT}/bin/piRNA_plot.R";

my $baseDir2 = "";
my $formatChoice = "fasta"; # Default prefered format for sequence files is FASTA where possible

################################
### Organise command options ###
################################

my $OptSuccess = GetOptions( 
            # Universal options 
            "stage=s"      => \$stage,                   ### Flag to specify which stage of the pipeline to initiate - currently "reaper", "filter", "align" OR "features"
            "directory=s"  => \$baseDir,                 ### Project directory
            "plotZeros"    => \$plotZeros,               ### Should the qc length plot consider the maximum 0-length count when plotting?
            "fastq"        => \$fastq,                   ### Where possible configure output in fastq rather than fasta format.
            "help"         => \$help,                    ### Get help
            "debug"        => \$debug,                   ### Print extra output to aid debugging

            # Stage one options - Reaper
            "geometry=s"      => \$geometry,             ### Geomtry string specifying which config file to read
            "reapConfig=s"    => \$configDir,            ### Path to Reaper config files
            "geoConfig=s"     => \$geoConfig,            ### Seqeunce_imp geometry config file 
            "reaper-extra=s"  => \$reaperExtra,          ### Passed to reaper as extra input variables to be pasted into the command line call
            "meta=s"          => \$metadata,             ### Annotation file of barcodes and adapters
            "basicName=s"     => \$base,                 ### Passed to Reaper as a base name for output - Alphanumeric or "_"
            "data=s"          => \$rawdata,              ### Path to the directory to the sample filesystem... Should contain a data/ sub-directory with the FASTQ file
            "perlUnq"         => \$perlUnq,              ### If flag included Perl version of uniquify will be used in preference to the tally version
            "paired"          => \$paired,               ### Specifies if paired end run. Will cause pipeline to stop after Reaper step as tally must be run in parallel

            # Stage two options - Filtering 
            "low=i"        => \$LCcutoff,                ### Low complexity cutoff threshold
            "maxSize=i"    => \$biggest,                 ### Largest reads to be selected for downstream analysis
            "minSize=i"    => \$smallest,                ### Smallest reads to be selected for downstream analysis
            "five=i"       => \$fiveTrim,                ### Bases to be removed from the 5' end of all reads
            "three=i"      => \$threeTrim,               ### Bases to be removed from the 3' end of all reads
            "reapData=s"   => \$reapData,                ### Path to the Reaper output files

            # Stage three options - Align
            "chunk=i"         => \$chunkmbs ,
            "mismatches=i"    => \$misMatches,              ### number of mismatches to be allowed by Bowtie
#            "strata"          => \$stratum,                ### Remove strata option - Strata set by default. NOTE: Due to Bowtie reporting biases etc. Strata always active and k always reater than 2. This will ensure column 7 always reports total reads returned by BOWTIE not a count of all exactly matching sites on the same strnad in the genome. 
            "maxHits=i"       => \$hitNo,                   ### Maximum number of hits to be reported by Bowtie - minimum 0f 2 to allow identification of unique hits
            "sam"             => \$bowtie_sam,              ### Specifies if SAM format output to be used by Bowtie

            # Stage three/stage four options -Align/Features
            "genome=s"        => \$genome ,                 ### Specify species for which analysis should be conducted
            "ensversion=i"    => \$ensversion,              ### Ensembl version to use - only required in some stage 4 instances
            "annotationDir=s" => \$annotDir,                ### Annotation base directory (should contain subdirectories for both the Ensembl and miRBase annotations)
            
            # Stage four options
            "feature=s"       => \$fType,                   ### A selected feature type to be quanitated - currently supported: "miRNA"
            
            # Stage four (feature = miRNA)
            "mirversion=i"       => \$mirversion,              ### Required in some instances of stage 4 if miRNA feature specified. Which miRBase version to use
            "overlap=i"          => \$overlap,                 ### Required in some instances of stage 4 if miRNA feature specified. Desired overlap between the sequence reads and feature of interest
            "separate_loci"      => \$no_merge,                ### Optional flag specifying whether miRNAs with same mature ID should have their expression values merged.
            "proportional"       => \$proportional,            ### Optional flag specifying whether to ascribe reads according to unique counts only (Only considers miRNA loci when dividing reads)
            "annot_conflict=s"   => \$annot_conflict,          ### Require a value. This specifies how overlapping miRNA loci from miRBase are handled.
            "collapse_method=s"  => \$collapse_method,         ### Specifies method to be used to merge miRNA counts.

            # Stage four (feature = repeat)
            "repversion=i"        => \$repversion,         ### Required in some instances of stage 4 if repeat feature specified. Specifies repeat annotation version to be selected. 
            "repMismatches=i"     => \$repMisMatches,       ### Required in some instances of stage 4 if repeat feature specified. Number of mismatches to allow between the sequences and canonical repeat
            "repMaxHits=i"        => \$repHitNo,            ### Required in some instances of stage 4 if repeat feature specified. Number of hits to allow between the sequences and the canonical repeat
            "repChunk=i"          => \$repChunkmbs,        ### Required in some instances of stage 4 if repeat feature specified. 

            # imp_commandline.pl synchronisation options
            "tag=s"           => \$sampleTag,               ### tag to indicate to user what sample is being run.
         
            # Paired sequence imp options
            "directory2=s"    => \$baseDir2               ### The second base directory containing all of the information relating to the second paired fastq sample.
         
         );


die "Unacceptable option specified, please check\n" if (!$OptSuccess);

print STDERR "\n************************************************************\n";

if ($sampleTag) {
   print STDERR "\n### NEW ANALYSIS:\n### Lane: $sampleTag ###";
}

print STDERR "\n############################\n## Begin the Imp Pipeline ##\n############################\n\n";


##############################
### HELP SUBROUTINE ##########
##############################

sub help {
   print <<EOH;
Usage:
   $progname [options]
Options:
--help            <FLAG>       Displays this page

___________________________________________________________________________

### UNIVERSAL OPTIONS ###

# Required:
--stage           <STRING>     The stage of the pipeline to initiate (currently "reaper", "filter", "align" OR "features"): REQUIRED
--directory       <STRING>     The base directory containing the data directory and metadata directory: REQUIRED

# Optional:
--plotZeros       <FLAG>       A flag specifying whether reads of length 0 should be considered when plotting QC plots: OPTIONAL
--paired          <FLAG>       Specifies whether paired files are being analysed. OPTIONAL - Default = "Unpaired"
--fastq           <FLAG>       Should the pipeline prefer output in fastq format where possible - Default = FASTA
___________________________________________________________________________

### STAGE ONE OPTIONS: reaper ###

# Required:
--geometry        <STRING>     The read geometry (currently "no_barcode", "5p_barcode", "3p_barcode", 
                               "5p_barcode_and_insert" and "3p_barcode_and_insert") : REQUIRED
--reapConfig      <STRING>     File path to the reaper geometry configuration directory : OPTIONAL - Default will be supplied config files
--geoConfig       <STRING>     File path to the sequence_imp.pl geometry file that contains information concerning 
                               which reaper config file is required for each specified geometry: OPTIONAL - Default will be supplied config files 
--meta            <STRING>     File path from the 'base directory' to the file containing sample adapter sequences and barcodes: REQUIRED
--data            <STRING>     File path from the 'base directory' to the sub-directory that contains a single FASTQ data file: REQUIRED

# Optional:
--basicName       <STRING>     Base name used as file prefix for results. WARNING: Alphanumeric or '_': OPTIONAL - Default = 'reap'
--reaper-extra    <STRING>     Any additional options to be passed directly to Reaper: OPTIONAL
--perlUnq         <FLAG>       Flag specifying whether perl should be used in the place of 'Tally' to collapse 
                               filtered FASTA files to unique reads: OPTIONAL - Default = 'Tally'
___________________________________________________________________________

### STAGE TWO OPTIONS: filter ###

# Required:
--reapData        <STRING>     The directory into which all reaper data has been placed relative to the 'base directory': REQUIRED

# Optional:
--low             <INTEGER>    Maximum low complexity score filter: OPTIONAL
--minSize         <INTEGER>    Minimum trimmed read length filter: OPTIONAL
--maxSize         <INTEGER>    Maximum trimmed read length filter: OPTIONAL
--five            <INTEGER>    Bases to be removed from the 5' end of all trimmed reads: OPTIONAL
--three           <INTEGER>    Bases to be removed from the 3' end of all reads: OPTIONAL
--directory2      <STRING>     The base directory for the second, paired sample to allow two sequence files to be processed together: OPTIONAL (REQUIRED if --paired specified)
___________________________________________________________________________

### STAGE THREE AND FOUR OPTIONS: align and features ###

## BOTH STAGES:
#Required:
--genome          <STRING>     The species for which the analysis will be performed: REQUIRED
--ensversion      <INTEGER>    Version number for the gene annotation to be used for the analyis. WARNING: Ensure the ensembl genome build and
                               miRBase builds are consistent: REQUIRED for 'align' and some 'features' (eg. --feature=repeat)
--annotationDir   <STRING>     The base directory in which all pipeline associated annotation is held (eg. Contains ENSEMBL and MIRBASE 
                               directories): REQUIRED

## STAGE THREE ONLY: align
# Optional:
--chunk           <INTEGER>    Value to be passed to the Bowtie chunkmbs option: OPTIONAL
--mismatches      <INTEGER>    Maximum number of mismatches to be allowed by Bowtie: OPTIONAL - Default = 0
--maxHits         <INTEGER>    Maximum number of hits to be reported by Bowtie (Reads that hit above this threshold discarded): OPTIONAL - Default = 10
--sam             <FLAG>       Specifies whether to run the pipeline with SAM/BAM compatible alignment format: OPTIONAL
NOTE: Bowtie always run with --strata flag set 

## STAGE FOUR ONLY: features
# Required:
--feature         <STRING>     A selected feature type to be quanitated - currently supported: "miRNA" and "repeat": REQUIRED

## STAGE FOUR ONLY: features: --feature=miRNA
# Required:
--mirversion      <INTEGER>    Required in some instances of stage four if miRNA feature specified. WARNING: Ensure the ensembl genome build and
                               miRBase builds are consistent: REQUIRED
--overlap         <INTEGER>    The nucleotide overlap between a read and a feature required for it to be included in the feature read count: REQUIRED 
--separate_loci   <FLAG>       Specifies whether the counts for miRNA loci with the same miRBase mature sequence ID should be merged: OPTIONAL - Default = merge
--proportional    <FLAG>       Specifies whether to divide reads between all alignments equally or divide multimappers between miRNA loci according to the number of uniquely mapping reads ascribed to each: OPTIONAL - Default = equal division
--annot_conflict  <STRING>     This determines how the pipeline will handle overlapping miRBase loci - currently supported: "ignore", "merge", "remove" : REQUIRED
--collapse_method <STRING>     Specified the method that will be used to merge miRNA counts unless 'seperate_loci' flag is set - currently supported: "mature_id" (for miRBase mature ID), "sequence" (for miRBase mature sequence) OPTIONAL - Default = mature_id

## STAGE FOUR ONLY: features: --feature=repeat
# Optional:
--repversion      <INTEGER>   The version number of the NCBI/repeat annotation to be used for the comparison: REQUIRED.
--repMaxHits      <INTEGER>   The maximum number of alignments allowed between a read and the canonical repeat sequences: OPTIONAL - Default 10   
--repMismatches   <INTEGER>   The maximum number of mismatches allowed between a read and the canonical repeat sequences or reference genome (used for normalisation): OPTIONAL - Default 0
--repChunk        <INTEGER>   Value to be passed to the Bowtie chunkmbs option: OPTIONAL

____________________________________________________________________________

##### PAIRED SEQUENCE IMP OPTIONS:
## STAGE: paired-tally

# Required:
--directory2      <STRING>    A second base directory containing the information relating to the second, paired fastq file: REQUIRED

____________________________________________________________________________
EOH
exit(0)
}

help() if $help;

#####################
### SANITY CHECKS ###
#####################

print STDERR "\n##################\n### DEBUG MODE ###\n##################\n\n" if $debug;

# Universal sanity checks

if(length($baseDir) == 0){
   die "Need a relative path to the directory of project files: --dir=Directory";
}

if (! -d $baseDir){
    die "Aha! $baseDir does not seem to be a directory please specify an appropriately constructed project directory";
}

### Stages acceptable to pipeline depend on whether the pipeline is run in standard or paired end mode BUT a stage is always required.

if ($paired){
   print STDERR "Running Sequence imp in paired end mode\n\n";
   die "Require a specified pipeline stage (--stage=reaper/filter)\n" if length($stage)==0;

   ############################## Alter check to allow paired filter stage to be called
   
   die "Paired sequence imp can only be initiated with the 'reaper' or 'filter' stages specified\n" if ($stage ne "reaper" && $stage ne "filter");
}else{
   die "Require a specified pipeline stage (--stage=reaper/filter/align/features)\n" if length($stage)==0;
}

print STDERR "Specified stage is $stage\n";


if ($stage eq "reaper"){
   
   print STDERR "Beginning stage one of the pipeline - REAPER call and quality control\n\n";
   
   # Stage one sanity checks

   if(length($rawdata) == 0){
      die "Need a directory for the raw data relative to the specified project directory: --data=Directory\n";
   }
   if(length($metadata) == 0){
      die "Need a file containing barcode and adapter sequences relative to the specified project directory: --meta=Directory/metadata.txt\n";
   }
   if(length($configDir) == 0){
      die "Need a config directory containing geometry specs for Reaper: --reapConfig=Directory\n";
   }
   if(length($geometry) == 0){
      die "Need a geometry specification: --geometry=Geometry_Spec\n";
   }
   
   if (! -d $configDir){
      die "Aha! $configDir does not seem to be a directory please specify a directory containing the geometry specification files\n";
   }

   if(length($geoConfig) == 0){
      die "Need a sequence imp geometry specification file --geoConfig=File\n";
   }

   die "The geometry configuration file seems to have some problems\n" if (! -s $geoConfig);
   die "The geometry configuration file does not seem to be a plain file\n" if (! -f $geoConfig);

   if (!length($perlUnq)==0){
      $tallyFfn ="";
      print STDERR "Will use Perl based uniquify system\n\n";
   }else{
      print STDERR "Will use C based uniquify system as set by default\n\n";
   }

   # Stage two sanity checks
}elsif($stage eq "filter"){

   if(length($reapData) == 0){
      die "Need a directory containing the Uniquify output: --reapData=Directory\n";
   }

   print STDERR "Beginning stage two of the pipeline - User defined read processing steps\n";

   if($paired){
   
      ########## Parameters that only need to be checked when the pipeline is filtering in paired mode (Performs simultaneous Tally)

      die "Require a second directory in 'paired-tally' mode\n" if (length ($baseDir2)==0);
      die "$baseDir2 does not seem to exist\n" if (! -d $baseDir2); 

   }

   # Stage three sanity checks
}elsif($stage eq "align"){

   print STDERR "Beginning stage three of the pipeline - Aligning reads to a genome\n";
   print STDERR "All alignments will be presented in SAM format\n" if $bowtie_sam;

   if (length($annotDir)==0){
      die "Require a path to the annotation base directory containing both Ensembl and miRBase subdirectories: --annotationDir=<PATH>\n";
   }

   if(!-d $annotDir){
      die "The --annotationDir directory specified does not exist, please check path\n";
   }

   if ($hitNo <1){
      die "Require at least a maxHits setting of 1 to allow for unique hits to be identified\n";            ### Reads that reach maximum hits will later be removed with the -m option to ensure that strand biases in reporting due to strand orientated reporting heirachy will be removed.
   }

   #if ($stratum > 0 && $hitNo < 2 ) {                                                                      ### Strata now set as default to account for Bowtie reporting biases 
   #   die "--strata only appropriate if --maxHits >1: See Bowite manual --strata -k respectively\n";
   #}  

   if($ensversion <= 0){
      die "Require a specified Ensembl version (a positive integer): --ensversion=<INTEGER>\n";
   }

   if (length($genome)==0){
      die "Require a genome to work with: --genome=<SPECIES>\n";
   }

#   if ($genome ne "mouse"){
#      die "$genome is not a recognised genome\n";                                                          ### Genome control depricated at this stage. Limited genome annotation available now from source. Genome sequence selection should be automated based on file structure.
#   }

}elsif($stage eq "features"){

   print STDERR "Beginning stage four of the pipeline - Analysis of genomic features\n";
      
      if (length($fType)==0){
         die "Require a specified feature type to work with: --feature=<FEATURE>\n";
      }
     

      if($fType eq "miRNA"){
         if($mirversion <= 0){
            die "Require a specified miRBase version (a positive integer): --mirversion=<INTEGER>\n";
         }
         if ($overlap ==0){
            die "Require a minimum overlap to be specifed for reads and features for the read to be included in feature depth calculations --overlap=<INTEGER>\n";
         }

         die "Can not specify a 'collapse_method' in the analysis config file if 'the seperate_loci' flag is used\n" if ($no_merge && $collapse_method ne "default");
         $collapse_method = "mature_id" if $collapse_method eq "default";
         die "--collapse_method=<STRING> can only take the options 'mature_id' and 'sequence'\n" if (($collapse_method ne "mature_id") && ($collapse_method ne "sequence"));
         
         if (length($annot_conflict)==0){
            die "Require a method to resolve overlapping miRNA loci --annot_conflict=<ignore/remove/merge>\n";
         }
         die "annot_conflict can only take values 'ignore', 'remove' or 'merge'\n" if ($annot_conflict ne "ignore" && $annot_conflict ne "remove" && $annot_conflict ne "merge");

      }elsif($fType eq "repeat"){
         if($repversion <= 0){
            die "Require a specified repeat version (a positive integer): --repversion=<INTEGER>\n";
         }
         if ($ensversion <= 0){
            die "Require a specified Ensembl version (a positive integer): --ensversion=<INTEGER>\n";
         }
         if ($repHitNo <1){
            die "Require at least a --repMaxHits setting of 1 to allow for unique hits to be identified\n";
         }
         if ($repMisMatches <0){
            die "Require a --repMismatches setting of 0 or more\n";
         }
         
      }else{
         die "Do not recognise feature type: Please try again (--feature).\n";
      }

      if (length($genome)==0){
         die "Require a genome to work with: --genome=<SPECIES>\n";
      }

      if (length($annotDir)==0){
         die "Require a path to the annotation base directory containing all of the relevant annotation directories: --annotationDir=<PATH>\n";
      }

      if(!-d $annotDir){
         die "The --annotationDir directory specified does not exist, please check path\n";
      }


}else{
   die "Require a specified pipeline stage (--stage=reaper/filter/align/features)\n";
}


###########################################
### READ AND CHECK FORMAT CONFIGURATION ###
###########################################


sub readformat {
  
   print STDERR "Reading format table\n";
   print STDERR "WARNING: &readformat should be passed a single file!\n" if @_ != 1;
   
   my $formatFfn = shift @_ ;

   die "Please check that file $formatFfn is compatible with the pipeline. Should be a tab delimited text file\n"
         if (! -e $formatFfn || -z _ || ! -f _ );
         
   my $header;
   my @header;
   my %columns = ();
   my $columnNo = 0;
   my %loadedTable = ();
   my %rows;
   
   open TABLE, "< $formatFfn" || die "Can't open $formatFfn: $!\n";
   
   $header = <TABLE>;
   chomp $header;
   @header = split "\t", $header, -1;
   $columns{$columnNo++} = $_++ for @header;
   
   while (<TABLE>) {
      chomp;
      my @gubbins = split "\t", $_,-1;
      die "The number of fields does not match the number of columns for row \'$gubbins[0]\'. Please check that the description file is tab seperated.\n" if ((scalar @header) != (scalar @gubbins)); 

      die "'rowname' column does not contain unique identifiers\n" if exists($rows{$gubbins[0]});
      $rows{$gubbins[0]}++;

      $loadedTable{$gubbins[0]} = {};
      for (my $i=1; $i < @header; $i++) {
         $loadedTable{$gubbins[0]}{$columns{$i}} = $gubbins[$i];
      }
   }
   close(TABLE);
   print STDERR "Reading of table complete\n";
   return(\%loadedTable);
}

die "Can't find the file format configuration file: $formatconfig" if (! -s $formatconfig);
my $formatTable = &readformat($formatconfig);

$formatChoice = "fastq" if $fastq;

#print "$formatChoice\n";
#print STDERR $formatTable->{fastq}{reaper_output}."\n";


#################
#################
### STAGE ONE ###
#################
#################

if ($stage eq "reaper"){

##################################################
### Checking geometry specification for Reaper ###
##################################################

   my $geoConfigLine;
   my @headerChoices = ();
   my %geoCols;
   my $geoId = 0;
   my %geomtable = ();

   die "Unable to open the sequence_imp geometry config file\n" unless open(GEOM, $geoConfig);

   $geoConfigLine = <GEOM>;
   chomp $geoConfigLine;
   @headerChoices = split "\t", $geoConfigLine;
   $geoCols{$geoId++} = $_++ for @headerChoices;
   
   while (<GEOM>) {
      chomp;
      my @geoLineChew = split "\t", $_,-1;
      $geomtable{$geoLineChew[0]} = {};
      for (my $i=1; $i<@headerChoices;$i++) {
         $geomtable{$geoLineChew[0]}{$geoCols{$i}} = $geoLineChew[$i];
      }
   }

   close(GEOM);
   
   print STDERR "Selecting the config file for a ",$geometry," Reaper run\n";
   die "Unrecognised geometry specified, please check\n" if (!exists $geomtable{$geometry});
   $configRoute = "$configDir/".$geomtable{$geometry}{"Config"};
   die "Please check that the ",$geomtable{$geometry}{"Config"}," file is present in the $configDir directory" unless -f $configRoute;

   #### Removed 120511: Updated to use aconfig file approach with geometries and associated data - see end
   ########### Alterations to geometry specifications
   ### Parse Config file with geomety names - Add additional columns for relevant Reaper config file
   ### Pass the config file to sequence_imp via a command line option
   ### Treat annotation slightly differently - deal with Ensembl version in mother_of_the_imp as this is more flexible. All geometry config is set in stone.

#   if($geometry eq "Illumina_3p_Barcode"){
#      $configRoute = "$configDir/Ill_3p_bar.config" ;
#      print STDERR "Selecting the config file for an Illumina, 3-prime-barcoded Reaper run\n";
#      die "Please check that the Ill_3p_bar.config file is present in the $configDir directory" unless -f $configRoute; 
#   } elsif ($geometry eq "Illumina_5p_Barcode") {
#      $configRoute = "$configDir/Ill_5p_bar.config" ;
#      print STDERR "Selecting the config file for an Illumina, 5-prime-barcoded Reaper run\n";
#      die "Please check that the Ill_5p_bar.config file is present in the $configDir directory" unless -f $configRoute;
#   }elsif($geometry eq  "Illumina_No_Barcode"){
#      $configRoute = "$configDir/Ill_no_bar.config";
#      print STDERR "Selecting the config file for an Illumina, non-barcoded Reaper run\n";
#      die "Please check that the Ill_no_bar.config file is present in the $configDir directory" unless -f $configRoute;
#   } else {
#      die "Unrecognised geometry specified, please check" ;
#   }
   
   ### Read specified Reaper config file and join into a Reaper command call
   
   if (! open CONFIG, $configRoute){
      die "Unable to open $configRoute, try checking permissions";
   }
   
   while(<CONFIG>){
       chomp;
       push @reaperConfig, $_;
   }
   
   close (CONFIG);
   
   my $reaperline = join "\t", @reaperConfig;
   #print "$reaperline\n";
   
   
   ##################################################################
   ### Checking input directories and creating output directories ###
   ##################################################################  ################ STAGE ONE ONLY
  
   print STDERR "\nChecking directories\n";

   my $metaPath = "$baseDir/$metadata";
   if (! -f $metaPath){
      die "Aha! $metaPath is not a plain file - please try again.";
   }
   
   my $dataDir = "$baseDir/$rawdata";
   if (! -d $dataDir){
      die "The specified data directory does not exist, please check";
   }
   
   my $qcDir = "$baseDir/QC/";
   if (-e $qcDir){
      print STDERR "WARNING: The QC directory already exists. Will use current directory ($qcDir). This may cause files to be overwritten.\n"; # The imp will no longer die if rerunning/overwiring - 041111
   }
      #   die "QC Directory exists at $baseDir, please tidy up";
   #}else{   
   #   mkdir $qcDir, 0755 or die "Cannot make QC directory: $!"; 
   #}
   die "Cannot make QC directory $qcDir: $!" if (system("mkdir -p -m 0755 $qcDir"));


   my $reaperDir = "$baseDir/REAPER/";
   if (-e $reaperDir){
      print STDERR "WARNING: The REAPER directory already exists. Will use current directory ($reaperDir). This may cause files to be overwritten.\n"; # The imp will no longer die if rerunning/overwiring - 041111
   }
      #   die "REAPER Directory exists at $baseDir, please tidy up"; 
   #}else{   
   #  mkdir $reaperDir, 0755 or die "Cannot make REAPER directory: $!"; 
   #}
   die "Cannot make REAPER directory $reaperDir: $!" if (system("mkdir -p -m 0755 $reaperDir"));  
   
   ##################################
   ### FINDING THE RAW INPUT FILE ###
   ################################## ############## STAGE ONE ONLY
   print "DataDir = $dataDir\n" if $debug;
   
   #print "[$dataDir]\n";
   my @datafiles = <$dataDir/*>;
   #print "have [@datafiles]\n";
   if (@datafiles != 1) {                                   # Expecting only a single input file in the input directory
      die "Expecting a single datafile in $dataDir";
   }
   
   print "DataFiles = @datafiles\n" if $debug;
   
   my $dataPath = "$datafiles[0]";
   #print "Hello [$dataPath] [$dataDir]\n";
   
   ################################
   ### REAPER CALL MUST GO HERE ###
   ################################

   my $formatCleanOpt = "";

   #### Arranging format according to user requirements - Paired reads must contain the read record number

   if ($paired){
      $formatCleanOpt = "-format-clean ".$formatTable->{$formatChoice}{reaper_output_paired};
   }else{
      $formatCleanOpt = "-format-clean ".$formatTable->{$formatChoice}{reaper_output};
   }
      
   my $reaperCall = "$reaperRoute -meta $metaPath -i $dataPath $formatCleanOpt -basename $reaperDir/$base $reaperline $reaperExtra";
   
   print STDERR "\n\nPreparing to trim and filter sequences: REAPER\n\n";
   print STDERR "-----------------------\nThe following command is being passed to Reaper:\n\n$reaperCall\n-----------------------\n\n";
   
   die "Reaper hasn't worked: $reaperCall\n" unless system($reaperCall)==0 ;
   
   
   #############################################################################
   ### Checking all the expected input files are present for QC and Uniquify ###
   #############################################################################
   print STDERR "\nChecking Reaper output\n";
   
   ### Collect the barcodes from the metadata file specified and used by Reaper
   
   my $barcodePattern;
   my $anyBarcode;

   if ($geomtable{$geometry}{"Barcodes"} eq "YES"){
      $barcodePattern = "barcode";
      $anyBarcode = 1;
   }elsif($geomtable{$geometry}{"Barcodes"} eq "NO"){
      $barcodePattern = "NO_PATTERN";
      $anyBarcode = 0;
   }else{
      die "Problem with 'Barcodes' specification in geometry specification file for",$geometry,"\n";
   }


   ### Removed 120511 - Replaced with new config file system and altered to operate with the new Reaper version   
#   # Identifying the column containing the barcode in the metadata file   ############# ALTER BARCODE COLLECTION FROM BARCODE FILE - Column name will be constant 1/0 will be assigned based on Barcode YES/NO column.
#   
#   if($geometry eq "Illumina_3p_Barcode"){
#      $barcodePattern = "^3p-bc\$" ;
#      $anyBarcode = 1;
#      print STDERR "Setting up a screen for the correct Reaper output for a 3 prime barcoded sample\n"; 
#   }elsif ($geometry eq "Illumina_5p_Barcode") {
#      $barcodePattern = "^5p-bc\$" ;
#      $anyBarcode = 1;
#      print STDERR "Setting up a screen for the correct Reaper output for a 5 prime barcoded sample\n";
#   }elsif($geometry eq "Illumina_No_Barcode"){
#      $barcodePattern = "NO_PATTERN"; 
#      $anyBarcode = 0;
#      print STDERR "The expected Reaper output will not be configured for barcodes\n";
#   }else{
#      die "Unrecognised geometry\n";
#   }
   
   # Checking all Reaper output files are present and correct
   
   my $file_result = reaper_file_check($anyBarcode,$barcodePattern, $base, $metaPath, $reaperDir, $debug);
      
   print STDERR $file_result ; 
      
   $R_qcZeros = length $plotZeros == 0 ? "--noZero=29" : "";
   $considerBarcode = "--multi=$anyBarcode";
   my $reapDebug = $debug ? "--debug" : "" ;
   my $R_reapDebug = $debug ? "--debug=29" : "" ;
   my $R_pairedOpt = $paired ? "--paired=29" : "";
   my $formatOpt = length($fastq) != 0 ? "--fastq" : "";

   ### YUP!

   if (! $paired){
  
      ###############################################
      ### Constructing subsiduary programme calls ###
      ###############################################
      
      my $tallyOpt = length($tallyFfn) != 0 ? "--tallyffn=$tallyFfn --inFormat=".$formatTable->{$formatChoice}{tally_input} : "";   

      #print "Parametres: $uniqMax\t $qcMax\n $uniqMin\t $qcMin\n $uniqLC \n" ;
      
      $uniqueCall = "$uniquifyRoute --input=$reaperDir $tallyOpt $formatOpt $reapDebug"; 
      
      print STDERR "\nPreparing to collapse duplicate sequences\n";
      print STDERR "$uniqueCall\n" if $debug;

      ### YUP!

      die "Uniquify hasn't worked: $uniqueCall\n" unless system($uniqueCall)==0 ;
      
      $qcplotCall = "R --slave --quiet --vanilla --args --stage=1 --dataDir=$reaperDir --outDir=$qcDir --basicID=$base --plotfunc=$kraken_plot $R_qcZeros $considerBarcode $R_reapDebug $R_pairedOpt <$qcRoute";
      
      print STDERR "\nPreparing to plot QC data for Reaper\n";
      print STDERR "$qcplotCall\n" if $debug;
      die "QC plots haven't been plotted: $qcplotCall\n" unless system($qcplotCall)==0 ;
      
      print STDERR "\n### COMPLETED: Fastq file has been trimmed. Please check all of the files produced. ###\n\n";
      
      exit(0);

   }else{
   
      print STDERR "\nSequence Imp is skipping Tally at this point to allow it to be inititated on paired files in parallel.\n";
     
      print STDERR "\nSummarising trinucleotide complexity for paired end files\n";

      my $simple_triCall = "$trinucl_light --input_directory=$reaperDir $formatOpt"; 

      die "Unable to summarise trinucleotide complexity for $reaperDir\n" if system($simple_triCall);


      ####################################### File check for trinucl output.

      $qcplotCall = "R --slave --quiet --vanilla --args --stage=1 --dataDir=$reaperDir --outDir=$qcDir --basicID=$base --plotfunc=$kraken_plot $R_qcZeros $considerBarcode $R_reapDebug $R_pairedOpt <$qcRoute";
      
      print STDERR "\nPreparing to plot QC data for Reaper\n";
      print STDERR "$qcplotCall\n" if $debug;
      die "QC plots haven't been plotted: $qcplotCall\n" unless system($qcplotCall)==0 ;
      
      print STDERR "\n### COMPLETED: Fastq file has been trimmed. Please check all of the files produced. ###\n\n";
      
      exit(0);

   }
}


#################
#################
### STAGE TWO ###
#################
#################

if ($stage eq "filter"){

   ### Organising the files and directories for filter stages (paired and unpaired)

   ### Input directory

   print STDERR "\nChecking directory structures\n";

   my $reapFfn = "$baseDir/$reapData";
   die "The data directory ($reapFfn) does not exist, please check\n" if (! -d $reapFfn);

   ### Output directories
   
   my $qcDir = "$baseDir/QC/";
   die "QC Directory cannot be accessed at $qcDir, please check file structure\n" if (! -d $qcDir);

   my $procDir = "$baseDir/PROCESSED/";
   print STDERR "WARNING: The PROCESSED directory already exists ($procDir). Will use current directory. This may cause files to be overwritten.\n" if (-e $procDir); # The imp will no longer die if rerunning/overwiring - 041111
   die "Cannot make PROCESSED directory $procDir: $!" if (system("mkdir -p -m 0755 $procDir"));
 
   ### Set up second directory set if paired

   my $reapFfn2;
   my $qcDir2;
   my $procDir2;

   if ($paired){
      $reapFfn2 = "$baseDir2/$reapData";
      die "The data directory ($reapFfn2) does not exist, please check\n" if (! -d $reapFfn2);
   
      $qcDir2 = "$baseDir2/QC/";
      die "QC Directory cannot be accessed at $qcDir2, please check file structure\n" if (! -d $qcDir2);

      $procDir2 = "$baseDir2/PROCESSED/";
      print STDERR "WARNING: The PROCESSED directory already exists ($procDir2). Will use current directory. This may cause files to be overwritten.\n" if (-e $procDir2); # The imp will no longer die if rerunning/overwiring - 041111
      die "Cannot make PROCESSED directory $procDir2: $!" if (system("mkdir -p -m 0755 $procDir2"));
   }

#      die "PROCESSED Directory exists at $baseDir, please tidy up";
#   }else{
#      mkdir $procDir, 0755 or die "Cannot make PROCESSED directory: $!";
#   }

   #########################################################################
   ### Identify base name based upon basename of the reaper output files ###
   #########################################################################
   print "Paired: $paired\n" if $debug;
   if (!$paired){
      #print "Hello\n";
      my $experiment = "";
      my %uniqueID = ();
      my @expts = ();
      #print "Hello\n";
      print STDERR "DIRECTORY: $reapFfn\n" if $debug;
      my @reaperfiles = <$reapFfn/*>;
      print "REAPER Files pre trim:",@reaperfiles,"\n" if $debug;
      foreach (@reaperfiles){
         s/$reapFfn\///;
      }
      print "REAPER Files:",@reaperfiles,"\n" if $debug;

      my @cleanFiles = grep /\.clean\.gz$/, @reaperfiles;

      print "CLEAN Files:",@cleanFiles,"\n" if $debug;

      for my $cleanNow (@cleanFiles){
         if ($cleanNow =~ /^([\w_]+)\.\w+\.clean\.gz$/){
            $uniqueID{$1} = 1;  
         }else{
            die "Can not identify the unique ID for reaper clean file $cleanNow\n";
         }
      } 

      print "EXPT IDs:",(keys %uniqueID),"\n" if $debug;

      if ((scalar(keys %uniqueID)) != 1){
         die "There aren't the expected number of experimental IDs in the $reapFfn directory\n";
      }else{
         print STDERR "The correct number of Experimental IDs have been identified files have been found\n" if $debug;
         @expts = keys %uniqueID; 
         $experiment = $expts[0];
         print "Experiment:$experiment\n" if $debug;
      }
      
      #############################
      ### Set up a trimmit call ###
      #############################

      if (length $biggest != 0){
         $trimmitMax = "--upper=$biggest" ;
         $trimQCMax = "--maxSize=$biggest";
      }
      if (length $smallest != 0){
         $trimmitMin = "--lower=$smallest";
         $trimQCMin = "--minSize=$smallest";
      }
      if ( length $LCcutoff != 0){
         $trimmitLC = "--tri=$LCcutoff";
      }
      if ( length $fiveTrim != 0){
         $trimmitFive = "--fiveP=$fiveTrim";
      }
      if ( length $threeTrim != 0){
         $trimmitThree = "--threeP=$threeTrim";
      }
      if (length $plotZeros == 0){
         $R_trimqcZeros = "--noZero=29";
      }

      my $R_trimmitDebug = $debug ? "--debug=29" : "" ;

      my $trimmitCall = "$trimmitFfn --inDir=$reapFfn --outDir=$procDir --tally=$tallyFfn $trimmitMax $trimmitMin $trimmitLC $trimmitFive $trimmitThree $R_trimmitDebug";
      
      my $trimSaved = "$procDir/total.saved.trimmit.tab"; 
      
      my $trimQCsaved = "--saved=$trimSaved";

      my $trimqcplotCall = "R --slave --quiet --vanilla --args --stage=2 --dataDir=$procDir --outDir=$qcDir --basicID=$experiment --plotfunc=$kraken_plot $trimQCsaved $R_trimqcZeros $trimQCMax $trimQCMin $R_trimmitDebug <$qcRoute";

      print STDERR "\nPreparing to begin filtering by user defined criteria\n";
      print STDERR "$trimmitCall\n$trimqcplotCall\n" if $debug;

      die "Timmit hasn't worked: $trimmitCall: $! \n" unless system($trimmitCall)==0 ;
      print STDERR "Filtering complete\n"; 
      die "It seems $trimSaved is not where expected. Please take a look.\n" if (! -s $trimSaved);
      
      print STDERR "\nPreparing to plot QC data\n";
      die "The Trimmit QC plot hasn't worked: $trimqcplotCall: $! \n" unless system($trimqcplotCall)==0 ;
      print STDERR "QC plots complete\n";

      print STDERR "\n### COMPLETED: Filtered and trimmed all of the reads to the desired specifications. ###\n\n";
      exit(0);

   }else{
      
      my $pairfilterMin  = "";
      my $pairfilterMax  = "";
      my $pairfilterLC   = "";
      my $pairfilterFive = "";

      die "Paired end filtering currently doesn't support 3' read trimming. Please set 'three' in the configuration file to NA.\n"if length $threeTrim != 0;
      if (length $smallest != 0){
         $pairfilterMin = "--lower=$smallest";
      }
      if (length $biggest != 0){
         $pairfilterMax = "--upper=$biggest" ;
      }
      if ( length $LCcutoff != 0){
         $pairfilterLC = "--tri=$LCcutoff";
      }
      if ( length $fiveTrim != 0){
         $pairfilterFive = "--fiveP=$fiveTrim";
      }

      print STDERR "\nBeginning to uniquify paired fastq files in parallel\n";
   
      my @inputsONE = <$reapFfn*\.clean\.gz>;
      my @inputsTWO = <$reapFfn2*\.clean\.gz>; 
   
      die "Too many trimmed sequence files (.clean.gz) in $reapFfn. Currently do not support multiplexed, paired end reads.\n" if ((scalar @inputsONE) != 1);
      die "Too many trimmed sequence files (.clean.gz) in $reapFfn2. Currently do not support multiplexed, paired end reads.\n" if ((scalar @inputsTWO) != 1);
   
      my $outputONE = "";
      my $intermedONE = "";
      my $outputTWO = "";
      my $intermedTWO = "";
      my $sstatONE = "";
      my $sstatTWO = "";
      my $readFormat = "";
      my $readOutFormat = "";
      my $baseNameONE = "";
      my $baseNameTWO = "";
      my $qcPlotONE = "";
      my $qcPlotTWO = "";

      if ($inputsONE[0] =~ /^$reapFfn(\S+)\.clean\.gz$/){
         if ($fastq){
            # 170217: Swap outputONE and outputTWO for a new temporary intermediate file - this will allow the subsequent indexing of reads before the final output is produced.
            $intermedONE = "$procDir/$1.intermed.fastq.gz";
            $outputONE = "$procDir/$1.tallied.fastq.gz";
         }else{
            $intermedONE = "$procDir/$1.intermed.fasta.gz";
            $outputONE = "$procDir/$1.tallied.fasta.gz";
         }
         $sstatONE = "$procDir/$1.proc.sumstat";
         my $baseNameTemp = $1;
         if ($baseNameTemp =~ /(\S+)\.\w+/){
            $baseNameONE = $1;
            $qcPlotONE = "$qcDir/$1_Processed_reads_qc.pdf";
         }
      }else{
         die "Not able to identify the base name for $inputsONE[0]";
      }
   
      if ($inputsTWO[0] =~ /^$reapFfn2(\S+)\.clean\.gz$/){
         if ($fastq){
            $intermedTWO = "$procDir2/$1.intermed.fastq.gz";
            $outputTWO = "$procDir2/$1.tallied.fastq.gz";
         }else{
            $intermedTWO = "$procDir2/$1.intermed.fasta.gz";
            $outputTWO = "$procDir2/$1.tallied.fasta.gz";
         }
         $sstatTWO = "$procDir2/$1.proc.sumstat";
         my $baseNameTemp = $1;
         if ($baseNameTemp =~ /(\S+)\.\w+/){
            $baseNameTWO = $1;
            $qcPlotTWO = "$qcDir2/$1_Processed_reads_qc.pdf";
         }
      }else{
         die "Not able to identify the base name for $inputsTWO[0]";
      }
  
      ### Organise output and input formats for tally:

      $readFormat = "--inFormat=".$formatTable->{$formatChoice}{tally_input_paired};
      $readOutFormat = "--outFormat=".$formatTable->{$formatChoice}{tally_output_paired};

      my $ptally_fastq_Opt = $fastq ? "--fastq" : "" ; 

      ##################################### Will need to incorporate maximum and minimum lengths, trinucleotide scores and 5' trimming.
      # Also need to organise sumstat output and links between paired directories.
   
      my $ptally_call = "$paired_tally_Ffn --tally=$tallyFfn $readFormat $readOutFormat --inONE=$inputsONE[0] --inTWO=$inputsTWO[0] --outONE=$intermedONE --outTWO=$intermedTWO $readFormat $readOutFormat --stat=$sstatONE $pairfilterMin $pairfilterMax $pairfilterLC $pairfilterFive $ptally_fastq_Opt";
   
      print $ptally_call."\n" if $debug;
      die "Paired end Tally has failed: $ptally_call\n" if system($ptally_call);
  
      die "The first Tally output file appears to be missing! Please check." if (! -s $intermedONE);
      die "The second Tally output file appears to be missing! Please check." if (! -s $intermedTWO);
      die "The Tally sstat output file appears to be missing! Please check." if (! -s $sstatONE);

      print STDERR "\nCopying sumstat summary file to both of the paired sample directories\n";
      die "Cannot copy the sumstat files for processed paired reads between sample directories: $!\n" if system ("cp $sstatONE $sstatTWO");

      ############################################
      ### Add a /1 or /2 index to paired files ###
      ############################################

      print STDERR "\nAdding indexes (/1 and /2) to paired reads\n";
      print STDERR "outputONE:$outputONE\n\noutputTWO:$outputTWO\n\n" if $debug;
      print STDERR "baseDir: $baseDir\nbaseDir2: $baseDir2\n" if $debug;

      # 160217:
      # At this stage insert the command to add the read index (/1 or /2) to the read IDs within each file 
      # $reformat_ids_output{$outputONE} = 
      # Index to be appended can be derived from the sample ID of each paired sample
      # $baseDir and $baseDir2 will contain the indexes. Order should not be assumed.
      # $procDir corresponds to $baseDir. $procDir2 corresponds to $baseDir2. So $outputONE is paired to $baseDir and $outputTWO corresponds with $baseDir2.
      # $formatChoice defines fasta or fastq.
      # Need to define the output file names for the renamed FASTQ files
      # Call to add_read_index.pl rename reads to meet paired convention 

      # baseDir and baseDir2 specified by the imp_commandline script - pairing happens upstream.
      # Is pairing random or ordered?
      # How is _1 and _2 specified. Is this ordered at the organise step and then random?
      # _1 and _2 is specified by the order of the files in the description file, but seperated to own lines in
      # metatable for trimming with _1 and _2 appended. Paired by pair column (see pairify).
      # Although files are paired, metatable is a hash, so order when passed to this script is not strictly controlled.
      # However, baseDirs will have _1 and _2 suffix recording file order in description file.

      # my $reapFfn = "$baseDir/$reapData";
      # $inputsONE[0]
      # $procDir
      # $baseNameONE
      # $reformat_ids_output{$outputONE} = 
     
      # $reapFfn2 = "$baseDir2/$reapData"
      # $inputsTWO[0]
      # $procDir2
      # $baseNameTWO
      # $reformat_ids_output{$outputTWO} = 

      # Retrieve the index from the corresponding baseDir
      
      sub retrieve_file_index {
         my $base_directory_name = shift;
         my $index_found = 0;
         
         if($base_directory_name =~ /_(\d+)[\/]*$/){
            $index_found = $1;
         }else{
            die "Unable to retrieve the file index from the base directory name\n";
         }
         
         die "Index must be 1 or 2\n" if (($index_found != 1)&&($index_found != 2));
         return($index_found);
      } 
    
      my $index_fastq_format = $fastq ? "fastq" : "fasta" ;

      # Pass the index, tallied file, output file and format to the add_read_index.pl script
      # add_read_index.pl --reads= --format= --index= --out=
      
      # First file
      my $first_index = &retrieve_file_index($baseDir);
      my $first_index_call = "$addIndexFfn --reads=$intermedONE --format=$index_fastq_format --index=$first_index --out=$outputONE";

      # Second file
      my $second_index = &retrieve_file_index($baseDir2);
      my $second_index_call = "$addIndexFfn --reads=$intermedTWO --format=$index_fastq_format --index=$second_index --out=$outputTWO";
      
      print STDERR "\nIndexing calls:\n$first_index_call\n$second_index_call\n" if $debug;
     
      die "Could not add the indexes to the first tallied file ($first_index_call)\n" if system($first_index_call);
      die "Could not add the indexes to the second tallied file ($second_index_call)\n" if system($second_index_call);

      die "The first reindexed tally file seems to be missing or empty. Please check ($outputONE).\n" if (! (-s $outputONE));
      die "The second reindexed tally file seems to be missing or empty. Please check ($outputTWO).\n" if (! (-s $outputTWO));

      # QC plot summarising the repairing of files

      my $PairedTQC_Call = "R --vanilla --slave --args --output=$qcPlotONE --left=$baseNameONE --right=$baseNameTWO --input=$sstatONE < $tally_sstat_Ffn";
      die "Not able to plot Tally sstat file: $PairedTQC_Call\n" if system($PairedTQC_Call);
      
      die "Cannot copy the sumstat QC plot files for processed paired reads between sample directories: $!\n" if system ("cp $qcPlotONE $qcPlotTWO");

      print STDERR "\n### COMPLETED: The paired sequence files have been filtered and collapsed to unique reads in parallel. ###\n\n";

      exit(0);
   
   }
}









######################
### STAGE THREE ######
######################

# boxBowtie
#              Will require an annotation base directory containing directories housing species-version/ directories for each Ensembl annotation
#              Will need to organise Bowtie genome selection such that coordination between genome and annotation is ensured

if ($stage eq "align"){                                    
   

   $versionDir = "$annotDir/ENSEMBL/$genome"."_$ensversion/";
   if(! -d $versionDir){
      die "Cannot find the directory for the annotation required for the genome and version. Please ensure the annotation has assembled as expected\n";
   }

# 150714: Annotation method change over.
#   my $APIDir = "$versionDir/API_extracted/";      
#   die "Cannot find the API_extracted directory for the annotation required for the genome and version requested. Please ensure it has been downloaded appropriately\n" unless -d $APIDir;
#
#   my @annotationFiles  = <$APIDir*>;   
#   print STDERR "Identifying Ensembl API extracted annotation files corresponding to desired annotation\n";
#   print STDERR "@annotationFiles\n" if $debug;
#
#   my @selFiles = grep /$genome-$ensversion.sirocco.chrom_len.tab.gz/, @annotationFiles; 
#   die "Too many chromosome length files found in $APIDir\n" if scalar @selFiles > 1;
#   die "Could not find the expected chromosome length file in $APIDir\n" unless scalar @selFiles == 1;
#   my $chrLenFile = $selFiles[0];
#   print STDERR "Selected chromosome length file: $chrLenFile\n" if $debug;
   
   my $genomeDir = "$versionDir/Genome/";      
   die "Cannot find the Genome directory for the annotation required for the genome and version requested ($genomeDir).\nPlease ensure the annotation has assembled as expected\n" unless -d $genomeDir;

   my @annotationFiles  = <$genomeDir*>;   
   print STDERR "Identifying Ensembl extracted annotation files corresponding to desired annotation\n";
   print STDERR "@annotationFiles\n" if $debug;
   
   ### Identifying desired chromosome length file
   
   my $chrLenFile = "$genomeDir/$genome"."_$ensversion.chrom_len.tab";                                        
   die "Cannot find the expected chromosome length file: $chrLenFile, please take a look.\n" unless -s $chrLenFile;
   print STDERR "Selected chromosome length file: $chrLenFile\n" if $debug;
 
   ### Identifying desired chromosome class file

   my $ChrTypeFile = "$genomeDir/$genome"."_$ensversion.chrom_class.tab";                                        
   die "Cannot find the expected chromosome class file: $ChrTypeFile, please take a look.\n" unless -s $ChrTypeFile;
   print STDERR "Selected chromosome class file: $ChrTypeFile\n" if $debug;

   ### Identifying the GRanges and Bowtie directories for the desired genomes - Bowtie index and GRanges RData Annotation files

   my $GRangesDir = "$versionDir/GRangesObjects/"; 
   die "Cannot find the GRanges directory for the annotation required for the genome and version requested. Please ensure it has been downloaded appropriately\n" unless -d $GRangesDir;
   
   my $BowtieIndDir = "$versionDir/BowIndex/";
   die "Cannot find the BowIndex directory for the annotation required for the genome and version requested. Please ensure it has been downloaded appropriately\n" unless -d $BowtieIndDir;

   my $bowtieBase;
   my @Bowtie_result = glob("$BowtieIndDir/*.rev.1.ebwt");
   die "expected only one file" unless @Bowtie_result == 1;
   if ($Bowtie_result[0] =~ /($BowtieIndDir.*)\.rev\.1\.ebwt/){
       $bowtieBase = $1;
   }

   my $GRangesClass = "$GRangesDir/$genome"."_$ensversion.annotation.classes.RData";
   die "Cannot locate the file describing the Repeat and Gene annotation sets : $GRangesClass, please check\n" unless -s $GRangesClass;

   my $reducedGRs = "$GRangesDir/$genome"."_$ensversion.reduced.GRanges.RData"; 
   die "Cannot locate the reduced GRanges annotation data for the requested genome and version, please check\n" unless -s $reducedGRs;

   ### Organising directory structure for alignment process - Checking and creating directories.
   
   my $procDir = "$baseDir/PROCESSED/"  ; 

   if (! -d $procDir){
      die "The PROCESSED directory is not present within $baseDir, please investigate";
   }

   my $bowtieDir = "$baseDir/BOWTIE/";
   
   print STDERR "Making a directory to store Bowtie results\n"; 
   print STDERR "$bowtieDir\n" if $debug;

   if (-e $bowtieDir){
      print STDERR "WARNING: The BOWTIE directory already exists ($bowtieDir). Will use current directory. This may cause files to be overwritten.\n"; # The imp will no longer die if rerunning/overwiring - 041111
   }
   die "Cannot make BOWTIE directory $bowtieDir: $!"if (system("mkdir -p -m 0755 $bowtieDir"));

#      die "BOWTIE Directory exists at $baseDir, please tidy up";
#   }else{   
#      mkdir $bowtieDir, 0755 or die "Cannot make BOWTIE directory: $!"; 
#   }

   my $qcDir = "$baseDir/QC/"  ; 

   if (! -d $qcDir){
      die " The QC directory is not present within $baseDir, please investigate\n";
   }
   
   ### Creating a Bowtie log

   $bowtieLog = "$bowtieDir/bowtie_log.txt";

   ####################################################################################
   ### Identify base name based upon basename of the trimmit processed output files ###
   ####################################################################################

   my $experimentBow = "";
   my %uniqueIDBow = ();
   my @exptsBow = ();

   print "DIRECTORY: $procDir\n" if $debug;
   my @procfiles = <$procDir/*>;
   print "PROCESSED Files pre trim:",@procfiles,"\n" if $debug;
   foreach (@procfiles){
      s/$procDir\///;
   }

   print "PROCESSED Files:",@procfiles,"\n" if $debug;

   my @trimmedFiles = grep /\.clean\.processed\.fa\.gz$/, @procfiles;

   print "FASTA Files:",@trimmedFiles,"\n" if $debug;

   for my $trimNow (@trimmedFiles){
      if ($trimNow =~ /^([\w_]+)\.\w+\.clean\.processed\.fa\.gz$/){
         $uniqueIDBow{$1} = 1;  
      }else{
         die "Can not identify the unique ID for trimmed file $trimNow\n";
      }
   } 

   print "Experiments:",(keys %uniqueIDBow),"\n\n" if $debug;

   if ((scalar(keys %uniqueIDBow)) != 1){
      die "There aren't the expected number of experimental IDs in the $procDir directory\n";
   }else{
      print STDERR "The correct number of Experimental IDs have been identified within the processed read directory\n" if $debug;
      @exptsBow = keys %uniqueIDBow; 
      $experimentBow = $exptsBow[0];
   }
#     exit;  
   #####################################
   ### CALLING BOWTIE ON A DIRECTORY ###
   #####################################

   my $bowtieDebug = $debug ? "--debug" : "" ;

   ### BOWTIE OPTIONS
   my $chunkOpt = $chunkmbs != 0 ? "--chunkmbs=$chunkmbs" : "";
   my $mismatchOpt = $misMatches != 0 ? "--mismatches=$misMatches" :"";
#   my $stratumOpt = "--strata" if $stratum != 0;                           # Now removed, see above.
   my $stratumOpt = "--strata";
#   my $hitNoOpt = "--number=$hitNo" if $hitNo !=0;                         # Hardcoded to ensure stratum always "ON" to avoid confusion and errors.
   my $hitNoOpt = "--number=$hitNo";
   my $samOpt = $bowtie_sam ? "--sam" : "";                                 # Run Bowtie steps in SAM format?
   my $expandOpt = $bowtie_sam ? "--opts=$depthFieldFfn": "";               # Required to add optional fields to SAM format for read depth etc.

   my $shipBowtieCall = "$shipBowtieFfn --in=$procDir --out=$bowtieDir --log=$bowtieLog --genome=$bowtieBase $chunkOpt $mismatchOpt $stratumOpt $hitNoOpt $samOpt $expandOpt $bowtieDebug";  


   ### BOXBOWTIE STEP 1 - SUMMARISING BOWTIE RESULTS.
#   my $maxHitsOpt = "--maxHits=$hitNo" if $hitNo !=0;                      # Now always 1 or more
   my $maxHitsOpt = "--maxHits=$hitNo";
   my $procTotFile = "$procDir/total.saved.trimmit.tab";
   die "Cannot locate the summary of processed reaads: $procTotFile" unless -e $procTotFile && -s $procTotFile; 
   my $procSumOpt = "--procSummary=$procTotFile";
   my $R_boxSamOpt = $bowtie_sam ? "--sam=29" : "";
   my $R_bowtieDebug = $debug ? "--debug=29" : "" ;

   my $boxBowtieONECall = "R --vanilla --slave --args --input=$bowtieDir --output=$bowtieDir --chrLen=$chrLenFile $procSumOpt $maxHitsOpt $R_boxSamOpt $R_bowtieDebug < $boxBowtieONEFfn";


   ### BOXBOWTIE STEP 2 - MAPPING BOWTIE RESULTS TO NON-REDUNDANT ANNOTATION 

   my $boxBowtieTWOCall = "R --vanilla --slave --args --input=$bowtieDir --output=$bowtieDir --annot=$reducedGRs $R_bowtieDebug < $boxBowtieTWOFfn";

   ### BOXBOWTIE QC - PLOTTING QC SUMMARY FOR BOWTIE STEP

#   my $qcMaxHitOpt = "--maxHits=$hitNo" if $hitNo !=0;                     # Now always 1 or more
   my $qcMaxHitOpt = "--maxHits=$hitNo";

   my $boxBowtieQCCall = "R --vanilla --slave --args --input=$bowtieDir --output=$qcDir --basicID=$experimentBow --grangesClass=$GRangesClass --chrTypes=$ChrTypeFile $procSumOpt $qcMaxHitOpt $R_bowtieDebug < $boxBowtieQCFFn";

   ### Coordinating systems  calls:

   print STDERR "Lined up shipBowtie call:\n$shipBowtieCall\n" if $debug;
   print STDERR "Lined up boxBowtie Step 1 call:\n$boxBowtieONECall\n" if $debug;
   print STDERR "Lined up boxBowtie Step 2 call:\n$boxBowtieTWOCall\n" if $debug;
   print STDERR "Lined up boxBowtie QC call:\n$boxBowtieQCCall\n" if $debug;

   print STDERR "\nSubmitting processed read files to Bowtie\n";

   die "shipBowtie hasn't worked: $shipBowtieCall: $! \n" unless system($shipBowtieCall)==0 ;
 
   print STDERR "\nCompleted Bowtie submissions\n\nPreparing to check Bowtie output\n";

   bowtie_file_check($bowtieDir, $procDir, $bowtie_sam);  # Check that all expected Bowtie output is present
   
   print STDERR "Checking complete\n\nBeginning GRanges conversion\n\n";

   die "GRanges conversion (boxBowtie_Step1) hasn't worked: $boxBowtieONECall: $! \n" unless system($boxBowtieONECall)==0;
   
   print STDERR "\nConversion complete and Bowtie results have been summarised to tables\n\nPreparing to compare with reference annotation\n";

   die "Comparisons to annotation have failed (boxBowtie_Step2): $boxBowtieTWOCall: $! \n" unless system($boxBowtieTWOCall)==0;

   print STDERR "Mapped the sequences to non-redundant annotation\n\nPreparing to plot QC results\n";

   die "QC plotting (boxBowtie_QC) hasn't worked: $boxBowtieQCCall: $! \n" unless system($boxBowtieQCCall)==0;

   print STDERR "\n### COMPLETED: Aligned all of the reads to the desired specifications. Please check the QC directory for summary plots ###\n\n";
    
   exit(0); 
}

#############################
### STAGE FOUR - FEATURES ###
#############################
            # Stage four options - Features   -> define checks dependent upon --features flag 
            # Will require genome - Need a config file for species - Relevant species should be determined at this level and required name passed to keep it central
            # Will require feature type - eg. miRNA
            # Need a better way to control the genome nomenclature
            # Will need to interpret, specify and check input directory
            # Will need to create, check and specify output directory
            # Include option to pass it own annotation files (correctly formatted)?
            # Will also require a version number. Will have to ensure mapping and annotation are version compatible?
            # Check annotation file - who's responsibility to ensure genomes match?
            # What directory level will annotation be specified at?

if ($stage eq "features"){
   my $featuresDebug = $debug ? "--debug" : "" ;
   my $R_featuresDebug = $debug ? "--debug=29" : "" ;

   if ($fType eq "miRNA"){

      my $annotationRData ;
      my $miRTableCall;   
      my $analysisDir;
      my $mapDir;

      my $R_methodOption = $no_merge ? "" : "--collapsetype=$collapse_method" ;

      my $R_mergingFlag = $no_merge ? "--nomerge=29" : "";
      my $R_propFlag = $proportional ? "--proportional=29" : "";

      # Set variables in a way that is feature dependent (eg. Select correct annotation file)
      
      $mapDir= "$baseDir/BOWTIE/";
   
      if (! -d $mapDir){
         die "Cannot find the BOWTIE/ directory in $baseDir, please check\n";
      }

      $annotationRData = "$annotDir/MIRBASE/$genome"."_$mirversion/$genome"."_$mirversion.miRNA_coordinates.RData";
      if(! -s $annotationRData){
         die "Cannot find the annotation file required genome and version ($annotationRData), please check.\n";
      }
      $analysisDir = "$baseDir/miRNA_ANALYSIS/";
     
      if (-e $analysisDir){
         print STDERR "WARNING: The miRNA_ANALYSIS directory already exists ($analysisDir). Will use current directory. This may cause files to be overwritten.\n"; # The imp will no longer die if rerunning/overwiring - 041111
      }

      ### HERE - Add collapse method option if nomerge FALSE
      die "Cannot make the miRNA analysis directory $analysisDir: $!\n" if (system("mkdir -p -m 0755 $analysisDir"));
   
      $miRTableCall = "R --slave --vanilla --args --inDir=$mapDir --outDir=$analysisDir --overlap=$overlap --annot=$annotationRData --treatOverlap=$annot_conflict $R_featuresDebug $R_propFlag $R_mergingFlag $R_methodOption < $miRTableFfn";
      
      print STDERR "\nPreparing to cross reference Bowtie results with miRBase miRNA annotation\n\n";
      print "Call to miR_table.R: $miRTableCall\n\n" if $debug;
      die "Mapping to miRBase annotation (miR_table.R) hasn't worked:$miRTableCall : $! \n" unless system($miRTableCall)==0;
      print STDERR "\n### COMPLETED: Successfully mapped reads to $fType features ###\n\n";
      exit;

   }elsif ($fType eq "repeat"){
      my %baseIDProc;
      my %repConsensi;
      my $bowRef;

      #### FILES AND DIRECTORIES

      my $preprocDir = "$baseDir/PROCESSED/";
      my $repeatAnnot = "$annotDir/REPEATS/$genome"."_$repversion/BowIndex/";
      my $repAnDir = "$baseDir/rep_ANALYSIS/";
      my $mappedOut = "$repAnDir/genome_mapped.txt";
      
      if (-e $repAnDir){
         print STDERR "WARNING: The rep_ANALYSIS directory already exists ($repAnDir). Will use current directory. This may cause files to be overwritten.\n"; # The imp will no longer die if rerunning/overwiring - 041111
      }

      die "Cannot make Repeat Analysis Directory $repAnDir: $!\n" if (system("mkdir -p -m 0755 $repAnDir"));
      
      my $QCrepDir = "$baseDir/QC/";

      die "Could not find a QC/ directory in $baseDir. Required for reporting results.\n" unless -d $QCrepDir;
      die "Could not identify the repeat BowIndex directory required for the species specified ($repeatAnnot). Please ensure the repeat data is available for the genome and version specified.\n" unless -d $repeatAnnot;

      #### GATHERING INPUT FILES - AUTOMATED: Will count how many consensi are available and will then run on all of them.

      print STDERR "Gathering input files\n";

      my @preprocdata = <$preprocDir/*.clean.processed.fa.gz>;
      my @repAvail = <$repeatAnnot/*rev.1.ebwt>;

      print "Preprocessed data files:",@preprocdata,"\n" if $debug;
      print "Available repeat/NCBI sequences:",@repAvail,"\n" if $debug;
      die "No repeat/NCBI sequences were found against which to perform repeat analysis\n" if scalar(@repAvail) == 0;

      #### CALCULATE THE BASE NAME FOR THE EXPERIMENT

      for my $preprocNow (@preprocdata){
         if ($preprocNow =~ /([\w_]+)\.\w+\.clean\.processed\.fa\.gz$/){
            $baseIDProc{$1} = 1;
         }else{
            die "Can not identify the unique ID for trimmed file $preprocNow\n";
         }
      }
      die "Could not identify a single lane associated 'base name' for the processed read files in the input directory\n"
                  if ((scalar(keys %baseIDProc)) != 1) ;
      my @generalID = keys %baseIDProc;

      #### PULL OUT ALL BASE NAMES FOR REPEAT BOWTIES

      for my $repeatNow (@repAvail){
         if ($repeatNow =~ /(.+)\.rev\.1\.ebwt$/){
            $repConsensi{$1} = 1;
         }else{
            die "Can not find the repeat name associated with $repeatNow\n";
         }
      }
      print STDERR "Identified ".scalar(keys(%repConsensi))." repeats to which the processed files will be compared\n";
   
      my $genAnnotDir = "$annotDir/ENSEMBL/$genome"."_$ensversion/BowIndex/";
      my $genPat = "$genAnnotDir*.rev.1.ebwt";

      # Find the genomic base name for mapping
      
      print STDERR "Identifying genomic reference sequence for generating scaling factor\n\n";
      
      my @indFileList = glob("$genPat");
      die "expected only one file\n" unless @indFileList == 1;
      print "Genomic Bowtie index:".$indFileList[0]."\n" if $debug;

      if ($indFileList[0] =~ /($genAnnotDir.*)\.rev\.1\.ebwt/){
         $bowRef = $1;
      }
      print "$bowRef\n" if $debug;

      # Open a file to record total mappings
      open(MAPPED, "> $mappedOut") || die "Could not open the file for recording genome mapping data\n";
      print MAPPED "Sample\tReads\n";

      print STDERR "Preparing to cycle through processed read files to identify reads that map to the genome\n";
      my $fileCount =1;
      
      for my $newInput (@preprocdata){
         my $mapTotal;
         my $normOut = "";
         my $sampleBass = "";
         my $fileBar = "";
         
         print STDERR "\nMapping file $fileCount\n";

         # Record the basename of the processed read file for naming output
         # Record barcode for summary file

         if ($newInput =~ /^$preprocDir\/(\S+)\.clean\.processed\.fa\.gz$/){
            $sampleBass = $1;
            $sampleBass =~ /\.(\w+)$/ && ($fileBar = $1);
            $normOut = "$repAnDir/$sampleBass.normalise.bowtie.temp.ebwt.gz";
         }

         # Bowtie call
         my $repMismatchOpt = $repMisMatches != 0 ? "-v $repMisMatches" : "";
         my $repChunkOpt = $repChunkmbs != 0 ? "--chunkmbs $repChunkmbs" : "";

         # Report only the first hit for each read. Only need to know if it maps somewhere for read count to be incorporated into the total
         my $normCall = "gzip -dcf $newInput | bowtie --time $repMismatchOpt $repChunkOpt --suppress 2,3,4,5,6,7,8 -f $bowRef - | gzip -c1 > $normOut";
         die "Genomic Bowtie hasn't worked\n" unless system($normCall)==0;
         print "$normCall\n" if $debug;

         # Open the mapped read file and count the totals
         open(GENBOW,"gzip -cd $normOut | cat |") || die "Could not open the genome aligned read file\n";
         while(<GENBOW>){
            chomp;
            /_x(\d+)/ && ($mapTotal += $1);
         }
         close GENBOW;

         # Record total read matches for each barcode
         print MAPPED "$fileBar\t$mapTotal\n";
         $fileCount++
      }
      close MAPPED;
      print STDERR "\nGenomic mapping complete\n";

      ### Check the summary file is present for downstream analysis
      die "Could not find genome mapping summary file\n" unless -e $mappedOut && -s $mappedOut;

      print STDERR "\nPreparing repeat analysis for all identified repeat sequences\n";
      print STDERR "Found ",scalar(keys(%repConsensi))," repeats. Will conduct mapping cycle for each\n";

      for my $currentRep (keys %repConsensi){

         ### Pull out repeat name and create a repeat specific file
         my $repeatType = "";
         $currentRep =~ /([^\/\.]+)$/ && ($repeatType = $1);
         
         print STDERR "\n### Beginning analysis for $repeatType ###\n";

         my $repeatOutDir = "$repAnDir/$repeatType/";
      
         if (-e $repeatOutDir){
            print STDERR "WARNING: A repeat specific directory already exists ($repeatOutDir). Will use current directory. This may cause files to be overwritten.\n"; # The imp will no longer die if rerunning/overwiring - 041111
         }

         die "Cannot make Repeat specific directory for $repeatType ($repeatOutDir): $!" if (system("mkdir -p -m 0755 $repeatOutDir"));

         ### Make a call to piRNA_bowtie.pl - bowtie against the canonical sequence
         print STDERR "\nMapping processed reads to canonical repeat sequence with Bowtie\n";
         my $piChunkOpt = $repChunkmbs != 0 ? "--chunkmbs=$repChunkmbs" : "";

         ###### Bowtie call passed to subsiduary script 
         ###### Pass file name for recording the alignment number per sample
        
         my $recordFile = "$repeatOutDir/$repeatType.repeat_alignments.txt";
         my $piBowCall = "$piBowtieFfn --in=$preprocDir --out=$repeatOutDir --maxHits=$repHitNo $piChunkOpt --mismatch=$repMisMatches --depthField=$depthFieldFfn --repeat=$currentRep --record=$recordFile $featuresDebug";
         print STDERR "Bowtie against repeat consensi: $piBowCall\n" if $debug;
         die "Bowtie against canonical repeat hasn't worked for $currentRep\n" unless system($piBowCall)==0;
         print STDERR "Mapping complete\n";

         # Skip this section if no alignments are found for any samples.
         
         open (ALIGN_REC, "< $recordFile") || die "Could not open the repeat alignment record: $recordFile\n";
         my $discarded_header = <ALIGN_REC>;
         my %aligned_counts;
         while (<ALIGN_REC>){
            chomp;
            my @line = split "\t", $_; 
            $aligned_counts{$line[0]} = $line[1];
         }
         close (ALIGN_REC) || die "Could not close the alignment record file: $recordFile\n";

         ### Only proceed some alignments are available
         my $no_alignments = 0;
         my $total_alignments = 0;

         foreach(keys %aligned_counts){
            $no_alignments++ if ($aligned_counts{$_}==0);
            $total_alignments+=$aligned_counts{$_};
         }

         print STDERR "$no_alignments read set had no alignments to $repeatType\n" if $no_alignments > 0;

         if($total_alignments > 0){

            ### boxBowtie call for creating GRanges object
            print STDERR "\nConverting mapped reads to GRanges objects\n";
            my $prTot = "$preprocDir/total.saved.trimmit.tab";
            die "Cannot locate the summary of processed reads: $prTot" unless -e $prTot && -s $prTot;
            print "Processed Summary File: $prTot\n" if $debug;

            my $chrlFile = "$currentRep.length.txt";
            die "Cannot locate the chromosome length file: $chrlFile" unless -e $chrlFile && -s $chrlFile;
            print "Chromosome length file: $chrlFile\n" if $debug;

            # This script will expect a range of output from $piBowCall - bases processing of samples on .bowtie.unique.output.sort.bam files present.
            # Therefore any samples with no alignments and hence no bam files from previous step will be lost.
            my $bbCall = "R --vanilla --slave --args --step=repeat --input=$repeatOutDir --output=$repeatOutDir --maxHits=$repHitNo --chrLen=$chrlFile --procSummary=$prTot --sam=29 $R_featuresDebug < $boxBowtieONEFfn";
            print STDERR "boxBowtie_Step1.R call: $bbCall\n" if $debug;
            die "boxBowtie_Step1.R has failed for $currentRep\n" unless system($bbCall)==0;
            
            print STDERR "GRanges conversion complete\n";   

            ### piRNA_coverage call for calculating a coverage object
            # Requires .coverage.RData, .pingpong.tab and .length_distrib.tab files to proceed. Must be present even if one barcode fails - looping based on .bowtie.conv.GR.RData files present.
            # Samples with no alignment would therefore be lost from output files.
            print STDERR "\nCalculating repeat coverage and read overlaps\n";
            my $piCovCall = "R --vanilla --slave --args --inputDir=$repeatOutDir --outputDir=$repeatOutDir --mapDepth=$mappedOut --repeatName=$repeatType $R_featuresDebug < $piCovFfn";
            print "$piCovCall\n" if $debug;
            print STDERR "Repeat coverage call: $piCovCall\n" if $debug;
            die "Could not calculate coverage for $currentRep\n" unless system($piCovCall)==0;
            print STDERR "Completed data collection\n";
         
         }else{
            print STDERR "Reads from this set of files do not align to $repeatType\n";
         }
      
         ### plot piRNA report QC
         # Ideally a record of the missing samples will be included in the first plot (0 counts) - Pass missing barcode info
         # If no samples have any reads then plot "No samples possessed alignments" in the file.
         print STDERR "\nPlotting repeat QC information\n";
         my $repPlotCall = "R --vanilla --slave --args --inputDir=$repeatOutDir --output=$QCrepDir --input=$repeatOutDir --genMap=$mappedOut --repName=$repeatType --basicID=$generalID[0] --alignRec=$recordFile $R_featuresDebug < $piPlotFfn";
         print "$repPlotCall\n" if $debug;
         print STDERR "Repeat coverage plots: $repPlotCall\n" if $debug;
         die "Could not plot repeat QC plots\n" unless system($repPlotCall)==0;
         print STDERR "QC plots complete\n";
      }
      print STDERR "\n### COMPLETED: Successfully mapped reads to $fType features ###\n\n";
      exit;
   }else{
      die "Do not recognise the specified --feature\n";
   }
}












#####################################
### Insert -annotation stage here ###
#####################################

# Will operate on ./sequence_imp.pl --annotation flag
# Will call imp_food.pl and download all of the required annotation from a particular species to a specified API_extracted/ base directory.
# Include the download and writing of Bowtie indexes.


###################
### SUBROUTINES ###
###################

# Checks that all of the expected Bowtie files have been produced.

sub bowtie_file_check {

   print STDERR "Beginning Bowtie output check to ensure all expected files are present\n";
   my ($bowStuff ,$procStuff, $samStuff) = @_ ;
   my @procObj;
   my @bowtObj;
   my @samObj;
   my @bamObj;
   my @tmpObj;
   my @indexObj;
   my @sortObj;
   my @barThingies;

   print STDERR "Processed file directory: $procStuff\n" if $debug;

   @procObj = <$procStuff/*\.clean\.processed\.fa\.gz>;
  
   print STDERR "Processed file directory includes: @procObj\n" if $debug;

   for my $file (@procObj){
      if ($file =~ /^$procStuff\/(\S+)\.clean\.processed\.fa\.gz$/){
         push @barThingies, $1;         
      }
   } 
   print STDERR "Looking for the following references in the Bowtie output files: @barThingies\n" if $debug;
   if ($samStuff){
      @tmpObj = <$bowStuff/*\.bowtie\.temp\.output\.sam\.gz>;
      @samObj = <$bowStuff/*\.bowtie\.unique\.output\.sam\.gz>;      
      @bamObj = <$bowStuff/*\.bowtie\.unique\.output\.bam>;
      @sortObj = <$bowStuff/*\.bowtie\.unique\.output\.sort\.bam>;
      @indexObj = <$bowStuff/*\.bowtie\.unique\.output\.sort\.bam\.bai>;
      for my $bar (@barThingies){
         die "Cannot find file $bar.bowtie.temp.output.sam.gz in $bowStuff" if ! grep /$bar/, @tmpObj;
         die "Cannot find file $bar.bowtie.unique.output.sam.gz in $bowStuff" if ! grep /$bar/, @samObj;
         die "Cannot find file $bar.bowtie.unique.output.bam in $bowStuff" if ! grep /$bar/, @bamObj;
         die "Cannot find file $bar.bowtie.unique.output.sort.bam in $bowStuff" if ! grep /$bar/, @sortObj;
         die "Cannot find file $bar.bowtie.unique.output.sort.bam.bai in $bowStuff" if ! grep /$bar/, @indexObj;
      }
   }else{
      @bowtObj = <$bowStuff/*\.bowtie\.output\.bwt\.gz> ;
      for my $bar (@barThingies){
         die "Cannot find file $bar.bowtie.output.bwt.gz in $bowStuff" if ! grep /$bar/, @bowtObj; 
      }   
   }
   #exit;
}



# Checks that the Reaper output is as expected

sub reaper_file_check {
   
   my ($isBar, $whatBar, $prefix, $barTable ,$reaperOut, $debuggle) = @_ ; 
   
   my $numberCycle = 0;
   my $column;
   my $matches;
   my @correctFiles;
   my @barcodes;
   
   if ($isBar) {

      if (! open BARS, $barTable){
          die "Oops, can't open the sample file for barcodes: try checking permissions";
      }
      chomp( my $firstline = <BARS>);
  
      #print "$firstline\n";
      
      my @firstBit = split "\t", $firstline;

      foreach (@firstBit){
         if (/$whatBar/){
            $column = $numberCycle;                   # Records which columns match the barcode header
            $numberCycle++;                           # Counts which column is being checked
            $matches++;                               # Records the number of matches to the barcode header
         }else{
            $numberCycle++;
         }
      }
   
   #print "Matches: $matches\n";
      die "Too many specified barcode columns in the sample metadata file" if ($matches != 1);

      ### Collect the barcodes in the correct columns

      while (<BARS>){
         chomp;
         my @arrs = split /\t/, $_, -1;
         push @barcodes, $arrs[$column];
      }

      close(BARS);

      print STDERR "Recognised barcodes for sample: @barcodes\n";

   ### Collect together a list of expected files in Reaper output

      push @correctFiles, "$prefix.total.report.input.nt";
      push @correctFiles, "$prefix.total.report.input.q";
      push @correctFiles, "$prefix.total.report.miss.nt";

      foreach my $code (@barcodes){
         push @correctFiles, "$prefix.$code.clean.gz";
         push @correctFiles, "$prefix.$code.report.clean.nt";
         push @correctFiles, "$prefix.$code.report.input.nt";
         push @correctFiles, "$prefix.$code.report.input.q";
   #   print"\n@correctFiles\n";
      }
   }elsif (!$isBar){
      print STDERR "There are no barcodes used within this lane, so screen for barcoded files skipped\n";
      push @correctFiles, "$prefix.lane.clean.gz";
      push @correctFiles, "$prefix.lane.report.clean.nt";
      push @correctFiles, "$prefix.lane.report.input.nt";
      push @correctFiles, "$prefix.lane.report.input.q";      
   }

   ### Ensure all the expected files are available in the directory

   @files = <$reaperOut/$prefix*>;
   foreach (@files){
      s/$reaperOut\///;
   }
   print STDERR "Expected files for downstream applications: @correctFiles\n" if $debuggle;
   print STDERR "Searching folder: $reaperOut\n" if $debuggle;
   print STDERR "REAPER directory contains: @files\n" if $debuggle;

   foreach my $oops (@correctFiles) {
   #     print "\n$oops\n";
         my $iggle = grep /$oops$/, @files;
         if($iggle == 0){
            die "The $oops file is missing. Please check your Reaper output...\n";
         }elsif ($iggle >= 2){
            die "There is a file name conflict amongst the Reaper output, file check failed\n";
         }else{
            print STDERR "$oops is present in the $reaperOut folder\n" if $debuggle;
         }
   #     print "$iggle\n";
   }


   ### Check the .tally file is present

   $tallyHo = grep /sumstat$/, @files;

   if ($tallyHo == 0){
            die "The Reaper 'sumstat' file is absent, please check the Reaper output\n";
   }elsif ($tallyHo > 1)  {
            die "There are too many Reaper 'sumstat' files, please check the Reaper output\n";
   }else{
            print STDERR "The correct number of sumstat files have been found\n" if $debuggle;
   }

   my $success = "All Reaper output files appear to be present\n";
   
   return $success;

}





