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

### This script will allow the initiation of batch jobs on multiple fastq files from the commandline.
### It includes organising the files into a directory structure, the movement of data files and the 
### creation of metadata files based on simple, tab delineated inputs.
### 'config' contains options for imp_commandline and sequence_imp. These are generally handled by
### parsing the config file and stringing the options together for each step. However, options required
### for multiple steps are handled seperately. In addition flags must be handled alone.


use Getopt::Long;
use strict;
use warnings;

BEGIN {
   my $path = $0;
   $path =~ s|[^/]*$||;
   if ($path =~ /\/?bin\/$/){
      $path =~ s/\/?bin\/$//;
   }else{
      $path =~ s/\/$//;
   }
   #print "$path\n";
   $ENV{SEQIMP_ROOT} = $path;
}

# Test to see if perl has been compiled to be compatible with threads.

my $threads_available=1;

if (!defined(eval qq{
   require threads;
   require Thread::Queue;
})) {
   print STDERR "\nWARNING: Perl module 'threads' is not available - disabling multiprocessor option\n";
   $threads_available = 0;
}

#exit(0);
print STDERR "\n###############################################\n### WELCOME TO THE SEQUENCE IMP COMMANDLINE ###\n###############################################\n";

print <<EOH;

Sequence Imp - Copyright (C) 2011 2012 2013 EMBL - European Bioinformatics Institute 
This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

Please contact us at: kraken\@ebi.ac.uk

For more info (including full manual) see: http://www.ebi.ac.uk/enright/kraken/

EOH

my $progname = 'imp_commandline.pl';
#my $seqimp = '/nfs/research2/enright/mat/Small_Stuff/miRNA_Seq/SIROCCO/pipeline/seqimp/bin/sequence_imp.pl';
my $seqimp = "$ENV{SEQIMP_ROOT}/bin/sequence_imp.pl";             # Ffn of sequence_imp.pl
my $defDir = "$ENV{SEQIMP_ROOT}/default_configuration_files/";    # Ffn of the default parametre files
my $sysChck = "$ENV{SEQIMP_ROOT}/bin/system_check.sh";           # Ffn of the system_check.sh script


###################
### Subroutines ###
###################

# Usage
sub example_comm {
  print <<EOH; 

Examples of usage:

* Recommendation:
Before using the pipeline for the first time we recommend running a quick check
to ensure all of the required software is installed and accessible on your
system.

seqimp/bin/imp_commandline.pl --system-check
_______________________________________________________________________________

* Initial pipeline step:
This should be run in the directory containing relevant data files and will
produce a new 'analysis_#' directory arranging the data to be used by
subsequent steps. Further options are available. Please note that the
description file must be supplied by the user. For details see documentation.

seqimp/bin/imp_commandline.pl --step=organise --description=description.txt --default-configuration=human_miRNA_analysis.txt
_______________________________________________________________________________

* Subsequent pipeline steps:
The example given will initiate the 'reaper' step of the pipeline for trimming
adapters and initial QC. For further details on available steps see
documentation. Note the name of the analysis directory will change, as a new
directory is created for each analysis, by the 'organise' step above.

seqimp/bin/imp_commandline.pl --step=reaper --description=description.txt --default-configuration=human_miRNA_analysis.txt --analysisDir=analysis_22077

Note the 'align' and 'features' steps require a single additional argument
specifying a directory containing the annotation information for the pipeline.
Annotation is available from the kraken website.
_______________________________________________________________________________

* Further help and documentation:

seqimp/bin/imp_commandline.pl --help

and read the documentation found within the seqimp directory.
_______________________________________________________________________________

* A list of the default configuration files supplied:

seqimp/bin/imp_commandline.pl --config-list
_______________________________________________________________________________

EOH
   exit(0);
}

# Help
sub help {
   print <<EOH;

Usage:
   $progname [options]
Options:
   ## UNIVERSAL:
   
   --help                  <FLAG>   OPTIONAL:   Displays this page
   
   --step                  <STRING> REQUIRED:   The pipeline step to initiate. 
                                                Can be 'organise','reaper',
                                                'filter','align' or 'features'. 
                                                Steps should be processed 
                                                sequentially.  
   
   --description           <STRING> REQUIRED:   Full file name of file 
                                                containing FASTQ/sample 
                                                descriptions.
   
   --user-configuration    <STRING> OPTIONAL/REQUIRED: A file containing the 
                                                settings for all of the options 
                                                required for the analysis run 
                                                (Should contain file path).
   
   --config-list           <FLAG>   OPTIONAL:   List the default configuration 
                                                files provided with the 
                                                pipeline.
   
   --default-configuration <STRING> OPTIONAL/REQUIRED: A default configuration 
                                                file containing the settings 
                                                for all of the options required 
                                                for the analysis run. (No file 
                                                path). 

   --system-check          <FLAG>   OPTIONAL:   Check that all external programs
                                                and libraries upon which the
                                                pipeline relies are installed
                                                and accessible.
   
   --paired                <FLAG>   OPTIONAL:   Specifying this flag
                                                prompts seqimp to perform a
                                                paired end analysis. Files
                                                must be provided in pairs in the
                                                description file and specific
                                                geometries must be supplied.

   ## STEP: organise
   
   --dataDir               <STRING> OPTIONAL:   The directory containing the 
                                                data files: DEFAULT: './'

   --outDir                <STRING> OPTIONAL:   The directory where you would 
                                                like analysis directories to be 
                                                built: DEFAULT: './' 

   --identifier            <STRING> OPTIONAL:   An identifier for the current 
                                                analysis: DEFAULT: 'analysis'

   --no-unique             <FLAG>   OPTIONAL:   Will not append process ID to
                                                analysis directory name

   ## STEP: reaper/filter/align/features

   --analysisDir           <STRING> REQUIRED:   The directory produced by the 
                                                'organise' step in which the 
                                                sample directories reside.

   --annotationDir         <STRING> REQUIRED:   A directory containing species 
                                                annotation data is required 
                                                for the 'align' and 'features' 
   
   --processors            <INT>    OPTIONAL:   Number of processors to utilise.
                                                DEFAULT: 1
______________________________________________________________________________

Configuration file options:

EOH

organise_help();
reaper_help();
filter_help();
align_help();
features_help();

example_comm();
exit(0);
}

sub organise_help{
   print <<EOH;
   --step=organise Configuration file options:

   reap2imp       <STRING>    OPTIONAL:   Full file name of configuration file
                                          that maps 'Sequence Imp' geometry 
                                          nomenclature to Reaper terms. 
                                          DEFAULT: Config files supplied with 
                                          the pipeline.

   geoConfig      <STRING>    OPTIONAL:   Full file name of configuration file 
                                          specifying 'description file' fields 
                                          required for each library geometry 
                                          and the corresponding Reaper 
                                          requirements, including the Reaper 
                                          config file associated with each 
                                          geometry. DEFAULT:  Config files 
                                          supplied with the pipeline.

EOH
}

sub reaper_help{
   print <<EOH;
   --step=reaper Configuration file options:
   
   fastq          <FLAG>       OPTIONAL:  Specifies FASTQ as the prefered 
                                          format for sequence files for the 
                                          pipeline. Please see the manual for 
                                          more details. DEFAULT: FASTA.
   
   plotZeros      <FLAG>       OPTIONAL:  A flag specifying whether reads of
                                          length 0 should be considered when 
                                          plotting QC plots. DEFAULT: y-axis
                                          set without considering the number 
                                          of reads of length 0.
   
   reapConfig     <STRING>     OPTIONAL:  Full path to the reaper geometry
                                          configuration directory. This 
                                          contains files corresponding to 
                                          each geometry containing Reaper 
                                          default parameters. DEFAULT: 
                                          Config files supplied with the 
                                          pipeline.

   perlUnq        <FLAG>       OPTIONAL:  Flag specifying whether perl should 
                                          be used in the place of 'Tally' to 
                                          collapse filtered FASTA files to 
                                          unique reads. DEFAULT: Will use 
                                          Tally.

EOH
}

sub filter_help{
   print <<EOH;
   --step=filter Configuration file options:
   
   fastq          <FLAG>       OPTIONAL:  Specifies FASTQ as the prefered 
                                          format for sequence files for the 
                                          pipeine. Please see the manual for 
                                          more details. DEFAULT: FASTA.

   plotZeros      <FLAG>       OPTIONAL:  A flag specifying whether reads of 
                                          length 0 should be considered when 
                                          plotting QC plots. DEFAULT: y-axis
                                          set without considering the number 
                                          of reads of length 0.

   low            <INTEGER>    OPTIONAL:  Maximum low complexity score filter.
                                          DEFAULT: No filter.

   minSize        <INTEGER>    OPTIONAL:  Minimum trimmed read length filter.
                                          DEFAULT: No filter. 

   maxSize        <INTEGER>    OPTIONAL:  Maximum trimmed read length filter.
                                          DEFAULT: No filter.

   five           <INTEGER>    OPTIONAL:  Bases to be removed from the 5' end 
                                          of all trimmed reads. DEFAULT: 0

   three          <INTEGER>    OPTIONAL:  Bases to be removed from the 3' end 
                                          of all reads. DEFAULT: 0

EOH
}

sub align_help{
   print <<EOH;
   --step=align Configuration file options:

   genome         <STRING>     REQUIRED:  The species for which the analysis 
                                          will be performed.

   ensversion     <INTEGER>    REQUIRED:  The version number specified for the
                                          compilation of the Ensembl gene 
                                          annotation to be used for run.

   chunk          <INTEGER>    OPTIONAL:  Value to be passed to the Bowtie 
                                          'chunkmbs' option. DEFAULT: Not
                                          passed to Bowtie unless set.

   mismatches     <INTEGER>    OPTIONAL:  Maximum number of mismatches to be 
                                          allowed by Bowtie. DEFAULT: 0

   maxHits        <INTEGER>    OPTIONAL:  Maximum number of hits to be reported 
                                          by Bowtie (Reads that hit more times  
                                          than this threshold are discarded)
                                          DEFAULT: 10

   sam            <FLAG>       OPTIONAL:  Specifies whether to run the pipeline 
                                          with SAM/BAM compatible alignment 
                                          format. DEFAULT: Standard Bowtie 
                                          format.

   NOTE: Bowtie always run with --strata flag set

EOH
}

sub features_help{
   print <<EOH;
   --step=features Configuration file options:

   genome         <STRING>     REQUIRED:  The species for which the analysis 
                                          will be performed.

   feature        <STRING>     REQUIRED:  A selected feature type to be 
                                          quanitated - currently supported: 
                                          "miRNA" and "repeat".
   
   feature miRNA:

   mirversion     <INTEGER>    REQUIRED/NOT REQUIRED:  miRNA annotation 
                                          version to which read alignments 
                                          will be compared. 
                                          Required only when "feature miRNA" 
                                          specified. WARNING: Ensure the ensembl 
                                          genome build and miRBase builds are 
                                          consistent.

   annot_conflict <STRING>     REQUIRED:  Instructs the pipeline how to handle
                                          overlapping miRBase miRNA loci -
                                          currenlty supported: "ignore", "remove"
                                          and "merge".

   overlap        <INTEGER>    REQUIRED/NOT REQUIRED:  The nucleotide overlap 
                                          between a read and a miRNA necessary 
                                          for it to be included in the miRNA 
                                          read count. Required only when "feature 
                                          miRNA" specified.

   proportional   <FLAG>       OPTIONAL:  Specifies whether the pipeline should 
                                          only consider alignments overlapping 
                                          miRNA loci when distributing reads 
                                          between miRNAs. Multi-mapping reads are
                                          distributed in a fashion that considers
                                          the distribution of uniquely mapping 
                                          reads. DEFAULT: Reads depths are split 
                                          equally between all of their genomic
                                          alignments.
  
   separate_loci  <FLAG>       OPTIONAL:  Specifies whether the counts for miRNA
                                          loci sharing a miRBase mature miRNA ID 
                                          should be merged. DEFAULT: Merge loci

   collapse_method <STRING>    OPTIONAL:  Specifies the method that will be used 
                                          to merge the counts of miRNAs as 
                                          determined by the separate_loci flag 
                                          above. The options include 'mature_id' to 
                                          merge based on the miRBase  mature ID 
                                          and 'sequence' if the miRNAs counts are 
                                          to be merged where mature miRNAs share an 
                                          identical sequence. This option can only 
                                          be provided if the separate_loci option 
                                          is not used. DEFAULT: mature_id
   
   feature repeat:

   ensversion     <INTEGER>    REQUIRED/NOT REQUIRED:  Ensembl version for which
                                          the analysis will be performed. 
                                          WARNING: Ensure the ensembl genome 
                                          build and miRBase builds are 
                                          consistent. Required only when 
                                          "feature  repeat" specified.
  
   repversion     <INTEGER>    REQUIRED/NOT REQUIRED: The version number of the
                                          ncbi bowtie index annotation to which 
                                          the samples will be compared. Required 
                                          only when "feature  repeat" specified.
                                          

   repMaxHits     <INTEGER>    OPTIONAL:  The maximum number of alignments 
                                          allowed between a read and the canonical 
                                          repeat sequences. Relevant only when 
                                          "feature repeat" specified. DEFAULT: 10

   repMismatches  <INTEGER>    OPTIONAL:  The maximum number of mismatches allowed 
                                          between a read and the canonical repeat 
                                          sequences or reference genome (used for 
                                          normalisation). Relevant only when 
                                          "feature repeat" specified. DEFAULT: 0 

   repChunk       <INTEGER>    OPTIONAL:  Value to be passed to the Bowtie 
                                          chunkmbs option. Relevant only when 
                                          "feature repeat" specified. DEFAULT: Not
                                          passed to Bowtie unless set.

EOH
}



####################
### Option check ###
####################

### UNIVERSAL OPTIONS
my $step = "";
my $help  =  0;
my $debug = 0;
my $metaffn = "";
my $configList = "";  
my $userConfig = "";
my $defaultConfig = "";
my $systemCheck = "";
my $processors = 1;
my $paired = "";

my $analysisConfig = "";
my $identifier = "analysis";
my $dirname  = "./";
my $outDir = "./";
my $analysisDir = "";
my $annotDirectory = "";
my $noUnique = 0;

#print "Number of arguments: ",scalar @ARGV,"\n";
#die "Ahahahahaha!\n\n";

example_comm() if (scalar(@ARGV)==0);

if
(! GetOptions(

   # Universal options
   "help"                    =>   \$help,
   "debug"                   =>   \$debug,
   "step=s"                  =>   \$step,
   "description=s"           =>   \$metaffn,
   "user-configuration=s"    =>   \$userConfig, 
   "default-configuration=s" =>   \$defaultConfig,
   "config-list"             =>   \$configList,
   "system-check"            =>   \$systemCheck,
   "paired"                  =>   \$paired,

   # Organise step requires this option
   "dataDir=s"       =>   \$dirname,
   "outDir=s"        =>   \$outDir,
   "identifier=s"    =>   \$identifier,
   "no-unique"       =>   \$noUnique,

   # Sequence_imp steps require these options
   "analysisDir=s"   =>   \$analysisDir,
   "processors=i"    =>   \$processors,

   # Specific sequence_imp steps require this option
   "annotationDir=s" =>   \$annotDirectory,
   )
)
   {  print STDERR "option processing failed\n";
      exit(1);
   }

if ($help) {
   help();
}

if ($configList){
   listConfigs($defDir);
}

if ($systemCheck){
   print STDERR "Checking all of the required software is installed correctly\n";
   die "Could not complete the system check\n" unless system("$sysChck -s complete -r") ==0;   
   print STDERR "All software installed and accessible\n";
   exit(0);
}

version($ENV{SEQIMP_ROOT},"Pipeline");

print STDERR "\n##################\n### DEBUG MODE ###\n##################\n\n" if $debug;

die "Require a 'step' (--step=<STEP>)\n" if (length($step)==0);
if($step eq "organize"){
   $step = "organise";
}
die "Require a 'description' file (--description=<FILE>)\n" if (length($metaffn)==0);

# Either one configuration file or the other must be supplied.
die "Require a 'configuration' file (--user-configuration=<FILE> OR --default-configuration=<FILE>)\n" if ((length($userConfig)==0) && (length($defaultConfig)==0));
die "Only require either --user-configuration=<FILE> OR --default-configuration=<FILE>"
                  if ((length($userConfig)!=0) && (length($defaultConfig)!=0));

$analysisConfig = $userConfig if length($userConfig)!=0;
$analysisConfig = "$defDir/$defaultConfig" if length($defaultConfig)!=0;
die "The configuration file doesn't exist or is empty: $analysisConfig\n" if (!-e $analysisConfig || !-s $analysisConfig);

# Limit the number of processors to be used

die "\nNumber of processors specified must be between 1 and 64.\n" if ($processors < 1 || $processors > 64);

if ( $step eq "organise" && $processors != 1 ){
   print STDERR "\nUnable to parallelise 'organise' step.\n";
}elsif ($processors > 1) {
   print STDERR "\nProcessing FASTQ files in parallel: $processors processors\n";
}

print STDERR "\nThis is the configuration file selected: $analysisConfig\n";
      
#############################################
### Parse the analysis configuration file ###
#############################################
      
print STDERR "\nParsing the analysis configuration file\n";
my ($analysisOpt, $pseudoCL ) = &parseconf($analysisConfig);

###################################
### Read in the geometry config ###
###################################

print STDERR "\nPreparing to read the geometry configuration file\n";
my $geoTable = &readtable($analysisOpt->{'organise'}{'geoConfig'});

#############################################
### Parse the metadata (description) file ###
#############################################

print STDERR "\nPreparing to load description file\n";
die "Can not access the description file\n" unless -e $metaffn && -r _ && -f _ ; 

my $metaTable = &readtable($metaffn);

&check_desc($metaTable, $geoTable, $paired);

die "Unacceptable sample names supplied in the description file: must only contain alphanumerics or underscores\n" if scalar(grep /\W/i, keys( %{$metaTable})) ;
      
#############################################################
### Check that stages are compatible with paired end mode ###
#############################################################

if ($paired){
   die "Paired imp_commandline.pl can only be initaited with the 'organise', 'reaper' and 'filter' steps\n" if ($step ne "organise" && $step ne "reaper" && $step ne "filter");
}


####################################################
### Reorganise the $metaTable for paired samples ###
####################################################

if ($paired){
   print STDERR "Converting the description file into a paired end, pipeline compatible format\n";
   $metaTable = &pairify($metaTable);
}

############################
## STEP: organise OPTIONS ##
############################

if ($step eq "organise"){

   print STDERR "\n####### Beginning step: ".uc($step)." #######\n\n";

   my $randomID = $$;
   my $reap2imp = "";
   
   if ($identifier eq "analysis"){
      print  STDERR "Using default 'identifier': $identifier\n";
   }else{
      print STDERR "Using user defined 'identifier': $identifier\n"
   }

   if ($noUnique){
      print STDERR "Process ID will not be appended to analysis directory name\n"
   }else{
      print  STDERR "Process number for analysis identification: $randomID\n" ;
   }

   print  STDERR "Using default 'dataDir': $dirname\n" if ($dirname eq "./");
   print  STDERR "Using default 'outDir': $outDir\n" if ($outDir eq "./");

   {

      ######################################################################################
      ### Read in the reap2imp config file for translation of description file to Reaper ###
      ######################################################################################
      
      print STDERR "\nPreparing to load Reaper interpreter\n";
      my $reap2impMap = &readtable($analysisOpt->{'organise'}{'reap2imp'});

      ########################################
      ### Identify all fastq files present ###
      ########################################

      # Check that all files in the metadata are present as expected
      die "Can not find the data directory $dirname\n" unless -d $dirname;
      
      my @fastqs = map {$metaTable -> { $_ }{"File"}} keys(%{$metaTable}); 

#      my @dirfastqs;
#      @dirfastqs = glob("$dirname/*.fastq* $dirname/*.fq*");
#      die "No fastq files (.fastq or .fq) found in specified input directory" if (scalar @dirfastqs) == 0 ;
      
      my %locations = ();

      for (@fastqs) {
         $locations{$_} = "$dirname/$_";
         die qq{$_ from the description file not found in $dirname} if (! -s "$dirname/$_");
         &check_fastq($locations{$_});
      }
      
      #foreach my $pattern (@fastqs){
      #   scalar (grep /\Q$pattern\E/, @dirfastqs) && (scalar (grep /$pattern/, @dirfastqs) == 1) || die "Could not find all of the files in the experimental description: $pattern\n";
      #   ($locations{$pattern}) = grep /\Q$pattern\E/, @dirfastqs;
      #print "The location for $pattern is ".$locations{$pattern}."\n";
      #}
      
      # Parse first 4 lines of each fastq to check format

      #######################################
      ### Make a fresh analysis directory ###
      #######################################
     
      my $analysisDirectory = "";
      
      if ($noUnique){
         $analysisDirectory = "$outDir/$identifier";
         die "A file or directory named with the identifier specified already exists\n" if (-e $analysisDirectory);
      }else{
         $analysisDirectory = "$outDir/$identifier"."_$randomID";
         die "The analysis directory $analysisDirectory already exists. Please initiate again to generate new process ID\n" if (-e $analysisDirectory);
      }

      print STDERR "\n***** Creating directory for current analysis: $analysisDirectory\n";

      die "The analysis directory already exists: $analysisDirectory\n" if (-e $analysisDirectory);
      (mkdir $analysisDirectory, 0755) || die "Unable to make $analysisDirectory directory: $!\n";
      
      for my $sample (keys %{$metaTable}){
      
         print STDERR "\nArranging Sequence Imp compatible data structure for $sample\n";

         ###########################################
         ### Assemble data in sample directories ###
         ###########################################
      
         my $sampleDir = "$analysisDirectory/$sample"; 
         die "The sample directory already exists: $sampleDir\n" if (-e $sampleDir);
         mkdir $sampleDir , 0755 || die "Could not make the sample directory for $sample: $!\n";
         my $dataDir = "$sampleDir/data/";
         my $dataFile = "";

         if ($noUnique){
            $dataFile = "$dataDir/$identifier"."_".$metaTable -> { $sample }{"File"}; 
         }else{
            $dataFile = "$dataDir/$identifier"."_$randomID"."_".$metaTable -> { $sample }{"File"}; 
         }

         die "The data directory already exists: $dataDir\n" if (-e $dataDir);
         mkdir $dataDir , 0755 ||"Could not make the data directory for $sample: $!\n";
         my $sampleData = $locations{$metaTable -> { $sample }{"File"}};
         #print "About to move $sampleData\n";
         die "Can not link the fastq file to data directory\n" unless system("ln $sampleData $dataFile") ==0;
         
         print STDERR "Creating Reaper compatible metadata file from description file\n";
         &sortMeta($sample, $sampleDir, $geoTable, $metaTable, $reap2impMap, $paired);
      }
   }
   exit(0);

   ###################################################################################################################################################################


   ###########################################
   ### STEPS: reaper/filter/align/features ###
   ###########################################

}elsif ($step eq 'reaper'|$step eq 'filter'|$step eq 'align'|$step eq 'features') {

   print STDERR "\n####### Beginning step: ".uc($step)." #######\n";

   my $annotOpt = "";

   print STDERR "\n----------------------\nRunning a system check:\n";
   die "System check failed: Closing pipeline.\nPlease initiate pipeline with --system-check to identify missing software.\n" unless system("$sysChck -s $step") ==0;   
   print STDERR "\nFor reference information for the dependencies please perform a system-check (--system-check) or refer to the manual\n\n";
   print STDERR "----------------------\n";

   #########################
   ### Check directories ###
   #########################

   if ($step eq 'align' | $step eq 'features'){
      die "Require an annotation directory for the 'align' and 'features' steps of the pipeline (--annotationDir)\n" if length($annotDirectory)==0;
      die "Can not find 'annotationDir'\n" if ! -d $annotDirectory;
      version($annotDirectory,"Annotation");
      $annotOpt = "--annotationDir=$annotDirectory";
   }

   die "Require an analysis directory (--analysisDir) to be specified\n" if length($analysisDir)==0;
   die "Can not find 'analysisDir'\n" if ! -d $analysisDir;

   my $retval = 0;

   {
   
      my $seqimpDebug = $debug ? "--debug" : "" ;
      my $pairOpt = $paired ? "--paired" : "";


      ##################################################
      ### Check for commands to be passed to seqimp  ###
      ##################################################

      print STDERR $pseudoCL->{$step},"\n" if $debug;

      #########################################################################
      ### Find all directories in the analysis directory: Check as expected ###
      #########################################################################

      ################################################################################ Address the lack of comparison here...

      my @cycleDirs = glob("$analysisDir*");
     
      my @expectedDirs = keys %{$metaTable};

      #################################################################
      ### Build sequence_imp.pl steps and cycle through directories ###
      #################################################################
  
      if ($threads_available){
      
         ### If threaded perl is available, use it to divide the jobs accross multiple processors.
         
         if ($step eq "filter" && $paired){
            # For paired end reads filter step has to be in performed in parallel
            $retval = &thread_call($threads_available, $processors, $debug, $step, $seqimp, $seqimpDebug, $analysisDir, $pseudoCL, $metaTable, $annotOpt, $pairOpt, \@expectedDirs, "Pair");
         }else{
            $retval = &thread_call($threads_available, $processors, $debug, $step, $seqimp, $seqimpDebug, $analysisDir, $pseudoCL, $metaTable, $annotOpt, $pairOpt, \@expectedDirs, "File");
         }

      }else{
     
         ### If threaded Perl is not available revert to old school looping
         
         if ($step eq "filter" && $paired){
            # For paired end reads filter step has to be in performed in parallel
            $retval = &nothreads_call($debug, $step, $seqimp, $seqimpDebug, $analysisDir, $pseudoCL, $metaTable, $annotOpt, $pairOpt, \@expectedDirs, "Pair");
         }else{ 
            $retval = &nothreads_call($debug, $step, $seqimp, $seqimpDebug, $analysisDir, $pseudoCL, $metaTable, $annotOpt, $pairOpt, \@expectedDirs, "File");
         }

      }
   }   
   exit($retval);
}else{
   die "Unrecognised step specified\n";
}


sub thread_call {

   # This subroutine organises threaded calls to sequence_imp.pl
   # Threaded calls are triggered by each sample directory (as determined by the description file)

   my $using_threads = $_[0];                # $threads_available
   my $using_processors = $_[1];             # $processors
   my $debugging = $_[2];                    # $debug
   my $step_thread = $_[3];                  # $step
   my $seqimp_thread =$_[4];                 # $seqimp
   my $seqimpDebug_thread = $_[5];           # $seqimpDebug
   my $analysisDir_thread = $_[6];           # $analysisDir
   my $pseudoCL = $_[7];                     # $pseudoCL - Smae as the original pseudoCL - passed a reference here.
   my $metaTable = $_[8];                    # $metatable - Same as the original metatable - passed a reference here.
   my $annotOpt_thread = $_[9];              # $annotOpt
   my $pairOpt_thread = $_[10];              # $pairOpt
   my @expectedDirs_thread = @{$_[11]};      # \@expectedDirs
   my $centric = $_[12];                     # 'File' || 'Pair'

   my $threaded_output = "";
   my $process = "";
   my @thread_farm;

   print STDERR "Threads are available using perl 'thread' module\n";
   
   my $DataQueue = Thread::Queue->new();

   for (my $threads=0; $threads<$using_processors; $threads++){
      print STDERR "\nCreated Thread ",$threads+1,"!\n";
      $thread_farm[$threads] =  threads->create(
         sub {
            while (my $D = $DataQueue->dequeue()) {
               invoke_syscall_thread($D,$using_processors, $debugging)
            }
         });
   }
   if ($centric eq "File"){
      foreach my $dir (@expectedDirs_thread){
         print_stage_info($dir);
         $process=check_dir($analysisDir_thread,$dir);
         $threaded_output = organise_output($using_processors, $analysisDir_thread, $dir, $step_thread, "");
         if ((length $threaded_output) > 0){
            if ($threaded_output =~ /2>>\s+(\S+)$/){
               print STDERR "Due to threading process, information being redirected to:\t$1\n";
            }else{
               die "Threaded output organisation has failed\n";
            }
         }
         my $thread_job = &make_call( $step_thread, $seqimp_thread, $seqimpDebug_thread, $dir, $process, $metaTable, $pseudoCL, $threaded_output, $annotOpt_thread, $pairOpt_thread, "");
         $DataQueue->enqueue("$thread_job");
      }
   }elsif ($centric eq "Pair"){
     
                              ################################# CHECK ORDER OF FASTQS MAINTAINED FOR LATER USE IN STRAND SPECIFIC MAPPINGS - NOT IMPORTANT FOR TALLY.


      my %pairingup;
      foreach my $fastqFile (keys %{$metaTable}){
         push @{$pairingup{$metaTable->{$fastqFile}{"Pair"}}}, $fastqFile;
      }
      
      #use Data::Dumper;
      #print Dumper(%pairingup);

      foreach my $thispair (keys %pairingup){
         die "Something has gone wrong with pairwise system call\n" if (scalar @{$pairingup{$thispair}} != 2);
         print_stage_info(@{$pairingup{$thispair}}[0]." and ".@{$pairingup{$thispair}}[1]);   
         my $first_dir_chosen = check_dir($analysisDir_thread, $pairingup{$thispair}[0]);
         my $second_dir_chosen = check_dir($analysisDir_thread, $pairingup{$thispair}[1]);
         $threaded_output = organise_output($using_processors, $analysisDir_thread, $pairingup{$thispair}[0],  $step_thread, $pairingup{$thispair}[1]);
        
         if ((length $threaded_output) > 0){
            if ($threaded_output =~ /2>>\s+(\S+)$/){
               print STDERR "Due to threading process, information being redirected to:\t$1\n";
            }else{
               die "Threaded output organisation has failed\n";
            }
         } 

         my $thread_job = &make_call( $step_thread, $seqimp_thread, $seqimpDebug_thread, @{$pairingup{$thispair}}[0]."_and_".@{$pairingup{$thispair}}[1] , $first_dir_chosen, $metaTable, $pseudoCL, $threaded_output, $annotOpt_thread, $pairOpt_thread, "--directory2=$second_dir_chosen");
         print STDERR "$thread_job\n" if $debugging; 
         $DataQueue->enqueue("$thread_job");
      }
   }else{
      die "What happened? No mode specified for preparing system calls!\n";
   }
   
   # Fill the end of the queue with blanks
   for (my $threads=0;$threads<$using_processors;$threads++){
      $DataQueue->enqueue(undef);
   }

   # Gather up the imps
   for (my $threads=0;$threads<$using_processors;$threads++){
      $thread_farm[$threads]->join();
   }

   for (my $threads=0;$threads<$using_processors;$threads++){
      my $err = $thread_farm[$threads]->error();
      if ($err) { 
         print STDERR "Failure during execution of threads!\n";
         return 14;
      }
   }
   return 0;
}

sub nothreads_call{
  
   # If imp_commandline does not have access to threaded perl this subroutine will run each imp call sequentially.

   my $debug_noth = $_[0];                 # $debug
   my $step_noth = $_[1];                  # $step
   my $seqimp_noth =$_[2];                 # $seqimp
   my $seqimpDebug_noth = $_[3];           # $seqimpDebug
   my $analysisDir_noth = $_[4];           # $analysisDir
   my $pseudoCL = $_[5];                   # $pseudoCL - Smae as the original pseudoCL - passed a reference here.
   my $metaTable = $_[6];                  # $metatable - Same as the original metatable - passed a reference here.
   my $annotOpt_noth = $_[7];              # $annotOpt
   my $pairOpt_noth = $_[8];              # $pairOpt
   my @expectedDirs_noth = @{$_[9]};      # \@expectedDirs
   my $centric_noth = $_[10];             # 'File' || "Pair"


   my $current_sample = "";
   my $retval = 0;

   print STDERR "\n\n *** As the current perl installation does not support threading, all steps of the pipeline will be performed sequentially. ***\n";   
   if ($centric_noth eq "File"){ 
      foreach my $dir_now (@expectedDirs_noth){
         
         #################################################################
         ### Build sequence_imp.pl steps and cycle through directories ###
         #################################################################
         
         print STDERR "\n\n******######******######******######******\n#######*** Beginning sequence imp pipeline for $dir_now ***#######\n******######******######******######******\n\n";
         $current_sample = check_dir($analysisDir_noth,$dir_now);


         my $nothread_job = &make_call( $step_noth, $seqimp_noth, $seqimpDebug_noth, $dir_now, $current_sample, $metaTable, $pseudoCL, "", $annotOpt_noth, $pairOpt_noth, "");
         print STDERR "$nothread_job\n\n" if $debug_noth;
         print STDERR "sequence_imp.pl hasn't worked: $nothread_job: $!\n" unless system($nothread_job)==0 ;
         $retval = 14;
      }
   }elsif($centric_noth eq "Pair"){
      my %pairingup_noth;
      foreach my $fastqFile (keys %{$metaTable}){
         push @{$pairingup_noth{$metaTable->{$fastqFile}{"Pair"}}}, $fastqFile;
      } 
      foreach my $thispair_noth (keys %pairingup_noth){
         
         print STDERR "\n\n******######******######******######******######******######******\n###*** Beginning sequence imp pipeline for ".@{$pairingup_noth{$thispair_noth}}[0]." and ".@{$pairingup_noth{$thispair_noth}}[1]." ***###\n******######******######******######******######******######******\n\n";
         
         my $first_dir_chosen_noth = check_dir($analysisDir_noth, $pairingup_noth{$thispair_noth}[0]);
         my $second_dir_chosen_noth = check_dir($analysisDir_noth, $pairingup_noth{$thispair_noth}[1]);
         my $nothread_job = &make_call( $step_noth, $seqimp_noth, $seqimpDebug_noth, @{$pairingup_noth{$thispair_noth}}[0]."_and_".@{$pairingup_noth{$thispair_noth}}[1], $first_dir_chosen_noth, $metaTable, $pseudoCL, "", $annotOpt_noth, $pairOpt_noth, "--directory2=$second_dir_chosen_noth");
         print STDERR "$nothread_job\n\n" if $debug_noth;
         print STDERR "sequence_imp.pl hasn't worked: $nothread_job: $!\n" unless system($nothread_job)==0 ;
         $retval = 14;
      }

   }else{
      print STDERR "What happened? No mode specified for preparing system calls in non-threaded system!\n";
      $retval = 14;
   }
   return $retval;
}


sub organise_output{
   
   my $procNo = shift;
   my $anaSel = shift;
   my $dirSel = shift;
   my $appl   = shift;
   my $dirSel2 = shift;

   my $directed_output = "";

   if ($procNo > 1) {
      if (length ($dirSel2)==0){
         my $errfile = "$anaSel/$dirSel/$dirSel"."_$appl"."_record.txt";
         $directed_output = "2>> $errfile";
      }else{
         my $errfile = "$anaSel/$dirSel"."_and_$dirSel2"."_$appl"."_record.txt";
         $directed_output = "2>> $errfile";
      }
   }
   
   return ($directed_output);
}

# Parse tab deliminated files
sub readtable {
  
   print STDERR "Reading table\n";
   print STDERR "WARNING: &readtable should be passed a single file!\n" if @_ != 1;

   my $fileFfn = shift @_ ;

   die "Please check that file $fileFfn is compatible with the pipeline. Should be a tab delimited text file\n"
         if (! -e $fileFfn || -z _ || ! -f _ );

   my $top;
   my @top;
   my %fileCols = ();
   my $fileId = 0;
   my %parsedTable = ();
   my %tableRows;

   open TABLE, "< $fileFfn" || die "Can't open $fileFfn: $!\n"; 
 
   $top = <TABLE>;
   chomp $top;
   @top = split "\t", $top, -1;
   $fileCols{$fileId++} = $_++ for @top;
 
   while (<TABLE>) {
      chomp;
      my @stuff = split "\t", $_,-1;
      die "The number of fields does not match the number of columns for row \'$stuff[0]\'. Please check that the description file is tab seperated.\n" if ((scalar @top) != (scalar @stuff)); 
      
      die "'Name' column does not contain unique identifier\n" if exists($tableRows{$stuff[0]});
      $tableRows{$stuff[0]}++;
      
      $parsedTable{$stuff[0]} = {};
      for (my $i=1; $i < @top; $i++) {
         $parsedTable{$stuff[0]}{$fileCols{$i}} = $stuff[$i];
      } 
   }
   close(TABLE);
   print STDERR "Reading of table complete\n";
   return(\%parsedTable);
}

# Check the parsed description file to make sure everything is correct and as expected
sub check_desc{
   print STDERR "Checking the description table\n";
   print STDERR "WARNING: &check_desc should be passed a hash and a pipeline pairing mode\n" if @_ != 3;
 
   my $desc_table = shift @_;
   my $geo_setup = shift @_;
   my $pairing = shift @_;

   my @expected_cols = ("File","Geometry","Barcodes","5p_ad","3p_ad","5p_seq_insert","3p_seq_insert");

   my @err_names = grep /\W/i, keys %{$desc_table};
   die ("Unacceptable 'Names' supplied in description file: ",(join ", ",@err_names),".Names may only contain alphanumerics or underscores.\n") if scalar @err_names;

   foreach my $curr_name (keys %{$desc_table}){
      my @err_columns = grep /\W/i, keys %{$desc_table->{$curr_name}};      
      die ("Unacceptable column names supplied in description file: ",(join ", ",@err_columns),".Names may only contain alphanumerics or underscores.\n") if scalar @err_columns;
      foreach my $curr_col (@expected_cols){
         die ("Could not find \'$curr_col\' in description table column names. Please check description file configuration.\n") if ! scalar(grep /$curr_col/, keys %{$desc_table->{$curr_name}});
      }
      foreach my $column_check (keys %{$desc_table->{$curr_name}}){
         die "There is an empty column in the description file for: $curr_name -> $column_check. Unused columns should contain a '-'.\n" if ! length $desc_table->{$curr_name}->{$column_check}; 
         die "There is whitespace found in the description file for: $curr_name -> $column_check\n" if $desc_table->{$curr_name}->{$column_check} =~ /\s+/i;
      }
      die "Unrecognised geometry specified in the description file ($curr_name). Please check.\n" unless exists($geo_setup->{$desc_table->{$curr_name}{Geometry}});
   
      # Check that the mode (paired/unpaired) in which the pipeline is being run is compatible with the geometry selected.
      if ($pairing){
         die "The 'Geometry' specified for $curr_name is incompatible with the pipeline 'Paired' mode.\n" if ($geo_setup->{$desc_table->{$curr_name}{Geometry}}{Paired} eq "NO");
      }elsif (! $pairing){
         die "The 'Geometry' specified for $curr_name is only compatible with the pipeline 'Paired' mode.\nIf intention is to perform a paired end analysis, please specify the '--paired' flag.\n" if ($geo_setup->{$desc_table->{$curr_name}{Geometry}}{Paired} eq "YES");
      }
      die "Multiple fastq files in a single line of the description file is only acceptable if pipeline is run in paired mode (specify the '--paired' flag).\n" if (!$pairing && $desc_table->{$curr_name}->{"File"} =~ /,/);
   }
}

# Parse the analysis config file
sub parseconf {

   # DEFAULTS:
   my $geoDefault =  "$ENV{SEQIMP_ROOT}/config/seqimp_config/geometry_config.txt";
   my $reap2impDefault = "$ENV{SEQIMP_ROOT}/config/seqimp_config/reap2imp_mappings.txt";

   print STDERR "Parsing the config file\n";
   my $config = "";
   my %Options = ();
   my $line = 0;
   my $stageOpt = "";

   $config = shift @_ ;

   open (CONF, "< $config") || die "Can not open the config file: $config: $!\n";
   
   while (<CONF>) {
      $line++;
      next if m/^#/;
      chomp;
      
      # Record the section/stage headers
      if (m/^@/) {
         $stageOpt = $_;
         $stageOpt =~ s/^@//;
         #print "$stageOpt\n";
         %{$Options{$stageOpt}} = ();
      }else{
      
         # Ensure that a header has been specifed
         die "Failed to parse config file. Config file stage unspecified (Must begin with @). Error at line: $line\n" if length($stageOpt) == 0;
         # Empty lines are not allowed
         die "Empty lines are not allowed in config file. Error at line: $line\n" if length($_) ==0;
         # Check line format
         die "Unrecognised format in config file. Variables and values must be seperated by a single tab. Spaces are not allowed. Error at line: $line\n" if (! ($_ =~ /^\S+\t\S+$/));
         
         my @currOpt = split "\t", $_;
         
         # Check option/config file format
         die "Failed to parse config file. Error at line: $line\n" if (scalar @currOpt != 2 || length $currOpt[0] == 0 || length $currOpt[1] ==0);
         #print "The valuse of ",$currOpt[0]," is ",$currOpt[1],"\n";
         
         # Skip all NA values
         #print "Skipping ".$currOpt[0]." because not applicable\n" if $currOpt[1] eq "NA";
         next if $currOpt[1] eq "NA";        # Use NA to specify which arguments should be ignored

         # Record all options in a hash   
         $Options{$stageOpt}->{$currOpt[0]} = $currOpt[1];
         #print "$stageOpt:$currOpt[0]:$currOpt[1]\n"
      }
   }
   
   close (CONF);
   
   # Check that all expected stages are accounted for
   
   my @expected_headers = ("organise", "reaper", "filter", "align+features", "align", "features");
   
   foreach (@expected_headers){
      if (! exists $Options{$_}){
         die "Missing an expected section break in the config file: $_\n";
      }
   }

   #die "Check all the expected headers are present in the decription file\n" unless 
   #               (exists $Options{'organise'} && exists $Options{'reaper+filter'} && 
   #                exists $Options{'reaper'} && exists $Options{'align+features'} && 
   #                exists $Options{'align'} && exists $Options{'features'} && exists $Options{'filter'});

   # Check that the required options for the 'organise' step is specified - set defaults.
   if (! exists $Options{'organise'}{'geoConfig'}){
      print STDERR "geoConfig not specified, using default: $geoDefault\n";
      $Options{'organise'}{'geoConfig'} = $geoDefault;
   }
   if (! exists $Options{'organise'}{'reap2imp'}){
      print STDERR "reap2imp not specified, using default: $reap2impDefault\n";
      $Options{'organise'}{'reap2imp'} = $reap2impDefault;
   }

   my %commandline = ();
  
   # Organise the options specified into a pseudo command line 
   for my $stage (keys %Options){
      my @optArray = ();
      for my $selectedOpt (keys %{$Options{$stage}}){
         my $collapsed = "";
         if ($Options{$stage}{$selectedOpt} eq "FLAG"){
            $collapsed = "--$selectedOpt";
            #print "Found a flag: $collapsed\n";
         }else{
            $collapsed = join "=", ($selectedOpt, $Options{$stage}{$selectedOpt});
            $collapsed = "--$collapsed";
         }
         push @optArray , $collapsed;
      }

      # Some organise options are a special cases and are required for reaper command lines too
      push @optArray, "--geoConfig=".$Options{'organise'}{'geoConfig'} if $stage eq 'reaper';
      $commandline{$stage} = join " ", @optArray;
   }

   print STDERR "Config file parsing complete\n";   

   return (\%Options, \%commandline);
   
}

# Oganise metadata files for Reaper based on fields in user supplied metadata and Reaper requirements 
sub sortMeta {
   
   print STDERR "\n";
  
   # Organise subroutine arguments
   my $sampleName = shift @_;
   my $sampleDirectory = shift @_;
   my $geometryInfo = shift @_;
   my $descriptionTable = shift @_;
   my $reapmap = shift @_;
   my $pairSel = shift @_;

   my $sampleGeo = $descriptionTable-> {$sampleName} -> {"Geometry"};
   my $metaDir = "$sampleDirectory/metadata/";

   # Create metadata directory
   die "The metadata directory already exists: $metaDir\n" if (-e $metaDir);
   mkdir $metaDir , 0755 ||"Could not make the metadata directory for $sampleName: $!\n";

   print STDERR "Defining metadata table for sample: $sampleName\n";
   print STDERR "Specified geometry: $sampleGeo\n";

   # Check that all expected fields are supplied by the user

   my $bar = $descriptionTable-> { $sampleName } -> {"Barcodes"}; 
   my $five_ad = $descriptionTable-> { $sampleName } -> {"5p_ad"};
   my $three_ad = $descriptionTable-> { $sampleName } -> {"3p_ad"};
   my $five_si =$descriptionTable-> { $sampleName } -> {"5p_seq_insert"};
   my $three_si = $descriptionTable-> { $sampleName } -> {"3p_seq_insert"};


   if ($geometryInfo-> {$sampleGeo}{"Barcodes"} eq "YES")  { 
                                                                  die "There is no 'Barcodes' info for $sampleName\n" if $bar =~ /-/;  
                                                                  print STDERR "Checking barcode information for $sampleName\n";
                                                            }else{
                                                                  die "\'Barcodes\' specification is incompatible with geometry selected for $sampleName. Please check your description file.\n" if $bar ne "-";
                                                                  $bar = "nobar"; 
                                                                  print STDERR "Barcode are not required for $sampleName\n";
                                                            }
   if ($geometryInfo-> {$sampleGeo}{"5p_ad"} eq "YES")     { 
                                                                  die "There is no '5p_ad' info for $sampleName\n" if $five_ad  =~ /-/; 
                                                                  print STDERR "Checking 5 prime adapter information for $sampleName\n";
                                                            }elsif ($geometryInfo-> {$sampleGeo}{"5p_ad"} eq "MAYBE"){
                                                                  print STDERR "5 prime adapter is optional for $sampleName\n";
                                                            }else{
                                                                  die "\'5p_ad\' specification is incompatible with geometry selected for $sampleName. Please check your description file.\n" if $five_ad ne "-";
                                                                  print STDERR "5 prime adapter is not required for $sampleName\n";
                                                            }
   if ($geometryInfo-> {$sampleGeo}{"3p_ad"} eq "YES")     { 
                                                                  die "There is no '3p_ad' info for $sampleName\n" if $three_ad  =~ /-/; 
                                                                  print STDERR "Checking 3 prime adapter information for $sampleName\n"; 
                                                            }else{
                                                                  die "\'3p_ad\' specification is incompatible with geometry selected for $sampleName. Please check your description file.\n" if $three_ad ne "-";
                                                                  print STDERR "3 prime adapter is not required for $sampleName\n";
                                                            }      
   if ($geometryInfo-> {$sampleGeo}{"5p_seq_insert"} eq "YES") { 
                                                                     die "There is no '5p_seq_insert' info for $sampleName\n" if $five_si  =~ /-/;
                                                                     print STDERR "Checking 5 prime sequence insert information for $sampleName\n" ;
                                                                }else{
                                                                     die "\'5p_seq_insert\' specification is incompatible with geometry selected for $sampleName. Please check your description file.\n" if $five_si ne "-";
                                                                     print STDERR "5 prime sequence insert is not required for $sampleName\n";
                                                                }
   if ($geometryInfo-> {$sampleGeo}{"3p_seq_insert"} eq "YES") { 
                                                                     die "There is no '3p_seq_insert' info for $sampleName\n" if $three_si  =~ /-/; 
                                                                     print STDERR "Checking 3 prime sequence insert information for $sampleName\n"; 
                                                                }else{ 
                                                                     die "\'3p_seq_insert\' specification is incompatible with geometry selected for $sampleName. Please check your description file.\n" if $three_si ne "-";
                                                                     print STDERR "3 prime sequence insert is not required for $sampleName\n";
                                                                }
   
   local $, = "\t";

   # print results to metadata file
   open (META, "> $metaDir/metadata.txt") || die "Can not open the metadata file: $metaDir/metadata.txt : $!\n";

   my $header = $geometryInfo-> {$sampleGeo}{"Columns"};
   my @header = split /,/, $header;
   my $headline = join "\t", @header;
   print META "$headline";

   my @barArray = split /,/, $bar;
   for my $curr (@barArray){
      print META "\n";
      my @currLine;
      for my $option (@header){
         # Use the reap2imp file to map column headers between the imp input requirements and Reaper
         if ($reapmap-> {$option}{"Imp"} eq "Barcode"){
            push @currLine, $curr;    
         }else{
            my $impHeader = $reapmap-> {$option}{"Imp"};
            my $field = $descriptionTable-> { $sampleName }{$impHeader};
            push @currLine, $field;
         }
      }
      print META @currLine;
   }
   close META;
   print STDERR "Metadata table for $sampleName completed.\n\n";
}


sub listConfigs {
   # List the default configuration files shipped with the pipeline   
   my $preConfigDir = shift @_;

   my $confCount = 1;
   print STDERR "\nThe following precomputed 'Configuration' files are available for your analysis:\n";
   
   my @preConfigFiles =  glob("$preConfigDir*.txt");
   
   for my $currPreConf (@preConfigFiles) {
      my $shortConfName = ""; 
      if ($currPreConf =~ /.*\/([^\/]+)$/){
         $shortConfName = $1;
      }
      print STDERR "\n$confCount. $shortConfName\n";
      open (PRECONF , "< $currPreConf") || die "Can not open the predefined config file: $currPreConf : $!\n";
      while (<PRECONF>){
         if (/^#/){
            chomp;
            print STDERR "$_\n";
         }else{
            next;
         }
      }
      close (PRECONF);
      $confCount++;
   }
   exit(0);
}

sub version {

   # Uses the VERSION files created upon the assembly of a release to report the pipeline version

   my $versDir = shift @_;
   my $versType = shift @_;
   die "Unable to determine the version for the ".lc($versType)."\n" if ! -s "$versDir/VERSION";
   open(VERSION, "< $versDir/VERSION") || die "Can't open the version records for the $versType\n";
   my $versInfo = <VERSION>;
   chomp $versInfo;
   close (VERSION);

   print STDERR "\n$versType version:\t$versInfo\n";
}

# Checks that the FASTQ files supplied are correct format eg. first line begins with @ and third line with +

sub check_fastq {
   my $selected_fastq = shift @_;
   open(INFASTQ, "gzip -cdf $selected_fastq | cat|") || die "Can't open FASTQ\n";
   my $lineONE = <INFASTQ>;
   my $lineTWO = <INFASTQ>;
   my $lineTHREE = <INFASTQ>;
   my $lineFOUR = <INFASTQ>;
   #print "$lineONE\t$lineTWO\t$lineTHREE\t$lineFOUR\n";
   close(INFASTQ);

   die "$selected_fastq does not appear to conform to FASTQ format" if ($lineONE !~ /^@/ || $lineTHREE !~ /^\+/);
}

# Handling the paired end description file (Expanding to one sample per lane while including the pairing information).

sub pairify {
    
   my %tempTable = ();
   my $descRef = shift @_;
   my $lineSel = 1;
  
   foreach my $pairname (keys %{$descRef}){
      print STDERR "\nOrganising paired sequencing reactions for $pairname\n";

      ### Currently only support a very simple geometry for paired end sequencing. If Barcodes or sequence inserts supplied in geometry file give up.
      my %supported = (File => 1, Geometry => 1, '5p_ad' => 1, '3p_ad' => 1); 
      
      my @notsupported = grep { !defined($supported{$_}) } keys %{$descRef->{$pairname}};
     
      for my $despised (@notsupported){
         if ($descRef->{$pairname}{$despised} ne '-'){
            die "Please check $pairname in the Description file. Currently only support 5p_ad and 3p_ad for Paired mode. Other sequence columns should contain a '-'\n";
         }
      }

      ### Seperate description file by fastq's but record pairing
      my @fileChoices = split ",", $descRef->{$pairname}{"File"};  
      die "In order to analyse a paired end experiment two sequence files are required per line in the description files (Comma seperated)\n" unless (scalar @fileChoices == 2);      
      #print "@fileChoices\n";
      
      my @FivePrimes  = ($descRef->{$pairname}{"5p_ad"});
      my @ThreePrimes = ($descRef->{$pairname}{"3p_ad"});

      ### Reverse compliment of adapter sequences for second sequence files
      push @FivePrimes, &revcomp($descRef->{$pairname}{"3p_ad"});
      push @ThreePrimes, &revcomp($descRef->{$pairname}{"5p_ad"});

      my $twin = 0;
      
      ### Geometry table should include paired geometries with paired flag requirement.

      while ($twin <= $#fileChoices){
         my $thisID = "$pairname"."_".($twin+1);
         
         print STDERR "\nArranging data for $fileChoices[$twin]\nCorresponding directory: $thisID\n5' Adapter sequence: $FivePrimes[$twin]\n3' Adapter sequence: $ThreePrimes[$twin]\n";

         $tempTable{$thisID} = {
            File              => $fileChoices[$twin],
            Geometry          => $descRef->{$pairname}{"Geometry"},
            Barcodes          => $descRef->{$pairname}{"Barcodes"},
            "5p_ad"           => $FivePrimes[$twin],
            "3p_ad"           => $ThreePrimes[$twin],
            "5p_seq_insert"   => $descRef->{$pairname}{"5p_seq_insert"},
            "3p_seq_insert"   => $descRef->{$pairname}{"3p_seq_insert"},
            "Pair"            => $lineSel
         };
         $twin++;
      }

      $lineSel++;
   }
   return(\%tempTable);
}

# Reverse complement sequences (in a fairly circumstance specific situation).

sub revcomp {
   my $sequence = shift @_;
   if ($sequence eq "-"){
      return($sequence);
   }elsif($sequence =~ /[^ATGC]+/i){
      die "Unrecognised base in adapters. Must be A, T, G or C.\n";
   }else{
      $sequence = uc $sequence;
      my %rc = ( A => 'T', C => 'G', G => 'C', T => 'A', U => 'A' );
      $sequence = join "" , (map { $rc{$_} } reverse split "", $sequence);
      return $sequence;
   }
}

# Set of subroutines used to organise threading of experiments

sub invoke_syscall_thread{
   my $syscall=$_[0];
   my $procNo = $_[1];
   my $dodebug = $_[2];
   print STDERR "$syscall\n\n" if $dodebug;
   if ($procNo > 1) {
      die "sequence_imp.pl hasn't worked: $syscall: $!\nFor diagnostic information please see the sample directory record files (See above).\n" unless system($syscall)==0;
   }
   else {
      die "sequence_imp.pl hasn't worked: $syscall: $!\n" unless system($syscall)==0;
   }
}


sub print_stage_info{
         print STDERR "\n\n####### Beginning sequence imp pipeline for $_[0] #######\n";
} 

sub check_dir{
   my $process="$_[0]/$_[1]/";
        -d $process || die "Could not identify the directories specified by the sample names in the experimental description file\n";
   return($process);
}

sub make_call {
   my $step_sub = $_[0];
   my $seqimp_sub = $_[1];
   my $seqimpDebug_sub = $_[2];
   my $dir_sub = $_[3];
   my $process_sub = $_[4];
   my $metaTable_sub = $_[5];
   my $pseudoCL_sub = $_[6];
   my $threaded_output_sub = $_[7];
   my $annotOpt_sub = $_[8];
   my $pairOpt_sub = $_[9];
   my $process_sub_2 =$_[10];

   my $job_build;
   
   if ($step_sub eq 'reaper'){
      $job_build = "$seqimp_sub --stage=$step_sub $seqimpDebug_sub --tag=$dir_sub --directory=$process_sub --geometry=".$metaTable_sub->{$dir_sub}{Geometry}." ".
                                          $pseudoCL_sub->{$step_sub}." --data=data/ --meta=metadata/metadata.txt --basicName=$dir_sub ".$pseudoCL_sub->{'reaper+filter'}." $pairOpt_sub $threaded_output_sub";
                  # print "$first_dir_chosen\t$second_dir_chosen\n"; 
   }elsif ($step_sub eq 'filter'){
      $job_build = "$seqimp_sub --stage=$step_sub $seqimpDebug_sub --tag=$dir_sub --directory=$process_sub --reapData=REAPER/ ".$pseudoCL_sub->{$step_sub}." ".$pseudoCL_sub->{'reaper+filter'}." $process_sub_2 $pairOpt_sub $threaded_output_sub";
   }elsif ($step_sub eq 'align'){
      $job_build= "$seqimp_sub --stage=$step_sub $seqimpDebug_sub --tag=$dir_sub --directory=$process_sub ".$pseudoCL_sub->{$step_sub}." ".$pseudoCL_sub->{'align+features'}." $annotOpt_sub $threaded_output_sub"; 
   }elsif ($step_sub eq 'features'){
      $job_build = "$seqimp_sub --stage=$step_sub $seqimpDebug_sub --tag=$dir_sub --directory=$process_sub ".$pseudoCL_sub->{$step_sub}." ".$pseudoCL_sub->{'align+features'}." $annotOpt_sub $threaded_output_sub";
   
      #### Fold back into the filter stage of the pipeline
      # Remove automatic paired tally step and incorporate into a second filter function.
      # CHeck that length and trinucleotide specification is consistent (eg. >= etc.).
      # Pass process_sub_2 as a "" : "--directory2 option" - Coordinate this via the --paired flag.
      # In sequence_imp die if 3' trimming specified. Tally can not do this.
      # Introduce --sumstat flag to tally and link between the paired end directories.
      # New tally has change in --record-format2 flag. Change this... and check.
      # Also pass $pairOpt_sub to sequence_imp to control the tally/filter principles used.

   }else{
      die "Unrecognised stage when building sequence_imp command lines\n\n";
   }
   return($job_build);
}

