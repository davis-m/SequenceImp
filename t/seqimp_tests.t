#!/usr/local/bin/perl

###########################################################################

#This script is part of the Kraken framework which aims to facilitate RNA 
#sequence analysis in a streamlined and efficient manner.
#Copyright (C) 2017 EMBL - European Bioinformatics Institute 

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


use strict;
use warnings;
use Test::More;
use Cwd 'abs_path';
use File::Spec;
use File::Basename;
use File::Path;
use IO::Select;
use Symbol 'gensym';
use IPC::Open3;
use Env;

my $SEQIMP = "imp_commandline.pl";

# Set up
my ($test_script,$test_directory) = fileparse(abs_path($0));
my $SEQIMP_ROOT = abs_path(File::Spec->catfile($test_directory, ".." ));
my $RESOURCES_PATH = File::Spec->catfile($test_directory,"resources");
my $TEST_DATA_PATH = File::Spec->catfile($SEQIMP_ROOT,"../../data-for-tests");

# Output
my $TEST_DIR = File::Spec->catfile($SEQIMP_ROOT, "test-run");
my $TEST_RUN =  File::Spec->catfile($TEST_DIR,"run");
my $TEST_MISC =  File::Spec->catfile($TEST_DIR,"misc");

# Add seqimp_annot bin to PATH
$ENV{PATH} = abs_path(File::Spec->catfile($SEQIMP_ROOT,"bin")).":".$ENV{PATH};

#print STDERR "$SEQIMP_ROOT\n$RESOURCES_PATH\n$TEST_DATA_PATH\n$TEST_DIR\n$TEST_RUN\n$TEST_MISC\n";

# Prepare test directories
sub prepare {
   mkpath($TEST_DIR);
   mkpath($TEST_RUN);
   mkpath($TEST_MISC);
}

die "Please remove previous test-run data: ($TEST_DIR)\n" if (-e $TEST_DIR);
&prepare();

# Cleanup test directories
sub cleanup {
   rmtree($TEST_DIR);
}

### PAIRED END TESTS
my $pairedID = "paired_test";
my $pairedDir = "$TEST_RUN/$pairedID/";
my $pairedDataIn = File::Spec->catfile($TEST_DATA_PATH,"paired");
my $pairedResources = File::Spec->catfile($RESOURCES_PATH,"paired");

# Paired organise
my $paired_organise_call = "$SEQIMP --step=organise --outDir=$TEST_RUN --description=$pairedResources/description.txt --user-configuration=$pairedResources/paired_mRNA_config.txt --dataDir=$pairedDataIn --paired --identifier=$pairedID --no-unique --debug";
runStep("Paired FASTQ organise",$paired_organise_call,"run");
#die "Could not run the paired end organise step\n" if system($paired_organise_call);

# Paired reaper
my $paired_reaper_call = "$SEQIMP --step=reaper --description=$pairedResources/description.txt --user-configuration=$pairedResources/paired_mRNA_config.txt --analysisDir=$pairedDir --paired --debug";
runStep("Paired FASTQ reaper",$paired_reaper_call,"run");
#die "Could not run the paired end reaper step\n" if system($paired_reaper_call);

# Paired filter
my $paired_filter_call = "$SEQIMP --step=filter --description=$pairedResources/description.txt --user-configuration=$pairedResources/paired_mRNA_config.txt --analysisDir=$pairedDir --paired --debug";
#print STDERR $paired_filter_call;
runStep("Paired FASTQ filter",$paired_filter_call,"run");
#die "Could not run the paired end filter step\n" if system($paired_filter_call);

# Paired end FASTQ tests
#
# Check that paired end files exist : .1 and .2

my $fastq_pair_1 = "$pairedDir/paired_end_AC_1/PROCESSED/paired_end_AC_1.lane.tallied.fastq.gz";
my $fastq_pair_2 = "$pairedDir/paired_end_AC_2/PROCESSED/paired_end_AC_2.lane.tallied.fastq.gz";
my $fastq_pair_3 = "$pairedDir/paired_end_GT_2/PROCESSED/paired_end_GT_2.lane.tallied.fastq.gz";

check_a_file_exists($fastq_pair_1, "Paired fastq (/1), first processed file exists (paired_end_AC_1)");
check_a_file_exists($fastq_pair_2, "Paired fastq (/2), second processed file exists (paired_end_AC_2)");
check_a_file_exists($fastq_pair_3, "Paired fastq second set (/2), second processed file exists (paired_end_GT_2)");

# Collapsing reads
# ID Check - first ID
# x2 begin with @
# Check first sequence matches duplicate
#
# Indexing
# ID Check - second ID
# /1 or /2

check_a_specific_line('^@read1_x2/1$',$fastq_pair_1,1,"Paired Fastq (/1), read collapsed and indexed");
check_a_specific_line('^ACTGCATACATCTGAAGAACAAAAACATCAACGTCTTTTGTCCAGCCTCTTTTTCTTCTGCTGTTCCACCTTTCTA$',$fastq_pair_1,2,"Paired Fastq (/1), collapsed read expected sequence");
check_a_specific_line('^@read2_x1/1$',$fastq_pair_1,5,"Paired Fastq (/1), read 2 ID in FASTQ format");
check_a_specific_line('^A',$fastq_pair_1,6,"Paired Fastq (/1), expected sequences (A)");

# File order
# Check first letter of read A (first file, first end) and T (second file, second end)
check_a_specific_line('^@read1_x2/2$',$fastq_pair_2,1,"Paired Fastq (/2), read collapsed and indexed");
check_a_specific_line('^CTGTCACTTTAGTACAACAAAATAGGATGTTAATCCATTACCACTGTTCTTCCATTCGAGAGGCAGCTGATGGTTA$',$fastq_pair_2,2,"Paired Fastq (/2), collapsed read expected sequence");

check_a_specific_line('^T',$fastq_pair_3,10,"Paired Fastq - second set (/2), file input -> expected file output");


### PAIRED END TEST FASTA
my $pairedFastaID = "paired_test_fasta";
my $pairedFastaDir = "$TEST_RUN/$pairedFastaID/";

# Paired organise
my $paired_fasta_organise_call = "$SEQIMP --step=organise --outDir=$TEST_RUN --description=$pairedResources/description.txt --user-configuration=$pairedResources/paired_mRNA_fasta_config.txt --dataDir=$pairedDataIn --paired --identifier=$pairedFastaID --no-unique --debug";
runStep("Paired FASTA organise",$paired_fasta_organise_call,"run");
#die "Could not run the paired end, FASTA format, organise step\n" if system($paired_fasta_organise_call);

# Paired reaper
my $paired_fasta_reaper_call = "$SEQIMP --step=reaper --description=$pairedResources/description.txt --user-configuration=$pairedResources/paired_mRNA_fasta_config.txt --analysisDir=$pairedFastaDir --paired --debug";
runStep("Paired FASTA reaper",$paired_fasta_reaper_call,"run");
#die "Could not run the paired end, FASTA format, reaper step\n" if system($paired_fasta_reaper_call);

# Paired filter
my $paired_fasta_filter_call = "$SEQIMP --step=filter --description=$pairedResources/description.txt --user-configuration=$pairedResources/paired_mRNA_fasta_config.txt --analysisDir=$pairedFastaDir --paired --debug";
#print STDERR $paired_fasta_filter_call;
runStep("Paired FASTA filter",$paired_fasta_filter_call,"run");
#die "Could not run the paired end, FASTA format, filter step\n" if system($paired_fasta_filter_call);

# Check that paired FASTA files produced
my $fasta_pair_1 = "$pairedFastaDir/paired_end_GT_1/PROCESSED/paired_end_GT_1.lane.tallied.fasta.gz";
my $fasta_pair_2 = "$pairedFastaDir/paired_end_GT_2/PROCESSED/paired_end_GT_2.lane.tallied.fasta.gz";

check_a_file_exists($fasta_pair_1, "Paired fasta (/1), first processed file exists (paired_end_GT_1)");
check_a_file_exists($fasta_pair_2, "Paired fasta (/2), first processed file exists (paired_end_GT_2)");

# Paired end FASTA tests
#
# ID Check - second ID
# /1 or /2
check_a_specific_line('>read3_x1/1',$fasta_pair_1,5,"Paired Fasta (/1), read identifier indexed as expected");
check_a_specific_line('>read5_x1/2',$fasta_pair_2,9,"Paired Fasta (/2), read identifier indexed as expected");
check_a_specific_line('^G',$fasta_pair_1,32,"Paired Fasta (/2), file input -> expected file output");
check_a_specific_line('^T',$fasta_pair_2,102,"Paired Fasta (/2), file input -> expected file output");

cleanup();
done_testing(22);

# Subroutines
# Check that a specified file exists
sub check_a_file_exists {

   my $expected_number = 2;
   die "check_a_file_exists: Needs $expected_number arguments\n" if (scalar(@_)!=$expected_number);

   my $the_file = shift;
   my $test_name = shift;

   ok(-e $the_file, "$test_name - Test file exists ($the_file)");
}


# Check a specific line in the file to make sure it is as expected
sub check_a_specific_line {

   my $expected_number = 4;
   die "check_a_specific_line: Needs $expected_number arguments\n" if  (scalar(@_)!=$expected_number);

   my $pattern = shift;
   my $file = shift;
   my $lineNumber = shift;
   my $test_name = shift;

   SKIP: {
      skip "Skipping - $file is missing or empty\n", 1 if (! -s $file);

      my $line_count = 0;
      open(LINE, "gzip -cdf $file |") || die "check_a_specific_line: Could not open $file\n";

      while(<LINE>){
         $line_count++;
         chomp;
         my $this_line = $_;
         if ($line_count == $lineNumber){
            ok($this_line=~/$pattern/, "$test_name - Test expected line ($pattern): $file line $lineNumber");
         }
      }

      close(LINE) || die "check_a_specific_line: Can not close the connection to $file\n";
   }
}

# Running each step
sub runStep {
   # Help from http://stackoverflow.com/questions/18373500/how-to-check-if-command-executed-with-ipcopen3-is-hung
   # And http://www.perlmonks.org/?node_id=150748

   my $expected_number = 3;
   die "runStep: Needs $expected_number arguments\n" if  (scalar(@_)!=$expected_number);

   my ($test_name, $command, $run_die) = @_;

   print STDERR "\nSubmitting command:\n$command\n";

   my ($proc_in,$proc_out,$proc_err);
   $proc_err = gensym;                                                                    # Vivified - STDERR can't use autovivified filehandle.

   my $pid = open3($proc_in, $proc_out, $proc_err, $command);                             # Open connections
   print STDERR "\nIN: $proc_in    OUT: $proc_out    ERR: $proc_err\nPID: $pid\n";        # List connection information

   close($proc_in) || die "STDIN feed would not close\n";                                 # Close standard in as input provided on the command line

   #my @outlines = <OUT>;
   #my @errlines = <ERR>;

   my $p_output = "";
   my $p_error = "";
   my $print_count = 0;
   print STDERR "\nKey:\nSTDERR read - .\nSTDOUT read - *\n";                             # Will track reading STDOUT and STDIN in case of hanging

   my $handles = IO::Select->new($proc_out, $proc_err);                                   # Prepare file handles to be tracked
   while (my @prepped_handles = $handles->can_read(150)) {                                # Wait for output (timeout)... If timeout is reached without output
                                                                                             # buffer clearing will cease and limit will eventually be reached
      for my $file_handle (@prepped_handles){
         my $buffer_length = sysread($file_handle, my $buffer, 4096);                     # Read from the prepped filehandles
         if(not defined $buffer_length){
            die "Error reading child output: $!\n";
         }elsif(length($buffer_length ==0)){                                              # Close a redundant file handle
            #print STDERR "\n\nRemoving file handle: $file_handle\n\n";
            $handles->remove($file_handle);
         }else{
            if($file_handle == $proc_out){                                                # Read filehandle and concatenate to collected output
               $print_count++;
               $p_output .= $buffer;

               if($print_count%10){                                                       # Print upon each read to track progress
                  print STDERR "*";
               }else{
                  my @time = localtime(time);
                  print STDERR "*\t$time[2]:$time[1]:$time[0]\n";
               }

            }elsif($file_handle == $proc_err){
               $print_count++;
               $p_error .= $buffer;

               if($print_count%10){
                  print STDERR ".";
               }else{
                  my @time = localtime(time);
                  print STDERR ".\t$time[2]:$time[1]:$time[0]\n";
               }

            }else{
               die "What is $file_handle\n";
            }
         }
      }
   }

   print STDERR "\nFinished collecting. Waiting for process to end\n";
   #print STDERR `ps`;

   waitpid($pid, 0);                                                                      # Wait for process to complete 

   my $exit_code = $?;

   print STDERR "\nProcess complete\n";
   #print STDERR `ps`;
   print STDERR "\n";

   close($proc_out) || die "STDOUT feed would not close\n";
   close($proc_err) || die "STDERR feed would not close\n";

   #print STDERR "out: $p_output\n";
   #print STDERR "err: $p_error\n";

   print "ERROR\n$p_error\n" if ($exit_code);                                             # If pocess fails print errors messages

   if($run_die eq "run"){
      is ($exit_code, 0, "$test_name exit code");
   }elsif($run_die eq "die"){
      isnt ($exit_code, 0, "$test_name exit code - non-zero");
   }else{
      die "runStep 'run_die' must be run/die\n";
   }

   return($exit_code);
}




