#!/usr/bin/env perl  

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Test::More;
use Test::Script;
use PiReT::Map;

$ENV{PERL5LIB} = "$Bin/../ext/lib/perl5:$ENV{PERL5LIB}";

#verify if trim_readmapping.sh is working
script_compiles('scripts/illumina_fastq_QC.pl');

my $command
    = 'perl scripts/illumina_fastq_QC.pl '
    . '-min_L 60 -n 5 -q 15 -lc 0.7 '
    . '-t 1 -prefix samp1 -d t/results/samp1/trimming_results '
    . '-p t/data/1.fastq t/data/2.fastq';
Map::executeCommand("mkdir -p t/results/samp1/trimming_results");
Map::executeCommand($command);

my $q_count = &count_lines('t/results/samp1/trimming_results/fastqCount.txt');
is( $q_count, 2, 'test for illumina_fastq_QC.pl' );

# verify if reads_mapping.pl compiles
script_compiles('scripts/reads_mapping.pl');
my $arg
    = '-P1 t/data/1.trimmed.fastq '
    . '-P2 t/data/2.trimmed.fastq -I t/data/euk_prok_index -K both '
    . '-E t/data/euk_test.fna -B t/data/prok_test.fna -W t/results '
    . '-S samp1';

my @arg = [
    '-P1 t/data/1.trimmed.fastq ',
    '-P2 t/data/2.trimmed.fastq',
    '-I t/data/euk_prok_index',
    '-K both',
    '-E t/data/euk_test.fna -B t/data/prok_test.fna -W t/results ',
    '-S samp1'
];

# verify if script runs without error
script_runs( [ 'scripts/reads_mapping.pl', $arg ], 'test reads_mapping.pl' );

my $map_command
    = 'perl scripts/reads_mapping.pl '
    .'-P1 t/data/1.trimmed.fastq '
    .'-P2 t/data/2.trimmed.fastq '
    .'-I t/data/euk_prok_index '
    .'-K both -E t/data/euk_test.fna '
    .'-B t/data/prok_test.fna '
    .'-W t/results -S samp1';
# Map::executeCommand("mkdir -p t/results/samp1/trimming_results");
Map::executeCommand($map_command);

my $pro_prk = &count_lines('t/results/samp1/mapping_results/paired_prok.sam');
is( $pro_prk, 8, 'test for reads_mapping.pl' );

# Function to count lines
sub count_lines {
    my $fn = shift;

    my $cnt;
    open( FH, $fn ) or die "Damn. $!";
    $cnt++ while <FH>;
    close FH;
    return $cnt;
}

done_testing();
