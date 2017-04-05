#!/usr/bin/env perl  

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use lib "$Bin/../ext/lib/perl5";
use Test::More;
use Test::Script;
use PiReT::Map;

$ENV{PERL5LIB} = "$Bin/../ext/lib/perl5:$ENV{PERL5LIB}";

#verify if illumina_fastq_QC is working
script_compiles('scripts/illumina_fastq_QC.pl');

my $command
    = 'perl scripts/illumina_fastq_QC.pl '
    . '-min_L 60 -n 5 -q 15 -lc 0.7 '
    . '-t 1 -prefix samp1 -d t/results/samp1/trimming_results '
    . "-p t/data/1.fastq t/data/2.fastq\n";
Map::executeCommand("mkdir -p t/results/samp1/trimming_results\n");
Map::executeCommand($command);

my $q_count = &count_lines('t/results/samp1/trimming_results/fastqCount.txt');
is( $q_count, 2, 'test for illumina_fastq_QC' );

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
