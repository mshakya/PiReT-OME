#!/usr/bin/env perl  

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use lib "$Bin/../ext/lib/perl5";
use Test::More;
use PiReT::Count;
use PiReT::Map;

# Verify module can be included via "use" pragma
BEGIN { use_ok('PiReT::Count') }

# Verify module can be included via "require" pragma
require_ok('PiReT::Count');

# Verify if counting function is working
Map::createHisatIndex(
    f1        => "t/data/euk_test.fna",
    numCPU    => 1,
    out_index => "t/results/euk_index"
);

Map::runMapping(
    r1             => "t/data/1.trimmed.fastq",
    r2             => "t/data/2.trimmed.fastq",
    hisat2options  => "--fast",
    IndexFile      => "t/results/euk_index",
    numCPU         => 1,
    splicesite     => "t/data/splice_sites_gff.txt",
    outsam         => "t/results/mapped.sam",
    mappingLogFile => "t/results/mapped.log"
);

Map::orderSAM(
    sam_file   => "t/results/mapped.sam",
    bam_file => "t/results/ordered_mapped.bam"
);

Count::stringtie(
    out_summary_gtf  => "t/results/out.gtf",
    gff              => "t/data/eukarya.gff",
    out_coverage_gtf => "t/results/coveaged.gtf",
    out_abun_tab     => "t/results/out.tab",
    sample           => "test1",
    in_bam           => "t/results/ordered_mapped.bam"
);

my $tab_cnt = count_lines("t/results/out.tab");
is( $tab_cnt, '488', "stringtie() IS test for tab file" );

my $gtf_cnt = count_lines("t/results/out.gtf");
is( $gtf_cnt, '8', "stringtie() IS test for gtf file" );

# test:4, # clean up!
Map::executeCommand("rm -rf t/results/*\n");
my $file_no = &files_in_dir("t/results/");
is($file_no, 0, 'test for file clearance');

done_testing();
# Function to count lines
sub count_lines {
    my $fn = shift;
    my $cnt;
    open( my $FH, '<', $fn ) or die "Damn. $!";
    $cnt++ while <$FH>;
    close $FH;
    return $cnt;
}

sub files_in_dir {
    my $dir = shift || '.';
    my $dh;
    opendir $dh, $dir;
    grep { -f } readdir $dh;
}

