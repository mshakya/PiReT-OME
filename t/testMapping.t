#!/usr/bin/env perl  

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use lib "$Bin/../ext/lib/perl5";
use Test::More;
use Cwd;

# use Test::File;

# Verify module can be included via "use" pragma
BEGIN { use_ok('PiReT::Map') }

# Verify module can be included via "require" pragma
require_ok('PiReT::Map');

# Verify if indexes are created correctly
Map::createHisatIndex(
    f1        => "t/data/euk_test.fna",
    numCPU    => 1,
    out_index => "t/results/euk_index"
);

my $in_cnt = count_lines("t/results/euk_index.6.ht2l");
is( $in_cnt, 1852, "createHisatIndex() IS test" );

# Verify if mapping produces an expected sam file
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

my $map_cnt = count_lines("t/results/mapped.sam");
is( $map_cnt, 64, "runMapping() IS test" );

#Verify if thr createFAI is working
Map::createFAI("t/data/euk_test.fna\n");
my $fai_euk = &count_lines("t/data/euk_test.fna.fai");
is( $fai_euk, 8, "createFAI() IS test with euk\n" );

# Verify if the parseFAI function is working as it should
my %fai_dic = Map::parseFAI("t/data/euk_test.fna.fai");
my $ind     = $fai_dic{"gi|347623741|ref|NT_175993.1|"};
is( $ind, '217846', "parseFAI() IS test" );

# Verify if order sam is working correctly
Map::orderSAM(
    sam_file   => "t/results/mapped.sam",
    bam_file => "t/results/ordered_mapped.bam"
);
my $bam_cnt = count_lines("t/results/ordered_mapped.bam");
is( $bam_cnt, '26', "orderSAM() IS test" );

# Verify if bwa index is working
my $full_ref_path = Cwd::abs_path("t/data/Sa_cervi.R64.fa");
my $full_out_path = Cwd::abs_path("t/results");
my $symlink = "ln -fs $full_ref_path $full_out_path/ref_bwa.fna";
Map::executeCommand($symlink);
Map::createBWAIndex(
    ref => "t/results/ref_bwa.fna",
    );

my $bwa_index_ann = count_lines("t/results/ref_bwa.fna.ann");
is( $bwa_index_ann, '35', "createBWAIndex() IS test (.ann file)" );
my $bwa_index_amb = count_lines("t/results/ref_bwa.fna.amb");
is( $bwa_index_amb, '1', "createBWAIndex() IS test (.amb file)" );

#Verify bwa mapping
Map::runBWAmem(
    ref => "t/results/ref_bwa.fna",
    r1 => "t/data/yeast1.fastq",
    r2 => "t/data/yeast2.fastq",
    outsam => "t/results/bwa_mapped.sam");

my $bwa_mapped_line = count_lines("t/results/bwa_mapped.sam");
is( $bwa_mapped_line, '25018', "runBWAmem() IS test" );

# Function to count lines
sub count_lines {
    my $fn = shift;
    my $cnt;
    open( my $FH, '<', $fn ) or die "Damn. $!";
    $cnt++ while <$FH>;
    close $FH;
    return $cnt;
}

done_testing();
