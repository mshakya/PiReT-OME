#!/usr/bin/env perl  

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use lib "$Bin/../ext/lib/perl5";
use Test::More;

# use Test::File;

# Verify module can be included via "use" pragma
BEGIN { use_ok('PiReT::Map') }

# Verify module can be included via "require" pragma
require_ok('PiReT::Map');

# Verify if indexes are created correctly
Map::createHisatIndex(
    f1        => "t/data/euk_test.fna",
    f2        => "t/data/prok_test.fna",
    numCPU    => 1,
    out_index => "t/results/euk_prok_index"
);

my $in_cnt = count_lines("t/results/euk_prok_index.6.ht2l");
is( $in_cnt, 2195, "createHisatIndex() IS test" );

# Verify if mapping produces an expected sam file
Map::runMapping(
    r1             => "t/data/1.trimmed.fastq",
    r2             => "t/data/2.trimmed.fastq",
    hisat2options  => "--fast",
    IndexFile      => "t/results/euk_prok_index",
    numCPU         => 1,
    splicesite     => "t/data/splice_sites_gff.txt",
    outsam         => "t/results/mapped.sam",
    mappingLogFile => "t/results/mapped.log"
);

my $map_cnt = count_lines("t/results/mapped.sam");
is( $map_cnt, 71, "runMapping() IS test" );

#Verify if thr createFAI is working
Map::createFAI("t/data/euk_test.fna\n");
my $fai_euk  = &count_lines("t/data/euk_test.fna.fai");
is( $fai_euk,  8, "createFAI() IS test with euk\n" );

# Verify if the parseFAI function is working as it should
my %fai_dic = Map::parseFAI("t/data/prok_test.fna.fai");
my $ind     = $fai_dic{"gi|50196905|ref|NC_007530.2|"};
print $ind
is( $ind, '5227419', "parseFAI() IS test" );


# Verify if parseMapEuk is working
Map::executeCommand("mkdir -p t/results/euk_map_test\n");
Map::parseMapEuk(
    mapDir  => "t/results/euk_map_test",
    samFile => "t/results/mapped.sam"
);
my $pp_euk = &count_lines("t/results/euk_map_test/paired_euk.sam");
is( $pp_euk, 14, "parseMapEuk() IS test" );

# Verify if parseMapBoth is working
Map::executeCommand("mkdir -p t/results/both_map_test\n");

my $pp_euk_both = &count_lines("t/results/both_map_test/paired_euk.sam");
is( $pp_euk_both, 14, "parseMapBoth() IS test" );


# Function to count lines
sub count_lines {
    my $fn = shift;

    my $cnt;
    open( my $FH,'<', $fn ) or die "Damn. $!";
    $cnt++ while <$FH>;
    close $FH;
    return $cnt;
}

done_testing();
