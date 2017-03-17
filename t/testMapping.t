#!/usr/bin/env perl  

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../lib";
# use Map::Mapping;
# use Test::Simple tests =>1;
use Test::More;
# use Test::File;

# Verify module can be included via "use" pragma
BEGIN { use_ok('Map::Mapping') };

# Verify module can be included via "require" pragma
require_ok( 'Map::Mapping' );

# Verify if mapping produces an expected sam file
Mapping::runMapping("data/1.trimmed.fastq", "data/2.trimmed.fastq",
                     "--fast", "data/euk_prok_index", 1,
                    "data/splice_sites_gff.txt", "results/mapped.sam",
                    "results/mapped.log");

my $map_cnt = count_lines("results/mapped.sam");
is($map_cnt, 50, "runMapping() IS test");

# Verify if reads that did not map are being parsed correctly
open( my $unmap1, ">", "results/unmap_1.fastq" ) or die "Damn. $!";
open( my $unmap2, ">", "results/unmap_2.fastq" ) or die "Damn. $!";
open(FH, "results/mapped.sam") or die "Damn. $!";
    while (<FH>) {
        chomp;
        next if (/^\@/);
        my $samline = $_;
        Mapping::parsePairedUnmapped($samline, $unmap1, $unmap2);
    }
close $unmap1;
close $unmap2;
close FH;

my $un_cnt1 = count_lines("results/unmap_1.fastq");
my $un_cnt2 = count_lines("results/unmap_2.fastq");

is($un_cnt1, 64, "parsePairedUnmapped() read 1 test");
is($un_cnt2, 64, "parsePairedUnmapped() read 2 test");





# Function to count lines
sub count_lines {
    my $fn = shift;

    my $cnt;
    open(FH, $fn) or die "Damn. $!";
    $cnt++ while <FH>;
    close FH;
    return $cnt;
}


done_testing();