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
BEGIN { use_ok('PiReT::Variant') }

# Verify module can be included via "require" pragma
require_ok('PiReT::Variant');

# Verify if bwa index is working
my $full_ref_path = Cwd::abs_path("t/data/Sa_cervi.R64.fa");
my $full_out_path = Cwd::abs_path("t/results");
my $symlink       = "ln -fs $full_ref_path $full_out_path/ref_bwa.fna";
Map::executeCommand($symlink);
Map::createBWAIndex( ref => "t/results/ref_bwa.fna", );

Map::runBWAmem(
    ref    => "t/results/ref_bwa.fna",
    r1     => "t/data/yeast1.fastq",
    r2     => "t/data/yeast2.fastq",
    outsam => "t/results/bwa_mapped.sam"
);

Map::orderSAM(
    sam_file => "t/results/bwa_mapped.sam",
    bam_file => "t/results/bwa_ordered_mapped.bam"
);

Variant::call_SNPs(
    ref       => "t/results/ref_bwa.fna",
    aln_file  => "t/results/bwa_ordered_mapped.bam",
    bc_file   => "t/results/bwa_bcfile.bcfile",
    vc_output => "t/results/bwa_vcf.vcf"
);

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
