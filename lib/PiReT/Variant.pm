#!/usr/bin/perl

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../";
use PiReT::Map;

=head1 NAME

PiReT::Variant

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

This calls variants from mapped reads

Perhaps a little code snippet.

    use PiReT::Variant;

    ...

=head1 SUBROUTINES/METHODS

=head2 function1

=cut


package Variant;
use Exporter;

################################################################################
sub call_SNPs {
    my %args             = @_;
    my $ref = $args{ref};
    my $aln_file = $args{aln_file};
    my $bc_file = $args{bc_file};
    my $vc_output = $args{vc_output};

    my $variant_call = "samtools mpileup -go "
                        ."$bc_file -f $ref $aln_file ";

    my $bcf_call = "bcftools call -vmO z -o $vc_output $bc_file";

    Map::executeCommand($variant_call);
    Map::executeCommand($bcf_call)
}
################################################################################

################################################################################
1;    # End of PiReT::Variant