#!/usr/bin/perl

=head1 NAME

PiReT::Map - The great new PiReT::Map!

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

Quick summary of what the module does.

Perhaps a little code snippet.

    use PiReT::Map;

    my $foo = PiReT::Map->new();
    ...

=head1 SUBROUTINES/METHODS

=head2 function1

=cut

package Map;
use Exporter;
@ISA    = ('Exporter');
@EXPORT = ('runMapping');

################################################################################
sub open_files {

    # sub routine to open files
    my $fh = shift;    # file handle
    my $fn = shift;    # full path to filename
    open( $fh, q{>}, $fn )
        or die "$! cannot open $fn\n";
}

################################################################################

sub createHisatIndex {
    my %args      = @_;
    my $fasta1    = $args{f1};
    my $fasta2    = $args{f2};
    my $numCPU    = $args{numCPU} || 1;
    my $out_index = $args{out_index};

    if ( not defined $fasta2 ) {
        $command
            = "hisat2-build -p $numCPU "
            . "--large-index -q "
            . "$fasta1 $out_index \n";
    }
    else {
        $command
            = "hisat2-build -p $numCPU "
            . "--large-index -q "
            . "$fasta1,$fasta2 $out_index \n";
    }
    &executeCommand($command);
}

################################################################################
sub runMapping {
    my %args               = @_;
    my $queryPairedFile_r1 = $args{r1};
    my $queryPairedFile_r2 = $args{r2};
    my $hisat2options      = $args{hisat2options};
    my $IndexFile          = $args{IndexFile};
    my $numCPU             = $args{numCPU} || 1;
    my $splicesite         = $args{splicesite};
    my $outsam             = $args{outsam};
    my $mappingLogFile     = $args{mappingLogFile};

    my @queryPairedFile1 = split /\,/, $queryPairedFile_r1;
    my @queryPairedFile2 = split /\,/, $queryPairedFile_r2;
    for ( my $i = 0; $i <= $#queryPairedFile1; $i++ ) {
        my $command;
        my $queryPairedFile1 = $queryPairedFile1[$i];
        my $queryPairedFile2 = $queryPairedFile2[$i];

        # if splice site file is present
        if ( -e $splicesite ) {
            $command
                = "hisat2 $hisat2options "
                . "--known-splicesite-infile $splicesite "
                . "-p $numCPU "
                . "-x $IndexFile "
                . "-1 $queryPairedFile1 "
                . "-2 $queryPairedFile2 "
                . "2>$mappingLogFile > "
                . "$outsam \n";
            &executeCommand($command);
        }
        else {
            $command
                = "hisat2 $hisat2options "
                . "-p $numCPU -x $IndexFile "
                . "-1 $queryPairedFile1 "
                . "-2 $queryPairedFile2 "
                . "2>$mappingLogFile > "
                . "$outsam \n";
            &executeCommand($command);
        }

    }
}
################################################################################
sub executeCommand {
    my $command = shift;
    print $command;
    `$command`;
}

################################################################################
# Function to reverse complement DNA
sub rev_comp {
    my $DNA = shift;

    my $revcom = reverse $DNA;
    $revcom =~ tr/ACGTacgt/TGCAtgca/;
    return $revcom;
}

################################################################################
sub orderSAM {
    my %args        = @_;
    my $mappedSAM   = $args{sam_file};
    my $orderedBAMs = $args{bam_file};

    my $order_samtools
        = "samtools view -Su $mappedSAM | samtools sort > $orderedBAMs";
    &executeCommand($order_samtools);
}

################################################################################
sub parseFAI {

    # makes hash from fai files
    my $fai = shift;
    my %seqln;
    open( GENOIN, $fai ) or die "$fai does not exist $!";

    # Process
    while (<GENOIN>) {
        chomp;
        my @line = split /\s+/, $_;
        for ( my $i = 1; $i <= $line[1]; $i++ ) {
            $seqln{ $line[0] } = $line[1];
        }
    }
    close GENOIN;
    return %seqln;
}

################################################################################
# Function to count lines
sub count_lines {
    my $fn = shift;
    my $cnt;
    open( FH, $fn ) or die "Damn. $!";
    $cnt++ while <FH>;
    close FH;
    return $cnt;
}

################################################################################
sub createFAI {
    my $fn = shift;
    $command = "samtools faidx $fn";
    &executeCommand($command);
}

=head1 AUTHOR

Migun Shakya, C<< <migun at lanl.gov> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-piret-map at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=PiReT-Map>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc PiReT::Map


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=PiReT-Map>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/PiReT-Map>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/PiReT-Map>

=item * Search CPAN

L<http://search.cpan.org/dist/PiReT-Map/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2017 Migun Shakya.

This program is free software; you can redistribute it and/or modify it
under the terms of the the Artistic License (2.0). You may obtain a
copy of the full license at:

L<http://www.perlfoundation.org/artistic_license_2_0>

Any use, modification, and distribution of the Standard or Modified
Versions is governed by this Artistic License. By using, modifying or
distributing the Package, you accept this license. Do not use, modify,
or distribute the Package, if you do not accept this license.

If your Modified Version has been derived from a Modified Version made
by someone other than you, you are nevertheless required to ensure that
your Modified Version complies with the requirements of this license.

This license does not grant you the right to use any trademark, service
mark, tradename, or logo of the Copyright Holder.

This license includes the non-exclusive, worldwide, free-of-charge
patent license to make, have made, use, offer to sell, sell, import and
otherwise transfer the Package with respect to any patent claims
licensable by the Copyright Holder that are necessarily infringed by the
Package. If you institute patent litigation (including a cross-claim or
counterclaim) against any party alleging that the Package constitutes
direct or contributory patent infringement, then this Artistic License
to you shall terminate on the date that such litigation is filed.

Disclaimer of Warranty: THE PACKAGE IS PROVIDED BY THE COPYRIGHT HOLDER
AND CONTRIBUTORS "AS IS' AND WITHOUT ANY EXPRESS OR IMPLIED WARRANTIES.
THE IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
PURPOSE, OR NON-INFRINGEMENT ARE DISCLAIMED TO THE EXTENT PERMITTED BY
YOUR LOCAL LAW. UNLESS REQUIRED BY LAW, NO COPYRIGHT HOLDER OR
CONTRIBUTOR WILL BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, OR
CONSEQUENTIAL DAMAGES ARISING IN ANY WAY OUT OF THE USE OF THE PACKAGE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


=cut

1;    # End of PiReT::Map
