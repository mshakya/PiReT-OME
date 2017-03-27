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
sub open_samfiles {

    # function to create and open necesary files
    my $kingdom = shift;
    my $mapDir  = shift;
    if ( $kingdom eq 'both' ) {
        &open_files();
        &open_files( PR_EUK,  "$mapDir/paired.eukarya_ref.sam" );
        &open_files( PR_PROK, "$mapDir/paired.prokaryote_ref.sam" );
        &open_files( NP_EUK,  "$mapDir/Notproperpaired.eukarya_ref.sam" );
        &open_files( NP_PROK, "$mapDir/Notproperpaired.prokaryote_ref.sam" );
        &open_files( FW_EUK,  "$mapDir/forward.eukarya_ref.sam" );
        &open_files( FW_PROK, "$mapDir/forward.prokaryote_ref.sam" );
        &open_files( BW_PROK, "$mapDir/backward.prokaryote_ref.sam" );
        &open_files( BW_EUK,  "$mapDir/backward.eukarya_ref.sam" );
    }
    elsif ( $kingdom eq 'eukayra' ) {
        &open_files( PR_EUK, "$mapDir/paired.eukarya_ref.sam" );
        &open_files( NP_EUK, "$mapDir/Notproperpaired.eukarya_ref.sam" );
        &open_files( BW_EUK, "$mapDir/backward.eukarya_ref.sam" )
            & open_files( FW_EUK, "$mapDir/forward.eukarya_ref.sam" );
    }
    elsif ( $kingdom eq 'prokaryotes' ) {
        &open_files( PR_PROK, "$mapDir/paired.prokaryote_ref.sam" );
        &open_files( NP_PROK, "$mapDir/Notproperpaired.prokaryote_ref.sam" );
        &open_files( FW_PROK, "$mapDir/forward.prokaryote_ref.sam" );
        &open_files( BW_PROK, "$mapDir/backward.prokaryote_ref.sam" );
    }
}
################################################################################
sub close_samfiles {

    # close all opened samfiles
    if ( $kingdom eq 'both' ) {
        close PR_PROK;
        close PR_EUK;
        close NP_EUK;
        close NP_PROK;
        close FW_EUK;
        close FW_PROK;
        close BW_PROK;
        close BW_EUK;
    }
    elsif ( $kingdom eq 'eukayra' ) {
        close PR_EUK;
        close NP_EUK;
        close FW_EUK;
        close BW_EUK;
    }
    elsif ( $kingdom eq 'prokaryotes' ) {
        close PR_PROK;
        close NP_PROK;
        close FW_PROK;
        close BW_PROK;
    }
}
################################################################################
sub open_files {

    # sub routine to open files
    my $fh = shift;    # file handle
    my $fn = shift;    # full path to filename
    open( $fh, ">$fn" )
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

sub sam2fastq {

    # takes a samfile line and converts to fastq
    # reverse complements if its reported mapped in samfile
    my $fh          = shift;
    my $pair_type   = shift;
    my (@samFields) = @_;

    my $header = @samFields[0];
    $header =~ s/\/\d$//;
    my $seq      = @samFields[9];
    my $qual     = @samFields[10];
    my $ref      = $samFields[2];
    my $map_qual = $samFields[4];

    if ( $pair_type eq 1 ) {
        print $fh "@"
            . $header
            . "/$pair_type\n"
            . $seq . "\n+\n"
            . $qual . "\n";
    }
    elsif ( ( $pair_type eq 2 ) and ( $ref eq '*' ) ) {
        print $fh "@"
            . $header
            . "/$pair_type\n"
            . $seq . "\n+\n"
            . $qual . "\n";
    }
    elsif ( ( $pair_type eq 2 ) and ( $ref ne '*' ) ) {
        if ( $map_qual eq 0 ) {
            print $fh "@"
                . $header
                . "/$pair_type\n"
                . $seq . "\n+\n"
                . $qual . "\n";
        }
        else {
            my $rc_seq = rev_comp($seq);
            print $fh "@"
                . $header
                . "/$pair_type\n"
                . $rc_seq . "\n+\n"
                . $qual . "\n";
        }
    }

}

################################################################################

sub parsePairedUnmapped {

    # reads (both) that are not mapped
    my $samline          = shift;
    my $unMappedMate1_fh = shift;
    my $unMappedMate2_fh = shift;

    my @samFields = split /\t/, $samline;
    if ( ( $samFields[1] & 4 ) and ( $samFields[1] & 8 ) ) {

        # write out the read pair 1
        if ( $samFields[1] & 64 ) {
            &sam2fastq( $unMappedMate1_fh, 1, @samFields );
        }

        # write out the read pair 2
        if ( $samFields[1] & 128 ) {
            &sam2fastq( $unMappedMate2_fh, 2, @samFields );
        }
    }
}

################################################################################

sub parseSingles {

    # Single reads that are unmapped
    my %args               = @_;
    my $samline            = $args{samline};
    my $unMappedForward_fh = $args{unMapFwd};
    my $unMappedReverse_fh = $args{unMapRev};
    my $MappedForward_fh   = $args{MapFwd};
    my $MappedReverse_fh   = $args{MapRev};

    my @samFields = split /\t/, $samline;

    if ( ( $samFields[1] & 4 ) and ( ( $samFields[1] & 8 ) eq 0 ) ) {
        if ( $samFields[1] & 64 ) {
            &sam2fastq( $unMappedForward_fh, 1, @samFields );
        }

        # if given segment is reverse read
        elsif ( $samFields[1] & 128 ) {
            &sam2fastq( $unMappedReverse_fh, 2, @samFields );
        }
    }

    if ( ( $samFields[1] & 8 ) and ( ( $samFields[1] & 4 ) eq 0 ) ) {
        if ( $samFields[1] & 64 ) {
            print $MappedForward_fh "$samline\n";
        }

        # if given segment is reverse read
        elsif ( $samFields[1] & 128 ) {
            print $MappedReverse_fh "$samline\n";
        }
    }

}

################################################################################

sub parsePairedMapped {

    # reads (both) that are correctly mapped (correct insertion size ~500)
    my %args         = @_;
    my $samline      = $args{samline};
    my $paired_fh    = $args{paired_out};
    my $nonproper_fh = $args{nonproper};

    my @samFields = split /\t/, $samline;
    if ( $samFields[1] & 2 ) {    # read mapped in proper pair
        if ( $samFields[1] & 32 ) {    # other mate is reverse complemented
            if ( $samFields[1] & 64 ) {    # first read
                print $paired_fh "$samline\n";
            }
        }
        elsif ( $samFields[1] & 16 ) {     # reverse complemented
            if ( $samFields[1] & 128 ) {    # second read
                print $paired_fh "$samline\n";
            }
        }
    }
    elsif (
        ( ( $samFields[1] & 2 ) eq 0 )      # read not mapped in proper pair
        and ( ( $samFields[1] & 4 ) eq 0 )  # first read mapped
        and ( ( $samFields[1] & 8 ) eq 0 )
        )                                   # second read mapped
    {
        print $nonproper_fh "$samline\n";

    }
}

################################################################################

sub parseMapping {

    ############################################################################
  #------------------------------FOR REFERENCE-------------------------------#
  # 1. QNAME Query template/pair NAME
  # 2. FLAG bitwise FLAG
  # 3. RNAME Reference sequence NAME
  # 4. POS 1-based leftmost POSition/coordinate of clipped sequence
  # 5. MAPQ MAPping Quality (Phred-scaled)
  # 6. CIGAR extended CIGAR string
  # 7. MRNM Mate Reference sequence NaMe (‘=’ if same as RNAME)
  # 8. MPOS 1-based Mate POSistion
  # 9. LEN inferred Template LENgth (insert size)
  # 10. SEQ query SEQuence on the same strand as the reference
  # 11. QUAL query QUALity (ASCII-33 gives the Phred base quality)
  # 12. OPT variable OPTional fields in the format TAG:VTYPE:VALUE
    ############################################################################
  #------------------------------FLAG SCORE----------------------------------#
  # Bit     Description
  # 1       template having multiple segments / pair in sequencing
  # 2       each segment properly aligned according to the aligner
  # 4       segment unmapped
  # 8       next segment in the template unmapped
  # 16      SEQ being reverse complemented
  # 32      SEQ of the next segment in the template being reverse complemented
  # 64      the first segment in the template
  # 128     the last segment in the template
  # 256     secondary alignment
  # 512     not passing filters, such as platform/vendor quality controls
  # 1024    PCR or optical duplicate
  # 2048    supplementary alignment
    ############################################################################
    my $workdir = shift;
    my $sample  = shift;
    &open_files( fh, "$workdir/$sample/mapping_results/mapped.sam" );
    while (<$fh>) {
        chomp;
        next if (/^\@/);
        my @samFields = split /\t/, $_;
        my $samline = $_;

        # when no quality information is stored
        if ( $samFields[10] eq "*" and !$outFasta ) {
            $samFields[10] = "f" x length( $samFields[9] );
        }

        &parsePairedUnmapped( $samFields, $unMappedMate1_fh,
            $unMappedMate2_fh );

        &parseNonPairedUnMapped( $samFields, $unMappedNonPairMate1_fh,
            $unMappedNonPairMate2_fh );

        &parsePairedMapped(
            $samFields,    $paired_fh,  $unpaired_fh,
            $testcoverage, $forward_fh, $reverse_fh
        );
    }
    close $fh;
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

sub sumMaps {
    my %args         = @_;
    my $fq           = $args{fastq};
    my $pro_paired   = $args{ppm};
    my $nopro_paired = $args{npm};
    my $fwd_map      = $args{fm};
    my $rev_map      = $args{rm};
    my $unmap        = $args{um};
    my $out_table    = $args{out};

    my $fq_reads   = &count_lines($fq);
    my $ppm_count  = &count_lines($pro_paired);
    my $npm_count  = &count_samreads($nopro_paired);
    my $fwd_count  = &count_lines($fwd_map);
    my $rev_count  = &count_lines($rev_map);
    my $un_count   = &count_lines($unmap);
    my %pair_reads = &parseMappedChromo($pro_paired);

    open( OUTSTAT, ">$out_table" ) or die "$! cannot open $out_table\n";
    print OUTSTAT "total_trimmed_reads\t" . int( $fq_reads / 4 ) . "\n";
    print OUTSTAT "total_unmapped_reads\t" . int( $un_count / 4 ) . "\n";
    print OUTSTAT "forward_only_mapped_reads\t" . int($fwd_count) . "\n";
    print OUTSTAT "reverse_only_mapped_reads\t" . int($rev_count) . "\n";
    print OUTSTAT "total_Mapped_reads\t"
        . int( $ppm_count + $fwd_count + $rev_count ) . "\n";
    print OUTSTAT "Not_proper_reads\t" . "$npm_count" . "\n";
    print OUTSTAT "proper_paired_reads\t" . int( $ppm_count / 2 ) . "\n";

    foreach ( sort keys %pair_reads ) {
        print OUTSTAT
            "proper_paried_reads_in_chromos:\t$_\t$pair_reads{$_}\n";
    }
    close OUTSTAT;
}

################################################################################
sub parseMappedChromo {
    my $fn         = shift;
    my $pair_reads = 0;
    my %pair_reads;
    open( FH, "$fn" ) or die "Damn. $!";
    while (< FH >) {
        $pair_reads++;
        my $samline = $_;
        my @samFields = split /\t/, $samline;
        $pair_reads{ $samFields[2] }++;
    }
    close FH;
    return %pair_reads;
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
sub count_samreads {
    my $fn = shift;
    my $npr;
    my $command = "awk -F'\\t' '{print \$1}'" . "< $fn | uniq | wc -l \n";
    $npr = executeCommand($command);
    return chomp($npr);
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
    &executeCommand($command)
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

1; # End of PiReT::Map