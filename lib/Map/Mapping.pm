#!/usr/bin/perl

package Mapping;
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
            $command = "hisat2 $hisat2options ",
                "-p $numCPU -x $IndexFile ",
                "-1 $queryPairedFile1 ",
                "-2 $queryPairedFile2 ",
                "2>$mappingLogFile > ",
                "$outsam \n";
            print $command;
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
    my $fh          = shift;
    my $pair_type   = shift;
    my (@samFields) = @_;

    my $header = @samFields[0];
    $header =~ s/\/\d$//;
    my $seq  = @samFields[9];
    my $qual = @samFields[10];

    print $fh "@" . $header . "/$pair_type\n" . $seq . "\n+\n" . $qual . "\n";
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

sub parseSingleUnMapped {

    # Single reads that are unmapped
    my %args               = @_;
    my $samline            = $args{samline};
    my $unMappedForward_fh = $args{unMapFwd};
    my $unMappedReverse_fh = $args{unMapRev};

    my @samFields = split /\t/, $samline;

    if ( ( $samFields[1] & 4 ) and ( ( $samFields[1] & 8) eq 0 ) ) {
        if ( $samFields[1] & 64 ) {
            &sam2fastq( $unMappedForward_fh, 1, @samFields );
        }
        # if given segment is reverse read
        elsif ( $samFields[1] & 128 ) {
                &sam2fastq( $unMappedReverse_fh, 2, @samFields );
            }
    }
    
    if ( ( $samFields[1] & 8 ) and ( ($samFields[1] & 4 eq 0 ) ) ) {
        if ( $samFields[1] & 64 ) {
            &sam2fastq( $unMappedForward_fh, 1, @samFields );
        }
        # if given segment is reverse read
        elsif ( $samFields[1] & 128 ) {
                &sam2fastq( $unMappedReverse_fh, 2, @samFields );
        }
    }

}
################################################################################

sub parseSingleMapped {
    my $samFields               = shift;
    my $unMappedNonPairMate1_fh = shift;
    my $unMappedNonPairMate2_fh = shift;

    if ( $samFields[1] & 4 ) {
        if ( $samFields[1] & 64 ) {
            &sam2fastq( $samFields, $unMappedNonPairMate1_fh );
            if ( $samFields[1] & 128 ) {
                &sam2fastq( $samFields, $unMappedNonPairMate2_fh );
            }
        }
    }
}

################################################################################

sub parsePairedMapped {

    # reads (both) that are correctly mapped (correct insertion size ~500)
    my %args      = @_;
    my $samline   = $args{samline};
    my $paired_fh = $args{paired_out};

    my @samFields = split /\t/, $samline;
    if ( $samFields[1] & 2 ) {

        #mapped within the correct insert size
        if ( $samFields[1] & 16 ) {
            if ( $samFields[1] & 64 ) {
                print $paired_fh "$samline\n";
            }
            elsif ( $samFields[1] & 128 ) {
                print $paired_fh "$samline\n";
            }
        }
        elsif ( $samFields[1] & 32 ) {
            if ( $samFields[1] & 64 ) {
                print $paired_fh "$samline\n";
            }
            elsif ( $samFields[1] & 128 ) {
                print $paired_fh "$samline\n";
            }
        }
    }
}

################################################################################

sub countParsedFiles {

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
