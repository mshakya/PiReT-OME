#!/usr/bin/perl

package Mapping;
use Exporter;
@ISA    = ('Exporter');
@EXPORT = ('runMapping');

sub runMapping {

    # inputs
    my $IndexFile          = shift;
    my $queryPairedFile_r1 = shift;
    my $queryPairedFile_r2 = shift;
    my $queryUnpairedFile  = shift;
    my $outputDir          = shift;    #output directory

    my $mappingLogFile         = "$outputDir/$prefix.mapping.log";
    my $statsfile              = "$outputDir/$prefix.stats.text";
    my $unalignedNonPairedFile = "$outputDir/unmapped.$prefix.unpaired.fastq";
    my $unalignedMate1File     = "$outputDir/unmapped.$prefix.1.fastq";
    my $unalignedMate2File     = "$outputDir/unmapped.$prefix.2.fastq";
    my $alignedNonPairedFile   = "$outputDir/unpaired.ref.fastq";
    my $refIdList              = "$outputDir/$prefix.refId.txt";
    my $numUnmappedPairedFile         = 0;    # non ref reads number
    my $numUnmappedUnpairedFile       = 0;    # non ref reads number
    my $numTotalReadsPaired           = 0;
    my $numTotalUnmappedReadsPaired   = 0;
    my $numTotalReadsUnpaired         = 0;
    my $numTotalUnmappedReadsUnpaired = 0;

    #sf
    my $numMapped = 0;
    my %eukarya_reads;
    my %prokaryote_reads;
    my %pair_prokaryote_reads;
    my %pair_eukarya_reads;
    my $eukarya_reads         = 0;
    my $prokaryote_reads      = 0;
    my $pair_prokaryote_reads = 0;
    my $pair_eukarya_reads    = 0;

    #outputs
    open( OUTNON, ">$workdir/$sample/mapping_results/paired.eukarya_ref.sam" )
        or die
        "$! cannot open $workdir/$sample/mapping_results/paired.eukarya_ref.sam\n";
    open( OUTALL,
        ">$workdir/$sample/mapping_results/paired.prokaryote_ref.sam" )
        or die
        "$! cannot open $workdir/$sample/mapping_results/paired.prokaryote_ref.sam\n";
    open( BADMAPEU,
        ">$workdir/$sample/mapping_results/Notproperpaired.eukarya_ref.sam" )
        or die
        "$! cannot open $workdir/$sample/mapping_results/Notproperpaired.eukarya_ref.sam\n";
    open( BADMAPPRO,
        ">$workdir/$sample/mapping_results/Notproperpaired.prokaryote_ref.sam"
        )
        or die
        "$! cannot open $workdir/$sample/mapping_results/Notproperpaired.prokaryote_ref.sam\n";
    open( OUTFW,
        ">$workdir/$sample/mapping_results/forward.prokaryote_ref.sam" )
        or die
        "$! cannot open $workdir/$sample/mapping_results/forward.prokaryote_ref.sam\n";
    open( OUTBW,
        ">$workdir/$sample/mapping_results/backward.prokaryote_ref.sam" )
        or die
        "$! cannot open $workdir/$sample/mapping_results/backward.prokaryote_ref.sam\n";
    open( OUTEUFW,
        ">$workdir/$sample/mapping_results/forward.eukarya_ref.sam" )
        or die
        "$! cannot open $workdir/$sample/mapping_results/forward.eukarya_ref.sam\n";
    open( OUTEUBW,
        ">$workdir/$sample/mapping_results/backward.eukarya_ref.sam" )
        or die
        "$! cannot open $workdir/$sample/mapping_results/backward.eukarya_ref.sam\n";

    #unlinked unPaired files
    unlink $unalignedNonPairedFile if ( -s $unalignedNonPairedFile );

    #print in screen the status
    print colored (
        "Running reads mapping to reference sequence ... and log file: $mappingLogFile",
        'yellow'
        ),
        "\n";

    # if reads are paired
    if ( $queryPairedFile_r1 && $queryPairedFile_r2 ) {
        open( my $unalignedMate1_fh, ">$unalignedMate1File" )
            or die "$! $unalignedMate1File";
        open( my $unalignedMate2_fh, ">$unalignedMate2File" )
            or die "$! $unalignedMate2File";
        open( my $unalignedNonPaired_fh, ">$unalignedNonPairedFile" )
            or die "$! $unalignedNonPairedFile";
        open( my $refId_fh, ">$refIdList" )
            or die "$refIdList $!\n"
            if ($outputHost);

        # spliting paired files
        my @queryPairedFile1 = split /\,/, $queryPairedFile_r1;
        my @queryPairedFile2 = split /\,/, $queryPairedFile_r2;

        # loop through query paired file
        for ( my $i = 0; $i <= $#queryPairedFile1; $i++ ) {
            my $command;
            my $queryPairedFile1 = $queryPairedFile1[$i];
            my $queryPairedFile2 = $queryPairedFile2[$i];

            # if splice site file is present
            if (-e "$workdir/differential_gene/eukarya/splice_sites_gff.txt" )
            {
                &executeCommand(
                    "hisat2 $Bowtie2Opts --known-splicesite-infile $workdir/differential_gene/eukarya/splice_sites_gff.txt -p $numCPU -x $IndexFile -1 $queryPairedFile1 -2 $queryPairedFile2  2>$mappingLogFile > $workdir/$sample/mapping_results/mapped.sam"
                );
            }
            else {
                $command
                    = "hisat2 $Bowtie2Opts -p $numCPU -x $IndexFile -1 $queryPairedFile1 -2 $queryPairedFile2  2>$mappingLogFile > $workdir/$sample/mapping_results/mapped.sam ";
                &executeCommand($command);
            }

            # open more files
            open( my $fh, "$workdir/$sample/mapping_results/mapped.sam" )
                or die "$! hisat2 $Bowtie2Opts failed\n";
            while (<$fh>) {
                chomp;
                next if (/^\@/);
                my @samFields = split /\t/, $_;
                my $samline = $_;
                if ( $samFields[10] eq "*" and !$outFasta ) {
                    $samFields[10] = "f" x length( $samFields[9] );
                }

                # bit operation [and] on the flag
                if ( ( $samFields[1] & 4 ) and ( $samFields[1] & 8 ) ) {

                    # both paired reads unmapped
                    if ( $samFields[1] & 64 ) {

                        # the read is the first read in a pair
                        $samFields[0] =~ s/\/\d$//;
                        print $unalignedMate1_fh "@"
                            . $samFields[0] . "/1\n"
                            . $samFields[9] . "\n+\n"
                            . $samFields[10] . "\n";
                    }
                    if ( $samFields[1] & 128 ) {

                        # the read is the second read in a pair
                        $samFields[0] =~ s/\/\d$//;
                        print $unalignedMate2_fh "@"
                            . $samFields[0] . "/2\n"
                            . $samFields[9] . "\n+\n"
                            . $samFields[10] . "\n";
                    }
                    $numUnmappedPairedFile++;
                }
                elsif ( $samFields[1] & 4 )    # query is unmapped
                {
                    if ( $samFields[1] & 64 ) {
                        $samFields[0] =~ s/\/\d$//;
                        if ($outFasta) {

  # print $unalignedNonPaired_fh  ">".$samFields[0]."/1\n".$samFields[9]."\n";
                        }
                        else {
                            print $unalignedNonPaired_fh "@"
                                . $samFields[0] . "/1\n"
                                . $samFields[9] . "\n+\n"
                                . $samFields[10] . "\n";
                        }
                    }
                    if ( $samFields[1] & 128 ) {
                        $samFields[0] =~ s/\/\d$//;
                        print $unalignedNonPaired_fh "@"
                            . $samFields[0] . "/2\n"
                            . $samFields[9] . "\n+\n"
                            . $samFields[10] . "\n";
                    }
                    $numUnmappedPairedFile++;
                }
                else    #mapped reads
                {

                    if ( $seqln{ $samFields[2] } ) {
                        if (   $samFields[1] == 99
                            || $samFields[1] == 147
                            || $samFields[1] == 83
                            || $samFields[1] == 163 )
                        {
                            print OUTALL "$samline\n";
                            $pair_prokaryote_reads++;
                            $pair_prokaryote_reads{ $samFields[2] }++;
                        }    #  print OUT "$samline\n";}
                        else { print BADMAPPRO "$samline\n"; }
                        if ( $samFields[1] == 99 || $samFields[1] == 147 ) {
                            if ( $testcoverage == 2 ) {
                                print OUTFW "$samline\n";
                            }
                        }
                        if ( $samFields[1] == 83 || $samFields[1] == 163 ) {
                            if ( $testcoverage == 2 ) {
                                print OUTBW "$samline\n";
                            }
                        }
                    }
                    else {
                        if (   $samFields[1] == 99
                            || $samFields[1] == 147
                            || $samFields[1] == 83
                            || $samFields[1] == 163 )
                        {
                            print OUTNON "$samline\n";
                            $pair_eukarya_reads++;
                            $pair_eukarya_reads{ $samFields[2] }++;
                        }
                        else { print BADMAPEU "$samline\n"; }
                        if ( $samFields[1] == 99 || $samFields[1] == 147 ) {
                            if ( $testcoverage == 2 ) {
                                print OUTEUFW "$samline\n";
                            }
                        }
                        if ( $samFields[1] == 83 || $samFields[1] == 163 ) {
                            if ( $testcoverage == 2 ) {
                                print OUTEUBW "$samline\n";
                            }
                        }
                    }

                }
            }
            close $fh;
            my @statReadspaired = &parseMappingLog($mappingLogFile);
            $numTotalReadsPaired         += $statReadspaired[0];
            $numTotalUnmappedReadsPaired += $statReadspaired[1];
            die
                "\nERROR: $queryPairedFile1 and $queryPairedFile2 Paired reads have different names.\n"
                if (
                `grep "paired reads have different names" $mappingLogFile`);
        }    # end foreach
        close $unalignedMate1_fh;
        close $unalignedMate2_fh;
        close $unalignedNonPaired_fh;
        close BADMAPEU;
        close BADMAPPRO;
    }

    close OUTFW;
    close OUTBW;
    close OUTEUFW;
    close OUTEUBW;
    close OUTNON;
    close OUTALL;

    my $trimmedreads  = $numTotalReadsUnpaired + $numTotalReadsPaired * 2;
    my $totalrawreads = $trimmedreads;

    if ( -e "$sampleDir/trimming_results/$prefix.stats.txt" ) {
        $totalrawreads = &parsetrimmingfile(
            "$sampleDir/trimming_results/$prefix.stats.txt");
    }

    my $unmappedreads
        = $numTotalUnmappedReadsUnpaired + $numTotalUnmappedReadsPaired * 2;
    my $mappedreads = $trimmedreads - $unmappedreads;
    open( OUTSTAT, ">$statsfile" ) or die "$! cannot open $statsfile\n";
    print OUTSTAT "total_reads\t$trimmedreads\n";
    print OUTSTAT "total_Unmapped_reads\t$unmappedreads\n";
    print OUTSTAT "total_Mapped_reads\t$mappedreads\n";
    print OUTSTAT "proper_paired_prokaryote_reads\t$pair_prokaryote_reads\n";
    print OUTSTAT "proper_paired_eukarya_reads\t$pair_eukarya_reads\n";

    foreach ( sort keys %pair_prokaryote_reads ) {
        print OUTSTAT
            "proper_paried_reads_in_prokaryote_chromos:\t$_\t$pair_prokaryote_reads{$_}\n";
    }
    foreach ( sort keys %pair_eukarya_reads ) {
        print OUTSTAT
            "proper_paried_reads_in_eukarya_chromos:\t$_\t$pair_eukarya_reads{$_}\n";
    }
    close OUTSTAT;
}

1;
