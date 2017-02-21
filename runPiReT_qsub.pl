#!/usr/bin/perl -W

use warnings;
use strict;
use Getopt::Long;
use File::Basename;
use Term::ANSIColor;
use Cwd;
use FindBin qw($Bin);
use POSIX qw(strftime);

$| = 1;    #?
$ENV{PATH} = "$Bin/bin/:$ENV{PATH}";
$ENV{PERL5LIB} = "$Bin/ext/lib/perl5:$ENV{PERL5LIB}"; 

#NOTE: need these paths to find qsub binaries
#TODO: change it so that its independent of the user system
#NOTE: also mention that qsub must be present
#NOTE: these paths should be set for any server with qsub capabalities
#$ENV{SGE_ROOT}         = "/cm/shared/apps/sge/2011.11p1"; #?
#$ENV{SGE_CELL}         = "default"; #?
#$ENV{SGE_CLUSTER_NAME} = "seqclust"; #?
#TODO: remove this later, once the this thing works
#foreach ( sort keys %ENV ) {
#    print "$_  =  $ENV{$_}\n";
#}

&checkDependedPrograms();

my $main_pid  = $$;
my $version   = "develop";
my $time      = time();
my $scriptDir = "$Bin/scripts";
my $jbBin = "$Bin/bin/JBrowse/bin"; 
my ($descriptfile,     $test,      $splice_file_out,
    $splice_file_in,   $pairfile,  $diffdir,
    $workdir,          $numCPU,    $eukarya_fasta,
    $prokaryote_fasta, $index_bt2, $gff_eukarya,
    $gff_prokaryote,   $coverage_fasta
);

#my $mapread='no';
$numCPU = 1;
$test   = 'prokaryote';
my $test_method      = 'both';
my $mptool           = 'bowtie2';
my $htseq            = 'gene';
my $p_cutoff         = 0.001;
my $memlim           = '10G';
my $rna_trimming_opt = 'yes';
my $rna_mapping_opt  = 'yes';
$gff_eukarya    = 'NONE';
$gff_prokaryote = 'NONE';
$eukarya_fasta = 'NONE';
$prokaryote_fasta = 'NONE';
$coverage_fasta = 'NONE';
$workdir = 'NONE';
$descriptfile = 'NONE';

#------------------------------------------------------------------------------#
GetOptions(
    'rna_mapping_opt=s'  => \$rna_mapping_opt,
    'rna_trimming_opt=s' => \$rna_trimming_opt,
    'E|exp=s'                 => \$descriptfile,
    'W|dir=s'                   => \$workdir,
    'P|cpu=i'                 => \$numCPU,             # bwa option
    'U|eukarya_fasta=s'       => \$eukarya_fasta,
    'R|prokaryote_fasta=s'    => \$prokaryote_fasta,
    'I|index_ref_bt2=s'       => \$index_bt2,
    'M|h_vmem=s'              => \$memlim,
    'G|gff_eukarya=s'         => \$gff_eukarya,
    'F|gff_prokaryote=s'      => \$gff_prokaryote,
    'C|gene_coverage_fasta=s' => \$coverage_fasta,
    'S|significant_pvalue=f'  => \$p_cutoff,
    'N|pair_comparison=s'     => \$pairfile,
    'A|geneopt=s'             => \$htseq
    , #count reads based on 'gene' or 'CDS' or 'tRNA' or 'mRNA' in annotation file, default ='gene';
      #  'BAM_ready=s' =>\$mapread,      #if mapping file are provided for samples by users
    'K|test_kingdom=s' => \$test,
    'T|test_method=s'  => \$test_method,
	'V|version'		 => sub{printVersion()},
    'help|?'         => sub {&Usage()}
);


#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
if ( $eukarya_fasta eq 'NONE' )    { $eukarya_fasta    = ""; } else {$eukarya_fasta=Cwd::abs_path($eukarya_fasta)};
if ( $prokaryote_fasta eq 'NONE' ) { $prokaryote_fasta = ""; } else {$prokaryote_fasta=Cwd::abs_path($prokaryote_fasta)};
if ( $gff_eukarya eq 'NONE' )      { $gff_eukarya      = ""; } else {$gff_eukarya=Cwd::abs_path($gff_eukarya)};
if ( $gff_prokaryote eq 'NONE' )   { $gff_prokaryote   = ""; } else {$gff_prokaryote=Cwd::abs_path($gff_prokaryote)};
if ( $coverage_fasta eq 'NONE' )   { $coverage_fasta   = ""; } else {$coverage_fasta=Cwd::abs_path($coverage_fasta)};

#------------------------------------------------------------------------------#

my $start_time_string = &getTmpNameByTime;

#------------------------------------------------------------------------------#

$scriptDir = Cwd::abs_path("$scriptDir");
$workdir   = Cwd::abs_path("$workdir");
$descriptfile = Cwd::abs_path("$descriptfile");
#------------------------------------------------------------------------------#

# creates a working directory
unless ( -d "$workdir" ) {
    mkdir "$workdir", 0777
        or die "failed: failed: can not make dir  $workdir $!";
}

#------------------------------------------------------------------------------#

# process log file
my $process_log_file = "$workdir/process.log";

# error log file
my $error_log_file   = "$workdir/error.log";

open( LOG, ">", $process_log_file )
    or die "failed: failed to write $process_log_file\n$!";
open( STDERR, '>&', STDOUT )
    or die "failed: failed: Can't redirect stderr: $!";
open( STDERR, '>', $error_log_file )
    or die "failed: failed: Can't redirect stderr: $!";

#------------------------------------------------------------------------------#

# if there is a waring or error write in in STDERR and lprint
$SIG{__WARN__} = sub { print STDERR @_; &lprint(@_) };
$SIG{__DIE__} = sub { print STDERR @_; &lprint(@_); exit 1 };

# &lprint("\nProject Start: $start_time_string\n");
# print LOG qx/ps -o args $$/;
# &lprint("Version: $version\n\n");

#NOTE: This looks like it reads the config.txt file, but for what??
# open( IN, "$workdir/config.txt" )
#     or die "failed: failed: not open $workdir/config.txt $!";
# my $linenumc = 0;
# while (<IN>) {
#     chomp;
#     $linenumc++;
#     my $line = $_;
#     if ( $linenumc < 28 ) {
#         &lprint("$line\n");
#     }
# }
# close IN;
# &lprint("[Checking Files]\nCheckingFiles=Always\n\n");
# &lprint(
#     "[Trimming and Mapping Reads]\nTrimmingMappingReads=$rna_trimming_opt\t$rna_mapping_opt\n\n"
# );
# &lprint(
#     "[De Novo Detection of small RNAs]\nDe Novo Detection of small RNAs =Always\n\n"
# );
# &lprint("[Differential Gene Analysis]\nDifferentialGeneAnalysis=Always\n\n");
#------------------------------------------------------------------------------#
unless ( $descriptfile
    && $workdir
    && $test
    && ( $gff_eukarya   || $gff_prokaryote )
    && ( $eukarya_fasta || $prokaryote_fasta )
    && $index_bt2 )
{
    &Usage
}

#------------------------------------------------------------------------------#
# create directory, throw error, if can't
unless ( -d "$workdir/sum_gene_count/" ) {
    mkdir "$workdir/sum_gene_count/", 0777
        or die
        "failed: failed: can not make dir  $workdir/sum_gene_count/ $!";
}
#------------------------------------------------------------------------------#
unless ( -d "$workdir/logdir/" ) {
    mkdir "$workdir/logdir/", 0777
        or die "failed: failed: can not make dir  $workdir/dir/ $!";
}
#------------------------------------------------------------------------------#
my %description;
my @colname;
my %expdescription;
my %allsample;
my @allsample;
my %pairs;
my %allrawreads1;
my %allrawreads2;

#------------------------------------------------------------------------------#

&lprint("[Checking Experimental Design File]\n Running\n\n");

#------------------------------------------------------------------------------#
# parse experimental design file
open( IN, "$descriptfile" )
    or die "failed: failed: can not open description file $descriptfile $!";
my $linenum = 0;
while (<IN>) {
    chomp;
    my @line = split /\s+/, $_;
    unless ( $#line >= 2 ) {
        &lprint("check format of the description file $_ \n");
        exit;
    }
    $linenum++;
    if ( $linenum == 1 ) {
        for ( my $i = 0; $i <= $#line; $i++ ) {
            $line[$i] =~ s/\s+//g;
            $colname[$i] = $line[$i];
        }
    }
    if ( $linenum > 1 ) {
        $allsample{ $line[0] } = 1;
        push @allsample, $line[0];
        for ( my $i = 3; $i <= $#line; $i++ ) {
            $line[$i] =~ s/\s+//g;
            $expdescription{ $line[0] }{ $colname[$i] } = $line[$i];
        }

        my @line1 = split /;/, $line[1];

#    foreach (@line1)
#  {
#    my @rawreadsfile=split /:/, $_;
#    foreach (@rawreadsfile)
#    {
#   if(&file_check($_)<0 ) { die "failed: failed: The reads file $_ doesn't exist or empty.\n";}
#    }
#  }
#
#

        for ( my $ill = 0; $ill <= $#line1; $ill++ ) {
            my @rawreadsfile = split /:/, $line1[$ill];

            push @{ $allrawreads1{ $line[0] } }, $rawreadsfile[0];
            push @{ $allrawreads2{ $line[0] } }, $rawreadsfile[1];

            for ( my $il = 0; $il <= $#rawreadsfile; $il++ ) {
                my $tmpreads = $rawreadsfile[$il];
                if ( &file_check($tmpreads) < 0 ) {
                    die "failed: The reads file $_ doesn't exist or empty.\n";
                }
                $rawreadsfile[$il] = $tmpreads;
            }
            $line1[$ill] = join ':', @rawreadsfile;
        }

        $description{ $line[0] }{group} = $line[2];
        push @{ $pairs{ $line[2] } }, $line[0];
        if ( $line[2] =~ /\s+/ ) {
            &lprint(
                "failed: please check format of experimental descpition \n group name can not have space as : $line[2]\n"
            );
            exit;
        }
        my $rawreadsline = $line1[0];
        for ( my $i = 1; $i <= $#line1; $i++ ) {
            $rawreadsline = $rawreadsline . '---ppp---' . $line1[$i];
        }
        $description{ $line[0] }{Rawreads_file} = $rawreadsline;
    }
}
close IN;
#------------------------------------------------------------------------------#

&lprint("[Checking Experimental Design File]\n Finished\n\n");

&lprint("[Checking if enough samples are present]\n Running\n\n");

#------------------------------------------------------------------------------#

# throw error if there are not enough samples
foreach ( keys %pairs ) {
    my $tmpgroup = $_;
    if ( ( ( $test_method eq 'both' ) || ( $test_method eq 'Deseq' ) )
        && $#{ $pairs{$tmpgroup} } < 2 )
    {
        &lprint(
            "failed: testing method using Deseq2 but the $tmpgroup has only $pairs{$tmpgroup} which is less than 3. please not using Deseq method or increase the number of dunplcates for $tmpgroup\n"
        );
        exit;
    }
}

#------------------------------------------------------------------------------#
&lprint("[Checking if enough samples are present]\n Finished\n\n");

#------------------------------------------------------------------------------#
my ( @eukaryagff, @prokaryotegffi, %allgff );

#------------------------------------------------------------------------------#
&lprint("[Creating additional directories]\n Running\n\n");
#------------------------------------------------------------------------------#
if ( -d "$workdir/Jbrowse/" ) {
    `rm -r $workdir/Jbrowse/`;
    mkdir "$workdir/Jbrowse/", 0777
        or die "can not make dir  $workdir/Jbrowse $!";
}
else {
    mkdir "$workdir/Jbrowse/", 0777
        or die "can not make dir  $workdir/Jbrowse $!";
}
mkdir "$workdir/Jbrowse/BigWig", 0777
    or die "can not make dir  $workdir/Jbrowse/BigWig $!";

unless ( -d "$workdir/differential_gene/" ) {
    mkdir "$workdir/differential_gene/", 0777
        or die
        "failed: failed: can not make dir  $workdir/differential_gene/ $!";
}

unless ( -d "$workdir/sum_gene_count/tmp_count/" ) {
    mkdir "$workdir/sum_gene_count/tmp_count/", 0777
        or die
        "failed: failed: can not make dir  $workdir/sum_gene_count/tmp_count/ $!";
}
unless ( -d "$workdir/sum_gene_count/read_count/" ) {
    mkdir "$workdir/sum_gene_count/read_count/", 0777
        or die
        "failed: failed: can not make dir  $workdir/sum_gene_count/read_count/ $!";
}
#------------------------------------------------------------------------------#
&lprint("[Creating additional directories]\n Finished\n\n");
#------------------------------------------------------------------------------#
my %allcontigs;

#------------------------------------------------------------------------------#
&lprint(
    "[Copying and creating indices of reference fasta files]\n Running\n\n");
#------------------------------------------------------------------------------#
if ($eukarya_fasta) {
    if ( &file_check($eukarya_fasta) < 0 ) {
        die
            "failed: failed: The eukarya_fast file $eukarya_fasta doesn't exist or empty.\n";
    }
    else {
        if ( -e "$workdir/eukarya.fa" ) { `rm $workdir/eukarya.fa`; }
        `ln -fs $eukarya_fasta $workdir/eukarya.fa`;
        `samtools faidx $workdir/eukarya.fa`;
        `$jbBin/prepare-refseqs.pl --trackLabel  DNA --seqType dna --key 'DNA+protein' --fasta  $workdir/eukarya.fa --out  $workdir/Jbrowse/`;
        my @contigs = &readfai("$workdir/eukarya.fa.fai");
        foreach (@contigs) { $allcontigs{$_}++; }
    }
}    #else {die "failed:  need eukarya sequence file\n";}
#------------------------------------------------------------------------------#
if ($prokaryote_fasta) {

    if ( &file_check($prokaryote_fasta) < 0 ) {
        die
            "failed: failed: The prokaryote_fast file $prokaryote_fasta doesn't exist or empty.\n";
    }
    else {
        if ( -e "$workdir/prokaryote.fa" ) { `rm $workdir/prokaryote.fa`; }
        `ln -fs $prokaryote_fasta $workdir/prokaryote.fa`;
		`samtools faidx $workdir/prokaryote.fa`;
        `$jbBin/prepare-refseqs.pl --trackLabel  DNA --seqType dna --key 'DNA+protein' --fasta  $workdir/prokaryote.fa --out  $workdir/Jbrowse/`;

        my @contigs = &readfai("$workdir/prokaryote.fa.fai");
        foreach (@contigs) { $allcontigs{$_}++; }
    }
}    #else  {die "failed: failed: need prokaryote sequence file\n";}

#------------------------------------------------------------------------------#
&lprint(
    "[Copying and creating indices of reference fasta files]\n Finished\n\n");
#------------------------------------------------------------------------------#

#------------------------------------------------------------------------------#
&lprint(
    "[Making sure if sequences in references are not named the same]\n Running\n\n"
);

#------------------------------------------------------------------------------#
foreach ( keys %allcontigs ) {
    if ( $allcontigs{$_} > 1 ) {
        &lprint(
            "contig name $_ is duplicated, please make sure the name of contigs is unique in reference fasta files (prokaryout and eukarya references)\n"
        );
        exit;
    }
}

#------------------------------------------------------------------------------#
&lprint(
    "[Making sure if sequences in references are not named the same]\n Finished\n\n"
);
#------------------------------------------------------------------------------#
my $pjcoverage = "$workdir/coverage.fa";
if ( -e "$workdir/coverage.fa.fai" ) { `rm $workdir/coverage.fa.fai`; }
if ($coverage_fasta) {
    if ( &file_check($coverage_fasta) < 0 ) {
        die
            "failed: The coverage_fasta file $coverage_fasta doesn't exist or empty.\n";
    }
    else {
        if ( -e $pjcoverage ) { `rm $pjcoverage`; }
        `ln -fs $coverage_fasta $pjcoverage`;
        `samtools faidx $pjcoverage`;
        my @contigs = &readfai("$workdir/prokaryote.fa.fai");
        foreach (@contigs) {
            unless ( $allcontigs{$_} ) {
                &lprint(
                    "coverage contig $_ is not part of mapping index references, please make sure it is included in either prokaryout or eukarya reference\n"
                );
                exit;
            }
        }
    }
}
else { $pjcoverage = 'NA' }

#------------------------------------------------------------------------------#
&lprint(
    "[Creating additional directories based on type of analysis]\n Running\n\n"
);
#------------------------------------------------------------------------------#

if ( $test eq 'both' || $test eq 'eukarya' ) {
    unless ( -d "$workdir/sum_gene_count/tmp_count/eukarya" ) {
        mkdir "$workdir/sum_gene_count/tmp_count/eukarya", 0777
            or die
            "failed: failed: can not make dir  $workdir/sum_gene_count/tmp_count/eukarya $!";
    }
    unless ( -d "$workdir/sum_gene_count/read_count/eukarya" ) {
        mkdir "$workdir/sum_gene_count/read_count/eukarya", 0777
            or die
            "failed: failed: can not make dir  $workdir/sum_gene_count/read_count/eukarya $!";
    }
    unless ( -d "$workdir/differential_gene/eukarya" ) {
        mkdir "$workdir/differential_gene/eukarya", 0777
            or die
            "failed: failed: can not make dir  $workdir/differential_gene/eukarya $!";
    }
    &lprint(
        "[Creating additional directories based on type of analysis]\n Finished\n\n"
    );

#------------------------------------------------------------------------------#
    &lprint("[parsing gff files]\n Running\n\n");
#------------------------------------------------------------------------------#

    my @tmpgff = split /,/, $gff_eukarya;

    foreach (@tmpgff) {
        my $tmpgff = $_;
        my @tmpeukarya = split '/', $tmpgff;
        $tmpeukarya[-1] =~ s/\.//g;
        $tmpeukarya[-1] =~ s/gff//g;
        if ( &file_check($tmpgff) < 0 ) {
            next &lprint(
                "The gff_eukarya file $tmpgff doesn't exist or empty.\n");
        }
        push @{ $allgff{eukarya} }, $tmpeukarya[-1];

       `$jbBin/flatfile-to-json.pl --gff $tmpgff --type CDS  --tracklabel CDS --out $workdir/Jbrowse`;
       `$jbBin/flatfile-to-json.pl --gff $tmpgff --type tRNA  --tracklabel tRNA --out $workdir/Jbrowse`;
       `$jbBin/flatfile-to-json.pl --gff $tmpgff --type exon  --tracklabel exon --out $workdir/Jbrowse`;
       `$jbBin/flatfile-to-json.pl --gff $tmpgff --type gene  --tracklabel gene --out $workdir/Jbrowse`;

        unless (
            -d "$workdir/sum_gene_count/tmp_count/eukarya/$tmpeukarya[-1]" )
        {
            mkdir "$workdir/sum_gene_count/tmp_count/eukarya/$tmpeukarya[-1]",
                0777
                or die
                "failed: can not make dir  $workdir/sum_gene_count/tmp_count/eukarya/$tmpeukarya[-1] $!";
        }
        unless (
            -d "$workdir/sum_gene_count/read_count/eukarya/$tmpeukarya[-1]" )
        {
            mkdir
                "$workdir/sum_gene_count/read_count/eukarya/$tmpeukarya[-1]",
                0777
                or die
                "failed: can not make dir  $workdir/sum_gene_count/read_count/eukarya/$tmpeukarya[-1] $!";
        }
        unless ( -d "$workdir/differential_gene/eukarya/$tmpeukarya[-1]" ) {
            mkdir "$workdir/differential_gene/eukarya/$tmpeukarya[-1]", 0777
                or die
                "failed: can not make dir  $workdir/differential_gene/eukarya/$tmpeukarya[-1] $!";
        }
        &lprint(
            "perl parse_eukarya_gfffile.pl $tmpgff $workdir/differential_gene/eukarya/$tmpeukarya[-1]/ $workdir/eukarya.fa.fai \n"
        );
        `perl $scriptDir/parse_eukarya_gfffile.pl $tmpgff $workdir/differential_gene/eukarya/$tmpeukarya[-1]/ $workdir/eukarya.fa.fai`;
        &lprint(
            "python $scriptDir/hisat2_extract_splice_sites.py $workdir/differential_gene/eukarya/$tmpeukarya[-1]/eukarya.gff > $workdir/differential_gene/eukarya/$tmpeukarya[-1]/splice_sites_gff.txt\n"
        );
        `python $scriptDir/hisat2_extract_splice_sites.py $workdir/differential_gene/eukarya/$tmpeukarya[-1]/eukarya.gtf > $workdir/differential_gene/eukarya/$tmpeukarya[-1]/splice_sites_gff.txt`;
        if (&file_check(
                "$workdir/differential_gene/eukarya/$tmpeukarya[-1]/splice_sites_gff.txt"
            ) > 0
            )
        {
            &lprint(
                "cat $workdir/differential_gene/eukarya/$tmpeukarya[-1]/splice_sites_gff.txt >> $workdir/differential_gene/eukarya/splice_sites_gff.txt\n"
            );
            `cat $workdir/differential_gene/eukarya/$tmpeukarya[-1]/splice_sites_gff.txt >> $workdir/differential_gene/eukarya/splice_sites_gff.txt`;
        }
    }
}

if ( $test eq 'both' || $test eq 'prokaryote' ) {
    unless ( -d "$workdir/sum_gene_count/tmp_count/prokaryote" ) {
        mkdir "$workdir/sum_gene_count/tmp_count/prokaryote", 0777
            or die
            "failed: can not make dir  $workdir/sum_gene_count/tmp_count/prokaryote $!";
    }
    unless ( -d "$workdir/sum_gene_count/read_count/prokaryote" ) {
        mkdir "$workdir/sum_gene_count/read_count/prokaryote", 0777
            or die
            "failed: can not make dir  $workdir/sum_gene_count/read_count/prokaryote $!";
    }
    unless ( -d "$workdir/differential_gene/prokaryote" ) {
        mkdir "$workdir/differential_gene/prokaryote", 0777
            or die
            "failed: can not make dir  $workdir/differential_gene/prokaryote $!";
    }

    my @tmpgff = split /,/, $gff_prokaryote;
    foreach (@tmpgff) {
        my $tmpgff = $_;
        my @tmpprokaryote = split '/', $tmpgff;
        $tmpprokaryote[-1] =~ s/\.//g;
        $tmpprokaryote[-1] =~ s/gff//g;
        if ( &file_check($tmpgff) < 0 ) {
            &lprint(
                "The gff_prokaryote file $tmpgff doesn't exist or empty.\n");
        }
        push @{ $allgff{prokaryote} }, $tmpprokaryote[-1];

        `$jbBin/flatfile-to-json.pl --gff $tmpgff --type CDS  --tracklabel CDS --out $workdir/Jbrowse`;
        `$jbBin/flatfile-to-json.pl --gff $tmpgff --type tRNA  --tracklabel tRNA --out $workdir/Jbrowse`;
        `$jbBin/flatfile-to-json.pl --gff $tmpgff --type exon  --tracklabel exon --out $workdir/Jbrowse`;
        `$jbBin/flatfile-to-json.pl --gff $tmpgff --type gene  --tracklabel gene --out $workdir/Jbrowse`;

        unless (
            -d "$workdir/sum_gene_count/tmp_count/prokaryote/$tmpprokaryote[-1]"
            )
        {
            mkdir
                "$workdir/sum_gene_count/tmp_count/prokaryote/$tmpprokaryote[-1]",
                0777
                or die
                "failed: can not make dir  $workdir/sum_gene_count/tmp_count/prokaryote/$tmpprokaryote[-1] $!";
        }
        unless (
            -d "$workdir/sum_gene_count/read_count/prokaryote/$tmpprokaryote[-1]"
            )
        {
            mkdir
                "$workdir/sum_gene_count/read_count/prokaryote/$tmpprokaryote[-1]",
                0777
                or die
                "failed: can not make dir  $workdir/sum_gene_count/read_count/prokaryote/$tmpprokaryote[-1] $!";
        }
        unless (
            -d "$workdir/differential_gene/prokaryote/$tmpprokaryote[-1]" )
        {
            mkdir "$workdir/differential_gene/prokaryote/$tmpprokaryote[-1]",
                0777
                or die
                "failed: can not make dir  $workdir/differential_gene/prokaryote/$tmpprokaryote[-1] $!";
        }
        &lprint(
            "perl $scriptDir/parse_prokaryote_gfffile.pl $tmpgff  $workdir/differential_gene/prokaryote/$tmpprokaryote[-1]/ $workdir/prokaryote.fa.fai \n"
        );
        `perl $scriptDir/parse_prokaryote_gfffile.pl $tmpgff $workdir/differential_gene/prokaryote/$tmpprokaryote[-1]/ $workdir/prokaryote.fa.fai`;
    }
}

my ( %pair1, %pair2, %usedexp );

if ($pairfile) {
    if ( &file_check($pairfile) < 0 ) {
        die
            "failed: The given pair file $pairfile  doesn't exist or empty.\n";
    }
    else {
        my $linenum = 0;
        open( IN, "$pairfile" )
            or die "failed:  can not open  file $pairfile $!";
        while (<IN>) {
            chomp;
            my @line = split /\t/, $_;
            unless ( $#line == 1 ) {
                &lprint(" failed: check format of the $pairfile file \n");
                exit;
            }
            $linenum++;
            if ( $linenum > 1 ) {
                if (   $#{ $pairs{ $line[0] } } >= 0
                    && $#{ $pairs{ $line[1] } } >= 0 )
                {
                    foreach ( @{ $pairs{ $line[0] } } ) { $usedexp{$_} = 1; }
                    foreach ( @{ $pairs{ $line[1] } } ) { $usedexp{$_} = 1; }
                    push @{ $pair1{ $line[0] } }, $line[1];
                }
                else {
                    &lprint(
                        "failed: $line[0] or $line[1] not defined in experimental descpition\n"
                    );
                    exit;
                }
            }
        }
    }
}
else {
    foreach ( sort keys %pairs ) {
        my $g1 = $_;
        $pair2{$g1} = 1;
        foreach ( @{ $pairs{$g1} } ) { $usedexp{$_} = 1; }
        foreach ( sort keys %pairs ) {
            unless ( $pair2{$_} ) { push @{ $pair1{$g1} }, $_; }
        }
    }
}

#my $command = "bowtie2-build -f $index_fasta $index_bt2";
my $checkIndexFile = join "", ( $index_bt2, '.5.ht2l' );
unless ( -s $checkIndexFile ) {

    #my $index_fasta1=join ',', ($prokaryote_fasta, $eukarya_fasta);
    # $index_fasta1=~s/:/,/g;
    &lprint( "Indexing the reference sequences", 'yellow', "\n" );

    #       &executeCommand($command);
    #     `bowtie2-build -f $index_fasta1 $index_bt2`;
    if ( $eukarya_fasta && $prokaryote_fasta ) {
        &lprint(
            "hisat2-build --large-index $eukarya_fasta,$prokaryote_fasta  $index_bt2\n"
        );
        `hisat2-build --large-index $eukarya_fasta,$prokaryote_fasta  $index_bt2`;
    }
    elsif ($eukarya_fasta) {
        &lprint(
            "hisat2-build --large-index $eukarya_fasta  $index_bt2\n"
        );
        `hisat2-build --large-index $eukarya_fasta  $index_bt2`;
    }
    elsif ($prokaryote_fasta) {
        &lprint(
            "hisat2-build --large-index $prokaryote_fasta $index_bt2\n"
        );
        `hisat2-build --large-index $prokaryote_fasta  $index_bt2`;
    }
    else { &lprint("failed: no INDEX files\n"); exit; }
}

if   ( -s $checkIndexFile ) { &lprint("done INDEX $index_bt2\n"); }
else                        { &lprint("failed: INDEX $index_bt2\n"); exit; }

&printRunTime($time);
&lprint("  Done Checking Files \n");

my $time1 = time();
&lprint("[Trimming and Mapping Reads]\n\tRunning \n\n");

foreach ( sort keys %description ) {
    my $sample   = $_;
    my $rawreads = $description{$sample}{Rawreads_file};
    my $jobname  = join '.', ( $sample, 'RNA_analysis' );

    if ( $rna_mapping_opt eq 'yes' ) {

        if ( $rna_trimming_opt eq 'yes' ) {

            $rawreads =~ s/---ppp---/ /g;
            $rawreads =~ s/:/ /g;

            my $outDir1 = join '/', ( $workdir, "$sample" );
            mkdir $outDir1 if ( !-e $outDir1 );
            if ( !-e $outDir1 ) { print "can not make dir $outDir1\n"; }

            my $troutDir = join '/', ( $outDir1, 'trimming_results' );
            mkdir $troutDir if ( !-e $troutDir );
            if ( !-e $troutDir ) { print "can not make dir $troutDir\n"; }
            #TODO: add a script to check for job id using `qstat -j` will tell you if the job is finished or not. 
            &lprint(
                "qsub -pe smp $numCPU   -l h_vmem=$memlim -v scriptDir=$scriptDir -v test=$test -v numCPU=$numCPU -v workdir=$workdir  -v sample=$sample -v rawreads=$rawreads  -v indexref=$index_bt2 -v descriptfile=$descriptfile  -o $workdir/logdir/$sample -N $jobname $scriptDir/trim_readmapping.sh \n\n"
            );

            `qsub -pe smp $numCPU -l h_vmem=$memlim -v scriptDir=$scriptDir  -v test=$test -v numCPU=$numCPU -v workdir=$workdir -v htseq=$htseq -v sample=$sample -v rawreads="$rawreads"  -v indexref=$index_bt2 -v  descriptfile=$descriptfile  -o $workdir/logdir/$sample -N $jobname $scriptDir/trim_readmapping.sh`;
        }
        else {

            my $outDir1 = join '/', ( $workdir, "$sample" );
            mkdir $outDir1 if ( !-e $outDir1 );
            if ( !-e $outDir1 ) { print "can not make dir $outDir1\n"; }

            my $troutDir = join '/', ( $outDir1, 'trimming_results' );
            mkdir $troutDir if ( !-e $troutDir );
            if ( !-e $troutDir ) { print "can not make dir $troutDir\n"; }

            &lprint("no need to trim rawreads\n");
            my $tmptrname1 = join '.', ( $sample, '1.trimmed.fastq' );
            my $tmptrname2 = join '.', ( $sample, '2.trimmed.fastq' );

            if ( $#{ $allrawreads1{$sample} } == 0 ) {
                `ln -sf $allrawreads1{$sample}[0] $troutDir/$tmptrname1`;
                `ln -sf $allrawreads2{$sample}[0] $troutDir/$tmptrname2`;
            }
            else {
                my $tmpreads1 = join ' ', @{ $allrawreads1{$sample} };
                my $tmpreads2 = join ' ', @{ $allrawreads2{$sample} };
                `cat $tmpreads1 > $troutDir/$tmptrname1`;
                `cat $tmpreads2 > $troutDir/$tmptrname2`;
            }

            &lprint(
                "qsub  -pe smp $numCPU   -l h_vmem=$memlim -v scriptDir=$scriptDir -v test=$test -v numCPU=$numCPU -v workdir=$workdir  -v sample=$sample -v rawreads=$rawreads  -v indexref=$index_bt2 -v descriptfile=$descriptfile  -o $workdir/logdir/$sample -N $jobname $scriptDir/readmapping.sh \n\n"
            );

            `qsub -pe smp $numCPU -l h_vmem=$memlim  -v scriptDir=$scriptDir  -v test=$test -v numCPU=$numCPU -v workdir=$workdir -v htseq=$htseq -v sample=$sample -v rawreads=$rawreads  -v indexref=$index_bt2 -v  descriptfile=$descriptfile  -o $workdir/logdir/$sample -N $jobname $scriptDir/readmapping.sh`;
        }

    }
    else {

        &lprint(
            "qsub -pe smp $numCPU   -l h_vmem=$memlim -v scriptDir=$scriptDir -v test=$test -v numCPU=$numCPU -v workdir=$workdir -v sample=$sample -v rawreads=$rawreads  -v indexref=$index_bt2 -v descriptfile=$descriptfile  -o $workdir/logdir/$sample -N $jobname $scriptDir/parse_BAMfile.sh \n\n"
        );

        `qsub -pe smp $numCPU -l h_vmem=$memlim  -v scriptDir=$scriptDir -v test=$test -v numCPU=$numCPU -v workdir=$workdir  -v sample=$sample -v rawreads=$rawreads  -v indexref=$index_bt2 -v  descriptfile=$descriptfile  -o $workdir/logdir/$sample -N $jobname $scriptDir/parse_BAMfile.sh`;

    }

}

my $alldone = keys(%allsample);
while ($alldone) {
    foreach ( sort keys %description ) {
        my $sample  = $_;
        my $tmpfile = "$workdir/$sample/mapping_results/$sample.stats.text";
        if ( &file_check($tmpfile) > 0 ) {
            $alldone--;
            print LOG"done samples : $alldone\n";
        }
        else { print "$tmpfile not done\n"; print LOG "$tmpfile not done\n"; }
    }
    if ( $alldone > 0 ) {
        print LOG"sample unfinished : $alldone\n";
        sleep(60);
        $alldone = keys(%allsample);
    }
    else {
        last;
    }
}

&printRunTime($time1);
&lprint("  Done Trimming and Mapping Reads \n");

my $time2 = time();
&lprint("[De Novo Detection of small RNAs]\n\tRunning\n\n");

foreach ( sort keys %description ) {
    my $sample = $_;
    if ($prokaryote_fasta) {
        &lprint(
            "qsub  -l h_vmem=$memlim  -o $workdir/$sample/mapping_results/prokaryote_rRNACoverageFold_plot.log -v scriptDir=$scriptDir -v sample=$sample  -v workdir=$workdir  $scriptDir/prokaryote_rRNACoverageFold_plot.sh\n"
        );
        `qsub -l h_vmem=$memlim -o $workdir/$sample/mapping_results/prokaryote_rRNACoverageFold_plot.log  -v scriptDir=$scriptDir  -v sample=$sample  -v workdir=$workdir  $scriptDir/prokaryote_rRNACoverageFold_plot.sh`;
    }

    if ($eukarya_fasta) {

#&lprint ( "qsub -l h_vmem=$memlim  -o $workdir/$sample/mapping_results/eukaryote_rRNACoverageFold_plot.log  -v scriptDir=$scriptDir  -v sample=$sample  -v workdir=$workdir  $Bin/eukaryote_rRNACoverageFold_plot.sh\n");
#`qsub -l h_vmem=$memlim -o $workdir/$sample/mapping_results/eukaryote_rRNACoverageFold_plot.log  -v scriptDir=$scriptDir -v sample=$sample  -v workdir=$workdir  $Bin/eukaryote_rRNACoverageFold_plot.sh`;
    }
}
$alldone = keys(%allsample);
print LOG "total $alldone samples\n";
while ($alldone) {
    foreach ( sort keys %description ) {
        my $sample   = $_;
        my $tmpfile1 = "$workdir/$sample/mapping_results/done.eukarya.txt";
        my $tmpfile2 = "$workdir/$sample/mapping_results/done.prokaryote.txt";

#my $tmpfile4="$workdir/$sample/mapping_results/backward.sorted.eukarya_ref.bedgraph";

        if ( $eukarya_fasta && $prokaryote_fasta ) {
            if ( &file_check($tmpfile1) > 0 && &file_check($tmpfile2) > 0 ) {
                $alldone--;
                print LOG"done samples : $alldone\n";
            }
            else {
                print LOG
                    "eukarya and  prokaryote rRNACoverageFold $sample not done\n";
            }
        }
        elsif ($prokaryote_fasta) {

            if ( &file_check($tmpfile2) > 0 ) {
                $alldone--;
                print LOG"done samples : $alldone\n";
            }
            else {
                print LOG "prokaryote rRNACoverageFold $sample not done\n";
            }

        }
        elsif ($eukarya_fasta) {

            if ( &file_check($tmpfile1) > 0 ) {
                $alldone--;
                print LOG"done samples : $alldone\n";
            }
            else { print LOG "eukarya rRNACoverageFold $sample not done\n"; }
        }
    }
    if ( $alldone > 0 ) {
        print LOG"sample unfinished : $alldone\n";
        sleep(60);
        $alldone = keys(%allsample);
    }
    else {
        last;
    }

}

my $allsample = join ',', @allsample;
my $makegffdone = 1;
if ($prokaryote_fasta) {

    if ( -e "$workdir/newprokaryotegffmade.txt" ) {
        `rm $workdir/newprokaryotegffmade.txt`;
    }

#&lprint ("qsub  -o $workdir/prokaryote_find_small_rna.log -v sample=$allsample  -v reffile=$prokaryote_fasta  -v workdir=$workdir -v gfffile=$gff_prokaryote $Bin/prokaryote_find_small_rna.sh\n");
##`qsub -o $workdir/prokaryote_find_small_rna.log  -v sample=$allsample -v reffile=$prokaryote_fasta  -v workdir=$workdir -v gfffile=$gff_prokaryote  $Bin/prokaryote_find_small_rna.sh`;
    &lprint(
        "perl $scriptDir/prokaryote_find_small_rna.pl $workdir $prokaryote_fasta  $gff_prokaryote $allsample\n"
    );
    `perl $scriptDir/prokaryote_find_small_rna.pl $workdir $prokaryote_fasta  $gff_prokaryote $allsample`;

    while ($makegffdone) {
        if ( -e "$workdir/newprokaryotegffmade.txt" ) {
            $makegffdone = 0;
            last;
        }
        else {
            sleep(60);
            print LOG "doing prokaryote_find_small_rna\n";
        }
    }
}
$makegffdone = 1;
if ( &file_check($eukarya_fasta) > 0 && &file_check($gff_eukarya) > 0 ) {

    if ( -e "$workdir/neweukaryotegffmade.txt" ) {
        `rm $workdir/neweukaryotegffmade.txt`;
    }

#print LOGFILE"qsub -o $workdir/prokaryote_find_small_rna.log  -v sample=$allsample   -v reffile=$eukarya_fasta -v workdir=$workdir -v gfffile=$gff_eukarya $Bin/eukaryote_find_small_rna.sh\n";
##`qsub  -o $workdir/prokaryote_find_small_rna.log -v sample=$allsample -v reffile=$eukarya_fasta  -v workdir=$workdir -v gfffile=$gff_eukarya $Bin/eukaryote_find_small_rna.sh`;
#&lprint ("perl $Bin/eukaryote_find_small_rna.pl $workdir $eukarya_fasta  $gff_eukarya $allsample\n");
#`perl $Bin/eukaryote_find_small_rna.pl $workdir $eukarya_fasta  $eukarya_fasta $allsample`;

    `pwd > $workdir/neweukaryotegffmade.txt`;

    while ($makegffdone) {
        if ( -e "$workdir/neweukaryotegffmade.txt" ) {
            $makegffdone = 0;
            last;
        }
        else {
            sleep(60);
            print LOG "doing eukaryote_find_small_rna\n";
        }
    }
}

&printRunTime($time2);
&lprint("Done De Novo Detection of small RNAs \n");

my $time3 = time();
&lprint("[Differential Gene Analysis]\n\tRunning\n\n");

foreach ( sort keys %description ) {
    my $sample = $_;

    &lprint(
        "qsub -l h_vmem=$memlim -o $workdir/$sample/mapping_results/htseq.log -v scriptDir=$scriptDir   -v test=$test  -v sample=$sample  -v workdir=$workdir $scriptDir/htseq-count.sh\n"
    );
    `qsub -l h_vmem=$memlim -o $workdir/$sample/mapping_results/htseq.log -v scriptDir=$scriptDir   -v test=$test  -v sample=$sample  -v workdir=$workdir $scriptDir/htseq-count.sh`;
}

$alldone = keys(%allsample);
print LOG "total $alldone samples\n";
while ($alldone) {
    foreach ( sort keys %description ) {
        my $sample = $_;
        my $tmpfile1
            = "$workdir/$sample/mapping_results/finish_htseq-count.txt";
        if ( &file_check($tmpfile1) > 0 ) {
            $alldone--;
            print LOG"done htseq-count samples : $alldone\n";
        }
        else { print LOG "htseq-count $sample not done\n"; }
    }
    if ( $alldone > 0 ) {
        print LOG"htseq-count sample unfinished : $alldone\n";
        sleep(60);
        $alldone = keys(%allsample);
    }
    else {
        last;
    }
}

foreach ( sort keys %description ) {
    my $sample = $_;

    my $jobname = join '.', ( $sample, 'RNA_analysis' );

    my $tmpname
        = $sample . $description{$sample}{group} . 'prokaryote_forward';
    `$jbBin/add-bw-track.pl --in $workdir/Jbrowse/trackList.json --out $workdir/Jbrowse/trackList.json --plot --label $tmpname  --bw_url Jbrowse/BigWig/$sample.forward.sorted.prokaryote_ref.bw --key $tmpname `;

    $tmpname = $sample . $description{$sample}{group} . 'prokaryote_backward';
    `$jbBin/add-bw-track.pl --in $workdir/Jbrowse/trackList.json --out $workdir/Jbrowse/trackList.json --plot --label $tmpname  --bw_url Jbrowse/BigWig/$sample.backward.sorted.prokaryote_ref.bw --key $tmpname `;

    $tmpname = $sample . $description{$sample}{group} . 'eukarya_forward';
    `$jbBin/add-bw-track.pl --in $workdir/Jbrowse/trackList.json --out $workdir/Jbrowse/trackList.json --plot --label $tmpname  --bw_url Jbrowse/BigWig/$sample.forward.sorted.eukarya_ref.bw --key $tmpname `;

    $tmpname = $sample . $description{$sample}{group} . 'eukarya_backward';
    `$jbBin/add-bw-track.pl --in $workdir/Jbrowse/trackList.json --out $workdir/Jbrowse/trackList.json --plot --label $tmpname  --bw_url Jbrowse/BigWig/$sample.backward.sorted.eukarya_ref.bw --key $tmpname `;

}

my %diffdir;

opendir( DIR, "$workdir/sum_gene_count/tmp_count/" ) or die $!;
while ( my $tmpdir = readdir(DIR) ) {
    next unless -d "$workdir/sum_gene_count/tmp_count/$tmpdir";
    next if ( $tmpdir =~ /^\./ );
    opendir( DIR1, "$workdir/sum_gene_count/tmp_count/$tmpdir" ) or die $!;
    while ( my $tmpdir1 = readdir(DIR1) ) {
        next unless -d "$workdir/sum_gene_count/tmp_count/$tmpdir/$tmpdir1";
        next if ( $tmpdir1 =~ /^\./ );
        print LOG"$tmpdir\t$tmpdir1\n";
        $diffdir{$tmpdir}{$tmpdir1} = 1;
    }
}

sleep(5);

my $allgoodsample = 0;

foreach ( sort keys %diffdir ) {
    my $kingdom = $_;

    unless ( $kingdom eq $test || $test eq 'both' ) { next; }
    foreach ( sort keys %{ $diffdir{$kingdom} } ) {
        my $strain = $_;
        print LOG"$kingdom\t$strain\n";
        sleep(5);
        my $rRNAdes
            = "$workdir/differential_gene/$kingdom/$strain/$kingdom.genedesc.rRNA.txt";
        my %rRNApos;

        if ( &file_check($rRNAdes) > 0 ) {

            open( GENOIN, "$rRNAdes" ) or die "failed: $rRNAdes $!";
            while (<GENOIN>) {

                chomp;

                my @line = split /\s+/, $_;
                for ( my $i = 0; $i <= $#line; $i++ ) {
                    $rRNApos{ $line[0] } = $line[-1];
                }
            }
            close GENOIN;
        }

        foreach ( sort keys %description ) {
            my $sample = $_;
            open( READIN,
                ">$workdir/sum_gene_count/read_count/$kingdom/$strain/$sample.$kingdom.name.htseq.locus_tag.txt"
                )
                or die
                "failed: $! $workdir/sum_gene_count/read_count/$kingdom/$strain/$sample.$kingdom.name.htseq.locus_tag.txt\n";
            open( IN,
                "$workdir/sum_gene_count/tmp_count/$kingdom/$strain/$sample.$kingdom.name.htseq.locus_tag.txt"
                )
                or die
                "failed:  can not open description file $workdir/sum_gene_count/tmp_count/$kingdom/$strain/$sample.$kingdom.name.htseq.locus_tag.txt $!";
            while (<IN>) {
                chomp;
                my $line = $_;
                if ( $line =~ m/^no_feature/ ) { last; }
                my @line = split /\s+/, $line;
                if ( $rRNApos{ $line[0] } )

                {
                }
                else {
                    print READIN "$line\n";
                }
            }

            close READIN;
            close IN;
        }

        my $Deseqdir = "$workdir/differential_gene/$kingdom/$strain/";

        unless ( -e "$workdir/differential_gene/$kingdom/$strain/" ) { next; }

        unless ( -d "$Deseqdir/Deseq/" ) {
            mkdir "$Deseqdir/Deseq/", 0777
                or die
                "failed: can not make dir  $workdir/differential_gene/$kingdom/$strain/Deseq/ $!";
        }
        unless ( -d "$Deseqdir/EdgeR/" ) {
            mkdir "$Deseqdir/EdgeR/", 0777
                or die
                "failed: can not make dir  $workdir/differential_gene/$kingdom/$strain//EdgeR/ $!";
        }
        unless ( -d "$Deseqdir/figures/" ) {
            mkdir "$Deseqdir/figures/", 0777
                or die
                "failed: can not make dir  $workdir/differential_gene/$kingdom/$strain//figures/ $!";
        }
        unless ( -d "$Deseqdir/significant_gene/" ) {
            mkdir "$Deseqdir/significant_gene/", 0777
                or die
                "failed: can not make dir  $workdir/differential_gene/$kingdom/$strain//significant_gene/ $!";
        }

        my %tablecounts;
        open( DESC, ">$Deseqdir/reads.table.txt" )
            or die "failed: $! $Deseqdir/reads.table.txt\n";
        print DESC "Sample";
        my %reads_sample;
        my %reads_desc;
        my $totalreadsmapped = 0;

        foreach ( sort keys %description ) {
            my $sample = $_;
            open( IN,
                "$workdir/sum_gene_count/read_count/$kingdom/$strain/$sample.$kingdom.name.htseq.locus_tag.txt"
                )
                or die
                "failed:  can not open description file $workdir/sum_gene_count/read_count/$kingdom/$strain/$sample.$kingdom.name.htseq.locus_tag.txt $!";
            while (<IN>) {
                chomp;
                my $line = $_;
                my @line = split /\s+/, $line;
                $tablecounts{ $line[0] }{$sample} = $line[1];
                $reads_sample{$sample}  += $line[1];
                $reads_desc{ $line[0] } += $line[1];
                $totalreadsmapped       += $line[1];
            }
            close IN;
        }

        my @goodsample;
        my $numsample = keys(%reads_sample);
        foreach ( sort keys %reads_sample ) {
            my $sample = $_;
            if ($reads_sample{$sample}
                >= 0.001 * $totalreadsmapped / $numsample
                && (   $totalreadsmapped / $numsample > 2
                    || $totalreadsmapped > 2000 )
                && $usedexp{$sample}
                )
            {
                push @goodsample, $sample;
                print DESC "\t$sample";
            }
        }
        print DESC "\n";
        if ( $#goodsample < 1 ) {
            &lprint("no good sample for $kingdom, $strain\n");
            next;
        }
        $allgoodsample++;
        foreach ( sort keys %tablecounts ) {
            my $tmpgene = $_;
            if ( $reads_desc{$tmpgene} <= 1 ) { next; }

            print DESC "$tmpgene";
            foreach (@goodsample) {
                print DESC "\t$tablecounts{$tmpgene}{$_}";
            }
            print DESC "\n";

        }
        close DESC;

        open( DESC, ">$Deseqdir/readcounts.expriment.txt" )
            or die "failed: $! $Deseqdir/readcounts.expriment.txt\n";

        print DESC "ID\tfiles\tgroup\n";

        my %goodgroup;

        foreach (@goodsample) {
            my $sample   = $_;
            my $tmpgroup = $sample;
            my @tmpgroup = split /_/, $tmpgroup;
            my $htseqfile
                = "$workdir/sum_gene_count/read_count/$kingdom/$strain/$sample.$kingdom"
                . '.name'
                . '.htseq'
                . '.locus_tag' . '.txt';
            print DESC "$sample\t$htseqfile\t$description{$sample}{group}";
            print DESC "\n";
            $goodgroup{ $description{$sample}{group} }++;
        }
        close DESC;

        &lprint("start differential gene finding using $test_method\n");

        chdir("$Deseqdir") or die "failed: $!";

        open( DESC, ">$Deseqdir/Deseq_EdgeRpairs.txt" )
            or die "failed: $! $Deseqdir/Deseq_EdgeRpairs.txt\n";

        print DESC "group1\tgroup2\n";

        foreach ( sort keys %pair1 ) {
            my $group1 = $_;
            next unless $goodgroup{$group1};
            foreach ( @{ $pair1{$group1} } ) {
                next unless $goodgroup{$_};
                print DESC "$group1\t$_\n";
            }
        }
        close DESC;

        if ( $test_method eq 'both' ) {

            &lprint("EdgeR and Deseq2 \t $Deseqdir\n");

            &lprint("Rscript $Bin/EdgeR.R  $p_cutoff\n");
            `Rscript $Bin/EdgeR.R  $p_cutoff`;

            &lprint("Rscript $Bin/Deseq.R  $p_cutoff\n");
            `Rscript $Bin/Deseq.R  $p_cutoff`;

            &lprint("Rscript $Bin/vennDiagram_Deseq_EdgeR.R\n");
            `Rscript $Bin/vennDiagram_Deseq_EdgeR.R`;

        }
        elsif ( $test_method eq 'EdgeR' ) {
            &lprint("EdgeR \t $Deseqdir \n");

            &lprint("Rscript $Bin/EdgeR.R  $p_cutoff\n");
            `Rscript $Bin/EdgeR.R  $p_cutoff`;
        }
        elsif ( $test_method eq 'Deseq' ) {

            &lprint("Deseq2 \t $Deseqdir \n");

            &lprint("Rscript $Bin/Deseq.R  $p_cutoff\n");
            `Rscript $Bin/Deseq.R  $p_cutoff`;
        }
        else {
            &lprint("failed: method $test_method is invalid\n");
            exit;
        }

        if ( $kingdom eq 'prokaryote' ) {
            &lprint(
                "perl $Bin/Differential_stats_prokaryote.pl $workdir $Deseqdir $descriptfile $workdir/differential_gene/$kingdom/$strain/prokaryote.NonrRNA.genedesc.txt  $p_cutoff @goodsample\n"
            );

            `perl $Bin/Differential_stats_prokaryote.pl $workdir $Deseqdir $descriptfile $workdir/differential_gene//$kingdom/$strain/prokaryote.NonrRNA.genedesc.txt $p_cutoff @goodsample`;

            `$jbBin/flatfile-to-json.pl --gff $workdir/sum_direction_prokaryote_ref.gff --type expressed_intergenic_region  --tracklabel expressed_intergenic_region  --out $workdir/Jbrowse`;

        }

        if ( $kingdom eq 'eukarya' ) {
            &lprint(
                "perl $Bin/Differential_stats_eukarya.pl $workdir $Deseqdir $descriptfile $workdir/differential_gene/$kingdom/$strain/eukarya.genedesc.txt  $p_cutoff @goodsample\n"
            );

            `perl $Bin/Differential_stats_eukarya.pl $workdir $Deseqdir $descriptfile $workdir/differential_gene/$kingdom/$strain/eukarya.genedesc.txt $p_cutoff @goodsample`;

            `$jbBin/flatfile-to-json.pl --gff $workdir/sum_direction_eukarya_ref.gff --type expressed_intergenic_region  --tracklabel expressed_intergenic_region  --out $workdir/Jbrowse`;

        }
        &lprint("done $kingdom\t$strain\n");

        &lprint("Rscript $Bin/scatt_plot.R  $p_cutoff\n");
        `Rscript $Bin/scatt_plot.R  $p_cutoff`;

        my @workdir = split '/', $workdir;
        my $jbrowsed;
        if ( $workdir[-1] ) {
            $jbrowsed = $workdir[-1];
        }
        elsif ( $workdir[-2] ) {
            $jbrowsed = $workdir[-2];
        }
        else {
            $jbrowsed = $workdir[-3];
        }

        `$jbBin/generate-names.pl --out  $workdir/Jbrowse`;

#if (-e "$Bin/../../edge_ui/JBrowse/data/$jbrowsed") {`unlink $Bin/../../edge_ui/JBrowse/data/$jbrowsed`;}
        #NOTE: I dont understand this? how is this url created
		`ln -s  $workdir/Jbrowse/ $Bin/bin/JBrowse/data/$jbrowsed`;
        &lprint(
            "Jbrowse link is at http://ergatis2.lanl.gov/jbrowse/?data=data/$jbrowsed\n"
        );

        sleep(5);
    }
}    #foreach (sort keys %diffdir)
if ( $allgoodsample > 0 ) {
    &printRunTime($time3);
    &lprint("  Done Differential Gene Analysis \n");
}
else {

    &lprint(
        "  Differential Gene Analysis failed because of no good samples \n");
    &printRunTime($time3);
    &lprint(
        "  failed to do Differential Gene Analysis failed because of no good samples \n"
    );
}

&lprint("Total");

&printRunTime($^T);

&lprint("All Done \n\n");

open( CURRENTLOGFILE, ">> $workdir/process_current.log" )
    or die "   $workdir/process_current.log $!";
print CURRENTLOGFILE "All Done\n";

close LOG;
close CURRENTLOGFILE;

sub Usage 
{
    print <<"END";

DESCRIPTION\n
	Pipeline that finds differentially expressed genes from raw fastq data.

 
Usage: perl $0 [options] -exp experimental_description_file.txt -d pipeline_test_both -prokaryote_fasta test_prok.fa -eukarya_fasta eukarya_test.fa -index_ref_bt2 indexfile -gff_prokaryote test_prokaryote.gff -gene_coverage_ref gene_coverage_reference.fa

  example: 
                       

OPTIONS
   Generic Program Information
		--help Print a usage message briefly summarizing these command-line options, then exit.

		-V, --version
			Print the version number to the standard output stream.  This version number should be included in all bug reports.

 
       	-W, --workdir
			Path to a directory where the whole analysis will be stored. Must have write permission.

		-F, --gff_prokaryote
			(Optional) prokaryote annotation file(s) in gff format. Multiple files can be seperated by commas. Needed for diffrential gene analysis.

        -G, --gff_eukarya
			(Optional) eukarya annotation file(s) in gff format. Multiple files seperated by commas. Needed for diffrential gene analysis.

		-U, --eukarya_fasta
			(Optional) comma separated list of reference sequences.

        -F, --prokaryote_fasta 
			(Optional) comma separated list of reference sequences. 

        -I, --index_ref_bt2
			hisat2 mapping index file. If the mapping index is already generated, one can pass that as well.

        -M, --h_vmem
			memory limit per node, string, default 20G

        -C, --gene_coverage_fasta
			(Optional) single fasta file  (for directional coverage analysis, sequnce must be prokaryote mapping reference sequence) 

        -K, --test_kingdom
			comparison of differential gene expression. `both` (for eukarya and prokaryote), `prokaryote`, or `eukarya`. 
			Default: prokaryote

        -T, --test_method
			R package/ method to use for differential gene expression analysis. `both` (for both EdgeR and Deseq2 method), EdgeR, or Deseq.
 			#TODO: change default to EdgeR
			Default: both (must have have at least 3 duplicates if using Deseq2)

        -P, --cpu
			number of cpu to be used (default 1)
		#TODO: check whats going on with this option
        -BAM_ready      #if mapping file are provided for samples by users (yes or no) default no
		
        -S, --significant_pvalue
			floating number cutoff to define significant differentially express genes.
			Default: 0.001

        -E, --exp
			tab delimited txt file descripting experimental design. Each row represents one sample.
                      Each coloumn is as:
                    (
                     ID:  uniq sample ID
                     Rawreads_file: absolute path, fastq format, pair reads seperated by colon, multiple datasets seprarate semicolon
                     group:    replicates group name for this project, each group must have uniqe name (without -pair_comparison option defined as below, all groups will be compared to each other in differential gene analysis)
                     experimental condiitons such as clock time, CFU etc:  one condition per name per coloumn, can be multiple colomuns, (for differentail gene analysis)
                     
                      Example: (the names and order of the first 3 colomums must be the exact the same as below) :  
                      ID		Rawreads_files/BAM_file                 group 
                      exp1      read1p1:read1p2;read1p1a:read1p2a       time0
                      exp2      read2p1:read2p2                         time0
                      .
                      .
                      .
                      exp21     read21p1:read21p2;read21p1a:read21p2a   timef
                      exp22     read22p1:read22p2;read22p1a:read22p2a   timef                              
                    )
        differential gene analysis  design text file: 
         -N, --pair_comparison
			tab delimited file describing pairwise comparison. If this file does NOT exist, all groups defined in the experimental design file will be compared to each other for  differential gene analysis. 
                               Each colomum is as:
                              (
                                group1:    first group name in pairwised comparison for this project (must be defined in the master design text file)
                                group2:    second group name in pairwised comparison for this project (must be defined in the master design text file)
                               Example: 
                                 group1         group2 
                                 time0          time1
                                 time0          time2
                                 time0          timef
                                 time2          timef 
                               )

#TODO: change this example to one in test
EXAMPLE:
	perl runPiReT_qsub.pl -geneopt gene -test_kingdom prokaryote  -significant_pvalue 0.001  -cpu 10 -exp /users/203270/scratch/momo_Rnaseq/Analysis_BTT_2015AUG/BTT_Experimetal_descriptions.txt -d ~/scratch/momo_Rnaseq/Analysis_BTT_2015AUG/  -prokaryote_fasta /users/203270//scratch/momo_Rnaseq/db/bowtie2/Bacillus_anthracis__Ames_Ancestor_uid58083.fa -eukarya_fasta /users/203270/scratch/momo_Rnaseq/db/bowtie2/cavPor3.fa -index_ref_bt2 /users/203270/scratch/momo_Rnaseq/db/bowtie2/Bacillus_anthracis__Ames_Ancestor_uid58083_CAVPor3i_hisat -gff_prokaryote /users/203270/scratch/momo_Rnaseq/db/Bacillus_anthracis__Ames_Ancestor_uid58083.gff -test_method EdgeR  -gene_coverage_fasta /users/203270/scratch/momo_Rnaseq/db/bowtie2/Bacillus_anthracis__Ames_Ancestor_uid58083.fa -pair_comparison ~/scratch/momo_Rnaseq/Analysis_BTT_2015AUG/pair_comparision.txt 

END
    exit;
}

sub file_check {

    #check file exist and non zero size
    my $file  = shift;
    my $exist = -1;
    if ( -e $file ) { $exist = 1 }
    if ( -z $file ) { $exist = -1 }
    return $exist;
}

sub getTmpNameByTime {
    my $now_string = strftime "%Y %b %e %H:%M:%S", localtime;
    return $now_string;
}

sub printRunTime {
    my $time        = shift;
    my $runTime     = time() - $time;
    my $time_string = sprintf(
        " Running time: %02d:%02d:%02d\n\n",
        int( $runTime / 3600 ),
        int( ( $runTime % 3600 ) / 60 ),
        int( $runTime % 60 )
    );
    &lprint($time_string);
}

sub lprint {
    my ($line) = @_;
    print LOG $line;
}

sub readfai {
    my $faifile = shift;
    my @tmpcontigs;
    open( IN, $faifile ) or die "failed: can not open $faifile $!";
    while (<IN>) {
        chomp;
        my $line = $_;
        my @line = split /\t+/, $line;
        push @tmpcontigs, $line[0];
    }
    close IN;
    return @tmpcontigs;
}

sub printVersion 
{
    print basename($0), " version: $version\n";
    exit;
}

sub checkDependedPrograms

    #TODO: Also check for appropriate version

{
    system("which bwa 1>/dev/null") == 0
        || die "\nbwa is not in your PATH\n $ENV{PATH}\n";
    system("which samtools 1>/dev/null") == 0
        || die "\nsamtools is not in your PATH\n $ENV{PATH}\n";
    system("which R 1>/dev/null") == 0
        || die "\nR is not in your PATH\n $ENV{PATH}\n";
    system("which hisat2-build 1>/dev/null") == 0
        || die "\nhisat2-build is not in your PATH\n $ENV{PATH}\n";
	system("which qsub 1>/dev/null") == 0
		|| die "qsub is not installed\n";   
# system("which parse_eukarya_gfffile.pl 1>/dev/null") == 0
    #     || die
    #     "\nparse_eukarya_gfffile.pl is not in your PATH\n $ENV{PATH}\n";
}
