#!/usr/bin/perl -W
use strict;

use Getopt::Long;
use File::Basename;
use Term::ANSIColor;

my ($descriptfile,
    $test,
    $splice_file_out,
    $splice_file_in,
    $pairfile,
    $diffdir,
    $workdir,
    $numCPU,
    $eukarya_fasta,
    $prokaryote_fasta,
    $index_bt2,
    $gff_eukarya,
    $gff_prokaryote,
    $coverage_fasta);

$gff_eukarya = 'unknown';
$gff_prokaryote = 'unknown';
my $mapread = 'no';
my  $rna_trimming_opt='yes';
$numCPU=1;
$test='prokaryote';
my $test_method='both';
my $mptool='bowtie2';
my $htseq='gene';
my $p_cutoff =0.001;
my $memlim='10G';
$|=1; #?

GetOptions(
            'exp=s' => \$descriptfile,
            'd=s' => \$workdir,
            'cpu=i' => \$numCPU, # bwa option
            'eukarya_fasta=s' => \$eukarya_fasta,
            'prokaryote_fasta=s' => \$prokaryote_fasta,
            'index_ref_bt2=s' => \$index_bt2,
            'h_vmem=s' => \$memlim, 
            'gff_eukarya=s' => \$gff_eukarya,
            'gff_prokaryote=s' => \$gff_prokaryote,
            'gene_coverage_fasta=s' => \$coverage_fasta,
            'significant_pvalue=f' => \$p_cutoff,
            'pair_comparison=s' => \$pairfile,
            'BAM_ready=s' =>\$mapread,      #if mapping file are provided for samples by users
            'trimming_reads' => $rna_trimming_opt, 
            'test_kingdom=s' => \$test,
            'test_method=s' => \$test_method,
            'help|?' => sub{&Usage()}
);

#?
my %allrawreads1;
my %allrawreads2;


unless ($descriptfile && $workdir && $test && ($gff_eukarya || $gff_prokaryote)  && ($eukarya_fasta || $prokaryote_fasta) && $index_bt2 )  {&Usage};

############checking permissions and creating required directories##############

# create working directory
unless (-d "$workdir") {mkdir "$workdir", 0777 or die "can not make dir  $workdir $!";}

# create Jbrowse directory
if (-d "$workdir/Jbrowse/") {
  `rm -r $workdir/Jbrowse/`; mkdir "$workdir/Jbrowse/", 0777 or die "can not make dir  $workdir/Jbrowse $!";
} 
else {
  mkdir "$workdir/Jbrowse/", 0777 or die "can not make dir  $workdir/Jbrowse $!";
}
# create BigWig directory
mkdir "$workdir/Jbrowse/BigWig", 0777 or die "can not make dir  $workdir/Jbrowse/BigWig $!";

#create sum_gene_count_directory and logdir
unless (-d "$workdir/sum_gene_count/") {mkdir "$workdir/sum_gene_count/", 0777 or die "can not make dir  $workdir/sum_gene_count/ $!";}
unless (-d "$workdir/logdir/") {mkdir "$workdir/logdir/", 0777 or die "can not make dir  $workdir/dir/ $!";}

################################################################################

########################opening logfile.txt#####################################
open (LOGFILE, ">$workdir/logfile.txt") or die "  can not open log file $workdir/logfile.txt $!";


#? what are these variables?
my %description;
my @colname;
my %expdescription;
my %allsample;
my @allsample;
my %pairs;


##########################parse description file################################
open (IN, "$descriptfile") or die " cannot open description file $descriptfile $!";
 my $linenum=0;
 while (<IN>)
{
         chomp; #?
        my @line=split /\s+/, $_; #?
       unless  ($#line >=2) { print LOGFILE"check format of the description file $_ \n"; exit;} #?
        $linenum++; #?
      if ($linenum ==1) #?
   {
         for (my $i=0; $i<=$#line; $i++)
    {
      $line[$i]=~ s/\s+//g;
      $colname[$i]=$line[$i];
    }
   }
        if ($linenum >1)
  {
    $allsample{$line[0]}=1;
     push @allsample, $line[0];
      for (my $i=3; $i<=$#line; $i++)
    {
     $line[$i]=~ s/\s+//g;
     $expdescription{$line[0]}{$colname[$i]}=$line[$i];
    }

      my @line1=split /;/, $line[1];
    foreach (@line1)
  {
    my @rawreadsfile=split /:/, $_;
    
     push @{$allrawreads1{$line[0]}}, $rawreadsfile[0];
    push @{$allrawreads2{$line[0]}}, $rawreadsfile[1];

      foreach (@rawreadsfile)
    {
   if(&file_check($_)<0 ) { die "The reads file $_ doesn't exist or empty.\n";}
    }
  }
  $description{$line[0]}{group}=$line[2];
 push @{$pairs{$line[2]}}, $line[0];
   if ($line[2] =~ /\s+/) {
   print LOGFILE"please check format of experimental descpition \n group name can not have space as : $line[2]\n";
   exit;
  }
   my $rawreadsline=$line1[0];
  # my @rawreadsline=split /:/, $rawreadsline;
   for (my $i=1; $i<=$#line1; $i++) {$rawreadsline=$rawreadsline.'---ppp---'.$line1[$i]         }
#    print LOGFILE"$rawreadsline:\n\n"; 
   $description{$line[0]}{Rawreads_file}=$rawreadsline;
  }
 }
close IN;
################################################################################

################################check if enough sample for DeSeq################
foreach (keys %pairs)
{
my $tmpgroup=$_;
if ( (($test_method eq 'both') || ($test_method eq 'Deseq') ) && $#{$pairs{$tmpgroup}}<2)
 {
   print LOGFILE"testing method using Deseq2 but the $tmpgroup has only $pairs{$tmpgroup} which is less than 3. please consider not using Deseq method or increase the number of replicates for $tmpgroup\n";
 exit;
 }
}
################################################################################



my (@eukaryagff, @prokaryotegffi, %allgff);
#####################create additional directories##############################
unless (-d "$workdir/differential_gene/") {mkdir "$workdir/differential_gene/", 0777 or die "can not make dir  $workdir/differential_gene/ $!";}
unless (-d "$workdir/sum_gene_count/tmp_count/") {mkdir "$workdir/sum_gene_count/tmp_count/", 0777 or die "can not make dir  $workdir/sum_gene_count/tmp_count/ $!";}
unless (-d "$workdir/sum_gene_count/read_count/") {mkdir "$workdir/sum_gene_count/read_count/", 0777 or die "can not make dir  $workdir/sum_gene_count/read_count/ $!";}

#unless (-d "$workdir/sum_gene_count/tmp_count/prokaryote/") {mkdir "$workdir/sum_gene_count/tmp_count/prokaryote/", 0777 or die "can not make dir  $workdir/sum_gene_count/tmp_count/prokaryote $!";}
#unless (-d "$workdir/sum_gene_count/read_count/eukarya/") {mkdir "$workdir/sum_gene_count/read_count/eukarya/", 0777 or die "can not make dir  $workdir/sum_gene_count/read_count/eukarya $!";}

#unless (-d "$workdir/sum_gene_count/read_count/prokaryote/") {mkdir "$workdir/sum_gene_count/read_count/prokaryote/", 0777 or die "can not make dir  $workdir/sum_gene_count/read_count/prokaryote $!";}
#unless (-d "$workdir/sum_gene_count/tmp_count/eukarya/") {mkdir "$workdir/sum_gene_count/tmp_count/eukarya/", 0777 or die "can not make dir  $workdir/sum_gene_count/tmp_count/eukarya $!";}
################################################################################

######################based on condition check if files are there###############
#and if present run the index
my %allcontigs;

if ($test eq 'both' || $test eq 'eukarya' )
{
if ($eukarya_fasta) 
 {
  if(&file_check($eukarya_fasta)<0  ) 
   {
     die "The eukarya_fast file $eukarya_fasta doesn't exist or empty.\n";
   } else {
   if (-e "$workdir/eukarya.fa") { `rm $workdir/eukarya.fa`;}
    #TODO: create a shortcut of eukarya fasta, but it might be better to copy it like i did for single processor version
      `ln -fs $eukarya_fasta $workdir/eukarya.fa`;
      `samtools faidx $workdir/eukarya.fa`;
    #TODO: install jbrowse and make sure this prepare-refseq.pl is in the path
      `/edge_dev/edge_ui/JBrowse/bin/prepare-refseqs.pl --trackLabel  DNA --seqType dna --key 'DNA+protein' --fasta  $workdir/eukarya.fa --out  $workdir/Jbrowse/`;
     my @contigs=&readfai("$workdir/eukarya.fa.fai");
    foreach (@contigs) {$allcontigs{$_}++;} #?
  } 
 }else {die "need eukarya sequence file\n";}
}
################################################################################
if ($test eq 'both' || $test eq 'prokaryote' )
{
if ($prokaryote_fasta) 
 {
   
   if( &file_check($prokaryote_fasta)<0 ) 
   { 
    die "The prokaryote_fast file $prokaryote_fasta doesn't exist or empty.\n";
   } else {
   if (-e "$workdir/prokaryote.fa") {`rm $workdir/prokaryote.fa`;}
    #TODO: see todos of eukaryotes
    `ln -fs $prokaryote_fasta $workdir/prokaryote.fa`;
    `samtools faidx $workdir/prokaryote.fa`;   
    `/edge_dev/edge_ui/JBrowse/bin/prepare-refseqs.pl --trackLabel  DNA --seqType dna --key 'DNA+protein' --fasta  $workdir/prokaryote.fa --out  $workdir/Jbrowse/`;
       my @contigs=&readfai("$workdir/prokaryote.fa.fai");
    foreach (@contigs) {$allcontigs{$_}++;}
   } 
 }else  {die "need prokaryote sequence file\n";}
} 
################################################################################
#?
foreach (keys %allcontigs) {
  if ($allcontigs{$_}>1) {
    print LOGFILE"contig name $_ is duplicated, please make sure the name of contigs is unique in reference fasta files (prokaryout and eukarya references)\n"; exit;} }

################################################################################

###########################coverage fasta#######################################
my $pjcoverage="$workdir/coverage.fa";
if (-e "$workdir/coverage.fa.fai") {`rm $workdir/coverage.fa.fai`;}
if ($coverage_fasta)
 {
   if(&file_check($coverage_fasta)<0 )
   { die  "The coverage_fasta file $coverage_fasta doesn't exist or empty.\n";
   } else {
   if (-e $pjcoverage) {`rm $pjcoverage`;}
    `ln -fs $coverage_fasta $pjcoverage`;
    `samtools faidx $pjcoverage`;
       my @contigs=&readfai("$workdir/prokaryote.fa.fai");
    foreach (@contigs) { unless ($allcontigs{$_}) {print LOGFILE"coverage contig $_ is not part of mapping index references, please make sure it is included in either prokaryout or eukarya reference\n"; exit;}}
   }
 } else {$pjcoverage='NA'}
###########################coverage fasta#######################################

if ($test eq 'both' || $test eq 'eukarya' )
{
      unless (-d "$workdir/sum_gene_count/tmp_count/eukarya") {mkdir "$workdir/sum_gene_count/tmp_count/eukarya", 0777 or die "can not make dir  $workdir/sum_gene_count/tmp_count/eukarya $!";}
     unless (-d "$workdir/sum_gene_count/read_count/eukarya") {mkdir "$workdir/sum_gene_count/read_count/eukarya", 0777 or die "can not make dir  $workdir/sum_gene_count/read_count/eukarya $!";}
     unless (-d "$workdir/differential_gene/eukarya") {mkdir "$workdir/differential_gene/eukarya", 0777 or die "can not make dir  $workdir/differential_gene/eukarya $!";}
   my @tmpgff=split /,/,$gff_eukarya;

  # if (-e "$workdir/differential_gene/eukarya/splice_sites_gff.txt") { `rm $workdir/differential_gene/eukarya/splice_sites_gff.txt`;} 

 foreach (@tmpgff) { 
     my $tmpgff=$_;
      my @tmpeukarya=split '/', $tmpgff; 
     $tmpeukarya[-1] =~ s/\.//g;
     $tmpeukarya[-1] =~ s/gff//g;  
     if(&file_check($tmpgff)<0 )
     { 
       next "The gff_eukarya file $tmpgff doesn't exist or empty.\n";
     }
     push @{$allgff{eukarya}}, $tmpeukarya[-1];

      `/edge_dev/edge_ui/JBrowse/bin/flatfile-to-json.pl --gff $tmpgff --type CDS  --tracklabel CDS --out $workdir/Jbrowse`;
       `/edge_dev/edge_ui/JBrowse/bin/flatfile-to-json.pl --gff $tmpgff --type tRNA  --tracklabel tRNA --out $workdir/Jbrowse`;
        `/edge_dev/edge_ui/JBrowse/bin/flatfile-to-json.pl --gff $tmpgff --type exon  --tracklabel exon --out $workdir/Jbrowse`;
         `/edge_dev/edge_ui/JBrowse/bin/flatfile-to-json.pl --gff $tmpgff --type gene  --tracklabel gene --out $workdir/Jbrowse`;

     unless (-d "$workdir/sum_gene_count/tmp_count/eukarya/$tmpeukarya[-1]") {mkdir "$workdir/sum_gene_count/tmp_count/eukarya/$tmpeukarya[-1]", 0777 or die "can not make dir  $workdir/sum_gene_count/tmp_count/eukarya/$tmpeukarya[-1] $!";}
     unless (-d "$workdir/sum_gene_count/read_count/eukarya/$tmpeukarya[-1]") {mkdir "$workdir/sum_gene_count/read_count/eukarya/$tmpeukarya[-1]", 0777 or die "can not make dir  $workdir/sum_gene_count/read_count/eukarya/$tmpeukarya[-1] $!";}
     unless (-d "$workdir/differential_gene/eukarya/$tmpeukarya[-1]") {mkdir "$workdir/differential_gene/eukarya/$tmpeukarya[-1]", 0777 or die "can not make dir  $workdir/differential_gene/eukarya/$tmpeukarya[-1] $!";}
    print LOGFILE"perl /users/203270/code/bin/parse_eukarya_gfffile.pl $tmpgff $workdir/differential_gene/eukarya/$tmpeukarya[-1]/ $workdir/eukarya.fa.fai \n";
    `perl /users/203270/code/bin/parse_eukarya_gfffile.pl $tmpgff $workdir/differential_gene/eukarya/$tmpeukarya[-1]/ $workdir/eukarya.fa.fai`;
    print LOGFILE"python /users/203270/scratch/hisat-master/extract_splice_sites.py $workdir/differential_gene/eukarya/$tmpeukarya[-1]/eukarya.gff > $workdir/differential_gene/eukarya/$tmpeukarya[-1]/splice_sites_gff.txt\n";
    `python /users/203270/scratch/hisat-master/extract_splice_sites.py $workdir/differential_gene/eukarya/$tmpeukarya[-1]/eukarya.gtf > $workdir/differential_gene/eukarya/$tmpeukarya[-1]/splice_sites_gff.txt`;
     if (&file_check("$workdir/differential_gene/eukarya/$tmpeukarya[-1]/splice_sites_gff.txt")>0) {
    print LOGFILE"cat $workdir/differential_gene/eukarya/$tmpeukarya[-1]/splice_sites_gff.txt >> $workdir/differential_gene/eukarya/splice_sites_gff.txt\n";
    `cat $workdir/differential_gene/eukarya/$tmpeukarya[-1]/splice_sites_gff.txt >> $workdir/differential_gene/eukarya/splice_sites_gff.txt`; 
   }
  }   
}


if ($test eq 'both' || $test eq 'prokaryote' )
{
     unless (-d "$workdir/sum_gene_count/tmp_count/prokaryote") {mkdir "$workdir/sum_gene_count/tmp_count/prokaryote", 0777 or die "can not make dir  $workdir/sum_gene_count/tmp_count/prokaryote $!";}
     unless (-d "$workdir/sum_gene_count/read_count/prokaryote") {mkdir "$workdir/sum_gene_count/read_count/prokaryote", 0777 or die "can not make dir  $workdir/sum_gene_count/read_count/prokaryote $!";}
     unless (-d "$workdir/differential_gene/prokaryote") {mkdir "$workdir/differential_gene/prokaryote", 0777 or die "can not make dir  $workdir/differential_gene/prokaryote $!";}

    my @tmpgff=split /,/, $gff_prokaryote; 
     foreach (@tmpgff) {
     my $tmpgff=$_;
          my @tmpprokaryote=split '/', $tmpgff;
       $tmpprokaryote[-1] =~ s/\.//g;
     $tmpprokaryote[-1] =~ s/gff//g;
       if(&file_check($tmpgff)<0 )      {
       next "The gff_prokaryote file $tmpgff doesn't exist or empty.\n";
     }
      push @{$allgff{prokaryote}},  $tmpprokaryote[-1];


      `/edge_dev/edge_ui/JBrowse/bin/flatfile-to-json.pl --gff $tmpgff --type CDS  --tracklabel CDS --out $workdir/Jbrowse`;
       `/edge_dev/edge_ui/JBrowse/bin/flatfile-to-json.pl --gff $tmpgff --type tRNA  --tracklabel tRNA --out $workdir/Jbrowse`;
        `/edge_dev/edge_ui/JBrowse/bin/flatfile-to-json.pl --gff $tmpgff --type exon  --tracklabel exon --out $workdir/Jbrowse`;
         `/edge_dev/edge_ui/JBrowse/bin/flatfile-to-json.pl --gff $tmpgff --type gene  --tracklabel gene --out $workdir/Jbrowse`; 


     unless (-d "$workdir/sum_gene_count/tmp_count/prokaryote/$tmpprokaryote[-1]") {mkdir "$workdir/sum_gene_count/tmp_count/prokaryote/$tmpprokaryote[-1]", 0777 or die "can not make dir  $workdir/sum_gene_count/tmp_count/prokaryote/$tmpprokaryote[-1] $!";}
     unless (-d "$workdir/sum_gene_count/read_count/prokaryote/$tmpprokaryote[-1]") {mkdir "$workdir/sum_gene_count/read_count/prokaryote/$tmpprokaryote[-1]", 0777 or die "can not make dir  $workdir/sum_gene_count/read_count/prokaryote/$tmpprokaryote[-1] $!";}
     unless (-d "$workdir/differential_gene/prokaryote/$tmpprokaryote[-1]") {mkdir "$workdir/differential_gene/prokaryote/$tmpprokaryote[-1]", 0777 or die "can not make dir  $workdir/differential_gene/prokaryote/$tmpprokaryote[-1] $!";}
     print LOGFILE"perl /users/203270/code/bin/parse_prokaryote_gfffile.pl $tmpgff  $workdir/differential_gene/prokaryote/$tmpprokaryote[-1]/ $workdir/prokaryote.fa.fai \n";
      print LOGFILE "test000\n";
     `perl /users/203270/code/bin/parse_prokaryote_gfffile.pl $tmpgff $workdir/differential_gene/prokaryote/$tmpprokaryote[-1]/ $workdir/prokaryote.fa.fai`;
      print STDERR "test001\n";
            print LOGFILE "test001\n";
     }
}


my (%pair1, %pair2, %usedexp);

if ($pairfile)
 {
 if(&file_check($pairfile)<0 ) 
   { die "The given pair file $pairfile  doesn't exist or empty.\n";
   } else {
   my $linenum=0;
   open (IN, "$pairfile") or die " can not open  file $pairfile $!";
 while (<IN>)
{
            chomp;
        my @line=split /\t/, $_;
       unless  ($#line ==1) { print LOGFILE"check format of the $pairfile file \n"; exit;}
        $linenum++;
        if ($linenum >1)
  {
     if ($#{$pairs{$line[0]}}>=0 && $#{$pairs{$line[1]}}>=0 ) {
     foreach(@{$pairs{$line[0]}}) {$usedexp{$_}=1;}
     foreach(@{$pairs{$line[1]}}) {$usedexp{$_}=1;}
     push @{$pair1{$line[0]}}, $line[1];
      } else { print LOGFILE"$line[0] or $line[1] not defined in experimental descpition\n"; exit
      }
   } 
 }
}
} else {
   foreach (sort keys %pairs) 
   {
    my $g1=$_;
    $pair2{$g1}=1;
    foreach(@{$pairs{$g1}}) {$usedexp{$_}=1;}
    foreach (sort keys %pairs) 
     {
     unless ( $pair2{$_} ) {push @{$pair1{$g1}}, $_;}
     }
   }
}

    #my $command = "bowtie2-build -f $index_fasta $index_bt2";
   my  $checkIndexFile = join "", ($index_bt2, '.5.bt2l');
   unless ( -s  $checkIndexFile)
   {
        #my $index_fasta1=join ',', ($prokaryote_fasta, $eukarya_fasta);
        # $index_fasta1=~s/:/,/g;
        print colored ("Indexing the reference sequences",'yellow'), "\n";
#       &executeCommand($command);
    #     `bowtie2-build -f $index_fasta1 $index_bt2`;
         if ($eukarya_fasta && $prokaryote_fasta) { print LOGFILE" /users/203270/scratch/hisat-master/hisat-build --large-index $eukarya_fasta,$prokaryote_fasta  $index_bt2";  `/users/203270/scratch/hisat-master/hisat-build --large-index $eukarya_fasta,$prokaryote_fasta  $index_bt2`;}
         elsif ($eukarya_fasta) {print LOGFILE" /users/203270/scratch/hisat-master/hisat-build --large-index $eukarya_fasta  $index_bt2"; `/users/203270/scratch/hisat-master/hisat-build --large-index $eukarya_fasta  $index_bt2`;}
         elsif ($prokaryote_fasta) {print LOGFILE" /users/203270/scratch/hisat-master/hisat-build --large-index $prokaryote_fasta $index_bt2"; `/users/203270/scratch/hisat-master/hisat-build --large-index $prokaryote_fasta  $index_bt2`;}
         else {print LOGFILE"no INDEX files";exit; }
   }

if ( -s  $checkIndexFile) {print LOGFILE"done INDEX $index_bt2\n"; } else {print LOGFILE"failed INDEX $index_bt2\n"; exit;}


foreach (sort keys %description ) {
my $sample=$_;
my $rawreads=$description{$sample}{Rawreads_file};
my $jobname=join '.', ($sample, 'RNA_analyis');

if ($mapread ne 'yes') 
{
#if(&file_check("$workdir/logdir/$sample")>0 ) { `rm $workdir/logdir/$sample`;}

###print LOGFILE"qsub -pe smp $numCPU  -v test=$test -v numCPU=$numCPU -v workdir=$workdir -v sample=$sample -v rawreads=$rawreads -v reference=$index_fasta -v indexref=$index_bt2 -v Gfffile=$Gfffile -v gff_reference_files=$coverage_fasta -o $workdir/logdir/$sample -N $jobname /users/203270/code/bin/Rnaanalysis.sh \n\n";

if ($rna_trimming_opt eq 'yes') {
print LOGFILE"qsub -pe smp $numCPU   -l h_vmem=$memlim -v test=$test -v numCPU=$numCPU -v workdir=$workdir  -v sample=$sample -v rawreads=$rawreads  -v indexref=$index_bt2 -v descriptfile=$descriptfile  -o $workdir/logdir/$sample -N $jobname /users/203270/code/bin/trim_readmapping.sh \n\n";

`qsub -pe smp $numCPU -l h_vmem=$memlim  -v test=$test -v numCPU=$numCPU -v workdir=$workdir -v htseq=$htseq -v sample=$sample -v rawreads=$rawreads  -v indexref=$index_bt2 -v  descriptfile=$descriptfile  -o $workdir/logdir/$sample -N $jobname /users/203270/code/bin/trim_readmapping.sh`;

} else {
  
     my $outDir1=join '/', ($workdir, "$sample");
  mkdir $outDir1 if ( ! -e $outDir1);
  if ( ! -e $outDir1) {print  "can not make dir $outDir1\n";}

     my $troutDir=join '/', ($outDir1, 'trimming_results');
    mkdir $troutDir if ( ! -e $troutDir);
    if ( ! -e $troutDir) {print  "can not make dir $troutDir\n";}


&lprint  ("no need to trim rawreads\n");
  my $tmptrname1=join '.', ($sample, '1.trimmed.fastq');
   my $tmptrname2=join '.', ($sample, '2.trimmed.fastq');

 if ($#{$allrawreads1{$sample}} ==0 )
  {
   `ln -sf $allrawreads1{$sample}[0] $troutDir/$tmptrname1`;
   `ln -sf $allrawreads2{$sample}[0] $troutDir/$tmptrname2`;
  } else {
     my $tmpreads1=join ' ', @{$allrawreads1{$sample}};
     my $tmpreads2=join ' ', @{$allrawreads2{$sample}};
    `cat $tmpreads1 > $troutDir/$tmptrname1`;
    `cat $tmpreads2 > $troutDir/$tmptrname2`;
  }

 print LOGFILE"qsub -pe smp $numCPU   -l h_vmem=$memlim -v test=$test -v numCPU=$numCPU -v workdir=$workdir  -v sample=$sample -v rawreads=$rawreads  -v indexref=$index_bt2 -v descriptfile=$descriptfile  -o $workdir/logdir/$sample -N $jobname /users/203270/code/bin/readmapping.sh \n\n";

`qsub -pe smp $numCPU -l h_vmem=$memlim  -v test=$test -v numCPU=$numCPU -v workdir=$workdir -v htseq=$htseq -v sample=$sample -v rawreads=$rawreads  -v indexref=$index_bt2 -v  descriptfile=$descriptfile  -o $workdir/logdir/$sample -N $jobname /users/203270/code/bin/readmapping.sh`;

 }

 } else {
  
  print LOGFILE"qsub -pe smp $numCPU   -l h_vmem=$memlim -v test=$test -v numCPU=$numCPU -v workdir=$workdir -v sample=$sample -v rawreads=$rawreads  -v indexref=$index_bt2 -v descriptfile=$descriptfile  -o $workdir/logdir/$sample -N $jobname /users/203270/code/bin/parse_BAMfile.sh \n\n";

`qsub -pe smp $numCPU -l h_vmem=$memlim  -v test=$test -v numCPU=$numCPU -v workdir=$workdir  -v sample=$sample -v rawreads=$rawreads  -v indexref=$index_bt2 -v  descriptfile=$descriptfile  -o $workdir/logdir/$sample -N $jobname /users/203270/code/bin/parse_BAMfile.sh`;
 
 }

}


  my $alldone=keys(%allsample);
  print LOGFILE "total $alldone samples\n";
   while ($alldone)
  {
    foreach (sort keys %description ) {
    my $sample=$_;
    my $tmpfile="$workdir/$sample/mapping_results/$sample.stats.text";
        if (&file_check($tmpfile)>0 ) {
        $alldone--;
       print LOGFILE"done samples : $alldone\n";
      } else {print "$tmpfile not done\n"; print LOGFILE "$tmpfile not done\n";}
    }
    if ($alldone>0)
   {
     print LOGFILE"sample unfinished : $alldone\n";
      sleep (60);
      $alldone=keys(%allsample);
     }else {
    last;
   }
 }
print STDERR "test1\n";

foreach (sort keys %description ) {
my $sample=$_;
 if ($prokaryote_fasta) {
print LOGFILE"qsub -o $workdir/logdir/$sample.rRNACoverageFold_prokaryote.log -v sample=$sample  -v workdir=$workdir  /users/203270/code/bin/rRNACoverageFold_prokaryote.sh\n";
`qsub -o $workdir/logdir/$sample.rRNACoverageFold_prokaryote.log  -v sample=$sample  -v workdir=$workdir  /users/203270/code/bin/rRNACoverageFold_prokaryote.sh`; 
  }

 if ($eukarya_fasta) {
print LOGFILE"qsub  -o $workdir/logdir/$sample.rRNACoverageFold_eukaryote.log -v sample=$sample  -v workdir=$workdir  /users/203270/code/bin/rRNACoverageFold_eukaryote.sh\n";
`qsub  -o $workdir/logdir/$sample.rRNACoverageFold_eukaryote.log  -v sample=$sample  -v workdir=$workdir  /users/203270/code/bin/rRNACoverageFold_eukaryote.sh`;
  }
}
  $alldone=keys(%allsample);
  print LOGFILE "total $alldone samples\n";
   while ($alldone)
  {
    foreach (sort keys %description ) {
    my $sample=$_;
    my $tmpfile1="$workdir/$sample/mapping_results/done.eukarya.txt";
    my $tmpfile2="$workdir/$sample/mapping_results/done.prokaryote.txt";
       
      if ($eukarya_fasta && $prokaryote_fasta)
   {
     if (&file_check($tmpfile1)>0 && &file_check($tmpfile2)>0  ) {
        $alldone--;
       print LOGFILE"done samples : $alldone\n";
      } else {print "eukarya and  prokaryote rRNACoverageFold $sample not done\n"; print LOGFILE "eukarya and  prokaryote rRNACoverageFold $sample not done\n";}
   } elsif  ( $prokaryote_fasta) {

        if ( &file_check($tmpfile2)>0  ) {
        $alldone--;
       print LOGFILE"done samples : $alldone\n";
      } else {print "prokaryote rRNACoverageFold $sample not done\n"; print LOGFILE "prokaryote rRNACoverageFold $sample not done\n";}

  }  elsif  ($eukarya_fasta) {

        if (&file_check($tmpfile1)>0  ) {
        $alldone--;
       print LOGFILE"done samples : $alldone\n";
      } else {print "eukarya rRNACoverageFold $sample not done\n"; print LOGFILE "eukarya rRNACoverageFold $sample not done\n";}
  }
}
    if ($alldone>0)
   {
     print LOGFILE"sample unfinished : $alldone\n";
      sleep (60);
      $alldone=keys(%allsample);
     }else {
    last;
   }
 
}

print LOGFILE"done rRNACoverageFold \n";

my $allsample=join ',', @allsample;
my $makegffdone=1;
 if ($prokaryote_fasta) {
if (-e "$workdir/newprokaryotegffmade.txt") {`rm $workdir/newprokaryotegffmade.txt`;}
#print LOGFILE"qsub  -o $workdir/$allsample.prokaryote_find_small_rna.log -v sample=$allsample  -v reffile=$prokaryote_fasta  -v workdir=$workdir -v gfffile=$gff_prokaryote /users/203270/code/bin/prokaryote_find_small_rna.sh\n";
#`qsub -o $workdir/$allsample.prokaryote_find_small_rna.log  -v sample=$allsample -v reffile=$prokaryote_fasta  -v workdir=$workdir -v gfffile=$gff_prokaryote  /users/203270/code/bin/prokaryote_find_small_rna.sh`;

print LOGFILE"perl /users/203270/code/bin/prokaryote_find_small_rna.pl $workdir $prokaryote_fasta  $gff_prokaryote $allsample\n";
`perl /users/203270/code/bin/prokaryote_find_small_rna.pl $workdir $prokaryote_fasta  $gff_prokaryote $allsample`;

   while ($makegffdone) 
  {
     if (-e "$workdir/newprokaryotegffmade.txt")
   {
      $makegffdone=0;
      last;
    } else {
    sleep (60);
    print LOGFILE "doing prokaryote_find_small_rna.sh\n";
    }
  }
 } 
 $makegffdone=1;
if(  &file_check($eukarya_fasta)>0 && &file_check($gff_eukarya)>0 )
 {
if (-e "$workdir/neweukaryotegffmade.txt") {`rm $workdir/neweukaryotegffmade.txt`;}
#print LOGFILE"qsub -o $workdir/prokaryote_find_small_rna.log  -v sample=$allsample   -v reffile=$eukarya_fasta -v workdir=$workdir -v gfffile=$gff_eukarya /users/203270/code/bin/eukaryote_find_small_rna.sh\n";
#`qsub  -o $workdir/prokaryote_find_small_rna.log -v sample=$allsample -v reffile=$eukarya_fasta  -v workdir=$workdir -v gfffile=$gff_eukarya /users/203270/code/bin/eukaryote_find_small_rna.sh`;
#print LOGFILE"perl /users/203270/code/bin/eukaryote_find_small_rna.pl $workdir $eukarya_fasta  $gff_eukarya $allsample\n";
#`perl /users/203270/code/bin/eukaryote_find_small_rna.pl $workdir $eukarya_fasta  $eukarya_fasta $allsample`;
`pwd > $workdir/neweukaryotegffmade.txt`;
    while ($makegffdone)
  {
     if (-e "$workdir/neweukaryotegffmade.txt")
   {
      $makegffdone=0;
      last;
    } else {
    sleep (60);
    print LOGFILE "doing eukaryote_find_small_rna.sh\n";
    }
  }
 }
foreach (sort keys %description ) {
my $sample=$_;

#print LOGFILE"perl /users/203270/code/bin/htseq-count.pl $workdir $sample $test\n";
#`perl /users/203270/code/bin/htseq-count.pl $workdir $sample $test`;

print LOGFILE"qsub -v test=$test  -v sample=$sample  -v workdir=$workdir /users/203270/code/bin/htseq-count.sh\n";
`qsub -o $workdir/$sample/mapping_results/htseq.log  -v test=$test  -v sample=$sample  -v workdir=$workdir /users/203270/code/bin/htseq-count.sh\n`;
}
foreach (sort keys %description ) {
my $sample=$_;

my $jobname=join '.', ($sample, 'RNA_analyis');

my $tmpname=$sample.$description{$sample}{group}.'prokaryote_forward';
 `/edge_dev/edge_ui/JBrowse/bin/add-bw-track.pl --in $workdir/Jbrowse/trackList.json --out $workdir/Jbrowse/trackList.json --plot --label $tmpname  --bw_url Jbrowse/BigWig/$sample.forward.sorted.prokaryote_ref.bw --key $tmpname `;

$tmpname=$sample.$description{$sample}{group}.'prokaryote_backward';
 `/edge_dev/edge_ui/JBrowse/bin/add-bw-track.pl --in $workdir/Jbrowse/trackList.json --out $workdir/Jbrowse/trackList.json --plot --label $tmpname  --bw_url Jbrowse/BigWig/$sample.backward.sorted.prokaryote_ref.bw --key $tmpname `;

 $tmpname=$sample.$description{$sample}{group}.'eukarya_forward';
  `/edge_dev/edge_ui/JBrowse/bin/add-bw-track.pl --in $workdir/Jbrowse/trackList.json --out $workdir/Jbrowse/trackList.json --plot --label $tmpname  --bw_url Jbrowse/BigWig/$sample.forward.sorted.eukarya_ref.bw --key $tmpname `;

 $tmpname=$sample.$description{$sample}{group}.'eukarya_backward';
  `/edge_dev/edge_ui/JBrowse/bin/add-bw-track.pl --in $workdir/Jbrowse/trackList.json --out $workdir/Jbrowse/trackList.json --plot --label $tmpname  --bw_url Jbrowse/BigWig/$sample.backward.sorted.eukarya_ref.bw --key $tmpname `;

 }


  $alldone=keys(%allsample);
  print LOGFILE "total $alldone samples\n";
   while ($alldone)
  {
    foreach (sort keys %description ) {
    my $sample=$_;
    my $tmpfile1="$workdir/$sample/mapping_results/finish_htseq-count.txt";
        if (&file_check($tmpfile1)>0 ) {
        $alldone--;
       print LOGFILE"done samples : $alldone\n";
      } else {print "$sample not done\n"; print LOGFILE "$sample not done\n";}
    }
    if ($alldone>0)
   {
     print LOGFILE"sample unfinished : $alldone\n";
      sleep (60);
      $alldone=keys(%allsample);
     }else {
    last;
   }
 }
my %diffdir;

 opendir(DIR, "$workdir/sum_gene_count/tmp_count/") or die $!;
 while (my $tmpdir = readdir(DIR))
  {
                 next unless -d "$workdir/sum_gene_count/tmp_count/$tmpdir";
        next if ($tmpdir =~ /^\./);
   opendir(DIR1, "$workdir/sum_gene_count/tmp_count/$tmpdir") or die $!; 
     while (my $tmpdir1 = readdir(DIR1)) 
      {
        next unless -d "$workdir/sum_gene_count/tmp_count/$tmpdir/$tmpdir1";
        next if ($tmpdir1 =~ /^\./);
        print LOGFILE"$tmpdir\t$tmpdir1\n";
        $diffdir{$tmpdir}{$tmpdir1}=1;
      }
  }

 sleep (5);

foreach (sort keys %diffdir)
{
 my $kingdom=$_;

unless ($kingdom eq $test || $test eq'both') {next;}
foreach (sort keys %{$diffdir{$kingdom}})
 {
   my $strain=$_;
 print LOGFILE"$kingdom\t$strain\n"; 
sleep(5);
my $rRNAdes="$workdir/differential_gene/$kingdom/$strain/$kingdom.genedesc.rRNA.txt";
my %rRNApos;

 if (&file_check($rRNAdes)>0) { 

open (GENOIN, "$rRNAdes") or die "$rRNAdes $!";
while (<GENOIN>) {

         chomp;

       my @line = split /\s+/,$_;
     for (my $i=0; $i<=$#line; $i++) {
      $rRNApos{$line[0]}=$line[-1];
      }
    }
close GENOIN;
 }

foreach (sort keys %description ) {
my $sample=$_;
open (READIN, ">$workdir/sum_gene_count/read_count/$kingdom/$strain/$sample.$kingdom.name.htseq.locus_tag.txt") or die "$! $workdir/sum_gene_count/read_count/$kingdom/$strain/$sample.$kingdom.name.htseq.locus_tag.txt\n"; 
     open (IN, "$workdir/sum_gene_count/tmp_count/$kingdom/$strain/$sample.$kingdom.name.htseq.locus_tag.txt") or die " can not open description file $workdir/sum_gene_count/tmp_count/$kingdom/$strain/$sample.$kingdom.name.htseq.locus_tag.txt $!";
        while (<IN>)
          {
                chomp;
                my $line=$_;
            if ($line =~m/^no_feature/) {last;}
           my @line=split /\s+/, $line;
         if ($rRNApos{$line[0]})

           {
           } else {
 #          print LOGFILE"$line\n";
           print READIN "$line\n";
           }
        }

     close READIN;
     close IN;
 }

my $Deseqdir = "$workdir/differential_gene/$kingdom/$strain/";

unless (-e "$workdir/differential_gene/$kingdom/$strain/") {next;}

unless (-d "$Deseqdir/Deseq/") {mkdir "$Deseqdir/Deseq/", 0777 or die "can not make dir  $workdir/differential_gene/$kingdom/$strain/Deseq/ $!";}
unless (-d "$Deseqdir/EdgeR/") {mkdir "$Deseqdir/EdgeR/", 0777 or die "can not make dir  $workdir/differential_gene/$kingdom/$strain//EdgeR/ $!";}
unless (-d "$Deseqdir/figures/") {mkdir "$Deseqdir/figures/", 0777 or die "can not make dir  $workdir/differential_gene/$kingdom/$strain//figures/ $!";}
unless (-d "$Deseqdir/significant_gene/") {mkdir "$Deseqdir/significant_gene/", 0777 or die "can not make dir  $workdir/differential_gene/$kingdom/$strain//significant_gene/ $!";}



my %tablecounts;
open (DESC, ">$Deseqdir/reads.table.txt") or die "$! $Deseqdir/reads.table.txt\n";
print DESC "Sample";
my %reads_sample;
my %reads_desc;
my $totalreadsmapped=0;

foreach (sort keys %description ) {
my $sample=$_;
      open (IN, "$workdir/sum_gene_count/read_count/$kingdom/$strain/$sample.$kingdom.name.htseq.locus_tag.txt") or die " can not open description file $workdir/sum_gene_count/read_count/$kingdom/$strain/$sample.$kingdom.name.htseq.locus_tag.txt $!";
        while (<IN>)
          {
                chomp;
                my $line=$_;
           my @line=split /\s+/, $line;
           $tablecounts{$line[0]}{$sample}=$line[1];
           $reads_sample{$sample}+=$line[1];
           $reads_desc{$line[0]}+=$line[1];
           $totalreadsmapped+=$line[1];
        }
      close IN;
 }

my @goodsample; 
my $numsample=keys(%reads_sample);
foreach (sort keys %reads_sample) 
 {
  my $sample=$_;
  if ($reads_sample{$sample}>=0.001*$totalreadsmapped/$numsample && ($totalreadsmapped/$numsample>2 || $totalreadsmapped>2000) && $usedexp{$sample} )
   {
    push @goodsample, $sample;
    print DESC "\t$sample";
   } 
 }
 print DESC "\n";
if ($#goodsample<1) {next;}

 foreach (sort keys %tablecounts) 
 {
  my $tmpgene=$_;
  if ($reads_desc{$tmpgene}<=1) {next;}
 
 print DESC "$tmpgene";
 foreach (@goodsample ) 
   {
   print DESC "\t$tablecounts{$tmpgene}{$_}";
   }
  print DESC "\n";
 
 }
close DESC;

open (DESC, ">$Deseqdir/readcounts.expriment.txt") or die "$! $Deseqdir/readcounts.expriment.txt\n";

 print DESC "ID\tfiles\tgroup\n";

my %goodgroup;

foreach (@goodsample ) {
my $sample=$_; 
my $tmpgroup=$sample;
my @tmpgroup=split /_/, $tmpgroup;
#my $htseqfile="$workdir/sum_gene_count/tmp_count/$sample".'.name'.'.htseq'.'.locus_tag'.'.txt';
my $htseqfile="$workdir/sum_gene_count/read_count/$kingdom/$strain/$sample.$kingdom".'.name'.'.htseq'.'.locus_tag'.'.txt';
print DESC "$sample\t$htseqfile\t$description{$sample}{group}";
 print DESC "\n";
$goodgroup{$description{$sample}{group}}++;
}
close DESC;


print LOGFILE"start differential gene finding using $test_method\n";

chdir("$Deseqdir") or die "$!";

open (DESC, ">$Deseqdir/Deseq_EdgeRpairs.txt") or die "$! $Deseqdir/Deseq_EdgeRpairs.txt\n";

 print DESC "group1\tgroup2\n";

foreach ( sort keys %pair1 ) {
my $group1=$_;
next unless  $goodgroup{$group1};
foreach (@{$pair1{$group1}}) 
  {
  next unless  $goodgroup{$_};
  print DESC "$group1\t$_\n";
  }
}
close DESC;

if ($test_method eq 'both') {

print LOGFILE"EdgeR and Deseq2 \t $Deseqdir\n";

print LOGFILE"Rscript /users/203270/code/bin/EdgeR.R  $p_cutoff\n";
`Rscript /users/203270/code/bin/EdgeR.R  $p_cutoff`;

print LOGFILE"Rscript /users/203270/code/bin/Deseq.R  $p_cutoff\n";
`Rscript /users/203270/code/bin/Deseq.R  $p_cutoff`;


} elsif ($test_method eq 'EdgeR') {
print LOGFILE"EdgeR \t $Deseqdir \n";

print LOGFILE"Rscript /users/203270/code/bin/EdgeR.R  $p_cutoff\n";
`Rscript /users/203270/code/bin/EdgeR.R  $p_cutoff`;
} elsif ($test_method eq 'Deseq') {

print LOGFILE"Deseq2 \t $Deseqdir \n";

print LOGFILE"Rscript /users/203270/code/bin/Deseq.R  $p_cutoff\n";
`Rscript /users/203270/code/bin/Deseq.R  $p_cutoff`;
} else {
print LOGFILE"method $test_method is invalid\n";exit;
}

if ($kingdom eq 'prokaryote') 
 {
print LOGFILE"perl /users/203270/code/bin/Differential_stats_prokaryote.pl $workdir $Deseqdir $descriptfile $workdir/differential_gene/$kingdom/$strain/prokaryote.NonrRNA.genedesc.txt  $p_cutoff @goodsample\n";

`perl /users/203270/code/bin/Differential_stats_prokaryote.pl $workdir $Deseqdir $descriptfile $workdir/differential_gene//$kingdom/$strain/prokaryote.NonrRNA.genedesc.txt $p_cutoff @goodsample`;


          `/edge_dev/edge_ui/JBrowse/bin/flatfile-to-json.pl --gff $workdir/sum_direction_prokaryote_ref.gff --type expressed_intergenic_region  --tracklabel expressed_intergenic_region  --out $workdir/Jbrowse`;
 

 }


if ($kingdom eq 'eukarya')
 {
print LOGFILE"perl /users/203270/code/bin/Differential_stats_eukarya.pl $workdir $Deseqdir $descriptfile $workdir/differential_gene/$kingdom/$strain/eukarya.genedesc.txt  $p_cutoff @goodsample\n";

`perl /users/203270/code/bin/Differential_stats_eukarya.pl $workdir $Deseqdir $descriptfile $workdir/differential_gene/$kingdom/$strain/eukarya.genedesc.txt $p_cutoff @goodsample`;

  `/edge_dev/edge_ui/JBrowse/bin/flatfile-to-json.pl --gff $workdir/sum_direction_eukarya_ref.gff --type expressed_intergenic_region  --tracklabel expressed_intergenic_region  --out $workdir/Jbrowse`;

 }
print LOGFILE"done $kingdom\t$strain\n";


my @workdir=split '/',$workdir;
my $jbrowsed;
if ($workdir[-1]) {
$jbrowsed=$workdir[-1];
} elsif ($workdir[-2]) {
$jbrowsed=$workdir[-2];
} else {
$jbrowsed=$workdir[-3];
}


  `/edge_dev/edge_ui/JBrowse/bin/generate-names.pl --out  $workdir/Jbrowse`;

  if (-e "/edge_dev/edge_ui/JBrowse/data/$jbrowsed") {`unlink /edge_dev/edge_ui/JBrowse/data/$jbrowsed`;}
  `ln -s  $workdir/Jbrowse/ /edge_dev/edge_ui/JBrowse/data/$jbrowsed`;
  print LOGFILE"Jbrowse link is at http://ergatis2.lanl.gov/jbrowse/?data=data/$jbrowsed\n";
  print "Jbrowse link is at http://ergatis2.lanl.gov/jbrowse/?data=data/$jbrowsed\n";

sleep(5);
}
} #foreach (sort keys %diffdir)
close LOGFILE;



sub readfai
{
my $faifile=shift;
my @tmpcontigs;
      open (IN, $faifile) or die "can not open $faifile $!";
        while (<IN>)
          {
                chomp;
                my $line=$_;
           my @line=split /\t+/, $line;
           push @tmpcontigs, $line[0]; 
        }
      close IN;
  return @tmpcontigs;
}

sub Usage
{
#perl ~/code/bin/rRNA_mapping_qsub.pl -test_kingdom prokaryote  -significant_pvalue 0.001  -cpu 10 -exp ~/scratch/momo_Rnaseq/Analysis_BTT/DeRNA_Experimetal_descriptions.txt -d ~/scratch/momo_Rnaseq/Analysis_BTT_new/ -prokaryote_fasta /users/203270/scratch/momo_Rnaseq/db/bowtie2/Bacillus_anthracis__Ames_Ancestor_uid58083.fa -eukarya_fasta /users/203270/scratch/momo_Rnaseq/db/bowtie2/cavPor3.fa  -index_ref_bt2 /users/203270/scratch/momo_Rnaseq/db/bowtie2/Bacillus_anthracis__Ames_Ancestor_uid58083_CAVPor3.fa  -gff_prokaryote /users/203270/scratch/momo_Rnaseq/db/Bacillus_anthracis__Ames_Ancestor_uid58083.gff  -gene_coverage_fasta /users/203270/scratch/momo_Rnaseq/db/bowtie2/Bacillus_anthracis__Ames_Ancestor_uid58083.fa

#perl ~/code/bin/rRNA_mapping_qsub.pl  -test_kingdom prokaryote  -significant_pvalue 0.001  -cpu 10 -exp ~/scratch/momo_Rnaseq/Analysis_BTT_new/BTT_Experimetal_descriptions.txt -d ~/scratch/momo_Rnaseq/Analysis_BTT_2015June/  -prokaryote_fasta /users/203270//scratch/momo_Rnaseq/db/bowtie2/Bacillus_anthracis__Ames_Ancestor_uid58083.fa -eukarya_fasta /users/203270/scratch/momo_Rnaseq/db/bowtie2/cavPor3.fa -index_ref_bt2 /users/203270/scratch/momo_Rnaseq/db/bowtie2/Bacillus_anthracis__Ames_Ancestor_uid58083_CAVPor3i_hisat -gff_prokaryote /users/203270/scratch/momo_Rnaseq/db/Bacillus_anthracis__Ames_Ancestor_uid58083.gff -test_method EdgeR  -gene_coverage_fasta /users/203270/scratch/momo_Rnaseq/db/bowtie2/Bacillus_anthracis__Ames_Ancestor_uid58083.fa -pair_comparison ~/scratch/momo_Rnaseq/Analysis_BTT_2015June/pair_comparision.txt


#~/code/bin/rRNA_mapping_qsub.pl    -significant_pvalue 0.001  -cpu 10 -exp /users/203270/scratch/EHC_bin/EHCExperimetal_descriptions.txt -d /users/203270/scratch/EHC_bin_hisat/ -test_kingdom both -eukarya_fasta /mnt/lustre/refdb/usrdb/human-GRCH38-2015-03-25/sequence/hs_ref_GRCh38.p2_combined.fa -prokaryote_fasta /users/203270/scratch/EHC_bin_hisat/YP_combined.fa -index_ref_bt2 /users/203270/scratch/EHC_bin_hisat/YP_human_bt2 -gff_eukarya /mnt/lustre/refdb/usrdb/human-GRCH38-2015-03-25/gff/ref_GRCh38.p2_top_level.gff3 -gff_prokaryote /users/203270/scratch/EHC_bin_hisat/YP_combined.gff -gene_coverage_fasta /users/203270/scratch/EHC_bin_hisat/YP_combined.fa

#perl ~/code/bin/rRNA_mapping_qsub.pl    -significant_pvalue 0.001  -cpu 10 -exp /users/203270/scratch/ATD_influ_shawn/ATDExperimetal_descriptions.txt -d /users/203270/scratch/ATD_influ_shawn/ -test_kingdom both -test_method EdgeR -eukarya_fasta /mnt/lustre/refdb/usrdb/human-GRCH38-2015-03-25/sequence/hs_ref_GRCh38.p2_combined.fa -prokaryote_fasta /users/203270/scratch/ATD_influ_shawn/refseq/IV_A_PR8_34.fna -index_ref_bt2 /users/203270/scratch/ATD_influ_shawn/IV_A_PR8_34_human_bt2 -gff_eukarya /mnt/lustre/refdb/usrdb/human-GRCH38-2015-03-25/gff/ref_GRCh38.p2_top_level.gff3 -gff_prokaryote /users/203270/scratch/ATD_influ_shawn/refseq/IV_A_PR8_34.gff 

#perl ~/code/bin/rRNA_mapping_qsub.pl -test_kingdom prokaryote  -significant_pvalue 0.001  -cpu 10 -exp /users/203270/scratch/momo_Rnaseq/Analysis_BTT_2015AUG/BTT_Experimetal_descriptions.txt -d ~/scratch/momo_Rnaseq/Analysis_BTT_2015AUG/  -prokaryote_fasta /users/203270//scratch/momo_Rnaseq/db/bowtie2/Bacillus_anthracis__Ames_Ancestor_uid58083.fa -eukarya_fasta /users/203270/scratch/momo_Rnaseq/db/bowtie2/cavPor3.fa -index_ref_bt2 /users/203270/scratch/momo_Rnaseq/db/bowtie2/Bacillus_anthracis__Ames_Ancestor_uid58083_CAVPor3i_hisat -gff_prokaryote /users/203270/scratch/momo_Rnaseq/db/Bacillus_anthracis__Ames_Ancestor_uid58083.gff -test_method EdgeR  -gene_coverage_fasta /users/203270/scratch/momo_Rnaseq/db/bowtie2/Bacillus_anthracis__Ames_Ancestor_uid58083.fa -pair_comparison ~/scratch/momo_Rnaseq/Analysis_BTT_2015AUG/pair_comparision.txt

#perl ~/code/bin/rRNA_mapping_qsub.pl  -test_kingdom both  -significant_pvalue 0.001  -cpu 10 -exp /users/203270/scratch/Influenza_momo/ref/expdesign.txt -d /users/203270/scratch/Influenza_momo/  -prokaryote_fasta /users/203270/scratch/Influenza_momo/ref/Influenza-A-California-07-2009.fna  -eukarya_fasta /users/203270/scratch/Influenza_momo/ref/canis.fna -index_ref_bt2  /users/203270/scratch/Influenza_momo/ref/canis_Influenza-A_hisat  -gff_eukarya /users/203270/scratch/Influenza_momo/ref/canis.gff -gff_prokaryote /users/203270/scratch/Influenza_momo/ref/Influenza-A-California-07-2009.gff  -test_method EdgeR  -gene_coverage_fasta /users/203270/scratch/Influenza_momo/ref/Influenza-A-California-07-2009.fna

# perl ~/code/bin/rRNA_mapping_qsub.pl  -test_kingdom prokaryote  -significant_pvalue 0.001  -cpu 10 -exp /users/203270/scratch/BIR_momo/ref/expdesign.txt -d /users/203270/scratch/BIR_momo/ -prokaryote_fasta /users/203270/scratch/BIR_momo/ref/Burk_E264-Ecoli_K-12.fna   -gene_coverage_fasta /users/203270/scratch/BIR_momo/ref/Burk_E264-Ecoli_K-12.fna  -gff_prokaryote /users/203270/scratch/BIR_momo/ref/Burk_E264.gff,/users/203270/scratch/BIR_momo/ref/Ecoli_K-12.gff -test_method EdgeR -pair_comparison /users/203270/scratch/BIR_momo/ref/group_comp.txt -index_ref_bt2 /users/203270/scratch/BIR_momo/ref/Burk_E264-Ecoli_K-12.fna

# perl ~/code/bin/rRNA_mapping_qsub.pl  -test_kingdom prokaryote  -significant_pvalue 0.001  -cpu 10 -exp /users/203270/scratch/OSI_rna1/OSIExperimetal_descriptions.txt -d /users/203270/scratch/OSI_rna1/ -prokaryote_fasta /users/203270/scratch/OSI_rna1/NC_004350.fna -gff_prokaryote /users/203270/scratch/OSI_rna1/differential_gene/prokaryote/NC_004350/prokaryote.gff -test_method both  -gene_coverage_fasta /users/203270/scratch/OSI_rna1/NC_004350.fna  -eukarya_fasta /users/203270/scratch/OSI_rna/host_Lactobacillus_casei.fa /users/203270/scratch/OSI_rna/Lac_casei-treptococcus_mutans_bt2 -index_ref_bt2   -significant_pvalue 0.001  -cpu 10

     print <<"END";
 Usage: perl $0 [options] -exp exp_descriptfile.txt -d workdir -prokaryote_fasta indexprokaryote.fa -eukarya_fasta indexeukarya.fa -index_ref_bt2 indexfile -gff_prokaryote prokaryote.gff -gene_coverage_ref gene_coverage_reference.fa

  example: 
                       
perl ~/code/bin/rRNA_mapping_qsub.pl  -test_kingdom prokaryote  -significant_pvalue 0.001  -cpu 10 -exp /users/203270/scratch/momo_Rnaseq/Analysis_BTT_2015AUG/BTT_Experimetal_descriptions.txt -d ~/scratch/momo_Rnaseq/Analysis_BTT_2015AUG/  -prokaryote_fasta /users/203270//scratch/momo_Rnaseq/db/bowtie2/Bacillus_anthracis__Ames_Ancestor_uid58083.fa -eukarya_fasta /users/203270/scratch/momo_Rnaseq/db/bowtie2/cavPor3.fa -index_ref_bt2 /users/203270/scratch/momo_Rnaseq/db/bowtie2/Bacillus_anthracis__Ames_Ancestor_uid58083_CAVPor3i_hisat -gff_prokaryote /users/203270/scratch/momo_Rnaseq/db/Bacillus_anthracis__Ames_Ancestor_uid58083.gff -test_method EdgeR  -gene_coverage_fasta /users/203270/scratch/momo_Rnaseq/db/bowtie2/Bacillus_anthracis__Ames_Ancestor_uid58083.fa -pair_comparison ~/scratch/momo_Rnaseq/Analysis_BTT_2015AUG/pair_comparision.txt 

 
        -d            string, working directory where the whole project will be under,absolute path,  must have permission to be created. 
        -gff_prokaryote string, absolute path, prokaryote annotation files in gff format, multiple files seperated by comma, needed for diffrential gene analysis (contigs must be in  mapping reference with the same names)   (optional)                  
        -gff_eukarya string, absolute path, eukarya annotation file in gff format, multiple files seperated by comma, needed for diffrential gene analysis (contigs must be in  mapping reference with the same names)  (optional) 
        -eukarya_fasta: comma-separated list of files with ref sequences,  absolute path, eukarya nucleotide sequence in fasta format (for making bowtie2 mapping index file)  (optional) 
        -prokaryote_fasta: comma-separated list of files with ref sequences, absolute path, prokaryote nucleotide sequence in fasta format, single file  (for making bowtie2 mapping index file)  (optional) 
        -index_ref_bt2:     string, absolute path, bowtie2 mapping #ndex file,  single file that already exists or you want generated from ref sequences for both eukarya and prokaryote fasta. (must have written permission).
        -h_vmem: memory limit per node, string, default 20G
        -gene_coverage_fasta: string, absolute path,  fasta format, single file  (for directional coverage analysis, sequnce  must be part of prokaryote mapping reference sequence)  (optional) 

        -test_kingdom         desired differential gene expression kingdom (both (for both eukarya and prokaryote), prokaryote, or eukarya (default prokaryote));
        -test_method          dessired differential gene expression method (both (for both EdgeR and Deseq2 method), EdgeR, or Deseq (default both)); must have have at least 3 duplicates if using Deseq2.
        -cpu          number of cpu to be used (default 1)
        -BAM_ready      #if mapping file are provided for samples by users (yes or no) default no
        -trimming_reads #is the raw reads needed to be trimmed
        -significant_pvalue: floating number cutoff to define significant differentially express genes, (default =0.001)

        master design text file:
        -exp         tab delimited txt file descripting experiments that each row represents one sample.
                      Each colomum is as:
                    (
                     ID:  uniq sample ID
                     Rawreads_file: absolute path, fastq format, pair reads seperated by colon, multiple datasets seprarate semicolon
                     group:    replicates group name for this project, each group must have uniqe name (without -pair_comparison option defined as below, all groups will be compared to each other in differential gene analysis)
                     experimental condisitons such as clock time, CFU etc:  one condition per name per colomumm, can be multiple colomums, (for differentail gene analysis)
                     
                      Example: (the names and order of the first 3 colomums must be the exact the same as below) :  
                      ID            Rawreads_files/BAM_file                      group 
                      exp1      read1p1:read1p2;read1p1a:read1p2a       time0
                      exp2      read2p1:read2p2                         time0
                      .
                      .
                      .
                      exp21     read21p1:read21p2;read21p1a:read21p2a   timef
                      exp22     read22p1:read22p2;read22p1a:read22p2a   timef                              
                    )
        differential gene analysis  design text file: 
         -pair_comparison       tab delimited txt file descripting pairwise comparison. If this file does NOT exist, all groups defined in the master design text file will be compared to each other in differential gene analysis, 
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
                                 
END
exit;
}

sub file_check 
{
    #check file exist and non zero size
    my $file=shift;
    my $exist=-1;
    if (-e $file) {$exist=1};
    if (-z $file) {$exist=-1};
    return $exist;
}


