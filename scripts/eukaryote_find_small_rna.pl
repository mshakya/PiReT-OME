#! /usr/bin/perl -w 
use File::Basename;

my $workdir=$ARGV[0];
my $reffile=$ARGV[1];
my $gfffile=$ARGV[2];
my $allsample=$ARGV[3];

my @sample=split /,/, $allsample;

my $headfile="$workdir/prokaryote.fa.fai"; # samtool head file

my $mincov=2;
my %seqln;
my %forwardcoverage;
my %backwardcoverage;
my %maximum;
my %minimum;
my $mixcutoff=0.3;
my %gffln;
open (GENOIN, "$headfile") or die "$headfile $!";
 while (<GENOIN>) {

         chomp;

       my @line = split /\s+/,$_;
     for (my $i=1; $i<=$line[1]; $i++) {
     $backwardcoverage{$line[0]}{$i}=0;
     $forwardcoverage{$line[0]}{$i}=0;
     $seqln{$line[0]}=$line[1];
      }
    }
close GENOIN;

foreach (@sample) 
{
my $sample=$_;
unless (  (-e "$workdir/$sample/mapping_results/forward.sorted.eukarya_ref.bedgraph") && (-e "$workdir/$sample/mapping_results/backward.sorted.eukarya_ref.bedgraph") ) {print "$sample doesn't have a eukarya_ref.bedgraph file\n"; next;}
my $forwardpileupfile="$workdir/$sample/mapping_results/forward.sorted.eukarya_ref.bedgraph";#=$ARGV[1]; # samtools pileup snp and indel call output 

my $backwardpileupfile="$workdir/$sample/mapping_results/backward.sorted.eukarya_ref.bedgraph"; #=$ARGV[2];

open (IN2,"$forwardpileupfile") or die ":!";
while(<IN2>)
{
   chomp;
   @array=split /\s+/,$_;
   $forwardcoverage{$array[0]}{$array[1]}+=$array[3];
}
close IN2;

open (IN2,"$backwardpileupfile") or die ":!";
while(<IN2>)
{
   chomp;
   @array=split /\s+/,$_;
   $backwardcoverage{$array[0]}{$array[1]}+=$array[3];
}
close IN2;


}

my @gfffile=split /,/, $gfffile;
foreach (@gfffile)
{
my $tmpfile=$_;
 unless (-e $tmpfile) {next;}
 open (my $fh, "$tmpfile") or die "can not open  $tmpfile $! \n";
        while (<$fh>)
  {
   chomp;
   if ( ( $_ =~ /\s+exon\s+/ ) || ( $_ =~ /\s+CDS\s+/ ) || ( $_ =~ /\s+gene\s+/ ) )
    {
      my @tmpln=split /\s+/, $_;
      if ($seqln{$tmpln[0]})
    {
           for (my $i=$tmpln[3]; $i<=$tmpln[4]; $i++)
      {
      $backwardcoverage{$tmpln[0]}{$i}=0;
      $forwardcoverage{$tmpln[0]}{$i}=0;
      }
    } else { print "$tmpln[0] is not in headfile \n";}
   }

  }
}

foreach  (sort keys %forwardcoverage)
{
 my $id=$_;
  foreach  (sort keys %{$forwardcoverage{$id}})
    {  my $pos=$_;
     my $tmpcov=$forwardcoverage{$id}{$pos};
     $forwardcoverage{$id}{$pos}=$tmpcov/($#sample+1);
    }
}

foreach  (sort keys %backwardcoverage)
{
 my $id=$_;
  foreach  (sort keys %{$backwardcoverage{$id}})
    {  my $pos=$_;
     my $tmpcov=$backwardcoverage{$id}{$pos};
     $backwardcoverage{$id}{$pos}=$tmpcov/($#sample+1);
    }
}


open (OUTG1,">$workdir/sum_direction_eukarya_ref.gff") or die " can not open description file $workdir/sum_direction_eukarya_ref.gff $!";
foreach  (sort keys %forwardcoverage)
{

       my %forward;
       my %backward;
        my %mixdir;

       my %forwardcov;
       my %backwardcov;


         my %deletion;

      my $forward=0;
       my $backward=0;
        my $mixdir=0;
        my $deletion=0;

                my $startif;
              my $endf;
             my $startb;
              my $endb;
             my $startm;
              my $endm;
             my $startd;
              my $endd;

      my $id=$_;
      my $tmpname=$_;
       my %tmpratio;
        my $adjratio;
        my $numadjratio=0;
       for my $pos1 (2..$seqln{$id})

  {
      #  if ($pos>3186) {exit;}
         my $pos=$pos1-1; 
         my $tmpratio;
         if ($forwardcoverage{$id}{$pos} >0 ||  $backwardcoverage{$id}{$pos} >0 )
     {
          $adjratio=($forwardcoverage{$id}{$pos1}+$backwardcoverage{$id}{$pos1})/($forwardcoverage{$id}{$pos}+$backwardcoverage{$id}{$pos});

         if ( $forwardcoverage{$id}{$pos} >= $backwardcoverage{$id}{$pos} )
         {
         $tmpratio=$backwardcoverage{$id}{$pos}/$forwardcoverage{$id}{$pos};
         } else {
         $tmpratio=$forwardcoverage{$id}{$pos}/$backwardcoverage{$id}{$pos};
         }
      } else {
     $adjratio=100000000;
     $tmpratio=1;
     }
         $tmpratio{$pos}=$tmpratio;
         if ($adjratio>5 || $adjratio<0.2) {$numadjratio++;} else {$numadjratio=0;}

             if ( $forwardcoverage{$id}{$pos}>$mincov || $backwardcoverage{$id}{$pos}>$mincov )
    {
             if ( $forwardcoverage{$id}{$pos} >= $backwardcoverage{$id}{$pos} )
        {
             $forward++;
             if ($forward==30) {$startf=$pos-$forward;}
                 if ($forward>=6 && $adjratio<5 && $adjratio>0.2 ) {
                $mixdir=0; $backward=0;$deletion=0;
               if ($forward>=30) {$endf=$pos;$forward{$startf}=$endf;$forwardcov{$startf}+=$forwardcoverage{$id}{$pos}; $backwardcov{$startf}+=$backwardcoverage{$id}{$pos};}
              } else  {
               $mixdir++;
               $startm=$pos-$mixdir;
                if ($mixdir>=30) { $backward=0; $forward=0;$deletion=0}
               if ($mixdir>=30) {$endm=$pos;$mixdir{$startm}=$endm;$forwardcov{$startm}+=$forwardcoverage{$id}{$pos}; $backwardcov{$startm}+=$backwardcoverage{$id}{$pos};}
              }

        } else {
               $backward++;
               if ($backward==30) {$startb=$pos-$backward;}
               if ($backward>=6 &&  $adjratio<5 && $adjratio>0.2) {
               $mixdir=0; $forward=0;$deletion=0;
               if ($backward>=30) {$endb=$pos;$backward{$startb}=$endb;$forwardcov{$startb}+=$forwardcoverage{$id}{$pos}; $backwardcov{$startb}+=$backwardcoverage{$id}{$pos};}
              } else {
               $mixdir++;
                if ($mixdir==30){$startm=$pos-$mixdir;}
               if ($mixdir>=6) {$forward=0; $backward=0;$deletion=0;}
               if ($mixdir>=30) {$endm=$pos;$mixdir{$startm}=$endm;$forwardcov{$startm}+=$forwardcoverage{$id}{$pos}; $backwardcov{$startm}+=$backwardcoverage{$id}{$pos};}
               }
             }


    } else  {
                $deletion++;
               if ($deletion>=9) {$backward=0;$forward=0;  $mixdir=0 } 
   }
  }


          foreach  (sort {$a<=>$b} keys %backward) {
           my $tmplength=$backward{$_}-$_+1;
           my $strand='-';
           my $SenseCoverge=$backwardcov{$_}/$tmplength;
           if ($SenseCoverge<$mincov) {next;}
           my $AntisenseCoverage;
           if ($forwardcov{$_}) { $AntisenseCoverage=$forwardcov{$_}/$tmplength;} else {$AntisenseCoverage=0;}
           my $Antisense2Sense_ratio=$AntisenseCoverage/$SenseCoverge;
           my $CombineCoverage=$SenseCoverge+$AntisenseCoverage;
           my $tmpid=join '_', ($tmpname,$_, $backward{$_});
         print OUTG1 "$id\tbased_on_coverage\tgene\t$_\t$backward{$_}\t.\t$strand\t0\tID=$tmpid;Length=$tmplength;SenseCoverge=$SenseCoverge;AntisenseCoverage=$AntisenseCoverage;CombineCoverage=$CombineCoverage;Direction=backward_strand;Antisense2Sense_ratio=$Antisense2Sense_ratio;gene_id=unknown$tmpid;locus_tag=unknown$tmpid\n";
                  print OUTG1 "$id\tbased_on_coverage\texpressed_intergenic_region\t$_\t$backward{$_}\t.\t$strand\t0\tID=$tmpid;Length=$tmplength;SenseCoverge=$SenseCoverge;AntisenseCoverage=$AntisenseCoverage;CombineCoverage=$CombineCoverage;Direction=backward_strand;Antisense2Sense_ratio=$Antisense2Sense_ratio;gene_id=unknown$tmpid;locus_tag=unknown$tmpid\n";
         }

                foreach  (sort {$a<=>$b} keys %forward) {
          my $tmplength=$forward{$_}-$_+1;
           my $strand='+';
           my $SenseCoverge=$forwardcov{$_}/$tmplength;
           if ($SenseCoverge<$mincov) {next;}
           my $AntisenseCoverage;
            if ($backwardcov{$_}){$AntisenseCoverage=$backwardcov{$_}/$tmplength;} else {$AntisenseCoverage=0;}
           my $Antisense2Sense_ratio=$AntisenseCoverage/$SenseCoverge;
           my $CombineCoverage=$SenseCoverge+$AntisenseCoverage;
           my $tmpid=join '_', ($tmpname,$_, $forward{$_});
           print OUTG1 "$id\tbased_on_coverage\tgene\t$_\t$forward{$_}\t.\t$strand\t0\tID=$tmpid;Length=$tmplength;SenseCoverge=$SenseCoverge;AntisenseCoverage=$AntisenseCoverage;CombineCoverage=$CombineCoverage;Direction=forward_strand;Antisense2Sense_ratio=$Antisense2Sense_ratio;gene_id=unknown$tmpid;locus_tag=unknown$tmpid\n";
         print OUTG1 "$id\tbased_on_coverage\texpressed_intergenic_region\t$_\t$forward{$_}\t.\t$strand\t0\tID=$tmpid;Length=$tmplength;SenseCoverge=$SenseCoverge;AntisenseCoverage=$AntisenseCoverage;CombineCoverage=$CombineCoverage;Direction=forward_strand;Antisense2Sense_ratio=$Antisense2Sense_ratio;gene_id=unknown$tmpid;locus_tag=unknown$tmpid\n";

         }

        foreach  (sort {$a<=>$b} keys %mixdir) {
           my $tmplength=$mixdir{$_}-$_+1;
           my $strand='?';
           my $ForwardCoverage=$forwardcov{$_}/$tmplength;
           my $BackCoverage=$backwardcov{$_}/$tmplength;
           my $CombineCoverage=$ForwardCoverage+$BackCoverage;
           my $Foorward2Back_ratio=$ForwardCoverage/$BackCoverage;
            if ($CombineCoverage<$mincov) {next;}
           my $tmpid=join '_', ($tmpname,$_, $mixdir{$_});
           print OUTG1 "$id\tbased_on_coverage\tgene\t$_\t$mixdir{$_}\t.\t$strand\t0\tID=$tmpid;Length=$tmplength;ForwardCoverge=$ForwardCoverage;BackwardCoverage=$BackCoverage;CombineCoverage=$CombineCoverage;Direction=undef;gene_id=unknown$tmpid;locus_tag=unknown$tmpid\n";
         print OUTG1 "$id\tbased_on_coverage\texpressed_intergenic_region\t$_\t$mixdir{$_}\t.\t$strand\t0\tID=$tmpid;Length=$tmplength;ForwardCoverge=$ForwardCoverage;BackwardCoverage=$BackCoverage;CombineCoverage=$CombineCoverage;Direction=undef;gene_id=unknown$tmpid;locus_tag=unknown$tmpid\n";
             }

   }

  close OUTG1;

  if (-e "$workdir/sum_direction_eukarya_ref.gff")
   {
       foreach (@gfffile)
        {
          my $tmpfile=$_;  
           `cat $workdir/sum_direction_eukarya_ref.gff >> $tmpfile`;
        }   
   opendir(DIR, "$workdir/differential_gene/eukarya/") or die $!;
   while (my $tmpdir = readdir(DIR))
   {
       next unless -d "$workdir/differential_gene/eukarya/$tmpdir";
       next if ($tmpdir =~ /^\./);
             my $tmpgff="$workdir/differential_gene/eukarya/$tmpdir/eukarya.gff";
       `cat $workdir/sum_direction_prokaryote_ref.gff >> $tmpgff`;
#           `/mnt/lustre/scratch/ergatis/apps/JBrowse/bin/flatfile-to-json.pl --gff $tmpgff --type CDS  --tracklabel CDS --out $workdir/Jbrowse`;
#      `/mnt/lustre/scratch/ergatis/apps/JBrowse/bin/flatfile-to-json.pl --gff $tmpgff --type tRNA  --tracklabel tRNA --out $workdir/Jbrowse`;
#       `/mnt/lustre/scratch/ergatis/apps/JBrowse/bin/flatfile-to-json.pl --gff $tmpgff --type exon  --tracklabel exon --out $workdir/Jbrowse`;
#        `/mnt/lustre/scratch/ergatis/apps/JBrowse/bin/flatfile-to-json.pl --gff $tmpgff --type gene  --tracklabel gene --out $workdir/Jbrowse`;

   }
  } 

open (OUTG1,">$workdir/neweukaryagffmade.txt") or die " can not open description file $workdir/neweukaryagffmade.txt$!";
print OUTG1 "done\n";
