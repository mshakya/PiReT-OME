#! /usr/bin/perl -w 
use File::Basename;

use FindBin qw($Bin);

$ENV{PATH} = "$Bin:$Bin/../:$ENV{PATH}";

my $sample=$ARGV[0];

my $workdir=$ARGV[1];

my $process_log_file="$workdir/process.log";
open ( LOG, ">>", $process_log_file) or die "failed: failed to write $process_log_file\n$!";

my $mincov=10;

 my $genedection=join '.', ('direction_prokaryote_ref', $sample,  'gff');
 open (OUTG1,">/$workdir/$sample/mapping_results/$genedection") or die " can not open description file /$workdir/$sample/mapping_results/$genedection $!";


my $headfile="$workdir/prokaryote.fa.fai"; # samtool head file


if (scalar(@ARGV)<2){              
   print STDERR "Usage:\nperl $0  samplename workdir replicate_file\n" ;
   exit;
}



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
      }

     $seqln{$line[0]}=$line[1];
    }
close GENOIN;

#`cp $workdir/prokaryote.fa.fai $workdir/coverage.fa.fai`; 

if(&file_check("$workdir/prokaryote.fa.fai")<0 ) {print LOG "failed: $workdir/prokaryote.fa.fai doesn't exist\n"; exit;}
if(&file_check("$workdir/$sample/mapping_results/forward.prokaryote_ref.bam")<0 &&  &file_check("$workdir/$sample/mapping_results/backward.prokaryote_ref.bam")<0 && &file_check("$workdir/$sample/mapping_results/forward.prokaryote_ref.sam")>0 &&  &file_check("$workdir/$sample/mapping_results/backward.prokaryote_ref.sam")>0)
 {
print LOG 
"          
samtools view -bt $workdir/prokaryote.fa.fai  $workdir/$sample/mapping_results/backward.prokaryote_ref.sam  > $workdir/$sample/mapping_results/backward.prokaryote_ref.bam\n
samtools view -bt $workdir/prokaryote.fa.fai  $workdir/$sample/mapping_results/forward.prokaryote_ref.sam  > $workdir/$sample/mapping_results/forward.prokaryote_ref.bam\n
";

`samtools view -bt $workdir/prokaryote.fa.fai  $workdir/$sample/mapping_results/backward.prokaryote_ref.sam  > $workdir/$sample/mapping_results/backward.prokaryote_ref.bam`;
`samtools view -bt $workdir/prokaryote.fa.fai  $workdir/$sample/mapping_results/forward.prokaryote_ref.sam  > $workdir/$sample/mapping_results/forward.prokaryote_ref.bam`;
#unlink "$workdir/$sample/mapping_results/backward.prokaryote_ref.sam";
#unlink "$workdir/$sample/mapping_results/forward.prokaryote_ref.sam";
}
if(&file_check("$workdir/$sample/mapping_results/sortedpos.forward.prokaryote_ref.bam")<0 &&  &file_check("$workdir/$sample/mapping_results/sortedpos.backward.prokaryote_ref.bam")<0 && &file_check("$workdir/$sample/mapping_results/forward.prokaryote_ref.bam")>0 &&  &file_check("$workdir/$sample/mapping_results/backward.prokaryote_ref.bam")>0)  
  {

print LOG
"          
samtools sort $workdir/$sample/mapping_results/backward.prokaryote_ref.bam -o $workdir/$sample/mapping_results/sortedpos.backward.prokaryote_ref.bam\n
samtools sort $workdir/$sample/mapping_results/forward.prokaryote_ref.bam -o $workdir/$sample/mapping_results/sortedpos.forward.prokaryote_ref.bam\n
";


`samtools sort $workdir/$sample/mapping_results/backward.prokaryote_ref.bam -o $workdir/$sample/mapping_results/sortedpos.backward.prokaryote_ref.bam`;
`samtools sort $workdir/$sample/mapping_results/forward.prokaryote_ref.bam -o $workdir/$sample/mapping_results/sortedpos.forward.prokaryote_ref.bam`;
#unlink "$workdir/$sample/mapping_results/backward.prokaryote_ref.bam";
#unlink "$workdir/$sample/mapping_results/forward.prokaryote_ref.bam";
}

print LOG
"          
samtools view -H $workdir/$sample/mapping_results/sortedpos.forward.prokaryote_ref.bam \n
";

        open GSZ, ">$workdir/$sample/mapping_results/sortedpos.forward.prokaryote_ref.genomeSize" or die $!;
        open REF, "samtools view -H $workdir/$sample/mapping_results/sortedpos.forward.prokaryote_ref.bam |" or die $!;
        while(<REF>){
            if( /SN:(\S+)\tLN:(\d+)/ ){
                print GSZ "$1\t$2\n";
            }
        }
        close REF;
        close GSZ;

print LOG
"          
samtools view -H $workdir/$sample/mapping_results/sortedpos.backward.prokaryote_ref.bam\n
";

        open GSZ, ">$workdir/$sample/mapping_results/sortedpos.backward.prokaryote_ref.genomeSize" or die $!;
        open REF, "samtools view -H $workdir/$sample/mapping_results/sortedpos.backward.prokaryote_ref.bam |" or die $!;
        while(<REF>){
            if( /SN:(\S+)\tLN:(\d+)/ ){
                print GSZ "$1\t$2\n";
            }
        }
        close REF;
        close GSZ;

print LOG
" 
genomeCoverageBed -split -bg -ibam $workdir/$sample/mapping_results/sortedpos.forward.prokaryote_ref.bam -g $workdir/$sample/mapping_results/sortedpos.forward.prokaryote_ref.genomeSize  > $workdir/$sample/mapping_results/forward.sorted.prokaryote_ref.bedgraph\n";
#wigToBigWig $workdir/$sample/mapping_results/forward.sorted.prokaryote_ref.bedgraph $workdir/$sample/mapping_results/sortedpos.forward.prokaryote_ref.genomeSize   $workdir/Jbrowse/BigWig/$sample.forward.sorted.prokaryote_ref.bw\n


`genomeCoverageBed -split -bg -ibam $workdir/$sample/mapping_results/sortedpos.forward.prokaryote_ref.bam -g $workdir/$sample/mapping_results/sortedpos.forward.prokaryote_ref.genomeSize  > $workdir/$sample/mapping_results/forward.sorted.prokaryote_ref.bedgraph`;
# `wigToBigWig $workdir/$sample/mapping_results/forward.sorted.prokaryote_ref.bedgraph $workdir/$sample/mapping_results/sortedpos.forward.prokaryote_ref.genomeSize   $workdir/Jbrowse/BigWig/$sample.forward.sorted.prokaryote_ref.bw`




print LOG
" 
genomeCoverageBed -split -bg -ibam $workdir/$sample/mapping_results/sortedpos.backward.prokaryote_ref.bam -g $workdir/$sample/mapping_results/sortedpos.backward.prokaryote_ref.genomeSize  > $workdir/$sample/mapping_results/backward.sorted.prokaryote_ref.bedgraph\n";
# wigToBigWig $workdir/$sample/mapping_results/backward.sorted.prokaryote_ref.bedgraph $workdir/$sample/mapping_results/sortedpos.backward.prokaryote_ref.genomeSize  $workdir/Jbrowse/BigWig/$sample.backward.sorted.prokaryote_ref.bw\n



`genomeCoverageBed -split -bg -ibam $workdir/$sample/mapping_results/sortedpos.backward.prokaryote_ref.bam -g $workdir/$sample/mapping_results/sortedpos.backward.prokaryote_ref.genomeSize  > $workdir/$sample/mapping_results/backward.sorted.prokaryote_ref.bedgraph`;
# `wigToBigWig $workdir/$sample/mapping_results/backward.sorted.prokaryote_ref.bedgraph $workdir/$sample/mapping_results/sortedpos.backward.prokaryote_ref.genomeSize  $workdir/Jbrowse/BigWig/$sample.backward.sorted.prokaryote_ref.bw`; 
  
 
my $forwardpileupfile="$workdir/$sample/mapping_results/forward.sorted.prokaryote_ref.bedgraph";#=$ARGV[1]; # samtools pileup snp and indel call output 

my $backwardpileupfile="$workdir/$sample/mapping_results/backward.sorted.prokaryote_ref.bedgraph"; #=$ARGV[2];

my (@array, $total_fold,$covered_base);
open (IN2,"$forwardpileupfile") or die ":!";
while(<IN2>)
{
#Sequence0000000001	6	N	1	^]T	h

   chomp;
   @array=split /\s+/,$_;
   $forwardcoverage{$array[0]}{$array[1]}=$array[3];
  # $total_fold += $array[3] ;
  # $covered_base++ ;
}
close IN2;



open (IN2,"$backwardpileupfile") or die ":!";
while(<IN2>)
{
#Sequence0000000001     6       N       1       ^]T     h

   chomp;
   @array=split /\s+/,$_;
   $backwardcoverage{$array[0]}{$array[1]}=$array[3];
  # $total_fold += $array[3] ;
  # $covered_base++ ;
}
close IN2;


chdir "$workdir/$sample/mapping_results/";

foreach  (sort keys %forwardcoverage)

{

     #  my $tmpratio;
       my %forward;
       my %backward;
        my %mixbdir;
        my %mixfdir;

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
      $tmpname=~s/\|/-/g;



       my $forwardfile=join '.', ('forward_prokaryote_ref', $sample, $tmpname, 'txt');
       my $backwardfile=join '.', ('backward_prokaryote_ref', $sample, $tmpname, 'txt'); 
       my $ratiofile=join '.', ('ratio_prokaryote_ref', $sample, $tmpname, 'txt');
       my $rscriptfile=join '.', ('Rscript_prokaryote_ref', $sample, $tmpname, 'txt');
       my $plotfile=join '.', ('prokaryote_ref_plot', $sample, $tmpname, 'jpg');
     open (OUT,">$forwardfile");
     open (OUT1,">$backwardfile");
           open (OUT2,">$ratiofile");
         
        my %tmpratio;


       for my $pos (1..$seqln{$id})

  {
         my $tmpratio;
         if ($forwardcoverage{$id}{$pos} >0 ||  $backwardcoverage{$id}{$pos} >0 )
     {
         if ( $forwardcoverage{$id}{$pos} >= $backwardcoverage{$id}{$pos} ) 
         {
         $tmpratio=$backwardcoverage{$id}{$pos}/$forwardcoverage{$id}{$pos};
         } else {
         $tmpratio=$forwardcoverage{$id}{$pos}/$backwardcoverage{$id}{$pos};
         } 
      } else {
     $tmpratio=1;
     }  
         $tmpratio{$pos}=$tmpratio;

        print OUT2 " $pos\t$tmpratio\t$forwardcoverage{$id}{$pos}\t$backwardcoverage{$id}{$pos}\n";  
      
             if ($forwardcoverage{$id}{$pos}>$mincov || $backwardcoverage{$id}{$pos}>$mincov)
    {
             if ( $forwardcoverage{$id}{$pos} >= $backwardcoverage{$id}{$pos} )
        {
        if ($tmpratio<=$mixcutoff) 
              { 
             $forward++;
               if ($forward==5) {$startf=$pos-4; $mixdir=0; $backward=0; $deletion=0 }
               if ($forward>=30) {$endf=$pos;$forward{$startf}=$endf;$forwardcov{$startf}+=$forwardcoverage{$id}{$pos}; $backwardcov{$startf}+=$backwardcoverage{$id}{$pos};}
               
              } elsif ($tmpratio<1) {
               $mixdir++;
               if ($mixdir==5) {$startm=$pos-4; $forward=0; $backward=0;$deletion=0 }
               if ($mixdir>=30) {$endm=$pos;$mixfdir{$startm}=$endm;$forwardcov{$startm}+=$forwardcoverage{$id}{$pos}; $backwardcov{$startm}+=$backwardcoverage{$id}{$pos};}
               

              }
          
	} else
            {
	            if ($tmpratio<=$mixcutoff) 
              {
               $backward++;
               if ($backward==5) {$startb=$pos-4; $mixdir=0; $forward=0; $deletion=0}
               if ($backward>=30) {$endb=$pos;$backward{$startb}=$endb;$forwardcov{$startb}+=$forwardcoverage{$id}{$pos}; $backwardcov{$startb}+=$backwardcoverage{$id}{$pos};}
               

              } elsif ($tmpratio<1) {
               $mixdir++;
               if ($mixdir==5) {$startm=$pos-4; $forward=0; $backward=0; $deletion=0}
               if ($mixdir>=30) {$endm=$pos;$mixbdir{$startm}=$endm;$forwardcov{$startm}+=$forwardcoverage{$id}{$pos}; $backwardcov{$startm}+=$backwardcoverage{$id}{$pos};}


              }

             } 
 
    } else  {
    $deletion++;
               if ($deletion==5) {$startd=$pos-4; $forward=0; $backward=0; $mixdir=0}
               if ($deletion>=30) {$endd=$pos;$deletion{$startd}=$endd;}
   }

          print OUT " $pos\t$forwardcoverage{$id}{$pos}\n";
          print OUT1 " $pos\t$backwardcoverage{$id}{$pos}\n";
  }
        my $genedeletion=join '.', ('deletion_prokaryote_ref', $sample, $tmpname, 'txt');
       open (OUTD1,">$genedeletion");

     
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
         }

        foreach  (sort {$a<=>$b} keys %mixfdir) {
           my $tmplength=$mixfdir{$_}-$_+1;
           my $strand='.';
           my $SenseCoverage=$forwardcov{$_}/$tmplength;
           my $AntisenseCoverage=$backwardcov{$_}/$tmplength;
           my $CombineCoverage=$SenseCoverage+$AntisenseCoverage;
           my $Antisense2Sense_ratio=$AntisenseCoverage/$SenseCoverage;
            if ($CombineCoverage<$mincov) {next;}
           my $tmpid=join '_', ($tmpname,$_, $mixfdir{$_});
         print OUTG1 "$id\tbased_on_coverage\tgene\t$_\t$mixfdir{$_}\t.\t$strand\t0\tID=$tmpid;Length=$tmplength;SenseCoverge=$SenseCoverge;AntisenseCoverage=$AntisenseCoverage;CombineCoverage=$CombineCoverage;Direction=forward_strand;Antisense2Sense_ratio=$Antisense2Sense_ratio;gene_id=unknown$tmpid;locus_tag=unknown$tmpid\n";

         }
         
                 foreach  (sort {$a<=>$b} keys %mixbdir) {
           my $tmplength=$mixbdir{$_}-$_+1;
           my $strand='.';
           my $AntisenseCoverage=$forwardcov{$_}/$tmplength;
           my $SenseCoverage=$backwardcov{$_}/$tmplength;
           my $Antisense2Sense_ratio=$AntisenseCoverage/$SenseCoverage;
           my $CombineCoverage=$SenseCoverage+$AntisenseCoverage;
           if ($CombineCoverage<$mincov) {next;}
           my $tmpid=join '_', ($tmpname,$_, $mixbdir{$_});
         print OUTG1 "$id\tbased_on_coverage\tgene\t$_\t$mixbdir{$_}\t.\t$strand\t0\tID=$tmpid;Length=$tmplength;SenseCoverge=$SenseCoverge;AntisenseCoverage=$AntisenseCoverage;CombineCoverage=$CombineCoverage;Direction=backward_strand;Antisense2Sense_ratio=$Antisense2Sense_ratio;gene_id=unknown$tmpid;locus_tag=unknown$tmpid\n";
         }

 

                foreach  (sort {$a<=>$b} keys %deletion) {
           my $tmplength=$deletion{$_}-$_+1;
         print OUTD1 "$id\tdeletion\t$_\t$deletion{$_}\t$tmplength\n";
         }


close OUT;
close OUT1;
close OUT2;
&coverageplot($forwardfile, $backwardfile, $ratiofile, $rscriptfile, $plotfile);

      my ($window_size, $step_size)=(500, 100);
        if ($seqln{$id} <= 1000) {
  ($window_size, $step_size)=(1,1);
  } elsif ($seqln{$id} <= 100000) {
  ($window_size, $step_size)=(100,10);
  } else {
   ($window_size, $step_size)=(1000,100);
 }


   
            my $wdforwardfile=join '.', ('window_plot',$window_size,$step_size,'forward_prokaryote_ref', $sample, $tmpname, 'txt');
       my $wdbackwardfile=join '.', ('window_plot',$window_size,$step_size,'backward_prokaryote_ref', $sample, $tmpname, 'txt');
       my $wdratiofile=join '.', ('window_plot',$window_size,$step_size,'ratio_prokaryote_ref', $sample, $tmpname, 'txt');
       my $wdrscriptfile=join '.', ('window_plot',$window_size,$step_size,'Rscript_prokaryote_ref', $sample, $tmpname, 'txt');
       my $wdplotfile=join '.', ('window_plot',$window_size,$step_size,'prokaryote_ref_plot', $sample, $tmpname, 'jpg');
     open (OUTWD,">$wdforwardfile");
     open (OUTWD1,">$wdbackwardfile");
           open (OUTWD2,">$wdratiofile");

    my $avg_pos = int($window_size/2);
    my ($bwpos_cov,$fwpos_cov,$rtpos_cov);
   my ($bwcov_sum,$fwcov_sum, $rtcov_sum)=(0,0,0);
   my $step=1;
   my ($bwwindow_sum,$fwwindow_sum,$rtwindow_sum)=(0,0,0);
   my ($bwstep_sum,$fwstep_sum,$rtstep_sum)=(0,0,0);
   my (@bwstep_sum,@fwstep_sum,@rtstep_sum);
   my (@bwstep_sum2,@fwstep_sum2,@restep_sum2);
   my ($bwstep_sum2,$step_sum2,$rtstep_sum2)=(0,0,0);


        for my $pos (1..$seqln{$id})
  {
            
         $bwcov_sum += $backwardcoverage{$id}{$pos};$fwcov_sum += $forwardcoverage{$id}{$pos};$rtcov_sum += $tmpratio{$pos};
         $bwstep_sum += $backwardcoverage{$id}{$pos};$fwstep_sum += $forwardcoverage{$id}{$pos};$rtstep_sum += $tmpratio{$pos}; 
      if (($pos % $step_size)==0)
      {
          push @fwstep_sum, $fwstep_sum; $fwstep_sum=0;
          push @bwstep_sum, $bwstep_sum; $bwstep_sum=0;
          push @rtstep_sum, $rtstep_sum; $rtstep_sum=0;
      }
      if ($pos == $window_size)
      {
          $step=1;
          $bwwindow_sum = $bwcov_sum;$fwwindow_sum = $fwcov_sum;$rtwindow_sum = $rtcov_sum;

                  print OUTWD $avg_pos,"\t",$fwwindow_sum/$window_size,"\n";
         print OUTWD1 $avg_pos,"\t",$bwwindow_sum/$window_size,"\n";
         print OUTWD2 $avg_pos,"\t",$rtwindow_sum/$window_size,"\n";

      }

      if ($pos > $window_size){
         $bwstep_sum2 += $backwardcoverage{$id}{$pos};$fwstep_sum2 += $forwardcoverage{$id}{$pos};$rtstep_sum2 += $tmpratio{$pos}; 
         if (($pos-$window_size)%$step_size == 0)
         {
                      push @fwstep_sum2, $fwstep_sum2; $fwstep_sum2=0;
          push @bwstep_sum2, $bwstep_sum2; $bwstep_sum2=0;
          push @rtstep_sum2, $rtstep_sum2; $rtstep_sum2=0;

         }
      }
      if ($pos == ($window_size+$step_size*$step))
      {
          my $bwprevious_step_sum = shift @bwstep_sum;
          my $bwafter_step_sum = shift @bwstep_sum2;
          $bwwindow_sum = $bwwindow_sum + $bwafter_step_sum - $bwprevious_step_sum;
        
              my $fwprevious_step_sum = shift @fwstep_sum;
          my $fwafter_step_sum = shift @fwstep_sum2;
          $fwwindow_sum = $fwwindow_sum + $fwafter_step_sum - $fwprevious_step_sum;


         my $rtprevious_step_sum = shift @rtstep_sum;
          my $rtafter_step_sum = shift @rtstep_sum2;
          $rtwindow_sum = $rtwindow_sum + $rtafter_step_sum - $rtprevious_step_sum;


          $avg_pos = $avg_pos + $step_size;
          $step++;

         print OUTWD $avg_pos,"\t",$fwwindow_sum/$window_size,"\n";
         print OUTWD1 $avg_pos,"\t",$bwwindow_sum/$window_size,"\n";
         print OUTWD2 $avg_pos,"\t",$rtwindow_sum/$window_size,"\n"; 

      }
  }

  close OUTWD;
close OUTWD1;
close OUTWD2;
&coverageplot($forwardfile, $backwardfile, $wdratiofile, $wdrscriptfile, $wdplotfile);
   
}

close OUTG1;

sub coverageplot
{

my $forwardfile=shift;
my $backwardfile=shift;
my $ratiofile=shift;
my $Rscript_plot=shift;
my $plotfile=shift;


open (Rscript, ">$Rscript_plot");

  print Rscript "

jpeg(filename=\"$plotfile\",width=1024,height=640,quality=100)

a1<-read.table(file=\"$forwardfile\")
b1<-read.table(file=\"$backwardfile\")

a<-subset(a1,a1\$V2>0)
b<-subset(b1,b1\$V2>0)

r1<-read.table(file=\"$ratiofile\")
r0<-subset(r1,r1\$V2>0)
r0a<-subset(r0,r0\$V3>=5)
r<-subset(r0a,r0a\$V4>=5)

a_coverage<-(length(a\$V2)/length(a1\$V2))*100
a_mean_fold<-mean(a1\$V2)
a_std_fold<-sd(a1\$V2)

b_coverage<-(length(b\$V2)/length(b1\$V2))*100
b_mean_fold<-mean(b1\$V2)
b_std_fold<-sd(b1\$V2)

r_mean<-mean(r1\$V2)
r_std<-sd(r1\$V2)

# init device
par(mar=c(0,5,5,2), mfrow = c(2,1), oma=c(5,0,2,0))
# setup plotting area
xlim <- range(c(a1\$V1,b1\$V1))
ylim <- range(c(a\$V2,1/b\$V2))

par( xpd=TRUE, cex.main=1.2, cex.lab=1.2, cex.axis=1.2)
plot(a\$V1,a\$V2, log=\"y\", type=\"p\",pch=\".\",col=\"blue\",cex=2,xlab=\"Position\",ylab=\"Coverage (fold)\",ylim=ylim,xlim=xlim,main=\"\")
par(new=T)
lines(b\$V1,1/b\$V2, type=\"p\",pch=\".\",col=\"red\")
par(new=F)
pa<-par(\'usr\')

leg.txt<-c(legend=\"forward\",paste(\"Coverage: \",format(a_coverage,digit=5),\"%\"),paste(\"Average fold: \",format(a_mean_fold,digit=4),\"±\",format(a_std_fold,digit=4)))
leg.txt1<-c(legend=\"backward\",paste(\"Coverage: \",format(b_coverage,digit=5),\"%\"),paste(\"Average fold: \",format(b_mean_fold,digit=4),\"±\",format(b_std_fold,digit=4)))
legend(\"bottomleft\",                       # x-y coordinates for location of the legend
       legend=c(\"forward\", \"backward\"),      # Legend labels
       col=c(\"blue\", \"red\"),   # Color of points or lines
      # pch=c(21,19,19),                 # Point type
       lty=c(1,3),                    # Line type
       lwd=c(1,1)                    # Line width
       )
legend(\"bottom\",leg.txt)
legend(\"bottomright\",leg.txt1)
#tmp<-dev.off()

rylim <- range(c(r\$V2))

par(mar=c(5,5,0,2))
plot(r\$V1,r\$V2,pch=19, col=\"black\",cex=0.2,xlab=\"Position\",ylab=\"Coverage ratio\",main=\"\",xlim=xlim,ylim=rylim)
pa<-par(\'usr\');
leg.txt<-c(paste(\"coverage_ratio\"),paste(\"Average ratio: \",format(r_mean,digit=4),\"±\",format(r_std,digit=4)))
legend(\"top\",leg.txt)
tmp<-dev.off()
";

#";

  system ("R --vanilla --slave --silent < $Rscript_plot  2>/dev/null");
#  unlink "Rscript$$";
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

`pwd > $workdir/$sample/mapping_results/done.prokaryote.txt`;
