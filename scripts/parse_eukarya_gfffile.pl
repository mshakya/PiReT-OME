#!/usr/bin/perl -W
use strict;

my $gfffile=$ARGV[0];
my $workdir=$ARGV[1];
my $failfile=$ARGV[2];

#my $workdir='/users/203270/scratch/momo_Rnaseq/Analysis_BTT/';
#my $gfffile='/users/203270/scratch/momo_Rnaseq/db/Bacillus_anthracis__Ames_Ancestor_uid58083_original.gff';
my $rRNArm;

$rRNArm=1;

my %genedesc;
my %deleline;
my $linecount=0;
my %printline;

my %seqhash;
 open (my $fai, "$failfile") or die "can not open  $failfile $! \n";
   my    $locuscount=0;
        while (<$fai>)
            {
            my @line=split /\t+/, $_;
            $seqhash{$line[0]}=1;
            }

close $fai;

open (DESC, ">$workdir/eukarya.gff") or die "$! $workdir/eukarya.gff\n";
open (GTF, ">$workdir/eukarya.gtf") or die "$! $workdir/eukarya.gtf\n";
 open (my $fh, "$gfffile") or die "can not open  $gfffile $! \n";
        while (<$fh>)
          {
          
           $linecount++;
                chomp;
            if ($_ =~ /^#/) {next;}
              if ( $_ =~ /\tgene\t/ ||  $_ =~ /\texon\t/ ) {
                        
               my @tmpln=split /\t/, $_;

                foreach (sort keys %seqhash) { my $tmp1=$_; $tmp1 =~ s/\|//g; my $tmp2=$tmpln[0]; $tmp2 =~ s/\|//g; if ( $tmp1 =~ /$tmp2/ || $tmp2 =~ /$tmp1/ ) { $tmpln[0]=$_;  } }

               my $tmpln=$tmpln[4]-$tmpln[3]+1;
	       my @tmpname=split /\;/, $tmpln[-1];
               my $tmpname;
               my $genename;
               my $transcript_id;
               my $gene_id;
               my $description;
               my $GeneID;
                  for( my $i=0; $i<= $#tmpname; $i++)
              {

                               unless ( $tmpname[$i]  =~ /=/ )
                   {
                print "please check gff file $gfffile format at line $linecount \n $tmpln[-1]\n";
                 $tmpname[$i]='unknownfeature'.$i.'='.$tmpname[$i];
                   print "$tmpname[$i]\n";
                   }


                if ( $tmpname[$i]  =~ /GeneID:/ )
                   {
                      my @genename=split 'GeneID:', $tmpname[$i];
                      my @GeneID=split ',', $genename[1];
                      $GeneID=$GeneID[0];
                   }
                 if($tmpname[$i]  =~ /gene_id=/)
                   {
                    $gene_id=$tmpname[$i]; 
                    $gene_id=~s/gene_id=//g;
                   } 
                 if($tmpname[$i]  =~ /ID=/) {
                   $tmpname=$tmpln[0].$tmpname[$i]; 
                  } else {
                $tmpname=$tmpln[0].$linecount;
                  }
                if ( $tmpname[$i]  =~ /product=/ )
                  {
                  $description=$tmpname[$i];
                  $description=~s/product=//g;
                  } 
                   if ( $tmpname[$i]  =~ /transcript_id=/ )
                  {
                  $transcript_id=$tmpname[$i];
                  $transcript_id=~s/transcript_id=//g;
                  }

             
               }
                  unless ($gene_id)
                {
                 if ($GeneID) { $gene_id=$GeneID;} else {$gene_id=$tmpname;}
                   $tmpln[-1]=join ';', (@tmpname, "gene_id=$gene_id");
                }

              my $tgfline;
                $genename=$gene_id;
               if ($genename)
              {
               $tmpname=$genename;
               unless ($transcript_id) {$transcript_id=$gene_id;}
                $tgfline=join '"',('gene_id ',$gene_id,';transcript_id ',$transcript_id,';')  
             }
               my  $tmpline=join "\t", @tmpln;
               print DESC "$tmpline\n";
               if ($tgfline) 
              {
              $tmpln[-1]=$tgfline; 
              $tmpline=join "\t", @tmpln;
              print GTF "$tmpline\n";
              }

                $tmpname=~s/"//g;
                $tmpname=~s/gene\=//g;  
               if ($description) {$genedesc{description}{$tmpname}=$description;}else{ $genedesc{description}{$tmpname}='unknown';}
               $genedesc{geneln}{$tmpname}+=$tmpln;
               $genedesc{gene_id}{$tmpname}=$tmpname;
             }
          } 
    close $fh;
    close DESC;
    close GTF;

#`python /users/203270/scratch/hisat-master/extract_splice_sites.py $workdir/eukarya.gff  > $workdir/splice_sites_gff.txt`;

open (DESC, ">$workdir/eukarya.genedesc.txt") or die "$! $workdir/eukarya.genedesc.txt\n";
foreach (sort { ($genedesc{gene_id}{$a} cmp $genedesc{gene_id}{$b})} keys %{$genedesc{gene_id}} )
 {
   print DESC "$_\t$genedesc{description}{$_}\t$genedesc{geneln}{$_}\n";
}
close DESC;


