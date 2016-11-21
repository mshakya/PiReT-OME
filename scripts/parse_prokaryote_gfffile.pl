#!/usr/bin/perl -W
use strict;

my $gfffile=$ARGV[0];
my $workdir=$ARGV[1];
my $failfile=$ARGV[2];

open (LOG, ">>$workdir/process.log") or die "failed: failed:  can not open log file $workdir/process.log $!";

#my $workdir='/users/203270/scratch/momo_Rnaseq/Analysis_BTT/';
#my $gfffile='/users/203270/scratch/momo_Rnaseq/db/Bacillus_anthracis__Ames_Ancestor_uid58083_original.gff';
my $rRNArm;

$rRNArm=1;

my %genedesc;
my %deleline;
my $linecount=0;
my %printline;
my %seqhash;
 open (my $fai, "$failfile") or die "failed: failed: can not open  $failfile $! \n";
   my    $locuscount=0;
        while (<$fai>)
            {
            my @line=split /\t+/, $_;
            $seqhash{$line[0]}=1; 
            }

close $fai;

my %locus;
 open (my $fh, "$gfffile") or die "failed: failed: can not open  $gfffile $! \n";
       $locuscount=0;
        while (<$fh>)
            {
              chomp;
               my $tmpline;
               $linecount++;
              my @tmpreads;
              my $produc='unknown';
              my $tmpname='comments';
              my  $tmplinef=$_;;#=join "\t", @tmpln;
              if ( $_ =~ /ID=/ ) {
               
               my @tmpln=split /\t/, $_;
               foreach (sort keys %seqhash) { my $tmp1=$_; $tmp1 =~ s/\|//g; my $tmp2=$tmpln[0]; $tmp2 =~ s/\|//g; if ( $tmp1 =~ /$tmp2/ || $tmp2 =~ /$tmp1/ ) { $tmpln[0]=$_;  } }
               $tmpline=join "\t", @tmpln;
#                print "$tmpln[0]\n\n";
                
               my $tmpln=$tmpln[4]-$tmpln[3]+1;
               my $tmppos=join '---', ($tmpln[3], $tmpln[4]);  
               my @tmpname=split /\;/, $tmpln[-1];
                for( my $i=0; $i<= $#tmpname; $i++)
               {
                unless ( $tmpname[$i]  =~ /=/ ) 
                   {
                print LOG "please check gff file $gfffile format at line $linecount \n $tmpln[-1]\n";
                 $tmpname[$i]='unknownfeature'.$i.'='.$tmpname[$i];
                   print LOG "$tmpname[$i]\n";
                   }
               }
               $tmpln[-1]=join ';', @tmpname;
               $tmpname=join "___", ($tmpln[0], $tmpln[4], $tmpln[3]);
               $genedesc{description}{$tmpname}='unknown';
               my $tmpID=$tmpname[0];
              $tmpID=~s/ID=//g;
              $tmpID=~s/\s+//g;
               $genedesc{geneln}{$tmpname}=$tmpln;

             if ( ($tmpline =~ /^\S+\t+\w+\t+gene/i ) ||  ($tmpline =~ /^\S+\t+\w+\t+CDS/i ) || ($tmpline =~ /^\S+\t+\w+\t+tRNA/i) ||  ($tmpline =~ /^\S+\t+\w+\t+mRNA/i) ||  ( $tmpline =~ /^\S+\t+\w+\t+(rRNA)/i) )
          { 
               $locuscount++;

              if ( $tmpline =~ /locus_tag=/)
            {
               my @tmpdes1=split /locus_tag=/, $tmpline;
               my @tmpdes=split /\;/, $tmpdes1[1];
              if ($tmpdes[0]) 
              {
               $genedesc{locus_tag}{$tmpname}=$tmpdes[0];
               $locus{$tmpdes[0]}=1
              }
             }  else {
               my $tmplocus=join '_', ('unknown',$tmpID, $locuscount);
               unless ($genedesc{locus_tag}{$tmpname})
              {
              $genedesc{locus_tag}{$tmpname}=$tmplocus;
               $tmpline = $tmpline.';locus_tag='.$tmplocus;
               $tmplinef=$tmplinef.';locus_tag='.$tmplocus;
              }
             }

               
              $genedesc{ID}{$tmpname}=$tmpID;
              $genedesc{END}{$tmpname}=$tmpln[4];
          }
           if ($rRNArm) 
          {
            if ( ( $tmpline =~ /(product=16S)/  ) ||  ( $tmpline =~ /^\S+\t+\w+\t+(rRNA)/)  ) {     
              $deleline{$tmpname}=1;
              }
          }

            if ( $tmpline =~ /product=/)
          {
               my @tmpdes1=split /product=/, $tmpline;
               my @tmpdes=split /\;/, $tmpdes1[1];
              if ($tmpdes[0])
              {
               $genedesc{description}{$tmpname}=$tmpdes[0];
              }
           }              
       
      
                   if ( $tmpline =~ /description=/)
            {
               my @tmpdes1=split /description=/, $tmpline;
               my @tmpdes=split /\;/, $tmpdes1[1];
              if ($tmpdes[0]) 
              {
               $genedesc{description}{$tmpname}=$tmpdes[0];
              }

           }    
            $printline{$linecount}{$tmpname}=$tmplinef;
          }
        }
    close $fh;
open (DESC, ">$workdir/prokaryote.NonrRNA.genedesc.txt") or die "failed: $! $workdir/prokaryote.NonrRNA.genedesc.txt\n";
open (DESCRRNA, ">$workdir/prokaryote.genedesc.rRNA.txt") or die "failed: $! $workdir/prokaryote.genedesc.rRNA.txt\n";
foreach (sort { ($genedesc{locus_tag}{$a} cmp $genedesc{locus_tag}{$b})} keys %{$genedesc{ID}} )
 {
  if ($deleline{$_}) {
 print DESCRRNA "$genedesc{locus_tag}{$_}\t$genedesc{description}{$_}\t$_---$genedesc{END}{$_}\n";
   } else {
   print DESC "$genedesc{locus_tag}{$_}\t$genedesc{description}{$_}\t$_---$genedesc{END}{$_}\n";
   }
}
close DESC;
close DESCRRNA;

open (DESC, ">$workdir/prokaryote.gff") or die "failed: $! $workdir/prokaryote.gff\n";
open (DESCNONRRNA, ">$workdir/prokaryote.NonrRNA.gff") or die "failed: $! $workdir/prokaryote.NonrRNA.gff\n";
open (DESCRRNA, ">$workdir/prokaryote.rRNA.gff") or die "failed: $! $workdir/prokaryote.rRNA.gff\n";
foreach (sort { $a <=> $b} keys %printline )
{
    my $tmpline=$_;
   foreach (sort  keys %{$printline{$tmpline}} )
 {

     my $tmpname=$_;
 if ($deleline{$tmpname})
  {
 print DESCRRNA "$printline{$tmpline}{$tmpname}\n"; 
  }
   else  {
 print DESCNONRRNA "$printline{$tmpline}{$tmpname}\n";
  }
  print DESC "$printline{$tmpline}{$tmpname}\n";
 }
}
close DESC;
close DESCNONRRNA;
close DESCRRNA;


