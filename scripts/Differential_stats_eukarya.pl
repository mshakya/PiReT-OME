#!/usr/bin/perl -W
use strict;
use FindBin qw($Bin);
use POSIX qw(strftime);

use lib "/lib:/lib64:/usr/lib:/usr/lib64:/usr/local/lib:/opt/apps/anaconda/lib";

$ENV{HTTP_PROXY}="http://proxyout.lanl.gov:8080";
$ENV{HTTPS_PROXY}="http://proxyout.lanl.gov:8080";

$|=1;
$ENV{PATH} = "$Bin:$Bin/hisat-0.1.5-beta/:$Bin/script/:$Bin/bin/:$ENV{PATH}";
my $main_pid = $$;

my $workdir=$ARGV[0];
my $diffgendir=$ARGV[1];#'/users/203270//scratch/momo_Rnaseq/Analysis_BTT/';
my $descriptionfile=$ARGV[3];
my $exprimentfile=$ARGV[2];
my $p_cutoff = $ARGV[4];# 0.001;
my @allsample= @ARGV[5..$#ARGV];

my %genedes;
my %geneln;
my %pathway;

my %allRPKM;
my %replicates;
my %averageRPKM;

open (LOGFILE, ">>$workdir/process.log") or die "can not create log file $workdir/process.log $!";


open (IN, "$exprimentfile") or die " can not open description file $exprimentfile $!";
 my $linenum=0;
 while (<IN>)
{
         chomp;
        my @line=split /\s+/, $_;
       unless  ($#line >=2) { print LOGFILE "failed: check format of the description file $_ \n"; exit;}
        $linenum++;
        if ($linenum >1)
  {
  $replicates{$line[0]}=$line[2];
  }
   if ($line[2] =~ /\s+/) {
   print LOGFILE " failed: please check format of experimental descpition \n group name can not have space as : $line[2]\n";
   exit;
  }
 }
close IN;
my %groupNUM;
foreach (@allsample) { $groupNUM{$replicates{$_}}++;}

open (my $fh, "$descriptionfile") or die "$! $descriptionfile  \n";

        while (<$fh>)
            {
                chomp;
               my $tmpline=$_;
               my @tmpline=split /\t/, $_;
                $geneln{$tmpline[0]}=$tmpline[2];
                $genedes{$tmpline[0]}=$tmpline[1];
            }
       close $fh;

# TODO: this is for pathway analysis, kegg_locus_orgnism is a large file, need to find a way to
# either package it or make user download it or gzip it read from gzip
# open (PATHIN, "./kegg_locus_orgnism.txt") or die "$! ./kegg_locus_orgnism.txt  \n";

#         while (<PATHIN>)
#             {
#                 chomp;
#                my @tmpline=split /\s+/, $_;
#                my $pathmap=join ':', ($tmpline[0],$tmpline[2]);
#                push @{$pathway{$tmpline[1]}}, $pathmap;
#             }
#        close PATHIN;

  
  # foreach (sort {$uniqdes{$b} <=> $uniqdes{$a}} keys %uniqdes ) {print "$_\t$uniqdes{$_}\n";}

 my %readcounts_gene; 
 my %genename;
my @samplename;
my %totalcounts_sample;

 open ( SUMSTAT, ">$diffgendir/sum_exp_stats.txt") or die "$! $diffgendir/sum_exp_stats.txt  \n";
 print  SUMSTAT "Sample_ID\tTotal_reads\tTotal_HighQuality_reads\tTotal_unmapped_reads\tTotal_eukarya_reads\tProper_paired_eukarya_reads\tTotal_prokaryote_reads\tProper_paired_prokaryote_reads\n";
 open ( COUNTIN, "$diffgendir/reads.table.txt") or die "$! $diffgendir/reads.table.txt  \n";
            my $linecount=0;
        while (<COUNTIN>)
         {
                chomp;
                $linecount++;
               my $tmpline=$_;
               my @tmpline=split /\s+/, $_;
               foreach (@tmpline) {$_=~s/\"//g;}
               if ($linecount==1) {@samplename=@tmpline; shift @samplename}# push @{$readcounts{$linecount}}, ('ID', @tmpline);}
               else
            {
 
                 $genename{$linecount}=$tmpline[0];
                for (my $i=1; $i<=$#tmpline; $i++ ) 
                {
                $readcounts_gene{$genename{$linecount}}{$samplename[$i-1]}=$tmpline[$i];
              #   $totalcounts_sample{$samplename[$i]}+=$tmpline[$i];
                }
            }
          }
       close COUNTIN;

my (%total_reads, %trim_reads,%map_reads, %nonmap_reads, %non_gff_ref_reads, %gff_ref_read, %proper_non_gff_ref_reads, %proper_gff_ref_read);

       foreach ( @allsample ) 
     {
        my $sample=$_;
        

           $trim_reads{$sample}=0;
           $nonmap_reads{$sample}=0;
           $map_reads{$sample}=0;
           $gff_ref_read{$sample}=0;
           $non_gff_ref_reads{$sample}=0;
           $proper_gff_ref_read{$sample}=0;
           $proper_non_gff_ref_reads{$sample}=0;


        my $trim_line=0;
        open ( COUNTIN, "$workdir/$sample/trimming_results/$sample.stats.txt") or die "$! $workdir/$sample/trimming_results/$sample.stats.txt  \n";
                while (<COUNTIN>)
         {
                chomp;
           my $tmpline=$_;
           my @tmpline=split /\s+/, $_;
            
            if ($tmpline=~m/^Reads \#:/) 
             {
              $trim_line++;
              if ($trim_line==1) { $total_reads{$sample}=$tmpline[2];}# print "$sample :: $tmpline :::: $tmpline[2]\n";}
           #   if ($trim_line==2) { $trim_reads{$sample}=$tmpline[2];}# print "$sample :: $tmpline :::: $tmpline[2]\n";}
             }
         } 
        close COUNTIN;

#         print "$workdir/$sample/mapping_results/$sample.stats.text\n";
        open ( COUNTIN, "$workdir/$sample/mapping_results/$sample.stats.text") or die "$! $workdir/$sample/mapping_results/$sample.stats.text  \n";
                while (<COUNTIN>)
         {
                chomp;
           my $tmpline=$_;
           my @tmpline=split /\s+/, $_;
            if ($tmpline=~m/^total_reads/) {$trim_reads{$sample}=$tmpline[1]; }
            if ($tmpline=~m/^total_Unmapped_reads/) {$nonmap_reads{$sample}=$tmpline[1]; }
            if ($tmpline=~m/^total_Mapped_reads/) {$map_reads{$sample}=$tmpline[1]; }
            if ($tmpline=~m/^prokaryote_reads/) {$gff_ref_read{$sample}=$tmpline[1]; }
            if ($tmpline=~m/^eukarya_reads/) {$non_gff_ref_reads{$sample}=$tmpline[1]; }
                        if ($tmpline=~m/^proper_paired_prokaryote_reads/) {$proper_gff_ref_read{$sample}=$tmpline[1]; }
            if ($tmpline=~m/^proper_paired_eukarya_reads/) {$proper_non_gff_ref_reads{$sample}=$tmpline[1]; }


           
         }
        close COUNTIN;
         $totalcounts_sample{$sample}=$proper_gff_ref_read{$sample}+$proper_non_gff_ref_reads{$sample};
        print  SUMSTAT "$sample\t$total_reads{$sample}\t$trim_reads{$sample}\t$nonmap_reads{$sample}\t$map_reads{$sample}\t$non_gff_ref_reads{$sample}\t$proper_non_gff_ref_reads{$sample}\t$gff_ref_read{$sample}\t$proper_gff_ref_read{$sample}\n";
     }
   close  SUMSTAT;
      
     my $RPKMfile="RPKM_all_gene.txt";
    open ( RPKMOUT, ">$diffgendir/$RPKMfile") or die "$! $diffgendir/$RPKMfile  \n";
    print RPKMOUT  "ID";
    for (my $i=0; $i<=$#samplename; $i++) {print RPKMOUT "\t$samplename[$i]";   } 
    print  RPKMOUT "\n";
    for (my $i=2; $i<=$linecount; $i++) 
   {
    unless ($geneln{$genename{$i}}) {print LOGFILE "$genename{$i}\n";next;}
   print RPKMOUT "$genename{$i}";
   for (my $k=0; $k<=$#samplename; $k++)  
      {
        my $counts=0; 
        if ($readcounts_gene{$genename{$i}}{$samplename[$k]} && $totalcounts_sample{$samplename[$k]})
      {
       $counts=$readcounts_gene{$genename{$i}}{$samplename[$k]}*1e9/($totalcounts_sample{$samplename[$k]}*$geneln{$genename{$i}});
       $allRPKM{$genename{$i}}{$samplename[$k]}=$counts;
       $averageRPKM{$genename{$i}}{$replicates{$samplename[$k]}}+=$counts;
      }
      print RPKMOUT "\t$counts";
      }
   print RPKMOUT "\n"; 
   }     
   close RPKMOUT;

  my (%sigEdgeRgene, %sigDeseqgene, @EdgeRfile, @Deseqfile);
  my $existedger=0;
  my $existdeseq=0;
    opendir(DIR, "$diffgendir/EdgeR") or die" canot open $diffgendir/EdgeR $!\n";
    while (my $edgfile = readdir(DIR))
   {
        if ($edgfile =~ m/_sig\.txt$/) 
       {
        my $tmpfile=$edgfile;
        $tmpfile=~s/_sig\.txt/\.txt/;



             if  (&file_check("$diffgendir/EdgeR/$edgfile")>0)
       {
         push @EdgeRfile, $edgfile;
         $existedger++;
         my $linecountsig=0;
      open ( COUNTIN, "$diffgendir/EdgeR/$edgfile ") or die "$! $diffgendir/EdgeR/$edgfile  \n";
            while (<COUNTIN>)
         {
                chomp;
               $linecountsig++;
               my @tmpline=split /\s+/, $_;
               foreach (@tmpline) {$_=~s/\"//g;}
               if ($linecountsig>1)
            { 
                 $sigEdgeRgene{$tmpfile}{$tmpline[0]}++;
            }
          }
       close COUNTIN;
        }
       }
   }
  close DIR;

    opendir(DIR, "$diffgendir/Deseq") or die" canot open $diffgendir/ $!\n";
    while (my $desfile = readdir(DIR))
   {
        if ($desfile =~ m/_sig\.txt$/)
       {
        my $tmpfile=$desfile;
         push @Deseqfile, $tmpfile;
        $tmpfile=~s/_sig\.txt/\.txt/;
 
      
       if  (&file_check("$diffgendir/Deseq/$desfile")>0) 
       {
                   push @Deseqfile, $desfile;
         $existdeseq++;
         my $linecountsig=0;
      open ( COUNTIN, "$diffgendir/Deseq/$desfile ") or die "$! $diffgendir/Deseq/$desfile  \n";
            while (<COUNTIN>)
         {
                chomp;
               $linecountsig++;
               my @tmpline=split /\s+/, $_;
               foreach (@tmpline) {$_=~s/\"//g;}
               if ($linecountsig>1)
            {
                 $sigDeseqgene{$tmpfile}{$tmpline[0]}++;
            }
          }
       close COUNTIN;
      }
    }
   }
  close DIR;



   unless (-d "$diffgendir/significant_gene/pathway") {mkdir "$diffgendir/significant_gene/pathway", 0777 or die "$diffgendir/significant_gene/pathway $!";} 

      my (@allpathwaysample, %allcommonname, %alldesdata, %alledgdata);
 
  if ($existdeseq >=1 && $existedger>=1)
 { 
    opendir(DIR, "$diffgendir/EdgeR") or die" canot open $diffgendir/EdgeR $!\n";

    while (my $edgfile = readdir(DIR))
   {

        next unless ($edgfile =~ m/\.txt$/);
        next if ($edgfile =~ m/_sig\.txt$/);
         my $desfile=$edgfile;
        $desfile=~s/EdgeR/Deseq/;
        if  (-f "$diffgendir/Deseq/$desfile")
       {
       #print "$edgfile\t$desfile\n";
     
      my $pathwaysample=$desfile;
     $pathwaysample=~s/Deseq\.//g;
     $pathwaysample=~s/_sig\.txt//g;   
     $pathwaysample=~s/\.txt//g;
      push @allpathwaysample, $pathwaysample;
   #start one compariosn at a time
     my (%allname, %commonname, %desname, %edgname, %sigdesname, %sigedgname, %desdata, %edgdata);

   #read over EdgeR
     my $linecountsig=0;
      my @edgcolname;
      open ( COUNTIN, "$diffgendir/EdgeR/$edgfile ") or die "$! $diffgendir/EdgeR/$edgfile  \n";  
            while (<COUNTIN>)
         {
                chomp;
                $linecountsig++;
               my $tmpline=$_;
               my @tmpline=split /\s+/, $_;
               foreach (@tmpline) {$_=~s/\"//g;}
               if ($linecountsig==1)
               {
                for (my $i=0; $i<=$#tmpline; $i++)
                {
                my $tmpcolname =join "_", ('EdgeR', $tmpline[$i]);
                $edgcolname[$i]=$tmpcolname;
                }
               }
               else
            {
                 $allname{$tmpline[0]}=$tmpline[-2];
                 $edgname{$tmpline[0]}=$tmpline[-2];
                for (my $i=1; $i<=$#tmpline; $i++ )
                {
                #print "$i\t$edgcolname[$i-1]\n";
              push  @{$alledgdata{$tmpline[0]}[$i-1]}, $tmpline[$i];
                $edgdata{$tmpline[0]}[$i-1]=$tmpline[$i];
                }
            }
          }
       close COUNTIN;

            my $heatmapfile=join "_", ('heatmap',$edgfile);
    open ( RPKMOUT, ">$diffgendir/significant_gene/$heatmapfile") or die "$! $diffgendir/significant_gene/$heatmapfile  \n";
    print RPKMOUT  "ID";
    for (my $i=0; $i<=$#samplename; $i++) {print RPKMOUT "\t$samplename[$i]";   }
    print  RPKMOUT "\n";
    foreach (sort {$edgname{$a} <=>$edgname{$b} } keys %edgname ) 
   {
    my $genename=$_; 
     unless ($geneln{$genename} ) {next;}
     unless ($sigEdgeRgene{$edgfile}{$genename} ) {next;} 
     $sigedgname{$genename}=$edgname{$genename};
     delete $edgname{$genename};
   print RPKMOUT "$genename";
   for (my $k=0; $k<=$#samplename; $k++)
      {
        my $counts=0; 
        if ($totalcounts_sample{$samplename[$k]} && $readcounts_gene{$genename}{$samplename[$k]} )
       {
       $counts=$allRPKM{$genename}{$samplename[$k]};#$readcounts_gene{$genename}{$samplename[$k]}*1e9/($totalcounts_sample{$samplename[$k]}*$geneln{$genename});
       }
      print RPKMOUT "\t$counts";
      }
   print RPKMOUT "\n";
   }
   close RPKMOUT;
    my $heatmapfig=join ".", ($heatmapfile, 'pdf');
   

     #read over Deseq
      $linecountsig=0;
      my @descolname;
      open ( COUNTIN, "$diffgendir/Deseq/$desfile ") or die "$! $diffgendir/Deseq/$desfile  \n";
            while (<COUNTIN>)
         {
                chomp;
                $linecountsig++;
               my $tmpline=$_;
               my @tmpline=split /\s+/, $_;
               #splice (@tmpline, 3,1);
               foreach (@tmpline) {$_=~s/\"//g;}
               if ($linecountsig==1) 
               {
                 splice (@tmpline, 3,1);
                for (my $i=0; $i<=$#tmpline; $i++)
                {
                my $tmpcolname =join "_", ('Deseq', $tmpline[$i]);
                $descolname[$i]=$tmpcolname;
                }
               }
               else
            {
                   splice (@tmpline, 4,1);
                 $allname{$tmpline[0]}=$tmpline[-2];
                 $desname{$tmpline[0]}=$tmpline[-2];
                for (my $i=1; $i<=$#tmpline; $i++ )
                {
                my $tmpl=$tmpline[$i];
                if ($i==2) {$tmpl=-$tmpline[$i];}
              push  @{$alldesdata{$tmpline[0]}[$i-1]}, $tmpl;
                $desdata{$tmpline[0]}[$i-1]=$tmpl;
                }
            }
          }
       close COUNTIN;

             $heatmapfile=join "_", ('heatmap',$desfile);
    open ( RPKMOUT, ">$diffgendir/significant_gene/$heatmapfile") or die "$! $diffgendir/significant_gene/$heatmapfile  \n";
    print RPKMOUT  "ID";
    for (my $i=0; $i<=$#samplename; $i++) {print RPKMOUT "\t$samplename[$i]";   }
    print  RPKMOUT "\n";
    foreach (sort {$desname{$a} <=>$desname{$b} } keys %desname )
   # foreach (sort {$desname{$a} cmp $desname{$b} } keys %desname )
   {
    my $genename=$_;
    unless ($geneln{$genename} ) {next;}
    unless ($sigDeseqgene{$desfile}{$genename} ) {next;}
     $sigdesname{$genename}=$desname{$genename};
     delete $desname{$genename};


   print RPKMOUT "$genename";
   for (my $k=0; $k<=$#samplename; $k++)
      {
      my $counts=0;
      if ($totalcounts_sample{$samplename[$k]} && $readcounts_gene{$genename}{$samplename[$k]} )
      {
      $counts=$allRPKM{$genename}{$samplename[$k]};#$readcounts_gene{$genename}{$samplename[$k]}*1e9/($totalcounts_sample{$samplename[$k]}*$geneln{$genename});
      }
      print RPKMOUT "\t$counts";
    }
   print RPKMOUT "\n";
   }
   close RPKMOUT;
     $heatmapfig=join ".", ($heatmapfile, 'pdf');


      #find Intercestion 

       my $sumfile=join "_", ('Sum_allgenes','EdgeR', $desfile);
     open ( SUMOUT, ">$diffgendir/$sumfile") or die "$! $diffgendir/$sumfile  \n";

        my $pathwayfile=join "_", ('EdgeR', $desfile);
               $pathwayfile=~s/\.txt//g;
unless (-d "$diffgendir/significant_gene/pathway/$pathwayfile/") {mkdir "$diffgendir/significant_gene/pathway/$pathwayfile/", 0777 or die "$diffgendir/significant_gene/pathway/$pathwayfile $!";}

     open ( PATHOUT, ">$diffgendir/significant_gene/pathway/$pathwayfile/pathway.stats.txt") or die "$! $diffgendir/significant_gene/pathway/$pathwayfile/pathway.stats.txt  \n";

    
             $heatmapfile=join "_", ('sum','heatmap','EdgeR',$desfile);
    open ( RPKMOUT, ">$diffgendir/significant_gene/$heatmapfile") or die "$! $diffgendir/significant_gene/$heatmapfile  \n";
           foreach (sort {$sigedgname{$a} <=>$sigedgname{$b} } keys %sigedgname )
      { 
       
       my $genename=$_;
       if ( $sigdesname{$genename}) 
          {
            $commonname{$genename}=$sigedgname{$genename}; 
            delete $sigdesname{$genename};
            delete $sigedgname{$genename};
          }
      
      }

           foreach (sort {$sigdesname{$a} <=>$sigdesname{$b} } keys %sigdesname )
      {

       my $genename=$_;
       if ( $sigedgname{$genename})
          {
            $commonname{$genename}=$sigdesname{$genename};
            delete $sigedgname{$genename};
            delete $sigdesname{$genename};
          }

      }

    my $RPKMscatfile=join "_", ('Scatter','RPKM','EdgeR',$desfile);
     my $tmpgroup=$desfile;
     $tmpgroup=~s/\.txt//g;
    $tmpgroup=~s/Deseq\.//g;
    my @tmpgroup=split /_/, $tmpgroup;
    open ( RPKMSCAT, ">$diffgendir/$RPKMscatfile") or die "$! $diffgendir/$RPKMscatfile  \n";
    print RPKMSCAT "gene\t$tmpgroup[0]\t$tmpgroup[1]\tsignificant\n";
    for (sort keys %averageRPKM) 
    {
    my $gene=$_;
    unless ($averageRPKM{$gene}{$tmpgroup[0]} && $groupNUM{$tmpgroup[0]} && $groupNUM{$tmpgroup[1]} && $averageRPKM{$gene}{$tmpgroup[1]}) {next;}
    my $avegroup0=$averageRPKM{$gene}{$tmpgroup[0]}/$groupNUM{$tmpgroup[0]};
    my $avegroup1=$averageRPKM{$gene}{$tmpgroup[1]}/$groupNUM{$tmpgroup[1]};  
    print RPKMSCAT "$gene\t$avegroup0\t$avegroup1\t";
    if ($commonname{$gene}) {print RPKMSCAT "yes\n";} else {print RPKMSCAT "no\n";}  
    }
    close RPKMSCAT;

    #chdir "$diffgendir";
    print LOGFILE "Rscript $Bin/scatt_plot.R  $RPKMscatfile $diffgendir\n";
    `Rscript $Bin/scatt_plot.R   $RPKMscatfile $diffgendir`;


         print RPKMOUT  "ID";
    for (my $i=0; $i<=$#samplename; $i++) {print RPKMOUT "\t$samplename[$i]";   }
    print  RPKMOUT "\n";
 
         print SUMOUT  "Gene\tDescription";
  foreach (@edgcolname )
      {
      print SUMOUT "\t$_";
      }
      foreach (@descolname )
      {
      print SUMOUT "\t$_";
      }
    print  SUMOUT "\n";

         print PATHOUT  "GeneID\t$edgcolname[2]\t$edgcolname[3]\tDescription\t$descolname[4]\tMap_ID\t$descolname[1]\tLOCUS_ID\t$edgcolname[0]\n";
   foreach (sort {$commonname{$a} <=>$commonname{$b} } keys %commonname )
   {
    my $genename=$_;
    unless ($geneln{$genename}) {next;}
   print RPKMOUT "$genename";
   for (my $k=0; $k<=$#samplename; $k++)
      {
       my $counts=0;
       if ($totalcounts_sample{$samplename[$k]} && $readcounts_gene{$genename}{$samplename[$k]} ) 
      {
      $counts=$allRPKM{$genename}{$samplename[$k]};#$readcounts_gene{$genename}{$samplename[$k]}*1e9/($totalcounts_sample{$samplename[$k]}*$geneln{$genename});
      }
      print RPKMOUT "\t$counts";
      }
   print RPKMOUT "\n";
   
  print SUMOUT "EdgeR-Deseq::$pathwaysample";
   print SUMOUT "::$genename\t$genedes{$genename}";
 
  foreach (@{$edgdata{$genename}} ) 
      {
      print SUMOUT "\t$_";
      } 
      foreach (@{$desdata{$genename}} )
      {
      print SUMOUT "\t$_";
      }
     print SUMOUT "\n";

   }

    foreach (sort {$commonname{$a} <=>$commonname{$b} } keys %commonname )
  {
     
     my $genename=$_;
     unless ($pathway{$genename}) {next;}

    @{$allcommonname{$genename}{Deslogfc}}= @{$alldesdata{$genename}[1]};
    @{$allcommonname{$genename}{Edgelogfc}}=@{$alledgdata{$genename}[0]};
    @{$allcommonname{$genename}{Despadj}}=@{$alldesdata{$genename}[4]};
    @{$allcommonname{$genename}{EdgeFDR}}=@{$alledgdata{$genename}[2]};

   foreach (@{$pathway{$genename}} )
  {
   my $mapid=$_;
   
   print PATHOUT "$pathwaysample";
   print PATHOUT "::$genename\t";
   print PATHOUT "$edgdata{$genename}[2]\t$edgdata{$genename}[3]\t$genedes{$genename}\t$desdata{$genename}[4]\t$mapid\t$desdata{$genename}[1]\t$genename\t$edgdata{$genename}[0]\n";
   }
  }
   close PATHOUT;

     foreach (sort {$sigedgname{$a} <=>$sigedgname{$b} } keys %sigedgname )
   {
    my $genename=$_;
   unless ($geneln{$genename} ) {next;}
   delete $edgname{$genename};
        delete $desname{$genename};
  print SUMOUT "EdgeR::$pathwaysample";
  print SUMOUT "::$genename\t$genedes{$genename}";

     foreach (@{$edgdata{$genename}} )
      {
      print SUMOUT "\t$_";
      }
      foreach (@{$desdata{$genename}} )
      {
      print SUMOUT "\t$_";
      }
     print SUMOUT "\n";

   }
 
    foreach (sort {$sigdesname{$a} <=>$sigdesname{$b} } keys %sigdesname )
   {
    my $genename=$_;
 unless ($geneln{$genename} ) {next;}
    delete $desname{$genename};
    delete $edgname{$genename};
  # print RPKMOUT "Deseq::$genename\t";
   for (my $k=0; $k<=$#samplename; $k++)
      
       {
       my $counts=0;
       if ($totalcounts_sample{$samplename[$k]} && $readcounts_gene{$genename}{$samplename[$k]} ) 
      {
       $counts=$allRPKM{$genename}{$samplename[$k]};#$readcounts_gene{$genename}{$samplename[$k]}*1e9/($totalcounts_sample{$samplename[$k]}*$geneln{$genename});
      }
      print RPKMOUT "$counts\t";
      }
   print RPKMOUT "\n";

   print SUMOUT "Deseq::$pathwaysample";
   print SUMOUT "::$genename\t$genedes{$genename}";

     foreach (@{$edgdata{$genename}} )
      {
      print SUMOUT "\t$_";
      }
      foreach (@{$desdata{$genename}} )
      {
      print SUMOUT "\t$_";
      }
     print SUMOUT "\n";
   }


        foreach (sort {$edgname{$a} <=>$edgname{$b} } keys %edgname )
   {
    my $genename=$_;
   unless ($geneln{$genename} ) {next;}
   delete $desname{$genename};
  print SUMOUT "NotSignificant::$pathwaysample";
  print SUMOUT "::$genename\t$genedes{$genename}";

     foreach (@{$edgdata{$genename}} )
      {
      print SUMOUT "\t$_";
      }
      foreach (@{$desdata{$genename}} )
      {
      print SUMOUT "\t$_";
      }
     print SUMOUT "\n";

   }

          foreach (sort {$desname{$a} <=>$desname{$b} } keys %desname )
   {
    my $genename=$_;
   unless ($geneln{$genename} ) {next;}
  print SUMOUT "NotSignificant::$pathwaysample";
  print SUMOUT "::$genename\t$genedes{$genename}";

     foreach (@{$edgdata{$genename}} )
      {
      print SUMOUT "\t$_";
      }
      foreach (@{$desdata{$genename}} )
      {
      print SUMOUT "\t$_";
      }
     print SUMOUT "\n";

   }


   close SUMOUT;   
   close RPKMOUT;
    $heatmapfig=join ".", ($heatmapfile, 'pdf');
    print LOGFILE "$Bin/heatmap_distinctZ_noClust_zeroRowAllow.py --in $diffgendir/significant_gene/$heatmapfile -s log --out $workdir/differential_gene/figures/$heatmapfig\n";
   `$Bin/heatmap_distinctZ_noClust_zeroRowAllow.py --in $diffgendir/significant_gene/$heatmapfile -s log --out $workdir/differential_gene/figures/$heatmapfig`;

    chdir "$diffgendir/significant_gene/pathway/$pathwayfile/";
  
 `R --vanilla --slave --silent  < $Bin/pathview.R  2>/dev/null`;
 #`perl ./json.pl`;

    } else { #end if $workdir/differential_gene/Deseq/$desfile
       print "$diffgendir/Deseq/$desfile does not exist\n";
   }

  } #end while (my $edgfile = readdir(DIR))

unless (-d "$diffgendir/significant_gene/pathway/all_comparison") {mkdir "$diffgendir/significant_gene/pathway/all_comparison", 0777 or die "$diffgendir/significant_gene/pathway/all_compari
son $!";}

chdir "$diffgendir/significant_gene/pathway/all_comparison";




open ( PATHOUT, ">$diffgendir/significant_gene/pathway/all_comparison/pathway.stats.txt") or die "$! $diffgendir/significant_gene/pathway/all_comparison//pathway.stats.txt  \n";


my $allpathwaysample=join ":", @allpathwaysample; 


 print PATHOUT  "GeneID";
 foreach (@allpathwaysample) 
 {
 my $tmpname=join "_", ($_, 'EdgeFDR');
 print PATHOUT "\t$tmpname";
 }
 print PATHOUT" \tDescription";
  foreach (@allpathwaysample) 
 {
 my $tmpname=join "_", ($_, 'Despadj');
 print PATHOUT"\t$tmpname";
 }
 print PATHOUT"\tMap_ID";
   foreach (@allpathwaysample)
 {
 my $tmpname=join "_", ($_, 'Deslogfc');
 print PATHOUT"\t$tmpname";
 }
 print PATHOUT" \tLOCUS_ID";
    foreach (@allpathwaysample)
 {
 my $tmpname=join "_", ($_, 'Edgelogfc');
 print PATHOUT"\t$tmpname";
 }
 print PATHOUT"\n";

   foreach (sort  keys %allcommonname )
{ 
  
my $genename=$_;

 unless ($pathway{$genename}) {next;}


      @{$allcommonname{$genename}{Deslogfc}}= @{$alldesdata{$genename}[1]};
    @{$allcommonname{$genename}{Edgelogfc}}=@{$alledgdata{$genename}[0]};
    @{$allcommonname{$genename}{Despadj}}=@{$alldesdata{$genename}[4]};
    @{$allcommonname{$genename}{EdgeFDR}}=@{$alledgdata{$genename}[2]};


   foreach (@{$pathway{$genename}} )
{
my $mapid=$_;
my $tmpname=join"::", ($allpathwaysample,$genename);
print PATHOUT "$tmpname";
foreach ( @{$allcommonname{$genename}{EdgeFDR}}) 
 {
  print PATHOUT"\t$_"; 
 }
 print PATHOUT"\t$genedes{$genename}";
foreach ( @{$allcommonname{$genename}{Despadj}})
 {
  print PATHOUT"\t$_";
 }
 print PATHOUT"\t$mapid";
 foreach (@{$allcommonname{$genename}{Deslogfc}})
 {
  print PATHOUT"\t$_";
 }
print PATHOUT"\t$genename";
  foreach (@{$allcommonname{$genename}{Edgelogfc}})
 {
  print PATHOUT"\t$_";
 }
 print PATHOUT"\n";
}
}
close PATHOUT;
`R --vanilla --slave --silent  < $Bin/pathview.R  2>/dev/null`;
#`perl ./json.pl`;
}  #if ($existdeseq >=1 && $existedger>=1)


 if ($existdeseq <1 && $existedger>=1)
 {
    opendir(DIR, "$diffgendir/EdgeR") or die" canot open $diffgendir/EdgeR $!\n";

    while (my $edgfile = readdir(DIR))
   {

        next unless ($edgfile =~ m/\.txt$/);
        next if ($edgfile =~ m/_sig\.txt$/);

      my $pathwaysample=$edgfile;
     #$pathwaysample=~s/EdgeR\.//g;
     $pathwaysample=~s/_sig\.txt//g;
     $pathwaysample=~s/\.txt//g;
      push @allpathwaysample, $pathwaysample;
   #start one compariosn at a time
     my (%allname,%commonname,  %edgname,%sigedgname,  %edgdata);

   #read over EdgeR
     my $linecountsig=0;
      my @edgcolname;
      open ( COUNTIN, "$diffgendir/EdgeR/$edgfile ") or die "$! $diffgendir/EdgeR/$edgfile  \n";
            while (<COUNTIN>)
         {
                chomp;
                $linecountsig++;
               my $tmpline=$_;
               my @tmpline=split /\s+/, $_;
               foreach (@tmpline) {$_=~s/\"//g;}
               if ($linecountsig==1)
               {
                for (my $i=0; $i<=$#tmpline; $i++)
                {
                my $tmpcolname =join "_", ('EdgeR', $tmpline[$i]);
                $edgcolname[$i]=$tmpcolname;
                }
               }
               else
            {
                 $allname{$tmpline[0]}=$tmpline[-2];
                 $edgname{$tmpline[0]}=$tmpline[-2];
                for (my $i=1; $i<=$#tmpline; $i++ )
                {
                #print "$i\t$edgcolname[$i-1]\n";
              push  @{$alledgdata{$tmpline[0]}[$i-1]}, $tmpline[$i];
                $edgdata{$tmpline[0]}[$i-1]=$tmpline[$i];
                }
            }
          }
       close COUNTIN;


      my $sumfile=join "_", ('Sum_allgenes', $edgfile);
     open ( SUMOUT, ">$diffgendir/$sumfile") or die "$! $diffgendir/$sumfile  \n";

#        my $pathwayfile=$edgfile;
#               $pathwayfile=~s/\.txt//g;
#unless (-d "$diffgendir/significant_gene/pathway/$pathwayfile/") {mkdir "$diffgendir/significant_gene/pathway/$pathwayfile/", 0777 or die "$diffgendir/significant_gene/pathway/$pathwayfile $!";}


            my $heatmapfile=join "_", ('heatmap',$edgfile);
    open ( RPKMOUT, ">$diffgendir/significant_gene/$heatmapfile") or die "$! $diffgendir/significant_gene/$heatmapfile  \n";
    print RPKMOUT  "ID";
    for (my $i=0; $i<=$#samplename; $i++) {print RPKMOUT "\t$samplename[$i]";   }
    print  RPKMOUT "\n";
    foreach (sort {$edgname{$a} <=>$edgname{$b} } keys %edgname )
   {
    my $genename=$_;
     unless ($geneln{$genename} ) {next;}
     unless ($sigEdgeRgene{$edgfile}{$genename} ) {next;}
     $sigedgname{$genename}=$edgname{$genename};
     $commonname{$genename}=$sigedgname{$genename};
     delete $edgname{$genename};
   print RPKMOUT "$genename";
   for (my $k=0; $k<=$#samplename; $k++)
      {
        my $counts=0;
        if ($totalcounts_sample{$samplename[$k]} && $readcounts_gene{$genename}{$samplename[$k]} )
       {
       $counts=$allRPKM{$genename}{$samplename[$k]};#$readcounts_gene{$genename}{$samplename[$k]}*1e9/($totalcounts_sample{$samplename[$k]}*$geneln{$genename});
       }
      print RPKMOUT "\t$counts";
      }
   print RPKMOUT "\n";
   }
   close RPKMOUT;
    my $heatmapfig=join ".", ($heatmapfile, 'pdf');
   `$Bin/heatmap_distinctZ_noClust_zeroRowAllow.py --in $diffgendir/significant_gene/$heatmapfile -s log --out $diffgendir/figures/$heatmapfig`;

      my $RPKMscatfile=join "_", ('Scatter','RPKM',$edgfile);
     my $tmpgroup=$edgfile;
     $tmpgroup=~s/\.txt//g;
    $tmpgroup=~s/EdgeR\.//g;
    my @tmpgroup=split /_/, $tmpgroup;
    print "@tmpgroup\n";
    open ( RPKMSCAT, ">$diffgendir/$RPKMscatfile") or die "$! $diffgendir/$RPKMscatfile  \n";
    print RPKMSCAT "gene\t$tmpgroup[0]\t$tmpgroup[1]\tsignificant\n";
    for (sort keys %averageRPKM)
    {
    my $gene=$_;
    unless ($averageRPKM{$gene}{$tmpgroup[0]} && $groupNUM{$tmpgroup[0]} && $groupNUM{$tmpgroup[1]} && $averageRPKM{$gene}{$tmpgroup[1]}) {next;}
    my $avegroup0=$averageRPKM{$gene}{$tmpgroup[0]}/$groupNUM{$tmpgroup[0]};
    my $avegroup1=$averageRPKM{$gene}{$tmpgroup[1]}/$groupNUM{$tmpgroup[1]};
    print RPKMSCAT "$gene\t$avegroup0\t$avegroup1\t";
    if ($commonname{$gene}) {print RPKMSCAT "yes\n";} else {print RPKMSCAT "no\n";}
    }
    close RPKMSCAT;

     # chdir "$diffgendir";
    print LOGFILE "Rscript $Bin/scatt_plot.R  $RPKMscatfile $diffgendir\n";
    `Rscript $Bin/scatt_plot.R   $RPKMscatfile $diffgendir`;


             $heatmapfile=join "_", ('sum',$edgfile);
    open ( RPKMOUT, ">$diffgendir/significant_gene/$heatmapfile") or die "$! $diffgendir/significant_gene/$heatmapfile  \n";
 

         print RPKMOUT  "ID";
    for (my $i=0; $i<=$#samplename; $i++) {print RPKMOUT "\t$samplename[$i]";   }
    print  RPKMOUT "\n";

         print SUMOUT  "Gene\tDescription";
  foreach (@edgcolname )
      {
      print SUMOUT "\t$_";
      }
    print  SUMOUT "\n";


 foreach (sort {$commonname{$a} <=>$commonname{$b} } keys %commonname )
   {
    my $genename=$_;
    unless ($geneln{$genename}) {next;}
   print RPKMOUT "$genename";
   for (my $k=0; $k<=$#samplename; $k++)
      {
       my $counts=0;
       if ($totalcounts_sample{$samplename[$k]} && $readcounts_gene{$genename}{$samplename[$k]} )
      {
      $counts=$allRPKM{$genename}{$samplename[$k]};#$readcounts_gene{$genename}{$samplename[$k]}*1e9/($totalcounts_sample{$samplename[$k]}*$geneln{$genename});
      }
      print RPKMOUT "\t$counts";
      }
   print RPKMOUT "\n";

  print SUMOUT "EdgeR::$pathwaysample";
   print SUMOUT "::$genename\t$genedes{$genename}";

  foreach (@{$edgdata{$genename}} )
      {
      print SUMOUT "\t$_";
      }
     print SUMOUT "\n";
   }

        foreach (sort {$edgname{$a} <=>$edgname{$b} } keys %edgname )
   {
    my $genename=$_;
   unless ($geneln{$genename} ) {next;}
  print SUMOUT "NotSignificant::$pathwaysample";
  print SUMOUT "::$genename\t$genedes{$genename}";

     foreach (@{$edgdata{$genename}} )
      {
      print SUMOUT "\t$_";
      }
     print SUMOUT "\n";
   }

  close SUMOUT;
   close RPKMOUT;

         my $pathwayfile=$edgfile;
               $pathwayfile=~s/\.txt//g;
unless (-d "$diffgendir/significant_gene/pathway/$pathwayfile/") {mkdir "$diffgendir/significant_gene/pathway/$pathwayfile/", 0777 or die "$diffgendir/significant_gene/pathway/$pathwayfile $!";}
 open ( PATHOUT, ">$diffgendir/significant_gene/pathway/$pathwayfile/pathway.stats.txt") or die "$! $diffgendir/significant_gene/pathway/$pathwayfile/pathway.stats.txt  \n";
 
  #print PATHOUT  "GeneID\t$edgcolname[2]\t$edgcolname[3]\tDescription\t$descolname[4]\tMap_ID\t$descolname[1]\tLOCUS_ID\t$edgcolname[0]\n";
  print PATHOUT  "GeneID\tDescription\tMap_ID\tLOCUS_ID\t$edgcolname[0]\n";  
   foreach (sort {$commonname{$a} <=>$commonname{$b} } keys %commonname )
  {

     my $genename=$_;
     unless ($pathway{$genename}) {next;}

  #  @{$allcommonname{$genename}{Deslogfc}}= @{$alldesdata{$genename}[1]};
    @{$allcommonname{$genename}{Edgelogfc}}=@{$alledgdata{$genename}[0]};
  #  @{$allcommonname{$genename}{Despadj}}=@{$alldesdata{$genename}[4]};
   # @{$allcommonname{$genename}{EdgeFDR}}=@{$alledgdata{$genename}[2]};

   foreach (@{$pathway{$genename}} )
  {
   my $mapid=$_;

   print PATHOUT "$pathwaysample";
   print PATHOUT "::$genename\t";
   #print PATHOUT "$edgdata{$genename}[2]\t$edgdata{$genename}[3]\t$genedes{$genename}\t$desdata{$genename}[4]\t$mapid\t$desdata{$genename}[1]\t$genename\t$edgdata{$genename}[0]\n";
   print PATHOUT "$genedes{$genename}\t$mapid\t$genename\t$edgdata{$genename}[0]\n";
   }
  }
   close PATHOUT;

   chdir "$diffgendir/significant_gene/pathway/$pathwayfile/";

 `R --vanilla --slave --silent  < $Bin/pathview.R  2>/dev/null`;
# `perl ./json.pl`;

  }#while (my $edgfile = readdir(DIR))

unless (-d "$diffgendir/significant_gene/pathway/all_comparison") {mkdir "$diffgendir/significant_gene/pathway/all_comparison", 0777 or die "$diffgendir/significant_gene/pathway/all_compari
son $!";}

chdir "$diffgendir/significant_gene/pathway/all_comparison";

open ( PATHOUT, ">$diffgendir/significant_gene/pathway/all_comparison/pathway.stats.txt") or die "$! $diffgendir/significant_gene/pathway/all_comparison//pathway.stats.txt  \n";

my $allpathwaysample=join ":", @allpathwaysample;


 print PATHOUT  "GeneID";

# foreach (@allpathwaysample)
# {
# my $tmpname=join "_", ($_, 'EdgeFDR');
# print PATHOUT "\t$tmpname";
# }
 print PATHOUT" \tDescription";
#  foreach (@allpathwaysample)
# {
# my $tmpname=join "_", ($_, 'Despadj');
# print PATHOUT"\t$tmpname";
# }
 print PATHOUT"\tMap_ID";
#   foreach (@allpathwaysample)
# {
# my $tmpname=join "_", ($_, 'Deslogfc');
# print PATHOUT"\t$tmpname";
# }
 print PATHOUT" \tLOCUS_ID";
    foreach (@allpathwaysample)
 {
 my $tmpname=join "_", ($_, 'Edgelogfc');
 print PATHOUT"\t$tmpname";
 }
 print PATHOUT"\n";

   foreach (sort  keys %allcommonname )
{

my $genename=$_;

 unless ($pathway{$genename}) {next;}

#      @{$allcommonname{$genename}{Deslogfc}}= @{$alldesdata{$genename}[1]};
    @{$allcommonname{$genename}{Edgelogfc}}=@{$alledgdata{$genename}[0]};
#    @{$allcommonname{$genename}{Despadj}}=@{$alldesdata{$genename}[4]};
#    @{$allcommonname{$genename}{EdgeFDR}}=@{$alledgdata{$genename}[2]};

   foreach (@{$pathway{$genename}} )
{
my $mapid=$_;
my $tmpname=join"::", ($allpathwaysample,$genename);
print PATHOUT "$tmpname";
#foreach ( @{$allcommonname{$genename}{EdgeFDR}})
# {
#  print PATHOUT"\t$_";
# }
 print PATHOUT"\t$genedes{$genename}";
#foreach ( @{$allcommonname{$genename}{Despadj}})
# {
#  print PATHOUT"\t$_";
# }
 print PATHOUT"\t$mapid";
# foreach (@{$allcommonname{$genename}{Deslogfc}})
# {
#  print PATHOUT"\t$_";
# }
print PATHOUT"\t$genename";
  foreach (@{$allcommonname{$genename}{Edgelogfc}})
 {
  print PATHOUT"\t$_";
 }
 print PATHOUT"\n";
}
}
close PATHOUT;
`R --vanilla --slave --silent  < $Bin/pathview.R  2>/dev/null`;
#`perl ./json.pl`;

}#  if ($existdeseq <1 && $existedger>=1)


 if ($existdeseq >=1 && $existedger < 1)
 {
    opendir(DIR, "$diffgendir/Deseq") or die" canot open $diffgendir/Deseq $!\n";

    while (my $desfile = readdir(DIR))
   {

        next unless ($desfile =~ m/\.txt$/);
        next if ($desfile =~ m/_sig\.txt$/);

      my $pathwaysample=$desfile;
     $pathwaysample=~s/_sig\.txt//g;
     $pathwaysample=~s/\.txt//g;
      push @allpathwaysample, $pathwaysample;
   #start one compariosn at a time
     my (%allname,%commonname, %desname, %sigdesname, %desdata);

   #read over Deseq
     my $linecountsig=0;
      my @descolname;
      open ( COUNTIN, "$diffgendir/Deseq/$desfile ") or die "$! $diffgendir/Deseq/$desfile  \n";
            while (<COUNTIN>)
         {
                chomp;
                $linecountsig++;
               my $tmpline=$_;
               my @tmpline=split /\s+/, $_;
               foreach (@tmpline) {$_=~s/\"//g;}
               if ($linecountsig==1)
               {
                for (my $i=0; $i<=$#tmpline; $i++)
                {
                my $tmpcolname =join "_", ('Deseq', $tmpline[$i]);
                $descolname[$i]=$tmpcolname;
                }
               }
               else
            {
                 $allname{$tmpline[0]}=$tmpline[-2];
                 $desname{$tmpline[0]}=$tmpline[-2];
                for (my $i=1; $i<=$#tmpline; $i++ )
                {
                #print "$i\t$edgcolname[$i-1]\n";
              push  @{$alldesdata{$tmpline[0]}[$i-1]}, $tmpline[$i];
                $desdata{$tmpline[0]}[$i-1]=$tmpline[$i];
                }
            }
          }
       close COUNTIN;


      my $sumfile=join "_", ('Sum_allgenes', $desfile);
     open ( SUMOUT, ">$diffgendir/$sumfile") or die "$! $diffgendir/$sumfile  \n";

#        my $pathwayfile=$edgfile;
#               $pathwayfile=~s/\.txt//g;
#unless (-d "$diffgendir/significant_gene/pathway/$pathwayfile/") {mkdir "$diffgendir/significant_gene/pathway/$pathwayfile/", 0777 or die "$diffgendir/significant_gene/pathway/$pathwayfile $!";}


            my $heatmapfile=join "_", ('heatmap',$desfile);
    open ( RPKMOUT, ">$diffgendir/significant_gene/$heatmapfile") or die "$! $diffgendir/significant_gene/$heatmapfile  \n";
    print RPKMOUT  "ID";
    for (my $i=0; $i<=$#samplename; $i++) {print RPKMOUT "\t$samplename[$i]";   }
    print  RPKMOUT "\n";
    foreach (sort {$desname{$a} <=>$desname{$b} } keys %desname )
   {
    my $genename=$_;
     unless ($geneln{$genename} ) {next;}
     unless ($sigDeseqgene{$desfile}{$genename} ) {next;}
     $sigdesname{$genename}=$desname{$genename};
     $commonname{$genename}=$sigdesname{$genename};
     delete $desname{$genename};
   print RPKMOUT "$genename";
   for (my $k=0; $k<=$#samplename; $k++)
      {
        my $counts=0;
        if ($totalcounts_sample{$samplename[$k]} && $readcounts_gene{$genename}{$samplename[$k]} )
       {
       $counts=$allRPKM{$genename}{$samplename[$k]};#$readcounts_gene{$genename}{$samplename[$k]}*1e9/($totalcounts_sample{$samplename[$k]}*$geneln{$genename});
       }
      print RPKMOUT "\t$counts";
      }
   print RPKMOUT "\n";
   }
   close RPKMOUT;
    my $heatmapfig=join ".", ($heatmapfile, 'pdf');
   `$Bin/heatmap_distinctZ_noClust_zeroRowAllow.py --in $diffgendir/significant_gene/$heatmapfile -s log --out $diffgendir/figures/$heatmapfig`;

      my $RPKMscatfile=join "_", ('Scatter','RPKM',$desfile);
     my $tmpgroup=$desfile;
     $tmpgroup=~s/\.txt//g;
    $tmpgroup=~s/Deseq\.//g;
    my @tmpgroup=split /_/, $tmpgroup;
    print "@tmpgroup\n";
    open ( RPKMSCAT, ">$diffgendir/$RPKMscatfile") or die "$! $diffgendir/$RPKMscatfile  \n";
    print RPKMSCAT "gene\t$tmpgroup[0]\t$tmpgroup[1]\tsignificant\n";
    for (sort keys %averageRPKM)
    {
    my $gene=$_;
    unless ($averageRPKM{$gene}{$tmpgroup[0]} && $groupNUM{$tmpgroup[0]} && $groupNUM{$tmpgroup[1]} && $averageRPKM{$gene}{$tmpgroup[1]}) {next;}
    my $avegroup0=$averageRPKM{$gene}{$tmpgroup[0]}/$groupNUM{$tmpgroup[0]};
    my $avegroup1=$averageRPKM{$gene}{$tmpgroup[1]}/$groupNUM{$tmpgroup[1]};
    print RPKMSCAT "$gene\t$avegroup0\t$avegroup1\t";
    if ($commonname{$gene}) {print RPKMSCAT "yes\n";} else {print RPKMSCAT "no\n";}
    }
    close RPKMSCAT;

     # chdir "$diffgendir";
    print LOGFILE "Rscript $Bin/scatt_plot.R  $RPKMscatfile $diffgendir\n";
    `Rscript $Bin/scatt_plot.R   $RPKMscatfile $diffgendir`;


                 $heatmapfile=join "_", ('sum',$desfile);
    open ( RPKMOUT, ">$diffgendir/significant_gene/$heatmapfile") or die "$! $diffgendir/significant_gene/$heatmapfile  \n";

         print RPKMOUT  "ID";
    for (my $i=0; $i<=$#samplename; $i++) {print RPKMOUT "\t$samplename[$i]";   }
    print  RPKMOUT "\n";

         print SUMOUT  "Gene\tDescription";
  foreach (@descolname )
      {
      print SUMOUT "\t$_";
      }
    print  SUMOUT "\n";


 foreach (sort {$commonname{$a} <=>$commonname{$b} } keys %commonname )
   {
    my $genename=$_;
    unless ($geneln{$genename}) {next;}
   print RPKMOUT "$genename";
   for (my $k=0; $k<=$#samplename; $k++)
      {
       my $counts=0;
       if ($totalcounts_sample{$samplename[$k]} && $readcounts_gene{$genename}{$samplename[$k]} )
      {
      $counts=$allRPKM{$genename}{$samplename[$k]};#$readcounts_gene{$genename}{$samplename[$k]}*1e9/($totalcounts_sample{$samplename[$k]}*$geneln{$genename});
      }
      print RPKMOUT "\t$counts";
      }
   print RPKMOUT "\n";

 print SUMOUT "Deseq::$pathwaysample";
   print SUMOUT "::$genename\t$genedes{$genename}";

      foreach (@{$desdata{$genename}} )
      {
      print SUMOUT "\t$_";
      }
     print SUMOUT "\n";

   }

        foreach (sort {$desname{$a} <=>$desname{$b} } keys %desname )
   {
    my $genename=$_;
   unless ($geneln{$genename} ) {next;}
  print SUMOUT "NotSignificant::$pathwaysample";
  print SUMOUT "::$genename\t$genedes{$genename}";

     foreach (@{$desdata{$genename}} )
      {
      print SUMOUT "\t$_";
      }
     print SUMOUT "\n";
   }

  close SUMOUT;
   close RPKMOUT;

         my $pathwayfile=$desfile;
               $pathwayfile=~s/\.txt//g;
unless (-d "$diffgendir/significant_gene/pathway/$pathwayfile/") {mkdir "$diffgendir/significant_gene/pathway/$pathwayfile/", 0777 or die "$diffgendir/significant_gene/pathway/$pathwayfile $!";}
 open ( PATHOUT, ">$diffgendir/significant_gene/pathway/$pathwayfile/pathway.stats.txt") or die "$! $diffgendir/significant_gene/pathway/$pathwayfile/pathway.stats.txt  \n";

  #print PATHOUT  "GeneID\t$edgcolname[2]\t$edgcolname[3]\tDescription\t$descolname[4]\tMap_ID\t$descolname[1]\tLOCUS_ID\t$edgcolname[0]\n";
  print PATHOUT  "GeneID\tDescription\tMap_ID\tLOCUS_ID\t$descolname[1]\n";
   foreach (sort {$commonname{$a} <=>$commonname{$b} } keys %commonname )
  {

     my $genename=$_;
     unless ($pathway{$genename}) {next;}

    @{$allcommonname{$genename}{Deslogfc}}= @{$alldesdata{$genename}[1]};
  #  @{$allcommonname{$genename}{Edgelogfc}}=@{$alledgdata{$genename}[0]};
  #  @{$allcommonname{$genename}{Despadj}}=@{$alldesdata{$genename}[4]};
   # @{$allcommonname{$genename}{EdgeFDR}}=@{$alledgdata{$genename}[2]};

   foreach (@{$pathway{$genename}} )
  {
   my $mapid=$_;

   print PATHOUT "$pathwaysample";
   print PATHOUT "::$genename\t";
   #print PATHOUT "$edgdata{$genename}[2]\t$edgdata{$genename}[3]\t$genedes{$genename}\t$desdata{$genename}[4]\t$mapid\t$desdata{$genename}[1]\t$genename\t$edgdata{$genename}[0]\n";
   print PATHOUT "$genedes{$genename}\t$mapid\t$genename\t$desdata{$genename}[1]\n";
   }
  }
   close PATHOUT;

   chdir "$diffgendir/significant_gene/pathway/$pathwayfile/";

 `R --vanilla --slave --silent  < $Bin/pathview.R  2>/dev/null`;
# `perl ./json.pl`;

  }#while (my $edgfile = readdir(DIR))


unless (-d "$diffgendir/significant_gene/pathway/all_comparison") {mkdir "$diffgendir/significant_gene/pathway/all_comparison", 0777 or die "$diffgendir/significant_gene/pathway/all_compari
son $!";}

chdir "$diffgendir/significant_gene/pathway/all_comparison";

open ( PATHOUT, ">$diffgendir/significant_gene/pathway/all_comparison/pathway.stats.txt") or die "$! $diffgendir/significant_gene/pathway/all_comparison//pathway.stats.txt  \n";

my $allpathwaysample=join ":", @allpathwaysample;

 print PATHOUT  "GeneID";

# foreach (@allpathwaysample)
# {
# my $tmpname=join "_", ($_, 'EdgeFDR');
# print PATHOUT "\t$tmpname";
# }
 print PATHOUT" \tDescription";
#  foreach (@allpathwaysample)
# {
# my $tmpname=join "_", ($_, 'Despadj');
# print PATHOUT"\t$tmpname";
# }
 print PATHOUT"\tMap_ID";
#   foreach (@allpathwaysample)
# {
# my $tmpname=join "_", ($_, 'Deslogfc');
# print PATHOUT"\t$tmpname";
# }
 print PATHOUT" \tLOCUS_ID";
    foreach (@allpathwaysample)
 {
 my $tmpname=join "_", ($_, 'Deslogfc');
 print PATHOUT"\t$tmpname";
 }
 print PATHOUT"\n";

   foreach (sort  keys %allcommonname )
{

my $genename=$_;
 unless ($pathway{$genename}) {next;}
      @{$allcommonname{$genename}{Deslogfc}}= @{$alldesdata{$genename}[1]};
#    @{$allcommonname{$genename}{Edgelogfc}}=@{$alledgdata{$genename}[0]};
#    @{$allcommonname{$genename}{Despadj}}=@{$alldesdata{$genename}[4]};
#    @{$allcommonname{$genename}{EdgeFDR}}=@{$alledgdata{$genename}[2]};

   foreach (@{$pathway{$genename}} )
{
my $mapid=$_;
my $tmpname=join"::", ($allpathwaysample,$genename);
print PATHOUT "$tmpname";
#foreach ( @{$allcommonname{$genename}{EdgeFDR}})
# {
#  print PATHOUT"\t$_";
# }
 print PATHOUT"\t$genedes{$genename}";
#foreach ( @{$allcommonname{$genename}{Despadj}})
# {
#  print PATHOUT"\t$_";
# }
 print PATHOUT"\t$mapid";
# foreach (@{$allcommonname{$genename}{Deslogfc}})
# {
#  print PATHOUT"\t$_";
# }
print PATHOUT"\t$genename";
  foreach (@{$allcommonname{$genename}{Deslogfc}})
 {
  print PATHOUT"\t$_";
 }
 print PATHOUT"\n";
}
}
close PATHOUT;
`R --vanilla --slave --silent  < $Bin/pathview.R  2>/dev/null`;
#`perl ./json.pl`;

}#  if ($existdeseq >=1 && $existedger < 1)
                                              

sub file_check 
{
    #check file exist and non zero size
    my $file=shift;
    my $exist=-1;
    if (-e $file) {$exist=1};
    if (-z $file) {$exist=-1};
    return $exist;
}







