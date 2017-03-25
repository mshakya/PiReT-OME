#!/usr/bin/perl

use strict;
use warnings;
use File::Fetch;

use FindBin qw($Bin);

$ENV{PATH} = "$Bin:$ENV{PATH}";



#my $url = 'http://rest.kegg.jp/list/organism';
#my $ff = File::Fetch->new(uri => $url);
#my $file = $ff->fetch() or die $ff->error;

unless (-d "$Bin/pathway") {mkdir "$Bin/pathway", 0777 or die "can not make dir  $Bin/pathway $!";}

unless (-d "$Bin/pathway/org_files") {mkdir "$Bin/pathway/org_files", 0777 or die "can not make dir  $Bin/pathway/org_files $!";}

#chdir '/users/203270/scratch/pathway/org_files';

my @organism;

open (ORGIN , "$Bin/pathway/organism") or die "$! $Bin/pathway/organism  \n";
while (<ORGIN>) {

 chomp;
        my @line=split /\t/, $_;
       $line[1]=~s/\s+//g;
       push @organism, $line[1];
#     my  $url_f = join "", ('http://rest.kegg.jp/link/pathway/', $line[1]);  
#     my $ff_f = File::Fetch->new(uri => $url_f);
#     my $file_f = $ff_f->fetch() or die $ff_f->error;
 
#     print "$url_f\n";

}
close ORGIN;

open (ORGOUT , ">$Bin/pathway/kegg_locus_orgnism.txt") or die "$! $Bin/pathway/kegg_locus_orgnism.txt  \n";
foreach (@organism) 

{
 my $filename=$_;
 my $filesize = -s $filename;
 if ($filesize <= 10) {next;}
open (ORGIN , "$filename") or die "$! $filename  \n";
   while (<ORGIN>) {
 chomp;
              
        my @line=split /\s+/, $_;
       my ($organism, $ID)=split /:/, $line[0];
       my ($dummy, $path)=split /:/, $line[1];
       $path=~s/$organism//;
       print ORGOUT "$organism\t$ID\t$path\n";
  }
close ORGIN;

}

close ORGOUT;
