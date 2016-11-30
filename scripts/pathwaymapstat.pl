#!/usr/bin/perl

use strict;
use warnings;

my $filename='/users/203270/code/bin/kegg_locus_orgnism.txt';
my %path;
open (ORGIN , "$filename") or die "$! $filename  \n";
   while (<ORGIN>) {
 chomp;
        my ($organism, $ID, $path)=split /\s+/, $_;
        my $tmpname=join ':', ($organism, $path);
        $path{$tmpname}{$ID}=1; 
  }
close ORGIN;
$filename='pathway.stats.txt';
my %expresspath;
open (ORGIN , "$filename") or die "$! $filename  \n";
   while (<ORGIN>) {
 chomp;
        my @line=split /\t/, $_;
        my @ID=split /::/, $line[0];
       if ($path{$line[5]}{$ID[-1]}) {
       push @{$expresspath{$line[5]}}, $ID[-1];
  #      print "good:$line[5]\t$ID[-1]\n";
   #    } else {
    #   print "bad:$line[5]\t$ID[-1]\n";
     }
  }
close ORGIN;

foreach (sort keys %expresspath)
{ 
   my $key=$_;
   my $expressid=join ',',  @{$expresspath{$key}};
  my $numid=keys(%{$path{$key}});
  my $numexpress=$#{$expresspath{$key}}+1;
   print "$key\t$numid\t$expressid\t$numexpress\n";
}
