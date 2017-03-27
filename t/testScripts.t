#!/usr/bin/env perl  

use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use Test::More;
use Test::Script;

script_compiles('scripts/reads_mapping.pl');
script_runs('scripts/reads_mapping.pl', '-P1 t/data/1.trimmed.fastq -P2 t/data/2.trimmed.fastq -I t/data/euk_prok_index -K both -E t/data/euk_test.fna -B t/data/prok_test.fna -W t/results -S samp1', 'test running reads_mapping.pl');
done_testing();