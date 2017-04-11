#!/usr/bin/env perl  
use strict;
use warnings;
use FindBin qw($Bin);
use lib "$Bin/../lib";
use lib "$Bin/../ext/lib/perl5";
use Test::More;

plan tests => 3;

BEGIN {
    use_ok( 'PiReT::Map' ) || print "Bail out!\n";
    use_ok( 'PiReT::Count' ) || print "Bail out!\n";
    use_ok( 'PiReT::Variant' ) || print "Bail out!\n";
}

# diag( "Testing PiReT::Map PiReT::Map::VERSION, Perl $], $^X" );
