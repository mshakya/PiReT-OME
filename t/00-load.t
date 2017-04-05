#!perl -T
use 5.006;
use strict;
use warnings;
use Test::More;

plan tests => 1;

BEGIN {
    use_ok( 'PiReT::Map' ) || print "Bail out!\n";
    # use_ok( 'PireT::Count' ) || print "Bail out!\n";
}

diag( "Testing PiReT::Map $PiReT::Map::VERSION, Perl $], $^X" );
