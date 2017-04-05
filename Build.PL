use strict;
use warnings;
use Module::Build;
use lib 'ext/lib/perl5';

my $v=1;

my $builder = Module::Build->new(
    dist_name      => 'PiReT',
    module_name    => 'PiReT',
    license             => 'perl',
    dist_abstract       => 'A collection of perl subrotuine for RNA seq analysis',
    dist_author         => 'Migun Shakya  <migun@lanl.gov>',
    dist_version        => $v,
    build_requires => {
        'Test::More' => '0.10',
        'Parallel::ForkManager' => '1.17',
        'String::Approx' => '3.27'
    },
);

$builder->create_build_script();