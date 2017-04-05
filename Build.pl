use strict;
use warnings;
use Module::Build;

my $v=1;

my $builder = Module::Build->new(
    dist_name      => 'PiReT',
    module_name    => 'PiReT',
    license        => 'perl',
    dist_abstract  => 'RNA seq',
    dist_author    => 'Migun Shakya <migun@lanl.gov>',
    dist_version        => $v,
    build_requires => { 'Test::More' => '0.10', },
);

$builder->create_build_script();
