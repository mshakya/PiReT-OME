use strict;
use warnings;
use Module::Build;

my $builder = Module::Build->new(
    module_name         => 'Map::Mapping',
    license             => 'perl',
    dist_abstract       => 'A collection of perl subrotuine for RNA seq analysis',
    dist_author         => 'Author Name <migun@lanl.gov>',
    build_requires => {
        'Test::More' => '0.10',
    },
);

$builder->create_build_script();