#!/usr/bin/env perl
# Mike Covington
# created: 2016-09-21
#
# Description:
#
use strict;
use warnings;
use autodie;
use feature 'say';
use Getopt::Long;
use List::Util qw(min max);

use FindBin;
use lib "$FindBin::Bin/../lib";
use gff;

# Default Options
my ( $gff_in_file, $bed_out_file );
my $cds_length = 500;
my $utr_length = 500;
my $transcript_delimiter = "_T";

my $options = GetOptions(
    "gff_in_file=s"          => \$gff_in_file,
    "bed_out_file=s"         => \$bed_out_file,
    "cds_length=i"           => \$cds_length,
    "utr_length=i"           => \$utr_length,
    "transcript_delimiter=s" => \$transcript_delimiter,
);

my $coding_regions = extract_cds_from_gff($gff_in_file);
convert_coding_regions_to_three_prime( $coding_regions, $cds_length,
    $utr_length );

my %genes;

for my $transcript_id ( keys %$coding_regions ) {
    my ($gene_id) = split /$transcript_delimiter/, $transcript_id;

    if ( exists $genes{$gene_id} ) {
        $genes{$gene_id}{'left_pos'} = min( $genes{$gene_id}{'left_pos'},
            $$coding_regions{$transcript_id}{'left_pos'} );
        $genes{$gene_id}{'right_pos'} = max( $genes{$gene_id}{'right_pos'},
            $$coding_regions{$transcript_id}{'right_pos'} );
    }
    else {
        $genes{$gene_id}{'left_pos'}
            = $$coding_regions{$transcript_id}{'left_pos'};
        $genes{$gene_id}{'right_pos'}
            = $$coding_regions{$transcript_id}{'right_pos'};
    }

    $genes{$gene_id}{'chr'}    = $$coding_regions{$transcript_id}{'chr'};
    $genes{$gene_id}{'strand'} = $$coding_regions{$transcript_id}{'strand'};
}

open my $bed_out_fh, ">", $bed_out_file;
for my $gene_id ( sort keys %genes ) {
    say $bed_out_fh join "\t",
        $genes{$gene_id}{'chr'},
        $genes{$gene_id}{'left_pos'},
        $genes{$gene_id}{'right_pos'},
        $gene_id,
        "",
        $genes{$gene_id}{'strand'};
}
close $bed_out_fh;
