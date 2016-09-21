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
use List::Util qw(min max);

use FindBin;
use lib "$FindBin::Bin/../lib";
use gff;

# TEMPORARY Defaults
my $gff_file = glob "~/git.repos/extract-utr/Zea_mays.AGPv3.21.gff3";
# my $gff_file = glob "~/git.repos/extract-utr/GRMZM2G176820.gff";
my $bed_file = glob "~/git.repos/extract-utr/out2.bed";
my $cds_length = 500;
my $utr_length = 500;
my $transcript_delimiter = "_";

my $coding_regions = extract_cds_from_gff($gff_file);
convert_coding_regions_to_three_prime( $coding_regions, $cds_length,
    $utr_length );

my %genes;

for my $transcript_id ( keys %$coding_regions ) {
    my ($gene_id) = split /_/, $transcript_id;

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

open my $bed_fh, ">", $bed_file;
for my $gene_id ( sort keys %genes ) {
    say $bed_fh join "\t",
        $genes{$gene_id}{'chr'},
        $genes{$gene_id}{'left_pos'},
        $genes{$gene_id}{'right_pos'},
        $gene_id,
        "",
        $genes{$gene_id}{'strand'};
}
close $bed_fh;
