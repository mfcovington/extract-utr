#!/usr/bin/env perl
# Mike Covington
# created: 2016-09-19
#
# Description:
#
use strict;
use warnings;
use autodie;
use feature 'say';
use List::Util 'sum';
use Number::RangeTracker;

use FindBin;
use lib "$FindBin::Bin/../lib";
use gff;

# TEMPORARY Defaults
# my $gff_file = "Zea_mays.AGPv3.21.gff3";
my $gff_file = glob "~/git.repos/extract-utr/GRMZM2G176820.gff";
my $alignment_file = "~/git.repos/extract-utr/accepted_hits.bam";
my $cds_length = 500;
my $utr_length = 500;

my $coding_regions = extract_cds_from_gff($gff_file);

use Data::Printer;
p $coding_regions;

convert_coding_regions_to_three_prime( $coding_regions, $cds_length,
    $utr_length );
my $chromosomal_ranges = merge_chromosomal_ranges($coding_regions);

p $coding_regions;

my %genes_per_chromosome;
my %gene_counts;
for my $gene ( sort keys %$coding_regions ) {
    push @{$genes_per_chromosome{$$coding_regions{$gene}{'chr'}}}, $gene;
    $gene_counts{$gene} = 0;
}

my $alignment_fh;
if ( $alignment_file =~ /.+\.sam$/i ) {
    open $alignment_fh, "<", $alignment_file;
}
elsif ( $alignment_file =~ /.+\.bam$/i ) {
    open $alignment_fh, "-|", "samtools view -h $alignment_file";
}
else {    # Should never happen because validate_options()
    die "File '$alignment_file' is not a .sam/.bam file\n";
}


while (<$alignment_fh>) {
    next if /^@/;

    my ( $seq_id, $aln_left, $cigar ) = (split)[ 2, 3, 5 ];
    next if $seq_id =~ /^\*$/;    # skip unmapped reads

    my $aln_right = $aln_left + get_alignment_length_from_cigar($cigar) - 1;

    next unless exists $$chromosomal_ranges{$seq_id};
    next
        unless $$chromosomal_ranges{$seq_id}->is_in_range($aln_left)
        or $$chromosomal_ranges{$seq_id}->is_in_range($aln_right);

    for my $gene ( @{$genes_per_chromosome{$seq_id}} ) {
        next
            unless $$coding_regions{$gene}{'range'}->is_in_range($aln_left)
            or $$coding_regions{$gene}{'range'}->is_in_range($aln_right);
        $gene_counts{$gene}++;
    }
}

close $alignment_fh;
exit;

sub get_alignment_length_from_cigar {
    my $cigar = shift;

    my $matches    = sum( $cigar =~ /(\d+)M/g ) // 0;
    my $insertions = sum( $cigar =~ /(\d+)I/g ) // 0;
    my $deletions  = sum( $cigar =~ /(\d+)D/g ) // 0;

    my @unsupported = $cigar =~ /(\d+[^DIM\d])/g;
    warn "Currently unsupported CIGAR scores detected: @unsupported\n"
        if @unsupported;

    return $matches + $deletions - $insertions;
}

sub merge_chromosomal_ranges {
    my $coding_regions = shift;
    my %chromosomal_ranges;

    for my $gene ( keys %$coding_regions ) {
        my $chr = $$coding_regions{$gene}{'chr'};

        $chromosomal_ranges{$chr} = Number::RangeTracker->new
            unless exists $chromosomal_ranges{$chr};
        $chromosomal_ranges{$chr}
            ->add( $$coding_regions{$gene}{'range'}->output );
    }

    for my $chr ( keys %chromosomal_ranges ) {
        $chromosomal_ranges{$chr}->collapse;
    }

    return \%chromosomal_ranges;
}
