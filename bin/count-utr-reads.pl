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
}

close $alignment_fh;
exit;

sub convert_coding_regions_to_three_prime {
    # \ {
    #     GRMZM2G176820_T01   {
    #         chr         1,
    #         left_pos    4752233,
    #         pos         [
    #             [0] [
    #                 [0] 4752233,
    #                 [1] 4752441
    #             ],
    #             [1] [
    #                 [0] 4753784,
    #                 [1] 4753847
    #             ],
    #             [2] [
    #                 [0] 4753949,
    #                 [1] 4754121
    #             ],
    #             [3] [
    #                 [0] 4754215,
    #                 [1] 4754422
    #             ]
    #         ],
    #         right_pos   4754422,
    #         strand      "-"
    #     }
    # }

    # WANT 500 CDS + 500 3'-UTR:
    # \ {
    #     GRMZM2G176820_T01   {
    #         chr         1,
    #         left_pos    4751733,
    #         pos         [
    #             [0] [
    #                 [0] 4751733,
    #                 [1] 4752232
    #             ],
    #             [1] [
    #                 [0] 4752233,
    #                 [1] 4752441
    #             ],
    #             [2] [
    #                 [0] 4753784,
    #                 [1] 4753847
    #             ],
    #             [3] [
    #                 [0] 4753949,
    #                 [1] 4754121
    #             ],
    #             [4] [
    #                 [0] 4754215,
    #                 [1] 4754268
    #             ]
    #         ],
    #         right_pos   4754268,
    #         strand      "-"
    #     }
    # }

    # IF '+' STRAND, EXPECT THIS:
    # \ {
    #     GRMZM2G176820_T01   {
    #         chr         1,
    #         left_pos    4752387,
    #         pos         [
    #             [0] [
    #                 [0] 4752387,
    #                 [1] 4752441
    #             ],
    #             [1] [
    #                 [0] 4753784,
    #                 [1] 4753847
    #             ],
    #             [2] [
    #                 [0] 4753949,
    #                 [1] 4754121
    #             ],
    #             [3] [
    #                 [0] 4754215,
    #                 [1] 4754422
    #             ],
    #             [4] [
    #                 [0] 4754423,
    #                 [1] 4754922
    #             ]
    #         ],
    #         right_pos   4754922,
    #         strand      "-"
    #     }
    # }

    my ( $coding_regions, $cds_length, $utr_length ) = @_;

    for my $gene ( keys %$coding_regions ) {
        my $strand = $$coding_regions{$gene}{'strand'};

        my @pos = @{ $$coding_regions{$gene}{'pos'} };
        @pos = reverse @pos if $strand eq '+';

        my $total_cds_length = 0;
        my @new_pos;

        for my $interval (@pos) {
            my $interval_length = @$interval[1] - @$interval[0] + 1;
            $total_cds_length += $interval_length;

            if ( $total_cds_length < $cds_length ) {
                push @new_pos, $interval;
            }
            else {
                my $excess_length = $total_cds_length - $cds_length;

                if ( $strand eq '+' ) {
                    my $new_left_pos = @$interval[0] + $excess_length;
                    push @new_pos, [ $new_left_pos, @$interval[1] ];
                    $$coding_regions{$gene}{'left_pos'} = $new_left_pos;
                }
                else {
                    my $new_right_pos = @$interval[1] - $excess_length;
                    push @new_pos, [ @$interval[0], $new_right_pos ];
                    $$coding_regions{$gene}{'right_pos'} = $new_right_pos;
                }
                last;
            }
        }

        my $utr;
        if ( $strand eq '+' ) {
            $utr = [
                $$coding_regions{$gene}{'right_pos'} + 1,
                $$coding_regions{$gene}{'right_pos'} + 500
            ];
            unshift @new_pos, $utr;
            @new_pos = reverse @new_pos;
            $$coding_regions{$gene}{'right_pos'} += 500;
        }
        else {
            $utr = [
                $$coding_regions{$gene}{'left_pos'} - 500,
                $$coding_regions{$gene}{'left_pos'} - 1
            ];
            unshift @new_pos, $utr;
            $$coding_regions{$gene}{'left_pos'} -= 500;
        }

        $$coding_regions{$gene}{'pos'} = \@new_pos;

        $$coding_regions{$gene}{'range'} = Number::RangeTracker->new;
        $$coding_regions{$gene}{'range'}->add(@new_pos);
    }
}

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
