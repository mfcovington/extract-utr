use strict;
use warnings;
use autodie;
use List::Util qw(min max);
use Number::RangeTracker;

sub build_coding_regions_hash {
    my ( $coding_regions, $gff_fh ) = @_;

    while ( my $feature = <$gff_fh> ) {
        chomp $feature;
        next if $feature =~ /^#/;
        my ( $chr, $type, $start, $end, $strand, $attributes ) =
          ( split /\t/, $feature )[ 0, 2 .. 4, 6, 8 ];
        next unless $type eq "CDS";

        my ($gene) = $attributes =~ /Parent=(?:mRNA:|transcript:)?([^;]+)/;
        $$coding_regions{$gene}{chr}    = $chr;
        $$coding_regions{$gene}{strand} = $strand;
        push @{ $$coding_regions{$gene}{pos} }, [ $start, $end ];
    }
}

sub check_gff_version {
    my ( $required_version, $gff_version ) = @_;
    die "Requires gff file to be version $required_version\n"
      unless $gff_version =~ /$required_version/;
}

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
                max( 0, $$coding_regions{$gene}{'left_pos'} - 500 ),
                max( 0, $$coding_regions{$gene}{'left_pos'} - 1 )
            ];
            unshift @new_pos, $utr;
            $$coding_regions{$gene}{'left_pos'} = $$utr[0];
        }

        $$coding_regions{$gene}{'pos'} = \@new_pos;

        $$coding_regions{$gene}{'range'} = Number::RangeTracker->new;
        $$coding_regions{$gene}{'range'}->add(@new_pos);
    }
}

sub extract_cds_from_gff {
    my $gff_file = shift;

    open my $gff_fh, "<", $gff_file;

    my $gff_header = <$gff_fh>;
    check_gff_version( 3, $gff_header );

    my %coding_regions;
    build_coding_regions_hash( \%coding_regions, $gff_fh );
    close $gff_fh;

    extract_start_end_of_full_cds( \%coding_regions );

    return \%coding_regions;
}

sub extract_start_end_of_full_cds {
    my $coding_regions = shift;

    for my $id ( keys %$coding_regions ) {
        $$coding_regions{$id}{left_pos}  = $$coding_regions{$id}{pos}[0][0];
        $$coding_regions{$id}{right_pos} = $$coding_regions{$id}{pos}[-1][1];
    }
}

sub find_utr_boundaries {
    my ( $cds_info, $utr_length, $fiveprime ) = @_;

    my $chr       = $$cds_info{chr};
    my $strand    = $$cds_info{strand};
    my $left_pos  = $$cds_info{left_pos};
    my $right_pos = $$cds_info{right_pos};

    my $utr_start;
    my $utr_end;

    # originally written for 3'
    # use opposite logic for 5'
    $strand =~ tr/+-/-+/ if $fiveprime;

    if ( $strand eq '+' ) {
        $utr_start = $right_pos + 1;
        $utr_end   = $right_pos + $utr_length;
    }
    elsif ( $strand eq '-' ) {
        $utr_start = $left_pos - $utr_length;
        $utr_end   = $left_pos - 1;
    }
    else { die "Problem with strand info\n" }

    return ( $utr_start, $utr_end );
}

1;
