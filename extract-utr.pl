#!/usr/bin/env perl
# Mike Covington
# created: 2013-12-19
#
# Description:
#
use strict;
use warnings;
use autodie;
use feature 'say';

my $gff_file = glob "~/git.repos/sample-files/annotation/ITAG2.3_gene_models.gff3";
my $cds_fa_file    = "~/git.repos/sample-files/fa/ITAG2.3_cds.fasta";
my $genome_fa_file = "~/git.repos/sample-files/fa/ITAG2.3_genomic.fasta";
my $output_fa_file = "out.fa";

my $utr_length  = 500;
my $gene_length = 500;

my $fa_width = 80;

my $coding_regions = extract_cds_from_gff($gff_file);

open my $output_fa_fh, ">", $output_fa_file;
for my $id ( keys %$coding_regions ) {
    my $chr    = $$coding_regions{$id}{chr};
    my $strand = $$coding_regions{$id}{strand};

    my ( $utr3_start, $utr3_end ) =
      find_utr_boundaries( $$coding_regions{$id}, $utr_length );

    my $utr_seq =
      extract_fa_seq( $genome_fa_file, $chr, $strand, $utr3_start, $utr3_end );

    my $gene_seq = extract_fa_seq( $cds_fa_file, $id );
    $gene_seq = substr $gene_seq, -$gene_length;

    my $combo_seq = "$gene_seq$utr_seq";
    output_fa( $id, $combo_seq, $output_fa_fh, $fa_width );
}
close $output_fa_fh;

exit;

sub extract_cds_from_gff {
    my $gff_file = shift;

    open my $gff_fh, "<", $gff_file;

    my $gff_header = <DATA>;
    check_gff_version( 3, $gff_header );

    my %coding_regions;
    build_coding_regions_hash( \%coding_regions, $gff_fh );
    close $gff_fh;

    extract_start_end_of_full_cds( \%coding_regions );

    return \%coding_regions;
}

sub check_gff_version {
    my ( $required_version, $gff_version ) = @_;
    die "Requires gff file to be version $required_version\n"
      unless $gff_version =~ /$required_version/;
}

sub build_coding_regions_hash {
    my ( $coding_regions, $gff_fh ) = @_;

    while ( my $feature = <$gff_fh> ) {
        next if $feature =~ /^#/;
        my ( $chr, $type, $start, $end, $strand, $attributes ) =
          ( split /\t/, $feature )[ 0, 2 .. 4, 6, 8 ];
        next unless $type eq "CDS";

        my ($gene) = $attributes =~ /Parent=(?:mRNA:)?([^;]+)/;
        $$coding_regions{$gene}{chr}    = $chr;
        $$coding_regions{$gene}{strand} = $strand;
        push @{ $$coding_regions{$gene}{pos} }, [ $start, $end ];
    }
}

sub extract_start_end_of_full_cds {
    my $coding_regions = shift;

    for my $id ( keys %$coding_regions ) {
        $$coding_regions{$id}{left_pos}  = $$coding_regions{$id}{pos}[0][0];
        $$coding_regions{$id}{right_pos} = $$coding_regions{$id}{pos}[-1][1];
    }
}

sub find_utr_boundaries {
    my ( $cds_info, $utr_length ) = @_;

    my $chr       = $$cds_info{chr};
    my $strand    = $$cds_info{strand};
    my $left_pos  = $$cds_info{left_pos};
    my $right_pos = $$cds_info{right_pos};

    my $utr3_start;
    my $utr3_end;

    if ( $strand eq '+' ) {
        $utr3_start = $right_pos + 1;
        $utr3_end   = $right_pos + $utr_length;
    }
    elsif ( $strand eq '-' ) {
        $utr3_start = $left_pos - $utr_length;
        $utr3_end   = $left_pos - 1;
    }
    else { die "Problem with strand info\n" }

    return ( $utr3_start, $utr3_end );
}

sub extract_fa_seq {
    my ( $fa_file, $seqid, $strand, $left_pos, $right_pos ) = @_;

    my $samtools_path = "~/installs/bin/samtools";
    my $faidx_cmd =
      defined $left_pos && defined $right_pos
      ? "$samtools_path faidx $fa_file $seqid:$left_pos-$right_pos"
      : "$samtools_path faidx $fa_file $seqid";

    my ( $fa_header, @fa_seq ) = `$faidx_cmd`;
    chomp @fa_seq;

    my $seq = join "", @fa_seq;
    $seq = reverse $seq
      if defined $strand
      && $strand eq '-';

    return $seq;
}

sub output_fa {
    my ( $seqid, $seq, $output_fa_fh, $fa_width ) = @_;

    $fa_width //= 80;
    my @fa_seq = unpack "(A$fa_width)*", $seq;

    say $output_fa_fh ">$seqid";
    say $output_fa_fh $_ for @fa_seq;
}
