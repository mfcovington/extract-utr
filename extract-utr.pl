#!/usr/bin/env perl
# Mike Covington
# created: 2013-12-19
#
# Description: Extract genic and extragenic sequence surrounding stop codon
#
use strict;
use warnings;
use autodie;
use feature 'say';
use Getopt::Long;

#TODO: Add README

# Defaults
my $gff_file = glob "~/git.repos/sample-files/annotation/ITAG2.3_gene_models.gff3";
my $cds_fa_file    = "~/git.repos/sample-files/fa/ITAG2.3_cds.fasta";
my $genome_fa_file = "~/git.repos/sample-files/fa/ITAG2.3_genomic.fasta";
my $output_fa_file = "out.fa";

my $utr_length  = 500;
my $gene_length = 500;

my $fa_width = 80;

my ( $fiveprime, $threeprime, $both );

my $options = GetOptions(
    "gff_file=s"       => \$gff_file,
    "cds_fa_file=s"    => \$cds_fa_file,
    "genome_fa_file=s" => \$genome_fa_file,
    "output_fa_file=s" => \$output_fa_file,
    "utr_length=i"     => \$utr_length,
    "gene_length=i"    => \$gene_length,
    "fa_width=i"       => \$fa_width,
    "fiveprime"        => \$fiveprime,
    "threeprime"       => \$threeprime,
    "both"             => \$both,
);

check_options( $fiveprime, $threeprime, $both );

$gff_file = glob $gff_file;

my $coding_regions = extract_cds_from_gff($gff_file);

open my $output_fa_fh, ">", $output_fa_file;
for my $id ( sort keys %$coding_regions ) {
    my $chr    = $$coding_regions{$id}{chr};
    my $strand = $$coding_regions{$id}{strand};

    $gene_length = -1 if $both;

    my $utr_seq = [];
    if ( $both ) {
        $$utr_seq[0] = get_utr_seq( $id, $coding_regions, $utr_length, 1,
            $genome_fa_file, $chr, $strand );
        $$utr_seq[1] = get_utr_seq( $id, $coding_regions, $utr_length, 0,
            $genome_fa_file, $chr, $strand );
    }
    else {
        $$utr_seq[0]
            = get_utr_seq( $id, $coding_regions, $utr_length, $fiveprime,
            $genome_fa_file, $chr, $strand );
    }

    my $combo_seq;
    if ( $gene_length == 0 ) {
        $combo_seq = $$utr_seq[0];
    }
    else {
        my $gene_seq = extract_fa_seq( $cds_fa_file, $id );
        $gene_seq
            = trim_seq( $gene_seq, $gene_length, $fiveprime, $threeprime,
            $both );

        $combo_seq = combine_seqs( $fiveprime, $threeprime, $both, $utr_seq,
            $gene_seq );
    }

    output_fa( $id, $combo_seq, $output_fa_fh, $fa_width );
}
close $output_fa_fh;

exit;

sub check_options {
    my ( $fiveprime, $threeprime, $both ) = @_;

    my $def_count = 0;
    map { $def_count++ if defined $_ } $fiveprime, $threeprime, $both;
    die "Specify '--fiveprime', '--threeprime', OR '--both'\n"
        unless $def_count == 1;
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

        my ($gene) = $attributes =~ /Parent=(?:mRNA:|transcript:)?([^;]+)/;
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

sub get_utr_seq {
    my ( $id, $coding_regions, $utr_length, $fiveprime, $genome_fa_file,
        $chr, $strand )
        = @_;
    my ( $utr_start, $utr_end ) =
      find_utr_boundaries( $$coding_regions{$id}, $utr_length, $fiveprime );
    return extract_fa_seq( $genome_fa_file, $chr, $strand, $utr_start,
        $utr_end );
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

sub extract_fa_seq {
    my ( $fa_file, $seqid, $strand, $left_pos, $right_pos ) = @_;

    my $faidx_cmd =
      defined $left_pos && defined $right_pos
      ? "samtools faidx $fa_file $seqid:$left_pos-$right_pos"
      : "samtools faidx $fa_file $seqid";

    my ( $fa_header, @fa_seq ) = `$faidx_cmd`;
    chomp @fa_seq;

    my $seq = join "", @fa_seq;

    if ( defined $strand && $strand eq '-' ) {
        $seq = reverse $seq;
        $seq =~ tr/ACGTacgt/TGCAtgca/;
    }

    return $seq;
}

sub trim_seq {
    my ( $gene_seq, $gene_length, $fiveprime, $threeprime ) = @_;

    if ( $gene_length == -1 ) {
        return $gene_seq;
    }
    elsif ( $fiveprime ) {
        return substr $gene_seq, 0, $gene_length;
    }
    elsif ( $threeprime ) {
        return substr $gene_seq, -$gene_length;
    }
}

sub combine_seqs {
    my ( $fiveprime, $threeprime, $both, $utr_seq, $gene_seq ) = @_;

    my $combo_seq;
    if ($fiveprime) {
        $combo_seq = $$utr_seq[0] . $gene_seq;
    }
    elsif ($threeprime) {
        $combo_seq = $gene_seq . $$utr_seq[0];
    }
    elsif ($both) {
        $combo_seq = $$utr_seq[0] . $gene_seq . $$utr_seq[1];
    }
    return $combo_seq;
}

sub output_fa {
    my ( $seqid, $seq, $output_fa_fh, $fa_width ) = @_;

    $fa_width //= 80;
    my @fa_seq = unpack "(A$fa_width)*", $seq;

    say $output_fa_fh ">$seqid";
    say $output_fa_fh $_ for @fa_seq;
}
