#!/usr/bin/env perl
# Mike Covington
# created: 2013-12-23
#
# Description:
#
use strict;
use warnings;
use autodie;
use Test::More tests => 3;
use FindBin qw($Bin);

my $base_extract_cmd = <<CMD;
$Bin/../bin/extract-utr.pl \\
  --gff_file $Bin/sample-files/ITAG2.3_gene_models.truncated.gff3 \\
  --cds_fa_file $Bin/sample-files/ITAG2.3_cds.truncated.fasta \\
  --genome_fa_file $Bin/sample-files/ITAG2.3_genomic.truncated.fasta \\
  --output_fa_file $Bin/got.fa \\
CMD
my $test_name;
my $extract_cmd;

subtest "5'-UTR extractions" => sub {
    plan tests => 3;

    $test_name   = "500-5prime+full-cds";
    $extract_cmd = "$base_extract_cmd --fiveprime --gene_length -1";
    compare_extracted_utr( $extract_cmd, $test_name );

    $test_name   = "500-5prime+500-cds";
    $extract_cmd = "$base_extract_cmd --fiveprime";
    compare_extracted_utr( $extract_cmd, $test_name );

    $test_name   = "500-5prime+no-cds";
    $extract_cmd = "$base_extract_cmd --fiveprime --gene_length 0";
    compare_extracted_utr( $extract_cmd, $test_name );
};

subtest "3'-UTR extractions" => sub {
    plan tests => 3;

    $test_name   = "full-cds+500-3prime";
    $extract_cmd = "$base_extract_cmd --gene_length -1 --threeprime";
    compare_extracted_utr( $extract_cmd, $test_name );

    $test_name   = "500-cds+500-3prime";
    $extract_cmd = "$base_extract_cmd --threeprime";
    compare_extracted_utr( $extract_cmd, $test_name );

    $test_name   = "no-cds+500-3prime";
    $extract_cmd = "$base_extract_cmd --gene_length 0 --threeprime";
    compare_extracted_utr( $extract_cmd, $test_name );
};

subtest "Special cases" => sub {
    plan tests => 2;

    $test_name   = "500-cds+500-3prime.60wide";
    $extract_cmd = "$base_extract_cmd --fa_width 60 --threeprime";
    compare_extracted_utr( $extract_cmd, $test_name );

    $test_name   = "500-5prime+full-cds+500-3prime";
    $extract_cmd = "$base_extract_cmd --both";
    compare_extracted_utr( $extract_cmd, $test_name );
};

sub compare_extracted_utr {
    my ( $extract_cmd, $test_name ) = @_;

    system($extract_cmd);

    my $expect_file = "$Bin/sample-files/expect.$test_name.fa";
    open my $expect_fh, "<", $expect_file;
    my @expected = <$expect_fh>;
    close $expect_fh;

    open my $got_fh, "<", "$Bin/got.fa";
    my @got = <$got_fh>;
    close $got_fh;

    is_deeply( \@got, \@expected, $test_name );

    unlink "$Bin/got.fa";
}
