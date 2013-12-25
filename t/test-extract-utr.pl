#!/usr/bin/env perl
# Mike Covington
# created: 2013-12-23
#
# Description:
#
use strict;
use warnings;
use autodie;
use Test::More tests => 1;

my $test_name = "500-cds+500-3prime";
my $extract_cmd = <<CMD;
../extract-utr.pl \\
  --gff_file sample-files/ITAG2.3_gene_models.truncated.gff3 \\
  --cds_fa_file sample-files/ITAG2.3_cds.truncated.fasta \\
  --genome_fa_file sample-files/ITAG2.3_genomic.truncated.fasta \\
  --output_fa_file got.fa \\
  --threeprime
CMD
compare_extracted_utr( $extract_cmd, "sample-files/expect.$test_name.fa", $test_name );

sub compare_extracted_utr {
    my ( $extract_cmd, $expect_file, $test_name ) = @_;

    system($extract_cmd);

    open my $expect_fh, "<", $expect_file;
    my @expected = <$expect_fh>;
    close $expect_fh;

    open my $got_fh, "<", "got.fa";
    my @got = <$got_fh>;
    close $got_fh;

    is_deeply( \@got, \@expected, $test_name );

}
