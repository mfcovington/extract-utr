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

use FindBin;
use lib "$FindBin::Bin/../lib";
use gff;

# TEMPORARY Defaults
# my $gff_file = "Zea_mays.AGPv3.21.gff3";
my $gff_file = glob "~/git.repos/extract-utr/GRMZM2G176820.gff";
my $alignment_file = "~/git.repos/extract-utr/accepted_hits.bam";

my $coding_regions = extract_cds_from_gff($gff_file);

use Data::Printer;
p $coding_regions;

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
