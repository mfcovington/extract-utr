#!/usr/bin/env perl
# Mike Covington
# created: 2016-10-21
#
# Description:
#
use strict;
use warnings;
use autodie;
use feature 'say';
use Getopt::Long;

my ( $input, $output, $delimiter );

my $options = GetOptions(
    "input=s"     => \$input,
    "output=s"    => \$output,
    "delimiter=s" => \$delimiter,
);

check_options( $input, $output, $delimiter );

open my $in_fh, "<", $input;
open my $out_fh, ">", $output;
while (<$in_fh>) {
    if (my ($id, $remainder) = $_ =~ /^>([^$delimiter]+)(.*)/) {
        say $out_fh ">$id - $id$remainder";
    }
    else {
        print $out_fh $_;
    }
}

close $in_fh;
close $out_fh;


exit;

sub check_options {
    my ( $input, $output, $delimiter ) = @_;

    die "Specify '--input FILE.'\n"  unless defined $input;
    die "Specify '--output FILE.'\n" unless defined $output;
    die "Specify the sequence ID the delimiter (i.e., '--delimiter .' " .
        "for the sequence IDs that look like 'ENSMUST00000178537.1').\n"
        unless defined $delimiter;
}
