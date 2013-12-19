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
use Data::Printer;

my %coding_regions;

my $gff_header = <DATA>;
check_gff_version(3, $gff_header);

while ( my $feature = <DATA> ) {
    next if $feature =~ /^#/;
    my ( $type, $start, $end, $strand, $attributes ) =
      ( split /\t/, $feature )[ 2 .. 4, 6, 8 ];
    next unless $type eq "CDS";

    my ( $gene ) = $attributes =~ /Parent=(?:mRNA:)?([^;]+)/;
    $coding_regions{$gene}{strand} = $strand;
    push @{ $coding_regions{$gene}{pos} }, [ $start, $end ];
}

p %coding_regions;
sub check_gff_version {
    my ( $required_version, $gff_version ) = @_;
    die "Requires gff file to be version $required_version\n"
      unless $gff_version =~ /$required_version/;
}

__DATA__
##gff-version 3
##feature-ontology http://song.cvs.sourceforge.net/*checkout*/song/ontology/sofa.obo?revision=1.93
##sequence-region SL2.40ch00 1 21805821
SL2.40ch00	ITAG_eugene	gene	16437	18189	.	+	.	Alias=Solyc00g005000;ID=gene:Solyc00g005000.2;Name=Solyc00g005000.2;from_BOGAS=1;length=1753
SL2.40ch00	ITAG_eugene	mRNA	16437	18189	.	+	.	ID=mRNA:Solyc00g005000.2.1;Name=Solyc00g005000.2.1;Note=Aspartic proteinase nepenthesin I (AHRD V1 **-- A9ZMF9_NEPAL)%3B contains Interpro domain(s)  IPR001461  Peptidase A1 ;Ontology_term=GO:0006508;Parent=gene:Solyc00g005000.2;from_BOGAS=1;interpro2go_term=GO:0006508;length=1753;nb_exon=2
SL2.40ch00	ITAG_eugene	exon	16437	17275	.	+	.	ID=exon:Solyc00g005000.2.1.1;Parent=mRNA:Solyc00g005000.2.1;from_BOGAS=1
SL2.40ch00	ITAG_eugene	five_prime_UTR	16437	16479	.	+	.	ID=five_prime_UTR:Solyc00g005000.2.1.0;Parent=mRNA:Solyc00g005000.2.1;from_BOGAS=1
SL2.40ch00	ITAG_eugene	CDS	16480	17275	.	+	0	ID=CDS:Solyc00g005000.2.1.1;Parent=mRNA:Solyc00g005000.2.1;from_BOGAS=1
SL2.40ch00	ITAG_eugene	intron	17276	17335	.	+	.	ID=intron:Solyc00g005000.2.1.1;Parent=mRNA:Solyc00g005000.2.1;from_BOGAS=1
SL2.40ch00	ITAG_eugene	exon	17336	18189	.	+	0	ID=exon:Solyc00g005000.2.1.2;Parent=mRNA:Solyc00g005000.2.1;from_BOGAS=1
SL2.40ch00	ITAG_eugene	CDS	17336	17940	.	+	2	ID=CDS:Solyc00g005000.2.1.2;Parent=mRNA:Solyc00g005000.2.1;from_BOGAS=1
SL2.40ch00	ITAG_eugene	three_prime_UTR	17941	18189	.	+	.	ID=three_prime_UTR:Solyc00g005000.2.1.0;Parent=mRNA:Solyc00g005000.2.1;from_BOGAS=1
###
SL2.40ch00	ITAG_eugene	gene	68062	68764	.	+	.	Alias=Solyc00g005020;ID=gene:Solyc00g005020.1;Name=Solyc00g005020.1;from_BOGAS=1;length=703
SL2.40ch00	ITAG_eugene	mRNA	68062	68764	.	+	.	ID=mRNA:Solyc00g005020.1.1;Name=Solyc00g005020.1.1;Note=Unknown Protein (AHRD V1);Parent=gene:Solyc00g005020.1;from_BOGAS=1;length=703;nb_exon=3;eugene_evidence_code=10F0H0E0IEG
SL2.40ch00	ITAG_eugene	exon	68062	68211	.	+	0	ID=exon:Solyc00g005020.1.1.1;Parent=mRNA:Solyc00g005020.1.1;from_BOGAS=1
SL2.40ch00	ITAG_eugene	CDS	68062	68211	.	+	0	ID=CDS:Solyc00g005020.1.1.1;Parent=mRNA:Solyc00g005020.1.1;from_BOGAS=1
SL2.40ch00	ITAG_eugene	intron	68212	68343	.	+	.	ID=intron:Solyc00g005020.1.1.1;Parent=mRNA:Solyc00g005020.1.1;from_BOGAS=1
SL2.40ch00	ITAG_eugene	exon	68344	68568	.	+	0	ID=exon:Solyc00g005020.1.1.2;Parent=mRNA:Solyc00g005020.1.1;from_BOGAS=1
SL2.40ch00	ITAG_eugene	CDS	68344	68568	.	+	0	ID=CDS:Solyc00g005020.1.1.2;Parent=mRNA:Solyc00g005020.1.1;from_BOGAS=1
SL2.40ch00	ITAG_eugene	intron	68569	68653	.	+	.	ID=intron:Solyc00g005020.1.1.2;Parent=mRNA:Solyc00g005020.1.1;from_BOGAS=1
SL2.40ch00	ITAG_eugene	exon	68654	68764	.	+	0	ID=exon:Solyc00g005020.1.1.3;Parent=mRNA:Solyc00g005020.1.1;from_BOGAS=1
SL2.40ch00	ITAG_eugene	CDS	68654	68764	.	+	0	ID=CDS:Solyc00g005020.1.1.3;Parent=mRNA:Solyc00g005020.1.1;from_BOGAS=1
###
SL2.40ch00	ITAG_eugene	gene	550920	551576	.	+	.	Alias=Solyc00g005040;ID=gene:Solyc00g005040.2;Name=Solyc00g005040.2;from_BOGAS=1;length=657
SL2.40ch00	ITAG_eugene	mRNA	550920	551576	.	+	.	ID=mRNA:Solyc00g005040.2.1;Name=Solyc00g005040.2.1;Note=Potassium channel (AHRD V1 ***- D0EM91_9ROSI)%3B contains Interpro domain(s)  IPR000595  Cyclic nucleotide-binding ;Parent=gene:Solyc00g005040.2;from_BOGAS=1;length=657;nb_exon=4
SL2.40ch00	ITAG_eugene	exon	550920	550945	.	+	.	ID=exon:Solyc00g005040.2.1.1;Parent=mRNA:Solyc00g005040.2.1;from_BOGAS=1
SL2.40ch00	ITAG_eugene	five_prime_UTR	550920	550945	.	+	.	ID=five_prime_UTR:Solyc00g005040.2.1.0;Parent=mRNA:Solyc00g005040.2.1;from_BOGAS=1
SL2.40ch00	ITAG_eugene	intron	550946	551033	.	+	.	ID=intron:Solyc00g005040.2.1.1;Parent=mRNA:Solyc00g005040.2.1;from_BOGAS=1
SL2.40ch00	ITAG_eugene	exon	551034	551132	.	+	.	ID=exon:Solyc00g005040.2.1.2;Parent=mRNA:Solyc00g005040.2.1;from_BOGAS=1
SL2.40ch00	ITAG_eugene	five_prime_UTR	551034	551042	.	+	.	ID=five_prime_UTR:Solyc00g005040.2.1.1;Parent=mRNA:Solyc00g005040.2.1;from_BOGAS=1
SL2.40ch00	ITAG_eugene	CDS	551043	551132	.	+	0	ID=CDS:Solyc00g005040.2.1.1;Parent=mRNA:Solyc00g005040.2.1;from_BOGAS=1
SL2.40ch00	ITAG_eugene	intron	551133	551217	.	+	.	ID=intron:Solyc00g005040.2.1.2;Parent=mRNA:Solyc00g005040.2.1;from_BOGAS=1
SL2.40ch00	ITAG_eugene	exon	551218	551250	.	+	0	ID=exon:Solyc00g005040.2.1.3;Parent=mRNA:Solyc00g005040.2.1;from_BOGAS=1
SL2.40ch00	ITAG_eugene	CDS	551218	551250	.	+	0	ID=CDS:Solyc00g005040.2.1.2;Parent=mRNA:Solyc00g005040.2.1;from_BOGAS=1
SL2.40ch00	ITAG_eugene	intron	551251	551342	.	+	.	ID=intron:Solyc00g005040.2.1.3;Parent=mRNA:Solyc00g005040.2.1;from_BOGAS=1
SL2.40ch00	ITAG_eugene	exon	551343	551576	.	+	0	ID=exon:Solyc00g005040.2.1.4;Parent=mRNA:Solyc00g005040.2.1;from_BOGAS=1
SL2.40ch00	ITAG_eugene	CDS	551343	551576	.	+	0	ID=CDS:Solyc00g005040.2.1.3;Parent=mRNA:Solyc00g005040.2.1;from_BOGAS=1
