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

my $gff_file = $ARGV[0];

my $cds_fa_file    = "~/git.repos/sample-files/fa/ITAG2.3_cds.fasta";
my $genome_fa_file = "~/git.repos/sample-files/fa/ITAG2.3_genomic.fasta";

my $utr_length  = 500;
my $gene_length = 500;

my $coding_regions = extract_cds_from_gff( $gff_file );

for my $id ( keys %$coding_regions ) {
    my $chr       = $$coding_regions{$id}{chr};
    my $strand    = $$coding_regions{$id}{strand};
    my $left_pos  = $$coding_regions{$id}{left_pos};
    my $right_pos = $$coding_regions{$id}{right_pos};

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
    else { die "Problem with strand info for $id\n" }

    say $id;
    say $utr3_start;
    say $utr3_end;
    system "~/installs/bin/samtools faidx $genome_fa_file $chr:$utr3_start-$utr3_end";
}

p $coding_regions;

exit;

sub extract_cds_from_gff {
    my $gff_file = shift;

    # open my $gff_fh, "<", $gff_file;

    my $gff_header = <DATA>;
    check_gff_version(3, $gff_header);

    my %coding_regions;
    # build_coding_regions_hash( \%coding_regions, $gff_fh );
    build_coding_regions_hash( \%coding_regions );
    # close $gff_fh;

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

    # while ( my $feature = <$gff_fh> ) {
    while ( my $feature = <DATA> ) {
        next if $feature =~ /^#/;
        my ( $chr, $type, $start, $end, $strand, $attributes ) =
          ( split /\t/, $feature )[ 0, 2 .. 4, 6, 8 ];
        next unless $type eq "CDS";

        my ( $gene ) = $attributes =~ /Parent=(?:mRNA:)?([^;]+)/;
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
###
SL2.40ch00	ITAG_eugene	gene	570544	575454	.	-	.	Alias=Solyc00g005050;ID=gene:Solyc00g005050.2;Name=Solyc00g005050.2;from_BOGAS=1;length=4911
SL2.40ch00	ITAG_eugene	mRNA	570544	575454	.	-	.	ID=mRNA:Solyc00g005050.2.1;Name=Solyc00g005050.2.1;Note=Arabinogalactan protein (AHRD V1 ***- B6SST2_MAIZE);Parent=gene:Solyc00g005050.2;from_BOGAS=1;length=4911;nb_exon=5
SL2.40ch00	ITAG_eugene	exon	570544	570897	.	-	.	ID=exon:Solyc00g005050.2.1.5;Parent=mRNA:Solyc00g005050.2.1;from_BOGAS=1
SL2.40ch00	ITAG_eugene	three_prime_UTR	570544	570803	.	-	.	ID=three_prime_UTR:Solyc00g005050.2.1.0;Parent=mRNA:Solyc00g005050.2.1;from_BOGAS=1
SL2.40ch00	ITAG_eugene	CDS	570804	570897	.	-	1	ID=CDS:Solyc00g005050.2.1.5;Parent=mRNA:Solyc00g005050.2.1;from_BOGAS=1
SL2.40ch00	ITAG_eugene	intron	570898	571002	.	-	.	ID=intron:Solyc00g005050.2.1.5;Parent=mRNA:Solyc00g005050.2.1;from_BOGAS=1
SL2.40ch00	ITAG_eugene	exon	571003	571208	.	-	1	ID=exon:Solyc00g005050.2.1.4;Parent=mRNA:Solyc00g005050.2.1;from_BOGAS=1
SL2.40ch00	ITAG_eugene	CDS	571003	571208	.	-	0	ID=CDS:Solyc00g005050.2.1.4;Parent=mRNA:Solyc00g005050.2.1;from_BOGAS=1
SL2.40ch00	ITAG_eugene	intron	571209	571284	.	-	.	ID=intron:Solyc00g005050.2.1.4;Parent=mRNA:Solyc00g005050.2.1;from_BOGAS=1
SL2.40ch00	ITAG_eugene	exon	571285	571356	.	-	0	ID=exon:Solyc00g005050.2.1.3;Parent=mRNA:Solyc00g005050.2.1;from_BOGAS=1
SL2.40ch00	ITAG_eugene	CDS	571285	571356	.	-	0	ID=CDS:Solyc00g005050.2.1.3;Parent=mRNA:Solyc00g005050.2.1;from_BOGAS=1
SL2.40ch00	ITAG_eugene	intron	571357	574466	.	-	.	ID=intron:Solyc00g005050.2.1.3;Parent=mRNA:Solyc00g005050.2.1;from_BOGAS=1
SL2.40ch00	ITAG_eugene	exon	574467	574572	.	-	0	ID=exon:Solyc00g005050.2.1.2;Parent=mRNA:Solyc00g005050.2.1;from_BOGAS=1
SL2.40ch00	ITAG_eugene	CDS	574467	574572	.	-	1	ID=CDS:Solyc00g005050.2.1.2;Parent=mRNA:Solyc00g005050.2.1;from_BOGAS=1
SL2.40ch00	ITAG_eugene	intron	574573	575209	.	-	.	ID=intron:Solyc00g005050.2.1.2;Parent=mRNA:Solyc00g005050.2.1;from_BOGAS=1
SL2.40ch00	ITAG_eugene	exon	575210	575454	.	-	1	ID=exon:Solyc00g005050.2.1.1;Parent=mRNA:Solyc00g005050.2.1;from_BOGAS=1
SL2.40ch00	ITAG_eugene	CDS	575210	575319	.	-	0	ID=CDS:Solyc00g005050.2.1.1;Parent=mRNA:Solyc00g005050.2.1;from_BOGAS=1
SL2.40ch00	ITAG_eugene	five_prime_UTR	575320	575454	.	-	.	ID=five_prime_UTR:Solyc00g005050.2.1.0;Parent=mRNA:Solyc00g005050.2.1;from_BOGAS=1
###
SL2.40ch00	ITAG_eugene	gene	723746	724018	.	-	.	Alias=Solyc00g005060;ID=gene:Solyc00g005060.1;Name=Solyc00g005060.1;from_BOGAS=1;length=273
SL2.40ch00	ITAG_eugene	mRNA	723746	724018	.	-	.	ID=mRNA:Solyc00g005060.1.1;Name=Solyc00g005060.1.1;Note=Unknown Protein (AHRD V1);Parent=gene:Solyc00g005060.1;from_BOGAS=1;length=273;nb_exon=1;eugene_evidence_code=10F0H0E0IEG
SL2.40ch00	ITAG_eugene	exon	723746	724018	.	-	0	ID=exon:Solyc00g005060.1.1.1;Parent=mRNA:Solyc00g005060.1.1;from_BOGAS=1
SL2.40ch00	ITAG_eugene	CDS	723746	724018	.	-	0	ID=CDS:Solyc00g005060.1.1.1;Parent=mRNA:Solyc00g005060.1.1;from_BOGAS=1
