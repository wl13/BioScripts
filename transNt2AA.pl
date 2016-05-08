#!/usr/bin/perl -w
#
#   transNt2AA.pl -- translate nucleotides to proteins
#
#   Author: Nowind
#   Created: 2010-09-29
#   Updated: 2012-08-20
#   Version: 1.0.0
#
#   Change logs:
#   Version 1.0.0 12/08/20: The initial version.

use strict;

use Bio::SeqIO;


######################## Main ########################

my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.0.0';
my $HEADER  = "# $CMDLINE\n# Version: $VERSION\n";

my ($input, $output) = @ARGV;


unless( $input ) {
    print <<EOF;

$0  -- translate nucleotides to proteins

Version: $VERSION

Usage:   perl $0 <input file> [output file]

EOF

    exit(1);
}

unless( $output ) { ($output = $input) =~ s/\.(fasta|fas)/_protein.fasta/; }

my $seqin  = Bio::SeqIO->new( -file   => "$input",
                              -format => 'fasta');
my $seqout = Bio::SeqIO->new( -file   => "> $output",
                              -format => 'fasta');

print "Translate nucleotide sequence to protein sequence...";
while (my $seqobj = $seqin->next_seq())
{
    my $trans = $seqobj->translate();
    $seqout->write_seq($trans);
}
print "\tdone!\n";
