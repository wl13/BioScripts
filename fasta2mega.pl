#!/usr/bin/perl -w
#
#   fasta2mega.pl -- convert fasta format to mega format.
#
#   Author: Nowind
#   Created: 2012-02-21
#   Updated: 2013-05-13
#   Version: 1.0.0
#
#   Change logs:
#   Version 1.0.0 13/05/13: The initial version.





use strict;

use Data::Dumper;
use Getopt::Long;


use MyPerl::FileIO qw(:all);

##################### Main ####################


my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.0.0';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";

my $SOURCE  = (scalar localtime()) . " Version: $VERSION";


my $word_wrap    = 60;
my $data_type    = 'DNA';
my $indel_symbol = '-';
my ($fasta_file, $output, $show_help);
GetOptions(
            "i|fasta=s"         => \$fasta_file,
            
            "output=s"          => \$output,
            
            "word-wrap=i"       => \$word_wrap,
            
            "data-type=s"       => \$data_type,
            "symbol=s"          => \$indel_symbol,
            
            "help|?"            => \$show_help,
           );

unless( $fasta_file && !$show_help ) {
    print <<EOF;

$0  -- sort sequences by ids in fasta format file.

Version: $VERSION

Usage:   perl $0 [options]

Options:
    -i, --fasta      <filename>
        sequence file to be sorted in fasta format, required

    -o, --output     <filename>
        output fasta file, default to STDOUT
    
    -d, --data-type  <string>
        specify data type of input fasta file, DNA or Protein, default: DNA
    -s, --symbol     <string>
        symbol for indels, default: "-"
    
    -w, --word-wrap  <int>
        maximum length of sequence to write per line, default 60
        
    -?, --help
        show this help message

EOF

    exit(1);
}




$|++;



if ($output) {
    open (STDOUT, "> $output") || die $!;
}


 
print STDERR "# $0 v$VERSION\n# " . (scalar localtime()) . "\n";


print STDERR ">> Start reading $fasta_file ... ";
my %SEQs = ();
my @ids  = parse_fasta_SEQs(\%SEQs, $fasta_file);
print STDERR "done!\n";



print STDERR ">> Start generating sorted file ... ";
print "#mega\n";
print "!Title $fasta_file;\n";
print "!Format DataType=$data_type indel=$indel_symbol;\n\n";
for my $id (@ids)
{
    my $formated_seq = format_fasta_SEQs($id, \$SEQs{$id}, $word_wrap);
       $formated_seq =~ s/\>/#/g;
    
    print "$formated_seq\n";
}
print STDERR "done!\n";


print STDERR "# " . (scalar localtime()) . "\n";


######################### Sub #########################



















