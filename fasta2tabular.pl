#!/usr/bin/perl -w
#
#   fasta2tabular.pl -- Convert sequence file from fasta format to tabular format.
#
#   Author: Nowind
#   Created: 2012-02-21
#   Updated: 2013-05-18
#   Version: 1.0.0
#
#   Change logs:
#   Version 1.0.0 13/05/18: The initial version.




use strict;

use Data::Dumper;
use Getopt::Long;


use MyPerl::FileIO qw(:all);

##################### Main ####################


###my $test_str = 'ACAT-GTC*C?CCAAANN??**---';
###print "$test_str\n";
###
###   $test_str =~ tr/ATGC\-N?*ZM/123456/;
###   
###print "$test_str\n";exit;


my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.0.0';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";

my $SOURCE  = (scalar localtime()) . " Version: $VERSION";


my ($fasta_file, $output, $numeric, $seperate, $show_help);
GetOptions(
            "i|fasta=s"         => \$fasta_file,
            "output=s"          => \$output,
            
            "numeric"           => \$numeric,
            
            "seperate"          => \$seperate,
            
            "help|?"            => \$show_help,
           );

unless( $fasta_file && !$show_help ) {
    print <<EOF;

$0  -- Convert sequence file from fasta format to tabular format.

Version: $VERSION

Usage:   perl $0 [options]

Options:
    -i, --fasta      <filename>
        sequence file to be sorted in fasta format, required
    -o, --output     <filename>
        output fasta file, default to STDOUT
    
    -n, --numeric    <filename>
        convert nucleotide bases to numbers according to the following
        conversions:
        A:1, T:2, G:3, C:4, -:5, N*:6
        *N stands for all unkown bases

    -s, --seperate
        seperate characters using comma
        
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
for my $id (@ids)
{
    if ($numeric) {
        $SEQs{$id} =~ tr/ATGC\-N?*ZM/123456/;
    }
    
    my $out_str = $seperate ? (join ',', (split //, $SEQs{$id})) : $SEQs{$id};
    
    print "$id,$out_str\n";
}
print STDERR "done!\n";


print STDERR "# " . (scalar localtime()) . "\n";


######################### Sub #########################



















