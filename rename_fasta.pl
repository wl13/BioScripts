#!/usr/bin/perl -w
#
#   rename_fasta.pl -- rename sequences ids in fasta format file.
#
#   Author: Nowind
#   Created: 2012-02-21
#   Updated: 2014-08-25
#   Version: 1.0.1
#
#   Change logs:
#   Version 1.0.0 13/05/30: The initial version.
#   Version 1.0.1 14/08/25: Add support for paired ids in ref_list.




use strict;

use Data::Dumper;
use Getopt::Long;


use MyPerl::FileIO qw(:all);

##################### Main ####################


my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.0.1';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";

my $SOURCE  = (scalar localtime()) . " Version: $VERSION";


my $word_wrap   = 0;
my ($fasta_file, $ref_list, $output, $show_help);
GetOptions(
            "i|fasta=s"         => \$fasta_file,
            
            "refer=s"           => \$ref_list,
            
            "output=s"          => \$output,
            
            "word-wrap=i"       => \$word_wrap,
            
            "help|?"            => \$show_help,
           );

unless( $fasta_file && $ref_list && !$show_help ) {
    print <<EOF;

$0  -- rename sequences ids in fasta format file.

Version: $VERSION

Usage:   perl $0 [options]

Options:
    -i, --fasta      <filename>
        sequence file to be sorted in fasta format, required

    -r, --refer      <filename>
        rename sequence name by related ids listed in this file

    -o, --output     <filename>
        output fasta file, default to STDOUT
    
    -w, --word-wrap  <int>
        maximum length of sequence to write per line, default each sequence
        per line
        
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



my @new_ids = ();
my %new_ids = ();
if ($ref_list) {
    print STDERR ">> Start parsing $ref_list ... ";
    my $fh = getInputFilehandle($ref_list);
    while (<$fh>)
    {
        next if (/#/ || /^\s+$/);
        
        my ($id, $new_id) = (split /\s+/);
        
        if ($new_id) {
            $new_ids{$id} = $new_id;
        }
        else {
            push @new_ids, $id;
        }
    }
    print STDERR "done!\n";
}


print STDERR ">> Start reading $fasta_file ... ";
my @SEQs = ();
my @ids  = parse_fasta_SEQs(\@SEQs, $fasta_file);
print STDERR "done!\n";




print STDERR ">> Start generating renamed file ... ";
for (my $i=0; $i<@ids; $i++)
{
    my $id = $ids[$i];
    
    unless($SEQs[$i]) {
        print STDERR "Error: $fasta_file: no records found for $id!\n"; exit(2);
    }
    
    if (@new_ids > 0) {
        print format_fasta_SEQs($new_ids[$i], \$SEQs[$i], $word_wrap);
    }
    else {
        print format_fasta_SEQs($new_ids{$id}, \$SEQs[$i], $word_wrap);
    }
}
print STDERR "done!\n";


print STDERR "# " . (scalar localtime()) . "\n";


######################### Sub #########################



















