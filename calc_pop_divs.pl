#!/usr/bin/perl -w
#
#   calc_pop_divs.pl -- Calculating population pi values.
#                          
#
#   Author: Nowind
#   Created: 2012-05-31
#   Updated: 2016-02-10
#   Version: 1.0.0
#
#   Change logs:
#   Version 1.0.0 16/02/10: The initial version.




use strict;

use Data::Dumper;
use Getopt::Long;
use File::Find::Rule;
use File::Basename;

use MyPerl::FileIO qw(:all);

################################# Main ###############################


my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.0.0';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";


my $min_info_perc = 50;
my ($query_file, $subject_file, $output);
GetOptions(
            "query=s"           => \$query_file,
            "subject=s"         => \$subject_file,
            "output=s"          => \$output,
            "min-info-perc=f"   => \$min_info_perc,
           );

unless( $query_file && $subject_file ) {
    print <<EOF;

$0  -- query fasta sequences

Version: $VERSION

Usage:   perl $0 [options]

Options:
    -q, --query    <filename>
        nucleotide differences between each pairs, required
    -s, --subject  <filename>
        informative sites between each pairs required, required
    -o, --output   <filename>
        output file, default to STDOUT
    
    -m, --min-info-perc <float>
        windows with informative sites less than this percentage will be
        discarded from calculation [default: 50]
    
EOF

    exit(1);
}

$|++;


print STDERR "# $0 v$VERSION\n# " . (scalar localtime()) . "\n";

if ($output) {
    open (STDOUT, "> $output") || die $!;
}


##
## read and parsing query file
##
print STDERR ">> Start parsing in $query_file ... ";
my %pair_diffs     = ();
my @query_pair_ids = ();
my $fh1   = getInputFilehandle($query_file);
while (<$fh1>)
{
    next if (/\#\#/ || /^\s+$/);
    
    my ($bin_id, $chrom, $bin_start, $bin_end, @pairs) = (split /\s+/);
    
    if (/^\#BIN_ID/) {
        for (my $i=0; $i<@pairs; $i++)
        {
            push @query_pair_ids, $pairs[$i];
        }
        next;
    }    

    for (my $i=0; $i<@pairs; $i++)
    {
        $pair_diffs{"$chrom\t$bin_start\t$bin_end"}->{$query_pair_ids[$i]} = $pairs[$i] * ($bin_end - $bin_start + 1);
    }
}
print STDERR "done!\n";


##
## retrieving query records in subject file
##
print STDERR ">> Start parsing $subject_file ... ";
print STDOUT "$HEADER##" . (scalar localtime()) . "\n";
print STDOUT "#BIN_ID\tCHROM\tBIN_START\tBIN_END\tNo_Of_Pairs\tPop_Divers\n";
my @sub_pair_ids   = ();
my $fh2   = getInputFilehandle($subject_file);
while (<$fh2>)
{
    next if (/\#\#/ || /^\s+$/);
    
    my ($bin_id, $chrom, $bin_start, $bin_end, @pairs) = (split /\s+/);
    
    if (/^\#BIN_ID/) {
        for (my $i=0; $i<@pairs; $i++)
        {
            push @sub_pair_ids, $pairs[$i];
        }
        next;
    }    
    
    my $bin_size = $bin_end - $bin_start + 1;
    
    my %pair_divs      = ();
       $pair_divs{sum} = 0;
       $pair_divs{num} = 0;
       
    for (my $i=0; $i<@pairs; $i++)
    {
        if (100 * $pairs[$i] / $bin_size >= $min_info_perc) {
            $pair_divs{sum} += $pair_diffs{"$chrom\t$bin_start\t$bin_end"}->{$sub_pair_ids[$i]} / $pairs[$i];
            $pair_divs{num} ++;
        } 
    }
    
    my $pop_divs = $pair_divs{num} > 0 ? $pair_divs{sum} / $pair_divs{num} : -1;
    
    print "$bin_id\t$chrom\t$bin_start\t$bin_end\t$pair_divs{num}\t$pop_divs\n";
}
print STDERR "done!\n";



print STDERR "# " . (scalar localtime()) . "\n";
