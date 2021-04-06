#!/usr/bin/perl -w
#
#   calc_pop_divs.pl -- Calculating population pi values.
#                          
#
#   Author: Nowind
#   Created: 2012-05-31
#   Updated: 2016-10-05
#   Version: 1.1.0
#
#   Change logs:
#   Version 1.0.0 16/02/10: The initial version.
#   Version 1.1.0 16/10/05: Updated: add option "--all-pairs" to output corrected diversities for all pairs.




use strict;

use Data::Dumper;
use Getopt::Long;
use File::Find::Rule;
use File::Basename;

use MyPerl::FileIO qw(:all);

################################# Main ###############################


my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.1.0';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";


my $min_info_perc = 50;
my ($query_file, $subject_file, $output, $out_all_pairs);
GetOptions(
            "query=s"           => \$query_file,
            "subject=s"         => \$subject_file,
            "output=s"          => \$output,
            "min-info-perc=f"   => \$min_info_perc,
            
            "all-pairs"         => \$out_all_pairs,
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
    
    -a, --all-pairs
        output corrected diversities for all pairs

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
unless ($out_all_pairs) {
    print STDOUT "#BIN_ID\tCHROM\tBIN_START\tBIN_END\tNo_Of_Pairs\tPop_Divers\n";
}
my @sub_pair_ids   = ();
my $fh2   = getInputFilehandle($subject_file);
while (<$fh2>)
{
    next if (/\#\#/ || /^\s+$/);
    
    my ($bin_id, $chrom, $bin_start, $bin_end, @pairs) = (split /\s+/);
    
    if (/^\#BIN_ID/) {
        if ($out_all_pairs) {
            print STDOUT;
        }
        
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
            my $divs = $pair_diffs{"$chrom\t$bin_start\t$bin_end"}->{$sub_pair_ids[$i]} / $pairs[$i];;
            
            $pair_divs{sum} += $divs;
            $pair_divs{num} ++;
            
            push @{$pair_divs{all}}, $divs;
        }
        else {
            push @{$pair_divs{all}}, -1;
        }
    }
    
    my $out_divs = "$pair_divs{num}\t";
    
    if ($out_all_pairs) {
        $out_divs = join "\t", @{$pair_divs{all}};
    }
    else {
        $out_divs .= $pair_divs{num} > 0 ? $pair_divs{sum} / $pair_divs{num} : -1;
    }
    
    print "$bin_id\t$chrom\t$bin_start\t$bin_end\t$out_divs\n";
}
print STDERR "done!\n";



print STDERR "# " . (scalar localtime()) . "\n";
