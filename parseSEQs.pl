#!/usr/bin/perl -w
#
#   parseSEQs.pl -- sequence related analysis
#
#   Author: Nowind
#   Created: 2012-02-21
#   Updated: 2015-12-15
#   Version: 1.2.4
#
#   Change logs:
#   Version 1.0.0 12/02/21: The initial version.
#   Version 1.0.1 12/05/29: Add stat of masked info.
#   Version 1.0.2 12/05/30: Error correct in stat repeat regions.
#   Version 1.1.0 12/06/21: Add support to gzipped or bzipped files.
#   Version 1.1.1 13/03/13: Change the way to import functions from MyPerl::FileIO.
#   Version 1.1.2 13/06/23: Change the output order to the same as the input fasta file.
#   Version 1.2.0 14/01/03: Add options to count codons.
#   Version 1.2.1 14/05/29: Add option "--rate" to calculate GC and repeat contents.
#   Version 1.2.2 14/06/17: Bug fixed while "--strict" option is specified; add option "--skip".
#   Version 1.2.3 15/12/07: Bug fixed: use informative length to calculate GC content and repeat content.
#   Version 1.2.4 15/12/15: Updated: add option "--sum-all" to output overall stats; add control
#                           for verbosity.


use strict;

use Data::Dumper;
use Getopt::Long;
use File::Find::Rule;
use File::Basename;


use MyPerl::FileIO qw(:all);

##################### Main ####################

my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.2.4';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";


my $verbosity = 1;
my ($fasta_file, @filter_strings, $output, $count_codon, $strict_mode,
    $calc_ratio, $summary_all);
GetOptions(
            "fasta=s"            => \$fasta_file,
            "output=s"           => \$output,

            "count-codon"        => \$count_codon,
            "K|strict"           => \$strict_mode,
            
            "rate"               => \$calc_ratio,
            
            "T|skip=s{,}"        => \@filter_strings,
            
            "A|sum-all"          => \$summary_all,
            
            "verbosity=i"        => \$verbosity,
           );
unless( $fasta_file ) {
    print <<EOF;

$0  -- parse sequences in fasta format

Version: $VERSION

Usage:   perl $0 [options]

Options:
    -F, --fasta  <filename>
        sequences file in fasta format, required
        
    -o, --output <filename>
        output file, default to STDOUT
        
    -c, --count-codon
        count codon triplets
    -K, --strict
        check and skip sequences not a multiple of 3 while count codons

    -r, --rate
        calculate gc and repeat contents instead of output numbers of all
        nucleotides

    -T, --skip  <strings>
        skip sequences with ids match these strings
    
    -A, --sum-all
        output summaries over all given sequences
        
    -v, --verbosity <number>
        change verbosity levels, 0 indicates no verbosity to STDERR
        [default: 1]
    
EOF

    exit(1);
}

$|++;



if ($output) {
    open (STDOUT, "> $output") || die $!;
}

if ($verbosity > 0) {
    print STDERR "# $0 v$VERSION\n# " . (scalar localtime()) . "\n";
}



my $filter_str = join '|', @filter_strings;

print STDERR ">> Start reading $fasta_file ... " if ($verbosity > 0);
my %SEQs = ();
my @ids  = parse_fasta_SEQs(\%SEQs, $fasta_file, "full");
print STDERR "done!\n" if ($verbosity > 0);



my @chars     = qw(A T G C N a t g c n);
my %codons    = ();
my %stats_all = ();



if ($verbosity > 0) {
    print "$HEADER##" . (scalar localtime()) . "\n";
}

if ($count_codon) {
    print "#Codon\tCount\n";
}
elsif ($calc_ratio) {
    print "#ID\tLength\tInformative_Length\tGC_Num\tGC_Ratio\tRepeat_Len\tRepeat_Ratio\n";
}
else {
    print "#ID\tLength\tCount:A\tT\tG\tC\tN\ta\tt\tg\tc\tn\n";
}
for my $id (@ids)
{
    if (@filter_strings > 0) {
        next if ($id =~ /($filter_str)/);
    }
    
    print STDERR "\r>> Start parsing sequences $id ... "  if $verbosity > 0;
    
    my $seq   = $SEQs{$id};
    my $len   = length $seq;
    
    $stats_all{len} += $len;
    
    my @count = ();
    
    if ($count_codon) {
        if ($strict_mode && ($len % 3)) {
            print STDERR "WARNING: $id sequence length not a multiple of 3!\n";
            next;
        }
        
        my $i = 0;
        while ($i+3 <= $len)
        {
            my $codon = substr($seq, $i, 3);
            $i += 3;
            
            $codons{$codon} ++ ;
        }
    }
    elsif ($calc_ratio) {
        my $info_len    = ($seq =~ tr/ATGCatgc/ATGCatgc/);
        
        my $gc_num      = ($seq =~ tr/GCgc/GCgc/);
        my $gc_rate     = $gc_num / $info_len;
        
        my $repeat_len  = ($seq =~ tr/atgc/atgc/);
        my $repeat_rate = $repeat_len / $info_len;
        
        print "$id\t$len\t$info_len\t$gc_num\t$gc_rate\t$repeat_len\t$repeat_rate\n";
        
        $stats_all{gc_num}   += $gc_num;
        $stats_all{repeat}   += $repeat_len;
        $stats_all{info_len} += $info_len;
    }
    else {
        for my $nt (@chars)
        {
            my $cnt = ($seq =~ s/$nt/$nt/g);
               $cnt = $cnt ? $cnt : 0;
            
            push @count, $cnt;
            
            $stats_all{$nt}   += $cnt;
        }
        
        my $count = join "\t", @count;
        
        print "$id\t$len\t$count\n";
    }
}


print STDERR "done!\n"  if $verbosity > 0;




if ($count_codon) {
    print STDERR ">> Start generating results ... " if $verbosity > 0;
    for my $codon (sort keys %codons)
    {
        my $count = $codons{$codon} ? $codons{$codon} : 0;
        print "$codon\t$count\n";
    }
    print STDERR "done!\n"  if $verbosity > 0;
}


##
## output summaries upon all sequences
##
if ($summary_all) {
    if ($calc_ratio) {
        my $gc_rate_all     = $stats_all{gc_num} / $stats_all{info_len};
        my $repeat_rate_all = $stats_all{repeat} / $stats_all{info_len};
        
        print "All\t$stats_all{len}\t$stats_all{info_len}\t$stats_all{gc_num}\t$gc_rate_all\t$stats_all{repeat}\t$repeat_rate_all\n";
    }
    else {
        my @counts_all = ();
        for my $nt (@chars)
        {
            push @counts_all, $stats_all{$nt};
        }
        
        my $counts_all = join "\t", @counts_all;
        
        print "All\t$stats_all{len}\t$counts_all\n";
    }
}


print STDERR "# " . (scalar localtime()) . "\n" if ($verbosity > 0);





######################### Sub #########################



=head2 

    About   : Determin ancestors of each loci due to comparison of different markers 
    Usage   : compare_markers(\%variants, $vcf_file);
    Args    : A hash to save all compare results
              Vcf file needs to be processed
    Returns : Null

=cut
sub count_base
{
    
}
