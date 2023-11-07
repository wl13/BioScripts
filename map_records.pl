#!/usr/bin/perl -w
#
#   map_records.pl -- mapping records from one file to another file
#                          
#
#   Author: Nowind
#   Created: 2012-05-31
#   Updated: 2023-10-30
#   Version: 1.2.0
#
#   Change logs:
#   Version 1.0.0 12/09/18: The initial version.
#   Version 1.1.0 12/11/19: Add option "--rows2" to specify rows in subject file.
#   Version 1.1.1 12/12/07: Query and subject rows now can have only a single value.
#   Version 1.1.2 12/12/29: Change the way to import functions from MyPerl::FileIO.
#   Version 1.1.3 13/01/23: Add option "--interval" and "--flanking" to map records in intervals.
#   Version 1.1.4 13/03/14: Change output header.
#   Version 1.1.5 13/03/15: Bug fixed while output file header.
#   Version 1.1.6 13/04/11: Bug fixed while reading query file from STDIN.
#   Version 1.1.7 16/05/31: Updated: Add option "--no-secondary"; remove mapped subject rows to reduce redundancy.
#   Version 1.1.8 16/06/13: Bug fixed while no remain fields present in subject file.
#   Version 1.1.9 16/02/16: Deprecated non-functional "-q" and "-s" options.
#   Version 1.2.0 23/10/30: Bug fixed: --no-dups prints empty results; add support for setting delimiters separately for query and subject files.

use strict;

use Data::Dumper;
use Getopt::Long;
use File::Find::Rule;
use File::Basename;

use MyPerl::FileIO qw(:all);

################################# Main ###############################


my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.2.0';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";


my $qdelimiter  = '\s+';
my $sdelimiter  = '\s+';
my $flank_len  = 0;
my ($query_file, $subject_file, $output, $interval,
    @query_rows, @subject_rows, $match_str, @sub_set, $out_first_only, $merge_sub_dups);
GetOptions(
            "query=s"          => \$query_file,
            "subject=s"        => \$subject_file,
            "output=s"         => \$output,
            "qdelimiter=s"     => \$qdelimiter,
            "sdelimiter=s"     => \$sdelimiter,
            "interval"         => \$interval,
            "flanking=i"       => \$flank_len,
            "match=s"          => \$match_str,
            "Q|rows1=i{,}"     => \@query_rows,
            "S|rows2=i{,}"     => \@subject_rows,
            "no-dups"          => \$out_first_only,
            "merge-dups"       => \$merge_sub_dups,
           );

unless( $query_file && $subject_file ) {
    print <<EOF;

$0  -- query fasta sequences

Version: $VERSION

Usage:   perl $0 [options]

Options:
    --query    <filename>
        query file, required
    --subject  <filename>
        subject file, required
    -o, --output   <filename>
        output file, default to STDOUT
    
    -Q, --rows1   <numbers>
        specify rows for comparing in query file [default: 0 1]
    -S, --rows2   <numbers>
        specify rows for comparing in subject file [default: 0 1]
    
    -i, --interval
        use this option  to map all records in the interval between
        query start and end, while this option is set, the -Q option
        should be specified with following rows in the order:
        
        start position, end position, [other keys...]
        
        while the -S option should be specified as:
        
        position, [other keys...]
    
    -f, --flanking  <interger>
        flanking region length to extend the interval, only valid while
        -i option is specified

    --qdelimiter <string>
    --sdelimiter <string>
        specify a delimiter while reading input files, such as ",",
        "\\t", multiple delimiters can be set such as ",|\\t"
        [default: "\\s+" for both query and subject]
    --match     <string>
        only considering lines matching a pattern, support perl
        regular expression
    
    -n, --no-dups
        if query record have multiple hits in subject file, only the first
        one would be written to output
    --merge-dups
        if query record have multiple hits in subject file, all hits would
        be written in a single record and seperate by semicolon
    
EOF

    exit(1);
}

$|++;


print STDERR "# $0 v$VERSION\n# " . (scalar localtime()) . "\n";

if ($output) {
    open (STDOUT, "> $output") || die $!;
}

unless(@query_rows >= 1){ @query_rows = (0, 1) };
unless(@subject_rows >= 1){ @subject_rows = (0, 1) };


##
## read and parsing query file
##
print STDERR ">> Read in $query_file ... ";
my %Query         = ();
my @Query_Records = ();
my $query_header = '';
my $fh1   = getInputFilehandle($query_file);
while (<$fh1>)
{
    push @Query_Records, $_;
    
    if (/^\#\w/) {
        chomp($query_header = $_);
    }    
    
    next if (/\#/ || /^\s+$/);
    next if ($match_str && !/$match_str/);
    
    if ($interval) {
        my ($start, $end, @cmp_rows) = (split /$qdelimiter/, $_)[@query_rows];
        
        for my $pos (($start-$flank_len)..($end+$flank_len))
        {
            my $cmp_rows = join "\t", ($pos, @cmp_rows);
            $Query{$cmp_rows}->{query} = 1;
        }
    }
    else {
        my $cmp_rows = join "\t", ((split /$qdelimiter/, $_)[@query_rows]);
        
        $Query{$cmp_rows}->{query} = 1;
    }
}
print STDERR "done!\n";

##
## retrieving query records in subject file
##
print STDERR ">> Search records in $subject_file ... ";
my $subject_header = '';
my $fh2   = getInputFilehandle($subject_file);
while (<$fh2>)
{  
    if (/\#\#/ || /^\s+$/) {
        next;
    }
    
    chomp;
    
    my @all_rows = (split /$sdelimiter/, $_);
    
    my $cmp_rows = join "\t", @all_rows[@subject_rows];
    
    @all_rows[@subject_rows] = ();
    my @remain_rows = grep defined, @all_rows;
    
    my $remain_rows = (@remain_rows > 0) ? join "\t", @remain_rows : "Found";
    
    if (/^\#\w/) {
        $subject_header = $remain_rows;
        next;
    }
    
    next unless( $Query{$cmp_rows}->{query} );
    
    push @{$Query{$cmp_rows}->{cmp}}, $remain_rows;
}
print STDERR "done!\n";


##
## generate results
##
print STDERR ">> Start generate results ... ";
print STDOUT "$HEADER##" . (scalar localtime()) . "\n";
foreach (@Query_Records)
{
    if (/\#\#/ || /^\s+$/) {
        print STDOUT; next;
    }
    elsif (/^\#\w/) {
        print STDOUT "$query_header\t$subject_header\n";
    }
    
    next if (/\#/ || /^\s+$/);
    next if ($match_str && !/$match_str/);
    
    chomp;
    
    if ($interval) {
        my ($start, $end, @cmp_rows) = (split /$qdelimiter/, $_)[@query_rows];
        
        for my $pos (($start-$flank_len)..($end+$flank_len))
        {
            my $cmp_rows = join "\t", ($pos, @cmp_rows);
            
            unless( $Query{$cmp_rows}->{cmp} ) {
                print STDOUT "$_\tN/A\n";
                next;
            }
            
            my @out_records = ();
            for my $record (@{$Query{$cmp_rows}->{cmp}})
            {
                
                if ($out_first_only) {
                    print STDOUT "$_\t$record\n";
                    last;
                }
                elsif ($merge_sub_dups) {
                    chomp($record);
                    push @out_records, $record;
                }
                else {
                    print STDOUT "$_\t$record\n";
                }
            }
            
            if ($merge_sub_dups) {
                my $merged_record = join ";", @out_records;
                
                print STDOUT "$_\t$merged_record\n";
            }
        }        
    }
    else {
        my $cmp_rows = join "\t", ((split /$qdelimiter/, $_)[@query_rows]);
        
        unless( $Query{$cmp_rows}->{cmp} ) {
            print STDOUT "$_\tN/A\n";
            next;
        }
        
        my @out_records = ();
        for my $record (@{$Query{$cmp_rows}->{cmp}})
        {
            
            if ($out_first_only) {
                print STDOUT "$_\t$record\n";
                last;
            }
            elsif ($merge_sub_dups) {
                chomp($record);
                push @out_records, $record;
            }
            else {
                print STDOUT "$_\t$record\n";
            }
        }
        
        if ($merge_sub_dups) {
            my $merged_record = join ";", @out_records;
            
            print STDOUT "$_\t$merged_record\n";
        }
    }

}
print STDERR "done!\n";


print STDERR "# " . (scalar localtime()) . "\n";
