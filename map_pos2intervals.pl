#!/usr/bin/perl -w
#
#   map_pos2intervals.pl -- Get the information of intervals where the query position reside.
#                          
#
#   Author: Nowind
#   Created: 2012-05-31
#   Updated: 2016-06-04
#   Version: 1.2.1
#
#   Change logs:
#   Version 1.0.0 15/08/18: The initial version.
#   Version 1.0.1 15/11/14: Bug fixed: no default value assigned for row numbers.
#   Version 1.1.0 16/05/25: Updated: Remove prerequisites of sample rows.
#   Version 1.2.0 16/05/28: Updated: Output all original records in query file; 
#                                    write all mapped subject records in a single line.
#   Version 1.2.1 16/06/04: Updated: update output file header lines.


use strict;

use Data::Dumper;
use Getopt::Long;
use File::Find::Rule;
use File::Basename;

use MyPerl::FileIO qw(:all);

################################# Main ###############################


my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.2.1';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";


my ($query_file, $subject_file, $output, @query_rows, @subject_rows);
GetOptions(
            "query=s"          => \$query_file,
            "subject=s"        => \$subject_file,
            "output=s"         => \$output,
            "Q|rows1=i{,}"     => \@query_rows,
            "S|rows2=i{,}"     => \@subject_rows,
           );

unless( $query_file && $subject_file ) {
    print <<EOF;

$0  -- Map query positions to regions, a more memory-efficient alternative of map_records.pl

Version: $VERSION

Usage:   perl $0 [options]

Options:
    -q, --query    <filename>
        query file, required
    -s, --subject  <filename>
        subject file, required
    -o, --output   <filename>
        output file, default to STDOUT
    
    -Q, --rows1   <numbers>
        specify chromosome, position and sample(optional) rows in query file;
        if the sample row was used, the subject file should specify the same
        field [default: 0 1]
    -S, --rows2   <numbers>
        specify rows of range (chromosome, start, end, [sample]) in subject
        file [default: 0 1 2]
    
EOF

    exit(1);
}

$|++;


print STDERR "# $0 v$VERSION\n# " . (scalar localtime()) . "\n";

if ($output) {
    open (STDOUT, "> $output") || die $!;
}


unless(@query_rows >= 1){ @query_rows = (0, 1) };
unless(@subject_rows >= 1){ @subject_rows = (0, 1, 2) };



##
## read and parsing query file
##
print STDERR ">> Start parsing $query_file ... ";
my %query_records  = ();
my @query_details  = ();
my $query_header   = 'Query_record';
my $fh1   = getInputFilehandle($query_file);
while (<$fh1>)
{
    if(/^\#\w+/) {
        chomp($query_header = $_);
    }
    
    next if (/\#/ || /^\s+$/);
    
    chomp(my $line = $_);
    
    my ($chrom, $pos, $sample) = (split /\s+/, $line)[@query_rows];
    
    $sample ||= "ALL";
    
    push @{$query_records{$sample}->{$chrom}}, $pos;
    
    push @query_details, $line;
}
print STDERR "done!\n";

##
## retrieving query records in subject file
##
print STDERR ">> Start searching records in $subject_file ... ";
my %subject_records = ();
my $fh2   = getInputFilehandle($subject_file);
while (<$fh2>)
{
    if (/\#/ || /^\s+$/) {
        next;
    }

    chomp(my $interval_info = $_);
    
    my ($chrom, $start, $end, $sample) = (split /\s+/, $interval_info)[@subject_rows];
    
    $sample ||= "ALL";
    
    next unless($query_records{$sample}->{$chrom});
    
    my @mapped_pos = grep {$_ >= $start && $_ <= $end} @{$query_records{$sample}->{$chrom}};
    
    for my $pos (@mapped_pos)
    {
        $interval_info =~ s/\s+/,/g;
        $subject_records{$sample}->{$chrom}->{$pos}->{$interval_info} ++;
    }
}
print STDERR "done!\n";



##
## writing results
##
print STDOUT "$HEADER##" . (scalar localtime()) . "\n";
print STDOUT "$query_header\tMapped_Intervals\n";
for (my $i=0; $i<@query_details; $i++)
{
    my ($chrom, $pos, $sample) = (split /\s+/, $query_details[$i])[@query_rows];
    
    $sample ||= "ALL";
    
    if ($subject_records{$sample}->{$chrom}->{$pos}) {
        my $mapped_records = join "|", (sort keys %{$subject_records{$sample}->{$chrom}->{$pos}});
        
        print STDOUT "$query_details[$i]\t$mapped_records\n";
    }
    else {
        print STDOUT "$query_details[$i]\tN/A\n";
    }
}


print STDERR "# " . (scalar localtime()) . "\n";
