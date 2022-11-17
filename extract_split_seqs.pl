#!/usr/bin/perl -w
#
#   extract_split_seqs.pl -- Query fasta sequences
#                          
#
#   Author: Nowind
#   Created: 2012-05-31
#   Updated: 2017-02-05
#   Version: 1.0.0
#
#   Change logs:
#   Version 1.0.0 17/02/05: The initial version.


=head1 NAME

extract_split_seqs.pl


=head1 SYNOPSIS

  extract_split_seqs.pl --help/?

=head1 DESCRIPTION

Fasta format file related processes.

=cut


use strict;

use Data::Dumper;
use Getopt::Long;
use File::Find::Rule;
use File::Basename;

use MyPerl::FileIO qw(:all);
use MyPerl::Convert qw(:all);

################################# Main ###############################


my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.0.0';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";
my $SOURCE  = (scalar localtime()) . " Version: $VERSION";

my $convert_pos = 'no';
my ($fasta_file, $output, $word_wrap, $query_file);
GetOptions(
            "fasta=s"          => \$fasta_file,
            "output=s"         => \$output,
            
            "query=s"          => \$query_file,
            
            "wordwrap=i"       => \$word_wrap,         ## Line feed for print
            
            "convert-pos=s"    => \$convert_pos,
           );

unless( $query_file && $fasta_file ) {
    print <<EOF;

$0  -- query fasta sequences

Version: $VERSION

Usage:   perl $0 [options]

Options:

    -f, --fasta  <filename>
        sequences file(s) in fasta format, required
    -q, --query  <filename>
        query file with sequence ids need to be extracted
    
    -o, --output <filename>
        output file, default to STDOUT
        
    -w, --wordwrap  <int>
        line feed for print, only valid while output in fasta format

    -c, --convert-pos <string>
        convert protein positions to cds positions "pro2cds" or
        vice versa "cds2pro"
        
EOF

    exit(1);
}

$|++;


print STDERR "# $0 v$VERSION\n# " . (scalar localtime()) . "\n";

if ($output) {
    open (STDOUT, "> $output") || die $!;
}

## read into all sequences
my %input_seqs = ();
print STDERR ">> Start reading sequences from $fasta_file ... ";
parse_fasta_SEQs(\%input_seqs, $fasta_file);
print STDERR "done!\n";


##
## query sequences
##
print STDERR ">> Start querying $query_file ... ";
extract_split_seqs($query_file);
print STDERR "done!\n";

print STDERR "# " . (scalar localtime()) . "\n";




######################### Sub #########################


=head2 extract_seqs

    About   : Query sequence data.
    Usage   : extract_seqs($fasta_file);
    Args    : Source fasta file.
    Returns : Null

=cut
sub extract_split_seqs
{
    my ($in) = @_;
    
    my $fh = getInputFilehandle($in);
    while (<$fh>)
    {
        next if (/\#/ || /^\s+$/);
    
        my ($query_id, $ranges) = (split /\s+/, $_)[0,1];
        
        if ($input_seqs{$query_id}) {
            my $seq  =  $input_seqs{$query_id};
               $seq  =~ s/\*$//;
            
            my @ranges = (split /\;/, $ranges);
               @ranges = sort {(split /\-/, $a)[0] <=> (split /\-/, $b)[0]} @ranges;
            
            for (my $i=0; $i+1<@ranges; $i++)
            {
                my ($prev_start, $prev_end) = (split /\-/, $ranges[$i]);
                my ($next_start, $next_end) = (split /\-/, $ranges[$i+1]);
                
                if ( $prev_end >= $next_start ) {
                    $ranges[$i] = "$prev_start\-$next_end";
                    splice(@ranges, $i+1, 1);
                }
            }
            
            my @out_ranges = ();
            
            my @seqs = ();
            for my $range (@ranges)
            {
                my ($start, $end) = (split /\-/, $range);
                
                if ($convert_pos eq 'cds2pro') {
                    $start = ($start+2)/3;
                    $end   = $end/3;
                }
                elsif ($convert_pos eq 'pro2cds') {
                    $start = $start*3-2;
                    $end   = $end*3;
                }
                
                push @out_ranges, "$start-$end";
                
                my $sub_seq = substr($seq, $start-1, $end-$start+1);
                push @seqs, $sub_seq;
            }
            
            unless (@seqs > 0) {
                print STDERR "$query_id\t$ranges\n";exit;
            }
            
            my $joined_seq   = join '', @seqs;
            my $out_ranges = join ';',  @out_ranges;
            
            print STDOUT ">$query_id $out_ranges\n$joined_seq\n";
        }
        else {
            print STDERR "no record found: $query_id\n";
        }
    }
}


