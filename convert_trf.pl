#!/usr/bin/perl -w
#
#   convert_trf.pl -- Convert results from TRF (Tandem Repeats Finder) to other format.
#
#   Author: Nowind
#   Created: 2011-09-18
#   Updated: 2016-05-25
#   Version: 1.0.0
#   References:
#       http://tandem.bu.edu/trf/trf.html
#       https://gist.github.com/lexnederbragt/3689ee2301493c34c8ab#file-trf2gff-py
#
#   Change logs:
#   Version 1.0.0 16/05/25: The initial version.


=head1 NAME

convert_trf.pl


=head1 SYNOPSIS

  convert_trf.pl --help/?

=head1 DESCRIPTION

Convert results from TRF (Tandem Repeats Finder) to other format like gff3 or bed.

=cut


use strict;

use Data::Dumper;
use Getopt::Long;
use File::Find::Rule;

use MyPerl::FileIO qw(:all);

################## Main ##################

my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.0.0';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";
my $SOURCE  = (scalar localtime()) . " Version: $VERSION";

my $out_format = 'gff';
my ($input, $output);
GetOptions(
            "input=s"           => \$input,
            "output=s"          => \$output,
            "format=s"          => \$out_format,
           );

unless( $input ) {
    print <<EOF;

$0  -- Convert results from TRF to other format like gff3 or bed.

Version: $VERSION

Usage:   perl $0 [options]
         
         convert_trf.pl -i trf.dat -o trf.gff
         
         
Options:
    -i,--input     <filename>
        input dat file generate from TRF, required
    -o,--output    <filename>
        output file, default to STDOUT

    -f,--format    <string>
        output format, could be set as "gff", "bed" or "tabular"
        [default: gff]
EOF

    exit(1);
}

$|++;


print STDERR "# $0 v$VERSION\n# " . (scalar localtime()) . "\n";

if ($output) {
    open (STDOUT, "> $output") || die $!;
}

print STDERR "Start parsing $input ... ";
parse_trf($input);
print STDERR "done!\n";


################# Sub #################

sub parse_trf
{
    my ($in) = shift;
    
    if ($out_format eq 'gff') {
        print STDOUT "##gff-version 3\n";
    }
    elsif ($out_format eq 'tabular') {
        print STDOUT "$HEADER##" . (scalar localtime()) . "\n";
        print STDOUT "#source\tsequence_id\ttype\tstart\tend\tperiod\tcopies\tconsensus_size\tperc_match\t" .
                     "perc_indels\talign_score\tperc_A\tperc_C\tperc_G\tperc_T\tentropy\tcons_seq\trepeat_seq\n";
    }
    
    
    my $fh = getInputFilehandle($in);
    my $sequence_id = '';
    my %counts      = ();
    while(<$fh>)
    {
        next if (/\#/ || /^\s+$/);
        
        if (/^Sequence\:/) {
            (undef, $sequence_id) = (split /\s+/);
            $counts{$sequence_id} = 0;
        }
        elsif (/^\d+/) {
            unless($sequence_id) {
                print STDERR "Error: No sequence id found!"; exit(2);
            }
            
            $counts{$sequence_id} ++;
            
            my ($start, $end, $period, $copies, $consensus_size, $perc_match,
                $perc_indels, $align_score, $perc_A, $perc_C, $perc_G, $perc_T,
                $entropy, $cons_seq, $repeat_seq) = (split /\s+/);
            
            
            if ($out_format eq 'gff') {
                my $tr_id = $sequence_id . "." . $counts{$sequence_id};
                
                print "$sequence_id\tTRF\ttandem_repeat\t$start\t$end\t$align_score\t+\t.\t" .
                        "ID=$tr_id;Name=($cons_seq)$copies;Target=$repeat_seq\n";
            }
            elsif ($out_format eq 'bed') {
                my $bed_start = $start - 1;
                print "$sequence_id\t$bed_start\t$end\t$repeat_seq([$cons_seq]x$copies)\n";
            }
            elsif ($out_format eq 'tabular') {
                print "TRF\t$sequence_id\ttandem_repeat\t$start\t$end\t$period\t$copies\t$consensus_size\t$perc_match\t" .
                      "$perc_indels\t$align_score\t$perc_A\t$perc_C\t$perc_G\t$perc_T\t$entropy\t$cons_seq\t$repeat_seq\n";
            }
        }
    }
}

