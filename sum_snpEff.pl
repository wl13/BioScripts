#!/usr/bin/perl -w
#
#   sum_snpEff.pl -- Reduce snpEff genes.txt by sum up relevant rows.
#
#   Author: Nowind
#   Created: 2011-09-18
#   Updated: 2017-05-27
#   Version: 1.0.1
#
#   Change logs:
#   Version 1.0.0 17/04/21: The initial version.
#   Version 1.0.1 17/05/27: Update: add support for snpEff v4.3o.



=head1 NAME

sum_snpEff.pl


=head1 SYNOPSIS

  sum_snpEff.pl --help/?

=head1 DESCRIPTION

Reduce snpEff genes.txt by sum up relevant rows.

=cut


use strict;

use Data::Dumper;
use Getopt::Long;
use File::Find::Rule;

use MyPerl::FileIO qw(:all);

################## Main ##################

my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.0.1';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";
my $SOURCE  = (scalar localtime()) . " Version: $VERSION";

my ($input, $output);
GetOptions(
            "input=s"           => \$input,
            "output=s"          => \$output,
           );

unless( $input ) {
    print <<EOF;

$0  -- Reduce snpEff genes.txt by sum up relevant rows.

Version: $VERSION

Usage:   perl $0 [options]
         
         sum_snpEff.pl -i genes.txt -o genes.sum.txt
         
Options:
    -i,--input     <filename>
        input genes.txt file generate by snpEff, required
        
    -o,--output    <filename>
        output file, default to STDOUT


EOF

    exit(1);
}

$|++;


print STDERR "# $0 v$VERSION\n# " . (scalar localtime()) . "\n";

if ($output) {
    open (STDOUT, "> $output") || die $!;
}

print STDOUT "$HEADER##" . (scalar localtime()) . "\n";

print STDERR "Start parsing $input ... ";
process_snpEff($input);
print STDERR "done!\n";



################# Sub #################

sub process_snpEff
{
    my ($in)     = shift;
    my $fh       = getInputFilehandle($in);
    my @effs     = ();
    my @out_effs = qw(frameshift stop_gain stop_lost start_lost
                      stop_var start_var inframe_var missense_var synonym_var
                      splice_var utr_var intron_var updown_var others);
    while(<$fh>)
    {
        if (/^\#Gene/) {
            my $out_effs = join "\t", @out_effs;
            print STDOUT "#GeneId\tGeneName\tBioType\t$out_effs\n";
            
            @effs = (split /\t/);
        }
        
        next if(/\#/ || /^\s+$/);
        
        my @fields = (split /\s+/);
        my %counts = ();
        
        for (my $i=0; $i<@fields; $i++)
        {
            next unless($effs[$i] =~ /Count/ || $effs[$i] =~ /variants_effect/);
            
            if ($effs[$i] =~ /frameshift/ || $effs[$i] =~ /exon_loss/) {
                $counts{frameshift}   += $fields[$i];
            }
            elsif ($effs[$i] =~ /stop_gain/) {
                $counts{stop_gain}    += $fields[$i];
            }
            elsif ($effs[$i] =~ /stop_lost/) {
                $counts{stop_lost}    += $fields[$i];
            }
            elsif ($effs[$i] =~ /start_lost/) {
                $counts{start_lost}   += $fields[$i];
            }
            elsif ($effs[$i] =~ /stop_retained/) {
                $counts{stop_var}     += $fields[$i];
            }
            elsif ($effs[$i] =~ /initiator_codon/) {
                $counts{start_var}    += $fields[$i];
            }
            elsif ($effs[$i] =~ /inframe/) {
                $counts{inframe_var}  += $fields[$i];
            }
            elsif ($effs[$i] =~ /missense/) {
                $counts{missense_var} += $fields[$i];
            }
            elsif ($effs[$i] =~ /synonymous/) {
                $counts{synonym_var}  += $fields[$i];
            }
            elsif ($effs[$i] =~ /splice/) {
                $counts{splice_var}   += $fields[$i];
            }
            elsif ($effs[$i] =~ /UTR/) {
                $counts{utr_var}      += $fields[$i];
            }
            elsif ($effs[$i] =~ /intron/) {
                $counts{intron_var}   += $fields[$i];
            }
            elsif ($effs[$i] =~ /upstream/ || $effs[$i] =~ /downstream/) {
                $counts{updown_var}   += $fields[$i];
            }
            else {
                $counts{others}       += $fields[$i];
            }
        }
        
        my @out_counts = ();
        for my $eff (@out_effs)
        {
            my $n = $counts{$eff} ? $counts{$eff} : 0;
            
            push @out_counts, $n;
        }
        
        my $out_counts = join "\t", @out_counts;
        
        print STDOUT "$fields[0]\t$fields[1]\t$fields[2]\t$out_counts\n";
    }
}

