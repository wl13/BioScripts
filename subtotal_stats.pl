#!/usr/bin/perl -w
#
#   subtotal_stats.pl -- Subtotal stats by different rows.
#                          
#
#   Author: Nowind
#   Created: 2012-05-31
#   Updated: 2017-05-27
#   Version: 2.0.4
#
#   Change logs:
#   Version 1.0.0 13/05/07: The initial version.
#   Version 1.0.1 13/05/09: Add option "--parts" to count values located in different
#                           intervals.
#   Version 1.0.2 13/05/13: Add option "--keys" to specify output keys.
#   Version 1.0.3 13/06/13: Add test for non-numeric values.
#   Version 1.0.4 14/06/10: Change "--parts" to "--percent" to use "percentile" function
#                           instead of "frequency_distribution_ref" function.
#   Version 1.1.0 14/07/02: Add function to stat lines with multiple values.
#   Version 1.2.0 14/12/18: Add function to calculate standard errors.
#   Version 1.2.1 14/12/19: Change input rows parsed in function count_multi2.
#   Version 2.0.0 14/12/23: Remove redundant functions, rewrite and rearrange most codes.
#   Version 2.0.1 15/10/20: Remove unused option "--rows"; bug fixed while processing multiple
#                           values; update explanation of some options.
#   Version 2.0.2 15/10/23: Bug fixed while no percentile value returned.
#   Version 2.0.3 15/11/12: Add median values in output results.
#   Version 2.0.4 17/05/27: Update: Add option "--skip-minus" to skip unwanted minus values.


use strict;

use Data::Dumper;
use Getopt::Long;
use File::Find::Rule;
use File::Basename;
use Statistics::Descriptive;
use Statistics::PointEstimation;
use Scalar::Util qw(looks_like_number);

use MyPerl::FileIO qw(:all);

######################## Main ########################


my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '2.0.4';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";


my %options = ();
my ($output, $multi_values, @percentiles, $skip_minus);
GetOptions(
            "input=s"               => \$options{input_file},
            
            "O|output=s"            => \$output,
            
            "D|out-orders=s{,}"     => \@{$options{out_orders}},
            
            "key-type=s"            => \$options{key_type},
            "value-type=s"          => \$options{value_type},
            
            "percent=f{,}"          => \@percentiles,
            
            "multi-values"          => \$multi_values,
            
            "skip-minus"            => \$skip_minus,
           );

unless( $options{input_file} ) {
    print <<EOF;

$0  -- Subtotal stats by different rows.

Version: $VERSION

Usage:   perl $0 [options]

Options:
    -i, --input   <filename>
        input file with first row used as key, remain rows as values, there
        could be only one key, but multiple values can be specified by using
        "--multi-values" option, required
        
    -O, --output  <filename>
        output filename, default to STDOUT
    
    --key-type   <string>
        manully specify the type of keys, "numeric" or "string", otherwise
        will be determined by automatic detection
    --value-type <string>
        manully specify the type of values, "numeric" or "string", otherwise
        will be determined by automatic detection, note when value type is
        string, only occurence will be counted
    
    -p, --percent <float>
        sort the data and returns the value that corresponds to those
        percentiles specified here, can have multiple values
    
    -m, --multi-values
        subtotal stats of lines with one master key with multiple values, the
        input file should start with an header line with content look like:
        
        #key    value_id1   value_id2   value_id3   ...
        ex1     1000        1000        1000        ...
        ex1     1500        1010        500         ...
        ex1     1300        2000        3000        ...
        ex2     1050        1000        1000        ...
        ...
    
        then the actual key field would be "value_id1:ex1", "value_id2:ex1",
        "value_id1:ex2", etc.

        
    -D, --out-orders <strings>
        specify keys manually, the output will be sorted according to the
        order specified here, this option is used to control the output
        order of key fields (or master key fields while "--multi-values" option
        is specified)
    
    -s, --skip-minus
        skip minus numbers in values
    
EOF

    exit(1);
}

$|++;

if ($output) {
    open (STDOUT, "> $output") || die $!;
}

 
print STDERR "# $0 v$VERSION\n# " . (scalar localtime()) . "\n";

print STDOUT "$HEADER##" . (scalar localtime()) . "\n";


count_subtotal(\%options);



print STDERR "# " . (scalar localtime()) . "\n";

######################### Sub #########################



=head2 count_subtotal

    About   : Count subtotals.
    Usage   : count_subtotal($file);
    Args    : Array reference to hold background blocks infos;
              File contains blocks infos.
    Returns : Null
    
=cut
sub count_subtotal
{
    my ($opts) = @_;
    
     
    printf STDERR ">> Start parsing $opts->{input_file} ... ";
    my @titles     = ();
    my %Stats      = ();
    my $key_type   = $options{key_type}   ? $options{key_type}   : 'string';
    my $value_type = $options{value_type} ? $options{value_type} : 'string';
    my $fh = getInputFilehandle($opts->{input_file});
    while (<$fh>)
    {
        if (/^\#/) {
            @titles = (split /\s+/);
        }
        
        next if (/^\#/ || /^\s+$/);
        
        my ($key, @values) = (split /\s+/);
        
        ##
        ## check type of key and values
        ##
        if(!$options{key_type} && looks_like_number($key)) {
            $key_type = 'numeric';
        }
        
        if(!$options{value_type} && looks_like_number($values[0])) {
            $value_type = 'numeric';
        }
        
        for (my $i=0; $i<@values; $i++)
        {
            next if ($skip_minus && $values[$i] < 0);
            
            if ($multi_values) {
                push @{$Stats{$key}->{$i}}, $values[$i];
            }
            else {
                push @{$Stats{$key}->{0}}, $values[$i];
            }
        }
    }
    print STDERR "done!\n";
    
    
    ##
    ## sort output by key values
    ##
    my @user_keys = ();
    if (@{$opts->{out_orders}} > 0) {
        @user_keys = @{$opts->{out_orders}};
    }
    else {
        if ($key_type eq 'string') {
            @user_keys = sort {$a cmp $b} keys %Stats;
        }
        else {
            @user_keys = sort {$a <=> $b} keys %Stats;
        }
    }
    
    
    ##
    ## generate stats
    ##
    if ($value_type eq 'numeric') {
        if (@percentiles > 0) {
            my $parts = join "\t", @percentiles;
            print "#ID\tCount\tSum\tMean\tMedian\tMin\tMax\tStdev\tStd_err\tPercentile:$parts\n";
        }
        else {
            print "#ID\tCount\tSum\tMean\tMedian\tMin\tMax\tStdev\tStd_err\n";
        }
    }
    else {
        print "#ID\tCount(All)\tCount(NoDup)\n";
    }
    
    
    print STDERR ">> Start generating results ... ";
    for my $sub_key (@user_keys)
    {
        unless($Stats{$sub_key}) {
            next;
        }
        
        for my $i (sort {$a <=> $b} keys %{$Stats{$sub_key}})
        {
            my $out_id    = $sub_key;
            my $out_stats = '';
            
            if ($multi_values) {
                my $column_id = (@titles > 0) ? $titles[$i+1] : ($i+1);
                
                $out_id = "$column_id:$sub_key";
            }
            
            if (($options{value_type} && $options{value_type} eq 'numeric') ||
               (!$options{value_type} && looks_like_number($Stats{$sub_key}->{$i}->[0]))) {
                my $stat = Statistics::Descriptive::Full->new();
                   $stat->add_data(@{$Stats{$sub_key}->{$i}});
                
                my $count  = $stat->count();
                my $sum    = $stat->sum();
                my $mean   = $stat->mean();
                my $median = $stat->median();
                my $min    = $stat->min();
                my $max    = $stat->max();
                my $var    = $stat->variance();
                my $stdev  = $stat->standard_deviation();
                my $stderr = $stdev / sqrt($count);
                
                if (@percentiles > 0) {
                    my @perc_counts = ();
                    for (@percentiles)
                    {
                        my $perc_count = defined($stat->percentile($_)) ? $stat->percentile($_) : '-';
                        push @perc_counts, $perc_count;
                    }
                    
                    my $perc_counts = join "\t", @perc_counts;
                    
                    $out_stats = "$count\t$sum\t$mean\t$median\t$min\t$max\t$stdev\t$stderr\t$perc_counts";
                }
                else {
                    $out_stats = "$count\t$sum\t$mean\t$median\t$min\t$max\t$stdev\t$stderr";
                }
            }
            else {
                ##
                ## only count frequency of each string
                ##
                my %dups = ();
                for my $j (@{$Stats{$sub_key}->{$i}})
                {
                    next if (exists $dups{$j});
                    $dups{$j}++;
                }
                
                my $cnt_all   = scalar @{$Stats{$sub_key}->{$i}};
                my $cnt_nodup = scalar (keys %dups);
                
                $out_stats = "$cnt_all\t$cnt_nodup";
            }
            
            print STDOUT "$out_id\t$out_stats\n";
        }
    }
    print STDERR "done!\n";
}
