#!/usr/bin/perl -w
#
#  Statistics.pm -- do some statistic jobs
#
#  Author: Nowind
#  Created: 2010-10-09
#  Updated: 2015-11-16
#  Version: 1.1.0
#
#  Change logs:
#  Version 1.0.0 12/10/18: The initial version.
#  Version 1.1.0 15/11/16: Updated: Add function chi_squared_test.




use strict;
use Carp qw< croak >;
use List::Util qw< sum >;
use Statistics::Distributions qw< chisqrprob >;


package MyPerl::Statistics;


sub new
{
    my ($pkg, %data)  = @_;
    
    bless {
        "numbers"       => $data{numbers},
        "observed"      => $data{observed},
        "expected"      => $data{expected},
        "correct"       => $data{yates_correction},
    }, $pkg;
}


sub count
{
    my $obj = shift;
    
    return (scalar @{$obj->{numbers}});
}

sub sum
{
    my $obj = shift;
    
    my $sum = 0;
    for my $num (@{$obj->{numbers}})
    {
        $sum += $num;
    }
    
    return $sum;
}

sub mean
{
    my $obj = shift;
    
    my $num = scalar @{$obj->{numbers}};
    my $sum = $obj->sum();
    
    my $mean = $sum / $num;
}


sub sqsum
{
    my $obj = shift;
    
    my $mean = $obj->mean();
    
    my $sqsum = 0;
    for my $num (@{$obj->{numbers}})
    {
        $sqsum += ( $num - $mean ) ** 2;
    }

    return $sqsum;
}

## sample standard deviation
sub stdev_s
{
    my $obj = shift;
    
    my $num = scalar @{$obj->{numbers}};
    
    return 0 unless( $num > 1 );
    
    my $sqsum = $obj->sqsum();
    
    $sqsum /= ( $num - 1 );      
    
    my $stdev = sqrt($sqsum);
    
    return $stdev;
}

sub stdev_p
{
    my $obj = shift;
    
    my $num = scalar @{$obj->{numbers}};
    
    my $sqsum = $obj->sqsum();
    
    $sqsum /= ( $num );      
    
    my $stdev = sqrt($sqsum);
    
    return $stdev;
}


=head2 chi_squared_test

    About   : Pearson's chi-square test.
    Source  : http://stackoverflow.com/questions/21204733/a-better-chi-square-test-for-perl
    Usage   : chi_squared_test("observed"   => \@observed,
                               "expected"   => \@expected,
                               "correct"    => $yates_correction,);
    Args    : Observed counts;
              Expected counts;
              Determine wether Yates' correction should be used.
    Returns : P-values.

=cut
sub chi_squared_test
{
    my $obj = shift;

    return -2 unless(@{$obj->{observed}} == @{$obj->{expected}});

    my $chi_squared = sum map { 
        ( $obj->{observed}->[$_] - $obj->{expected}->[$_] )**2 / $obj->{expected}->[$_]
    } 0 .. $#{$obj->{observed}};
    
    ## Yates's correction for continuity
    if ($obj->{correct}) {
        $chi_squared = sum map { 
            ( abs($obj->{observed}->[$_] - $obj->{expected}->[$_]) - 0.5 )**2 / $obj->{expected}->[$_]
        } 0 .. $#{$obj->{observed}};
    }
    
    my $degrees_of_freedom = @{$obj->{observed}} - 1;
    
    my $probability        = chisqrprob(
        $degrees_of_freedom,
        $chi_squared
    );

    return $probability;
}


1;