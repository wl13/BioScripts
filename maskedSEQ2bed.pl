#!/usr/bin/perl -w
#
#   maskedSEQ2bed.pl -- sequence related analysis
#
#   Author: Nowind
#   Created: 2012-02-21
#   Updated: 2014-06-04
#   Version: 1.0.1
#
#   Change logs:
#   Version 1.0.0 12/05/24: The initial version.
#   Version 1.0.1 14/06/04: Add some comments.

use strict;

use Data::Dumper;


##################### Main ####################

my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.0.1';
my $HEADER  = "# $CMDLINE\n# Version: $VERSION\n";


unless( @ARGV == 1 ) {
    print <<EOF;

$0 -- sequence related analysis

Usage:
perl $0 <fasta file>

EOF

    exit(0);
}

$|++;

print STDERR "# $0 v$VERSION\n# " . (scalar localtime()) . "\n";

print STDERR ">> Read in $ARGV[0]...";
my %SEQs = ();
parseFasta(\%SEQs, $ARGV[0]);
print STDERR "done!\n";

#print "$HEADER# " . (scalar localtime()) . "\n";

## BED -- BED type, a general purpose format for representing   
## genomic interval data, useful for masks and other interval 
## outputs. Please note that the bed format is 0-based (most 
## other formats are 1-based). 

print STDERR ">> Start analysis sequences...";
for my $id (sort keys %SEQs)
{
    my $seq   = $SEQs{$id};
    
    #my @nts = split //, $seq;
    #
    #for my $nt (@nts)
    #{
    #    
    #}
    
    ## the matching process is really slow and a new solution is urged
    while ( $seq =~ m/(N|[a-z])+/g ) 
    {
        my $end   = (pos $seq);     ## BED ends are one-based
        my $len   = length $&;     
        my $start = $end - $len;    ## BED starts are zero-based
        
        print STDOUT "$id\t$start\t$end\n";
    }
}
print STDERR "done!\n";


print STDERR "# " . (scalar localtime()) . "\n";

####################### Sub ###########################

#**************************************************
# Save all the sequences in the fasta
# file into a hash or an array
#**************************************************
sub parseFasta
{
    my ($rf_seq, $in, $desc) = @_;
    
    my $match_str = '^(.*?)\s+';
    
    if ($desc) {
        if ($desc eq 'full') {
            $match_str = '^(.*?)\n' 
        }
        elsif ($desc eq 'gbk') {
            $match_str = 'gb\|(.*?)\.';
        }        
    }

    open (F, "< $in") or die $!;               
    my $fas = do { local $/; <F> };
    close F;
    
    my @fas = split /\>/, $fas;
    
    shift @fas;
    
    for my $str (@fas)
    {
        $str =~ /$match_str/;
        my $id = $1;
        
        $str =~ s/.*\n//;
        #$str =~ s/\>//;
        $str =~ s/\s+//g;
        
        if (ref($rf_seq) eq 'HASH') {
            $rf_seq->{$id} = $str;
        }
        elsif (ref($rf_seq) eq 'ARRAY') {
            push @{$rf_seq}, $str;
        }
        
    }
}
