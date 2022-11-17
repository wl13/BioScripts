#!/usr/bin/perl -w
#
#   downsample_vcf.pl -- generate certain number of loci by random.
#
#   Author: Nowind
#   Created: 2011-10-19
#   Updated: 2014-03-11
#   Version: 1.0.1
#
#   Change logs:
#   Version 1.0.0 14/03/05: The initial version.
#   Version 1.0.1 14/03/11: Bug fixed in generate suffix of file names; add option
#                           "--pos-only".





use strict;

use Data::Dumper;
use Getopt::Long;
use File::Find::Rule;
use Statistics::Descriptive;
use Statistics::PointEstimation;
use Data::Random qw(:all);

use MyPerl::FileIO qw(:all);
use MyPerl::Statistics;

################### Main #################

my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.0.1';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";
my $SOURCE  = (scalar localtime()) . " Version: $VERSION";

my %opts = ();
   $opts{times}   = 1;
   $opts{number}  = 100;
my ($out_pos_only);
GetOptions(
            "input=s"          => \$opts{input},
            "output=s"         => \$opts{prefix},
            "times=i"          => \$opts{times},
            "size=i"           => \$opts{size},
            
            "pos-only"         => \$out_pos_only,
           );

unless( $opts{input} ) {
    print <<EOF;

$0  -- random select SNP markers and count directions 

Version: $VERSION

Usage:   perl $0 [options]

Options:
    -i, --input  <filename>
        file contains chromosome lengths
    -o, --output <string>
        output file prefix, default: rand[001.vcf ...]

    -t, --times
        random times [default: 1]
    -s, --size
        numbers of loci [default: 100]

    -p, --pos-only
        only output chromosome and positions
        
EOF

    exit(0);
}

$|++;

print STDERR "# $0 v$VERSION\n# " . (scalar localtime()) . "\n";

unless ($opts{prefix}) {
    $opts{prefix} = 'rand';
}

print STDERR ">> Start parsing $opts{input} ... ";
my $vcf_header  = '';
my $sample_line = '';
my @records_all = ();
parse_vcf($opts{input});
print STDERR "\tdone!\n";


gen_random_vars($opts{times}, $opts{size});


print STDERR "# " . (scalar localtime()) . "\n";

######################### Sub #########################

sub parse_vcf
{
    my ($in) = shift;
    
    my $fh = getInputFilehandle($in);
    while (<$fh>)
    {
        if (/\#\#/) {
            $vcf_header .= $_;
        }
        elsif (/^\#CHROM/) {
            if ($out_pos_only) {
                $sample_line = "#CHROM\tPOS\n";
            }
            else {
                $sample_line = $_;
            } 
        }
        
        next if (/\#/ || /^\s+$/);
        
        if ($out_pos_only) {
            my ($chrom, $pos) = (split /\s+/, $_)[0,1];
            push @records_all, "$chrom\t$pos";
        }
        else {
            chomp;
            push @records_all, $_;
        }
        
    }
}


sub gen_random_vars
{
    my ($rand_times, $rand_size) = @_;
    
    my %rand_values = ();
    
    my $i = 0;
    my $tag = '0' x (length($rand_times));
    while (++$i <= $rand_times)
    {
        print STDERR "\r>> Start random process ... duplicate $i\/$rand_times";
        
        open (my $fh, "> $opts{prefix}.$tag.vcf") || die $!;
        
        print {$fh} "$vcf_header";
        print {$fh} "##source=$SOURCE $CMDLINE\n";
        print {$fh} "$sample_line";
        
        my @rand_vars = rand_set( set => [@records_all], size => $rand_size, shuffle => 0 );
        
        print {$fh} (join "\n", @rand_vars);
        
        $tag++;
    }
    print STDERR "\tdone!\n";
}



