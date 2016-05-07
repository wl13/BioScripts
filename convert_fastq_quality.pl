#!/usr/bin/perl -w
#
#   convert_fastq_qaulity.pl -- Convert Illumina quality (+64) to Sanger quality (+33) or reverse.
#
#   Author: Nowind
#   Created: 2012-02-21
#   Updated: 2016-03-13
#   Version: 1.0.0
#
#   Change logs:
#   Version 1.0.0 16/03/13: The initial version.



use strict;

use Data::Dumper;
use Getopt::Long;

use MyPerl::FileIO qw(:all);

##################### Main ####################


my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.0.0';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";

my $SOURCE  = (scalar localtime()) . " Version: $VERSION";

my $convert_type = 'illumina2sanger';
my ($fastq_file, $output);
GetOptions(
            "input=s"          => \$fastq_file,

            "output=s"         => \$output,
            
            "convert-type"     => \$convert_type,
           );

my $show_help = ($CMDLINE =~ /\-help/) ? 0 : 1;

unless( $fastq_file && $show_help ) {
    print <<EOF;

$0  -- Extract all reads with the 

Version: $VERSION

Usage:   perl $0 [options]

Options:

    -i, --input  <filename>
        fastq file, default from STDIN, required

    -o, --output <filename>
        output file, default to STDOUT

    -c, --convert-type
        illumina2sanger: Convert Illumina quality (+64) to Sanger quality (+33)
        sanger2illumina: Convert Sanger quality (+33) to Illumina quality (+64)
        default: illumina2sanger
    
EOF

    exit(1);
}




$|++;



print STDERR "# $0 v$VERSION\n# " . (scalar localtime()) . "\n";




print STDERR ">> Start converting $fastq_file ... ";
convert_quality($fastq_file);
print STDERR "done!\n";



print STDERR "# " . (scalar localtime()) . "\n";


######################### Sub #########################


=head2 illumina2sanger

    About   : Convert Illumina quality (+64) to Sanger quality (+33) or reverse.
    Usage   : convert_quality($fastq_file)
    Args    : Fastq file.
    Returns : Null
    Reference: http://wiki.bioinformatics.ucdavis.edu/index.php/IllQ2SanQ.pl
               https://github.com/jstjohn/KentLib/tree/master/examples/qseqToFastq

=cut
sub convert_quality
{
    my ($in) = @_;
    
    my $fh = getInputFilehandle($in);
    while (<$fh>)
    {
        chomp(my $seq_id = $_);
        chomp(my $seq    = <$fh>);
        chomp(my $fill   = <$fh>);
        chomp(my $qual   = <$fh>);
        
        my @prev_quals = split(//,$qual);
        my @curr_quals = ();
        
        for (my $i=0; $i<=$#prev_quals; $i++) {
            if ($convert_type eq 'illumina2sanger') {
                ## Convert a phred64 string into a phred64 string assuming illumina's minimum 'B' score thing
                my $char = $prev_quals[$i] eq 'B' ? '!' : chr(ord($prev_quals[$i]) - 31);
                push @curr_quals, $char;
            }
            else {
                ## Convert a phred33 string into a phred64 string, doing illumina's minimum 'B' score thing
                my $char = $prev_quals[$i] eq '!' ? 'B' : chr(ord($prev_quals[$i]) + 31);
                push @curr_quals, $char;
            }
        }
        
        my $curr_qual = join '', @curr_quals;
        
        if (length($qual) != length($curr_qual)) {
            print STDERR "$seq_id\n$seq\n$fill\n$qual\n$curr_qual\n";
            exit;
        }
        
        print STDOUT "$seq_id\n$seq\n$fill\n$curr_qual\n";
    }
}

