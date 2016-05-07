#!/usr/bin/perl -w
#
#   sam2fastq.pl -- Convert sam format to fastq format.
#
#   Author: Nowind
#   Created: 2012-02-21
#   Updated: 2016-03-02
#   Version: 1.1.0
#
#   Change logs:
#   Version 1.0.0 14/01/21: The initial version.
#   Version 1.1.0 16/03/02: Updated: remove duplicated records with same id.



use strict;

use Data::Dumper;
use Getopt::Long;

use MyPerl::FileIO qw(:all);

##################### Main ####################


my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.1.0';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";

my $SOURCE  = (scalar localtime()) . " Version: $VERSION";

my ($input, $output, $no_rc, $use_rg_id);
GetOptions(
            "I|input=s"        => \$input,

            "output=s"         => \$output,
            
            "no-rc"            => \$no_rc,
            "use-rg"           => \$use_rg_id,
           );

my $show_help = ($CMDLINE =~ /\-help/) ? 0 : 1;

unless( $input && $show_help ) {
    print <<EOF;

$0  -- Extract all reads with the 

Version: $VERSION

Usage:   perl $0 [options]

Options:

    -i, --input  <filename>
        input file of query positions, required

    -o, --output <filename>
        output prefix, required

    -n, --no-rc
        do not reverse complement sequence with negtive strand
    -u, --use-rg
        add read group id to extracted records
    
EOF

    exit(1);
}




$|++;



print STDERR "# $0 v$VERSION\n# " . (scalar localtime()) . "\n";




print STDERR ">> Start converting $input ... ";
sam2fastq($input, $output);
print STDERR "done!\n";



print STDERR "# " . (scalar localtime()) . "\n";


######################### Sub #########################


=head2 sam2fastq

    About   : Convert sam file to fastq files.
    Usage   : sam2fastq($sam_file, $out_prefix);
    Args    : File in SAM format;
              Prefix of output files.
    Returns : Null

=cut
sub sam2fastq
{
    my ($in, $out_prefix) = @_;

    my %read_ids = ();

    open (FQ1, "> $out_prefix" . "_1.fq") || die $!;
    open (FQ2, "> $out_prefix" . "_2.fq") || die $!;
    
    my $fh = getInputFilehandle($input);
    while (<$fh>)
    {
        next if (/^@/ || /^\s+$/); ## skip header
        
        chomp(my $record = $_);
        
        my ($QNAME, $FLAG, $RNAME, $POS, $MAPQ, $CIGAR, 
            $MRNM, $NPOS, $TLEN, $SEQ, $QUAL, @OPT) = (split /\s+/, $record);

        my $rg_id = ($record =~ /RG:Z:(.*?)\s+/);
            
        $QNAME =~ s/\/\d$//;
        
        my $read_id = $use_rg_id ? "$rg_id:$QNAME" : $QNAME;
        
        ## re-reverse bases and qualities of reads with negative strand
        ## flag set
        if (($FLAG & 16) && !$no_rc) {
            $SEQ =~ tr/ATGCatgc/TACGtacg/;
            $SEQ = reverse $SEQ;
            
            $QUAL = reverse $QUAL;
        }
        
        my $pair_id = 0;
        if ($FLAG & 64) {      ## first in pair
            $pair_id = 1;
            
            next if (exists($read_ids{"$read_id:$pair_id"}));
            
            print FQ1 "\@$read_id/1\n$SEQ\n\+\n$QUAL\n";
        }
        elsif ($FLAG & 128) {  ## second in pair
            $pair_id = 2;
            
            next if (exists($read_ids{"$read_id:$pair_id"}));
            
            print FQ2 "\@$read_id/2\n$SEQ\n\+\n$QUAL\n";
        }
        
        $read_ids{"$read_id:$pair_id"}++;
    }

    close FQ1;
    close FQ2;
}


