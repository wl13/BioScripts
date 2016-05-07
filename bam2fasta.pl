#!/usr/bin/perl -w
#
#   bam2fasta.pl -- Convert aligned sequences from a SAM format file
#
#
#   Author: Nowind
#   Created: 2012-02-21
#   Updated: 2016-02-25
#   Version: 1.0.0
#
#   Change logs:
#   Version 1.0.0 16/02/25: The initial version.




use strict;

use Data::Dumper;
use Getopt::Long;

use MyPerl::FileIO qw(:all);

##################### Main ####################


my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.0.0';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";

my $SOURCE  = (scalar localtime()) . " Version: $VERSION";

my $out_format  = 'fasta';
my (@bam_files, $no_rc, $output, $extend_size, $use_rg_id, $samtools_opts,
    $min_seq_len, $max_clipping, $min_insert_size, $max_insert_size);
GetOptions(
            "bam=s{,}"         => \@bam_files,
            
            "output=s"         => \$output,
            "format=s"         => \$out_format,
            "no-rc"            => \$no_rc,
            
            "samtools=s"       => \$samtools_opts,
            "min-insert=i"     => \$min_insert_size,
            "max-insert=i"     => \$max_insert_size,
            "min-len=i"        => \$min_seq_len,
            "max-clipping=i"   => \$max_clipping,
            "use-rg"           => \$use_rg_id,
           );

my $show_help = ($CMDLINE =~ /\-help/) ? 0 : 1;

unless( (@bam_files >= 1) && $show_help ) {
    print <<EOF;

$0  -- Extract all reads with the 

Version: $VERSION

Usage:   perl $0 [options]

Input Options:

    -b, --bam    <filename>
        bam file(s), required

Output Options:

    -o, --output <filename>
        output filename if output in sam format, or output filename
        prefix if output if fastq format
    -f, --format <string>
        output format, default output in sam format, can be set to
        fastq format
        
    -n, --no-rc
        do not reverse complement sequence with negtive strand
    -u, --use-rg
        add read group id to extracted records
        
Filtering Options:
    -s, --samtools <string>
        directly pass samtools view options to this script, e.g.
        "-f 4 -F 8"

    --min-len      <int>
        minium sequence length

    --max-clipping <int>
        maximum allowed clipping length, include both soft and hard
        clipping bases
    
    --min-insert   <int>
    --max-insert   <int>
        screen out records with insert size wihtin this range

EOF

    exit(1);
}




$|++;



print STDERR "# $0 v$VERSION\n# " . (scalar localtime()) . "\n";

if ($output) {
    open (STDOUT, "> $output") || die $!;
}


bam2fasta();


print STDERR "# " . (scalar localtime()) . "\n";


######################### Sub #########################


=head2 bam2fasta

    About   : Convert SAM records to fasta/fastq sequences
    Usage   : bam2fasta();
    Args    : Null
    Returns : Null

=cut
sub bam2fasta
{
    for (my $i=0; $i<@bam_files; $i++)
    {
        my $pipe_str = "samtools view $bam_files[$i] |";
        
        if ($samtools_opts) {
            $pipe_str = "samtools view $samtools_opts $bam_files[$i] |";
        }
        
        open (my $fh, $pipe_str) || die $!;
        while (<$fh>)
        {
            next if (/^@/ || /^\s+$/); ## skip header
            
            chomp(my $record = $_);
            
            my ($QNAME, $FLAG, $RNAME, $POS, $MAPQ, $CIGAR, 
                $MRNM, $NPOS, $TLEN, $SEQ, $QUAL, @OPT) = (split /\s+/, $record);
            
            
            ##
            ## check if this record should be update
            ##
            my $rg_id = ($record =~ /RG:Z:(.*?)\s+/);
                
            $QNAME =~ s/\/\d$//;
            
            my $read_id = $use_rg_id ? "$rg_id:$QNAME" : $QNAME;
            
            next if ($min_seq_len && length($SEQ) < $min_seq_len);
            
            if (defined $max_clipping) {
                my $soft_clipped = ($CIGAR =~ /(\d+)S/) ? $1 : 0;
                my $hard_clipped = ($CIGAR =~ /(\d+)H/) ? $1 : 0;
                
                next if ($soft_clipped + $hard_clipped > $max_clipping);
            }
            
            
            if (($FLAG & 16) && !$no_rc) {      ## reverse strand
                $SEQ =~ tr/ATGCatgc/TACGtacg/;
                $SEQ = reverse $SEQ;
                
                $QUAL = reverse $QUAL;
            }
            
            my $pair_id = 0;
            if ($FLAG & 64) {      ## first in pair
                $pair_id = 1;
            }
            elsif ($FLAG & 128) {  ## second in pair
                $pair_id = 2;
            }
            
            if ($out_format eq 'fastq') {
                print "\@$read_id/$pair_id\n$SEQ\n\+\n$QUAL\n";
            }
            else {
                print ">$read_id/$pair_id\n$SEQ\n";
            }
            
        }
    }
}


