#!/usr/bin/perl -w
#
#   extract_sam_seqs.pl -- extract aligned sequences and related reference sequence
#                          from a SAM format file
#
#   Author: Nowind
#   Created: 2012-02-21
#   Updated: 2013-01-03
#   Version: 1.0.1
#
#   Change logs:
#   Version 1.0.0 12/12/31: The initial version.
#   Version 1.0.1 13/01/03: Change output sequence id.





use strict;

use Data::Dumper;
use Getopt::Long;

use MyPerl::FileIO qw(:all);

##################### Main ####################


my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.0.1';
my $HEADER  = "# $CMDLINE\n# Version: $VERSION\n";

my $SOURCE  = (scalar localtime()) . " Version: $VERSION";

my ($input, $ref_file, $output);
GetOptions(
            "input=s"          => \$input,
            "refer=s"          => \$ref_file,
            "output=s"         => \$output,
           );

my $show_help = ($CMDLINE =~ /\-help/) ? 0 : 1;

unless( $input && $ref_file && $show_help ) {
    print <<EOF;

$0  -- extract aligned sequences and related reference sequence from a SAM format file

Version: $VERSION

Usage:   perl $0 [options]

Options:
    -i, --input   <filename>
        input file in SAM format, use "-" to indicate STDIN,
        required
    -r, --refer   <filename>
        reference file, in fasta format, required
    -o, --output  <filename>
        output filename, output extracted sequences in fasta
        format, default to STDOUT

EOF

    exit(1);
}




$|++;

if ($output) {
    open (STDOUT, "> $output") || die $!; 
}

 
print STDERR "# $0 v$VERSION\n# " . (scalar localtime()) . "\n";


print STDERR ">> Read in reference sequences $ref_file ... ";
my %SEQs = ();
parse_fasta_SEQs(\%SEQs, $ref_file);
print STDERR "done!\n";

print STDERR ">> Start extract sequences from $input ... ";
extract_sam_seqs(\%SEQs, $input);
print STDERR "done!\n";

print STDERR "# " . (scalar localtime()) . "\n";


######################### Sub #########################


=head2 extract_seqs

    About   : Extract sequences from SAM format file
    Usage   : extract_seqs(\%Reference_SEQs, $sam_file);
    Args    : Hash of reference sequences
              File in SAM format
    Returns : Null

=cut
sub extract_sam_seqs
{
    my ($rh_ref_seqs, $in) = @_;
    
    my @names = ();
    my %SEQs  = ();
    
    my $fh = getInputFilehandle($in);
    while (<$fh>)
    {
        next if (/^@/ || /^\s+$/); ## skip header
        
        my ($QNAME, $FLGA, $RNAME, $POS, $MAPQ, $CIGAR, 
            $MRNM, $NPOS, $TLEN, $SEQ, $QUAL, @OPT) = (split /\s+/);
        
        push @{$SEQs{$RNAME}->{$POS}}, "$QNAME,$SEQ";
    }
    
    for my $chrom (sort keys %SEQs)
    {
        my @POS = sort {$a <=> $b} (keys %{$SEQs{$chrom}});
        
        my $max_start = $POS[0];
        my $max_end   = $POS[-1] + (length $SEQs{$chrom}->{$POS[-1]}->[0]) - 1;
        my $max_range = $max_end - $max_start + 1;
        
        my $ref_seq = substr($rh_ref_seqs->{$chrom}, $max_start-1, $max_range);
        my $ref_id  = "$chrom $max_start-$max_end|$max_range";
        
        print format_fasta_SEQs($ref_id, \$ref_seq);
        
        for my $start (@POS)
        {
            for my $record (@{$SEQs{$chrom}->{$start}})
            {
                my ($id, $seq) = (split /\,/, $record);
                
                my $len = length $seq;
                my $end = $start + $len - 1;
                
                print format_fasta_SEQs("$id $start-$end|$len", \$seq);                
            }

        }
    }
    
    
}


















