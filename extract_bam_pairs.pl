#!/usr/bin/perl -w
#
#   extract_bam_pairs.pl -- extract aligned sequences and related reference sequence
#                          from a SAM format file
#
#   Author: Nowind
#   Created: 2012-02-21
#   Updated: 2015-01-21
#   Version: 1.1.2
#
#   Change logs:
#   Version 1.0.0 14/02/23: The initial version.
#   Version 1.0.1 14/02/25: Add support for multiple bam files.
#   Version 1.0.2 14/03/02: Add option "--extend".
#   Version 1.1.0 14/05/29: Change input file format; add option "--rows"; add options
#                           for filtering bam records.
#   Version 1.1.1 14/06/01: Add option "--use-rg" to add readgroup id to bam records
#                           instead of bam file index.
#   Version 1.1.2 15/01/21: Add option "--replace-file" to replace extracted sam records
#                           which matches those in the replace file.



use strict;

use Data::Dumper;
use Getopt::Long;

use MyPerl::FileIO qw(:all);

##################### Main ####################


my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.1.2';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";

my $SOURCE  = (scalar localtime()) . " Version: $VERSION";

my $out_format  = 'sam';
my ($input, @bam_files, $no_rc, @rows, $output, $patch_file,
    $extend_size, $use_rg_id, $samtools_opts,
    $min_seq_len, $max_clipping, $min_insert_size, $max_insert_size);
GetOptions(
            "I|input=s"        => \$input,
            "R|rows=i{,}"      => \@rows,
            
            "extend=i"         => \$extend_size,
            
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
            
            "patches=s"        => \$patch_file,
           );

my $show_help = ($CMDLINE =~ /\-help/) ? 0 : 1;

unless( $input && (@bam_files > 0) && $show_help ) {
    print <<EOF;

$0  -- Extract all reads with the 

Version: $VERSION

Usage:   perl $0 [options]

Input Options:

    -I, --input  <filename>
        input file of query positions, required
        
    -R, --rows   <numbers>
        specify the row fields of chromosome, start position and end
        position (0-based), in the query block file [default: 0 1 2]
    
    -p, --patches <filename>
        update sam records with same QNAME and FLAG found in this file
    
    -e, --extend <int>
        extend regions of this size to retrieve for read pairs

    -b, --bam    <filename>
        bam file(s), at least one bam file should be specified


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


unless(@rows > 0){ @rows = qw(0 1 2) };

if ($out_format eq 'fastq') {
    open (FQ1, "> $output" . "_1.fq") || die $!;
    open (FQ2, "> $output" . "_2.fq") || die $!;
}
else {
    if ($output) {
        open (STDOUT, "> $output") || die $!;
    }
}



my %new_sam_records = ();
if ($patch_file) {
    print STDERR ">> Patch data found, start parsing records in $patch_file ... ";
    my $patch_fh = getInputFilehandle($patch_file);
    while (<$patch_fh>)
    {
        next if (/^@/ || /^\s+$/); ## skip header
        
        chomp;
        my ($QNAME, $FLAG) = (split /\s+/)[0,1];
        
        $new_sam_records{"$QNAME\t$FLAG"} = $_;
    }
    print STDERR "done!\n";
    
    ###print STDERR Dumper(%new_sam_records);exit;
}



print STDERR ">> Read in query ids in $input ... ";
my %query_ids  = ();
my %counts_all = ();
my $fh = getInputFilehandle($input);
while (<$fh>)
{
    next if (/^#/ || /^\s+$/); ## skip header
    
    my ($chrom, $start, $end) = (split /\s+/)[@rows];
    
    my %reads = ();
    
    search_read_id(\%reads, $chrom, $start, $end);
    extract_pairs(\%reads, $chrom, $start, $end);
    
    
    for my $id (sort keys %{$reads{pairs}})
    {
        $counts_all{reads} ++;
        
        next unless(keys %{$reads{pairs}->{$id}} == 2);
        
        $counts_all{paired} ++;
        
        if ($out_format eq 'fastq') {
            print FQ1 "\@" . "$id/1\n" . $reads{pairs}->{$id}->{1} . "\n";
            print FQ2 "\@" . "$id/2\n" . $reads{pairs}->{$id}->{2} . "\n";
        }
        else {
            print $reads{pairs}->{$id}->{1} . "\n";
            print $reads{pairs}->{$id}->{2} . "\n";
        }
    }
}
print STDERR "done!\n";

my $find_perc = ($counts_all{reads} > 0) ? ($counts_all{paired} / $counts_all{reads} * 100) : 0;
   $find_perc = sprintf("%.2f", $find_perc);
  
print STDERR <<EOF;
# $counts_all{paired} pairs were found during traversal out of $counts_all{reads} total pairs ($find_perc%)
EOF

if ($out_format eq 'fastq') {
    close FQ1;
    close FQ2;
}


print STDERR "# " . (scalar localtime()) . "\n";


######################### Sub #########################


=head2 search_bam_records

    About   : Extract sequences from SAM format file
    Usage   : extract_seqs(\%Reference_SEQs, $sam_file);
    Args    : Hash of reference sequences
              File in SAM format
    Returns : Null

=cut
sub search_read_id
{
    my ($rh_reads, $chrom, $start, $end) = @_;
    
    for (my $i=0; $i<@bam_files; $i++)
    {
        my $pipe_str = "samtools view $bam_files[$i] $chrom:$start-$end |";
        
        if ($samtools_opts) {
            $pipe_str = "samtools view $samtools_opts $bam_files[$i] $chrom:$start-$end |";
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
            if ($new_sam_records{"$QNAME\t$FLAG"}) {
                $record = $new_sam_records{"$QNAME\t$FLAG"};
                
               ($QNAME, $FLAG, $RNAME, $POS, $MAPQ, $CIGAR, 
                $MRNM, $NPOS, $TLEN, $SEQ, $QUAL, @OPT) = (split /\s+/, $record);
            }
            
            
            my $rg_id = ($record =~ /RG:Z:(.*?)\s+/);
            
            next if ($min_seq_len && length($SEQ) < $min_seq_len);
            
            next if ($min_insert_size && abs($TLEN) < $min_insert_size);
            next if ($max_insert_size && abs($TLEN) > $max_insert_size);
            
            if (defined $max_clipping) {
                my $soft_clipped = ($CIGAR =~ /(\d+)S/) ? $1 : 0;
                my $hard_clipped = ($CIGAR =~ /(\d+)H/) ? $1 : 0;
                
                next if ($soft_clipped + $hard_clipped > $max_clipping);
            }
            
            $QNAME =~ s/\/\d$//;
            
            my $read_id = $use_rg_id ? "$rg_id:$QNAME" : $QNAME;
            
            $rh_reads->{id}->{$read_id} = 1;
        }
    }
    

}

=head2 search_bam_records

    About   : Extract sequences from SAM format file
    Usage   : extract_seqs(\%Reference_SEQs, $sam_file);
    Args    : Hash of reference sequences
              File in SAM format
    Returns : Null

=cut
sub extract_pairs
{
    my ($rh_reads, $chrom, $start, $end) = @_;
    
    if ($extend_size) {
        $start -= $extend_size;
        $end   += $extend_size;
        
        $start = 0 if $start < 0;
    }
    
    for (my $i=0; $i<@bam_files; $i++)
    {
        my $pipe_str = "samtools view $bam_files[$i] $chrom:$start-$end |";
        
        if ($samtools_opts) {
            $pipe_str = "samtools view $samtools_opts $bam_files[$i] $chrom:$start-$end |";
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
            if ($new_sam_records{"$QNAME\t$FLAG"}) {
                $record = $new_sam_records{"$QNAME\t$FLAG"};
                
               ($QNAME, $FLAG, $RNAME, $POS, $MAPQ, $CIGAR, 
                $MRNM, $NPOS, $TLEN, $SEQ, $QUAL, @OPT) = (split /\s+/, $record);
            }
            
            
            my $rg_id = ($record =~ /RG:Z:(.*?)\s+/);
                
            $QNAME =~ s/\/\d$//;
            
            my $read_id = $use_rg_id ? "$rg_id:$QNAME" : $QNAME;
            
            next unless($rh_reads->{id}->{$read_id});
            
            next if ($min_seq_len && length($SEQ) < $min_seq_len);
            
            if (defined $max_clipping) {
                my $soft_clipped = ($CIGAR =~ /(\d+)S/) ? $1 : 0;
                my $hard_clipped = ($CIGAR =~ /(\d+)H/) ? $1 : 0;
                
                next if ($soft_clipped + $hard_clipped > $max_clipping);
            }
            
            
            if ($rh_reads->{id}->{$read_id}) {
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
                    $rh_reads->{pairs}->{$read_id}->{$pair_id} = "$SEQ\n\+\n$QUAL";
                }
                else {
                    $rh_reads->{pairs}->{$read_id}->{$pair_id} = $use_rg_id ? "$rg_id:$record" : $record;
                }
            }
        }
    }
}


