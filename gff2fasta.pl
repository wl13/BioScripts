#!/usr/bin/perl -w
#
#   gff2fasta.pl -- Extract sequences from gff3 file to fasta file
#                          
#
#   Author: Nowind
#   Created: 2012-05-31
#   Updated: 2016-05-06
#   Version: 1.1.1
#
#   Change logs:
#   Version 1.0.0 13/11/13: The initial version.
#   Version 1.1.0 15/02/12: Skip loci where chromosome sequence is not given; add
#                           and change several options.
#   Version 1.1.1 16/05/06: Bug fixed: use Parent id instead of ID field in CDS.






=head1 NAME

gff2fasta.pl


=head1 SYNOPSIS

  gff2fasta.pl --help/?

=head1 DESCRIPTION

Extract sequences from gff3 file to fasta file

=cut




use strict;

use Data::Dumper;
use Getopt::Long;
use File::Find::Rule;
use File::Basename;

use MyPerl::FileIO qw(:all);
use MyPerl::Convert;
use MyPerl::Compare;

######################### Main #########################

my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.1.1';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";
my $SOURCE  = (scalar localtime()) . " Version: $VERSION";


my $out_features = "mRNA";
my $out_format   = "fasta";
my ($gff_file, $output, $ref_seq, $out_details, $word_wrap, $show_help);
GetOptions(
            "i|gff=s"            => \$gff_file,
            "s|seqs=s"           => \$ref_seq,
            
            "output=s"           => \$output,

            "F|format=s"         => \$out_format,
            
            "R|features=s"       => \$out_features,
            
            "wordwrap=i"         => \$word_wrap, 
            
            "details"            => \$out_details,
            
            "help|?"             => \$show_help,
           );

unless( !$show_help && $gff_file && $ref_seq ) {
    print <<EOF;

$0  -- Extract sequences from gff3 file to fasta file.

Version: $VERSION

Usage:   perl $0 [options]

Options:
    -i, --gff     <filename>
        annotation file in gff3 format, required
    -s, --seqs    <filename>
        reference sequences in fasta format, required
        
    -o, --output  <dirname>
        output filename, default to STDOUT
    
    -F, --format  <string>
        output file format, can be set to fasta or tabular,
        default: fasta
    
    -R, --features
        output features, mRNA or cds, default: mRNA

    -w, --wordwrap  <int>
        line feed for print
    
    -d, --details
        output verbose info in header line
    
    -?, --help
        show this help message
    
EOF

    exit(1);
}

$|++;


print STDERR "# $0 v$VERSION\n# " . (scalar localtime()) . "\n";


if ($output) {
    open (STDOUT, "> $output") || die $!;
}

$out_features = lc($out_features);

print STDERR ">> Start reading $ref_seq ... ";
my %SEQs = ();
my @ids  = parse_fasta_SEQs(\%SEQs, $ref_seq);
print STDERR "done!\n";




print STDERR ">> Start converting $gff_file ... ";
convert_gff2seqs($gff_file);
print STDERR "done!\n";

print STDERR "# " . (scalar localtime()) . "\n";

######################### Sub #########################



=head2 convert_gff2seqs

    About   : Convert annotation file from  gff format to table-delimited format.
    Usage   : convert_gff2tables($gff_file, \%representative_trans);
    Args    : Annotation file in gff3 format.
    Returns : Null
=cut
sub convert_gff2seqs
{
    my ($in) = @_;
    
    my %gff_features = ();
    
    my $fh = getInputFilehandle($in);    
    while (<$fh>)
    {
        chomp;
        my ($chrom, $source, $feature, $start, $end,
            $score, $strand, $frame, $attribute) = (split /\t/);
        
        next unless($SEQs{$chrom});
        
        my $ID = '';
        
        if ($feature eq 'CDS' && $attribute =~ /Parent\=(.*?)(,|;|$)/) {  ## CDS
            $ID = $1;
        }
        elsif ($attribute =~ /^ID\=(.*?)(,|;|$)/) {                       ## mRNA
            $ID = $1;
        }
        
        unless( $ID ) {
            print STDERR "Error: No ID info found for line $.: $_\n"; exit;
        }
        
        if ($feature eq 'mRNA') {
            $gff_features{gene}->{$chrom}->{$ID} = [$start, $end, $strand];
        }
        elsif ($feature eq 'CDS') {
            push @{$gff_features{cds}->{$chrom}->{$ID}}, [$start, $end, $strand];
        }
    }

    if ($out_features eq 'mrna') {
        for my $chrom (sort keys %{$gff_features{gene}})
        {
            my @ids = sort { $gff_features{gene}->{$chrom}->{$a}->[0] <=>
                             $gff_features{gene}->{$chrom}->{$b}->[0] }
                      keys %{$gff_features{gene}->{$chrom}};
            
            for my $id (@ids)
            {
                my ($start, $end, $strand) = @{$gff_features{gene}->{$chrom}->{$id}};
                
                my $strand_info = ($strand eq "+") ? "FORWARD" : "REVERSE";
                my $seq_length  = $end - $start + 1;
                
                my $out_id = $out_details ? "$id\t$chrom($start..$end)\t$strand_info\tLENGTH=$seq_length" : $id;
                
                my $seq = substr($SEQs{$chrom}, $start-1, $end-$start+1);
                
                if ($strand eq '-') {
                    $seq =~ tr/ATGC/TACG/;
                    $seq = reverse $seq;
                }
                
                if ($out_format eq 'tabular') {
                    print STDOUT "$out_id\t$seq\n";
                }
                else {
                    print STDOUT format_fasta_SEQs($out_id, \$seq, $word_wrap);
                }
                
            }
        }
    }
    elsif ($out_features eq 'cds') {
        for my $chrom (sort keys %{$gff_features{cds}})
        {
            my %cds_seqs = ();
            for my $id (keys %{$gff_features{cds}->{$chrom}})
            {
                my @cds_parts = @{$gff_features{cds}->{$chrom}->{$id}};
                
                if ($cds_parts[0]->[-1] eq '+') {
                    @cds_parts = sort { $a->[0] <=> $b->[0] } @cds_parts;
                    
                    $cds_seqs{$id}->{start}  = $cds_parts[0]->[0];
                    $cds_seqs{$id}->{strand} = '+';
                    
                    for my $part (@cds_parts)
                    {
                        my ($start, $end, $strand) = @{$part};
                        my $seq = substr($SEQs{$chrom}, $start-1, $end-$start+1);
                        
                        push @{$cds_seqs{$id}->{seq}}, $seq;
                        push @{$cds_seqs{$id}->{pos}}, "$start..$end";
                    }
                }
                elsif ($cds_parts[0]->[-1] eq '-') {
                    @cds_parts = sort { $b->[0] <=> $a->[0] } @cds_parts;
                    
                    $cds_seqs{$id}->{start}  = $cds_parts[-1]->[0];
                    $cds_seqs{$id}->{strand} = '-';
                    
                    for my $part (@cds_parts)
                    {
                        my ($start, $end, $strand) = @{$part};
                        
                        my $seq = substr($SEQs{$chrom}, $start-1, $end-$start+1);
                           $seq =~ tr/ATGC/TACG/;
                           $seq = reverse $seq;
                        
                        push @{$cds_seqs{$id}->{seq}}, $seq;
                        push @{$cds_seqs{$id}->{pos}}, "$start..$end";
                    }
                }
            }
            
            my @cds_ids = sort { $cds_seqs{$a}->{start} <=>
                                 $cds_seqs{$b}->{start} } keys %cds_seqs;
            
            for my $id (@cds_ids)
            {
                my $cds_seq = join '',  @{$cds_seqs{$id}->{seq}};
                my $cds_pos = join ',', @{$cds_seqs{$id}->{pos}};
                my $strand  = $cds_seqs{$id}->{strand};
                
                my $strand_info = ($strand eq "+") ? "FORWARD" : "REVERSE";
                my $seq_length  = length($cds_seq);
                
                my $out_id = $out_details ? "$id\t$chrom($cds_pos)\t$strand_info\tLENGTH=$seq_length" : $id;
                
                if ($out_format eq 'tabular') {
                    print STDOUT "$out_id\t$cds_seq\n";
                }
                else {
                    print STDOUT format_fasta_SEQs($out_id, \$cds_seq, $word_wrap);
                }
            }
        }
    }
}



