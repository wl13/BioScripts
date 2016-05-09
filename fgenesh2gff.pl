#!/usr/bin/perl -w
#
#   fgenesh2gff.pl -- Convert results from fgenesh to gff3 or fasta format.
#
#   Author: Nowind
#   Created: 2011-09-18
#   Updated: 2016-05-09
#   Version: 2.1.0
#
#   Change logs:
#   Version 1.0.0 11/09/18: The initial version.
#   Version 2.0.0 15/05/06: Rewrite all codes.
#   Version 2.0.1 15/05/07: Bug fixed: failed to correctly split each field due to begining blanks.
#   Version 2.1.0 15/05/09: Bug fixed: try to choose the correct start and end for CDS;
#                           Update: add support for direct output sequences in fasta format.


=head1 NAME

fgenesh2gff.pl


=head1 SYNOPSIS

  fgenesh2gff.pl --help/?

=head1 DESCRIPTION

Convert results from fgenesh to gff3 or fasta format.

=cut


use strict;

use Data::Dumper;
use Getopt::Long;
use File::Find::Rule;

use MyPerl::FileIO qw(:all);

################## Main ##################

my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '2.1.0';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";
my $SOURCE  = (scalar localtime()) . " Version: $VERSION";

my ($input, $output, $write_seq, $word_wrap);
GetOptions(
            "input=s"           => \$input,
            "output=s"          => \$output,
            "S|write-seq=s"     => \$write_seq,
            "L|wordwrap=i"      => \$word_wrap,
           );

unless( $input ) {
    print <<EOF;

$0  -- Convert results from fgenesh to gff3 or fasta format.

      *Note: if the original results contain sequences with no reliable
             predictions, please remove those records or replace those
             lines with " no reliable predictions " with "//"

Version: $VERSION

Usage:   perl $0 [options]
         
         genesh2gff.pl -i fgenesh.txt -o fgenesh.gff
         
         sed 's/ no reliable predictions /\\/\\//' fgenesh.txt | \
            genesh2gff.pl -i - -o fgenesh.gff
         
Options:
    -i,--input     <filename>
        input file generate from fgenesh, required
    -o,--output    <filename>
        output file, gene ids will be renamed as "sequence_id.fgeneshN",
        where N stands for predicted gene number, default to STDOUT

    -S,--write-seq <string>
        output sequences in fgenesh results in fasta format rather than
        generate gff results, can be set to
        
            cds: whole coding sequences, sequence titile start
                 with FGENESH:[mRNA]
            exon: exon sequences, sequence titile start with
                 FGENESH:[exon]
            protein: protein sequences

    -L,--wordwrap  <int>
        line feed for print, only valid while output in fasta format

EOF

    exit(1);
}

$|++;


print STDERR "# $0 v$VERSION\n# " . (scalar localtime()) . "\n";

if ($output) {
    open (STDOUT, "> $output") || die $!;
}


if ($write_seq && $write_seq eq 'cds') {
    $write_seq = 'mRNA';
}

unless ($write_seq) {
    print STDOUT "##gff-version 3\n";
}

print STDERR "Start parsing $input ... ";
my $fh   = getInputFilehandle($input);
my @predicts = do{local $/ = '//'; <$fh>};
parse_fgenesh($_) for @predicts;
print STDERR "done!\n";



################# Sub #################

sub parse_fgenesh
{
    my ($str) = shift;
    
    ## get sequence id
    $str =~ /Seq\s+name\:\s+(.*?)\s+/;
    my $contig_id = $1;
    
    if ($str =~ / no reliable predictions /) {
        print STDERR "Error: Found sequences with no reliable predictions, please remove those "
                   . "records or replace those lines to \"//\" before converting!\n";
        exit(2);
    }
    
    return -1 unless $contig_id;
    
    ## parse gene structures infos and sequences
    my ($predicted_infos, @predicted_seqs) = (split /\>/, $str);
    
    ##
    ## parse predicted sequences
    ##
    my %predicted_seqs = ();
    for (my $i=0; $i<@predicted_seqs; $i++)
    {
        my ($seq_title, @seqs) = (split /\n/, $predicted_seqs[$i]);
        
        my $seq    = join '', @seqs;
           $seq    = uc($seq);
           $seq    =~ s/\///g;
        my $out_id = '';
        my $type   = '';
        
        if ($seq_title =~ /FGENESH:\[exon\]/) {  ## exon
            ## FGENESH:[exon] Gene:  1 Exon:  1  Pos:  3947  -   4027    81 bp., chain -
            my ($gene_id, $exon_id, $start, $end, $len, $strand) =
               (split /\s+/, $seq_title)[2,4,6,8,9,12];

            $predicted_seqs{$gene_id}->{exon}->{$exon_id}->{pos} = [$start, $end];
            $predicted_seqs{$gene_id}->{exon}->{$exon_id}->{seq} = $seq;
            
            $type   = 'exon';
            $out_id = "$contig_id\.fgenesh$gene_id:$exon_id";
        }
        else {                                    ## mRNA or protein
            ## FGENESH:[mRNA]   1  12 exon (s)   3947  -  15342  1269 bp, chain -
            ## FGENESH:   1  12 exon (s)   3947  -  15342   422 aa, chain -
            my ($gene_id, $exon_num, $start, $end, $len, $strand) =
               (split /\s+/, $seq_title)[1,2,5,7,8,11];
            
            $type = ($seq_title =~ /mRNA/) ? 'mRNA' : 'protein';
            
            $predicted_seqs{$gene_id}->{$type}->{pos}    = [$start, $end];
            $predicted_seqs{$gene_id}->{$type}->{strand} = $strand;
            $predicted_seqs{$gene_id}->{$type}->{seq}    = $seq;
            $predicted_seqs{$gene_id}->{$type}->{exon}   = $exon_num;
            
            $out_id = "$contig_id\.fgenesh$gene_id";
        }
        
        if ($write_seq && uc($write_seq) eq uc($type)) {
            print STDOUT format_fasta_SEQs($out_id, \$seq, $word_wrap);
        }
    }
    
    return 0 if $write_seq;
    
    ##
    ## parse predicted features
    ##
    my @predicted_infos    = split /\n/, $predicted_infos;
    my %predicted_features = ();
    for (my $i=0; $i<@predicted_infos; $i++)
    {
        ## http://linux1.softberry.com/berry.phtml?topic=fgenesh&group=help&subgroup=gfind
        ## Fgenesh output:
        ##  G               - predicted gene number, starting from start of sequence; 
        ##  Str             - DNA strand (+ for direct or - for complementary); 
        ##  Feature         - type of coding sequence:
        ##             CDSf - First (Starting with Start codon),
        ##             CDSi - internal (internal exon),
        ##             CDSl - last coding segment (ending with stop codon); 
        ##  TSS             - Position of transcription start (TATA-box position and score); 
        ##  Start and End   - Position of the Feature; 
        ##  Weight          - Log likelihood*10 score for the feature; 
        ##  ORF             - start/end positions where the first complete codon starts and the last codon ends.
        
        if ($predicted_infos[$i] =~ /(PolA|TSS)/) {
            $predicted_infos[$i] =~ s/^\s+//;
            my ($G, $Str, $Feature, $Start, $Score) = (split /\s+/, $predicted_infos[$i]);
            
            push @{$predicted_features{$G}->{$Feature}}, "$contig_id\tFgenesh\t$Feature\t$Start\t.\t$Score\t$Str\t.\t" .
                                                         "ID=$contig_id\.fgenesh$G;Parent=$contig_id\.fgenesh$G";
        }
        elsif ($predicted_infos[$i] =~ /CDS/) {
               $predicted_infos[$i] =~ s/^\s+//;
            
            my ($G, $Str, $exon_id, $Feature, $Start, $End, $Score, $ORF_Start, $ORF_End)
            = (split /\s+/, $predicted_infos[$i])[0..4,6..8,10];
            
            my $exon_num = $predicted_seqs{$G}->{mRNA}->{exon};
            
            my $cds_start = $Start;
            my $cds_end   = $End;
            
            ## correct the position of start codon
            if ($exon_id == 1) {
                $cds_start = $ORF_Start;
            }
            
            ## correct the position of stop codon
            if ($exon_id == $exon_num) {
                $cds_end   = $ORF_End;
            }
            
            push @{$predicted_features{$G}->{CDS}}, "$contig_id\tFgenesh\tCDS\t$cds_start\t$cds_end\t$Score\t$Str\t.\t" .
                                                    "ID=$contig_id\.fgenesh$G:cds\_$exon_id;Parent=$contig_id\.fgenesh$G";
        }
    }
    
    ##
    ## generate results
    ##
    my @mRNA_infos   = grep {$_ =~ /\[mRNA\]/} @predicted_seqs;
    for (my $i=0; $i<@mRNA_infos; $i++)
    {
        my ($gene_id, $exon_num, $start, $end, $len, $strand) =
           (split /\s+/, $mRNA_infos[$i])[1,2,5,7,8,11];
        
        my $mRNA_feature   = "$contig_id\tFgenesh\tmRNA\t$start\t$end\t.\t$strand\t.\t" .
                             "ID=$contig_id\.fgenesh$gene_id;Name=$contig_id\.fgenesh$gene_id;Parent=$contig_id\.fgenesh$gene_id";
    
        my $cds_features   = join "\n", @{$predicted_features{$gene_id}->{CDS}};
        my $TSS_features   = join "\n", @{$predicted_features{$gene_id}->{TSS}}  if $predicted_features{$gene_id}->{TSS};
        my $PolA_features  = join "\n", @{$predicted_features{$gene_id}->{PolA}} if $predicted_features{$gene_id}->{PolA};
        
        if ($strand eq '+' && $predicted_features{$gene_id}->{TSS}) {
            print STDOUT "$TSS_features\n";
        }
        
        if ($strand eq '-' && $predicted_features{$gene_id}->{PolA}) {
            print STDOUT "$PolA_features\n";
        }
        
        print STDOUT "$mRNA_feature\n$cds_features\n";
        
        if ($strand eq '-' && $predicted_features{$gene_id}->{TSS}) {
            print STDOUT "$TSS_features\n";
        }
        
        if ($strand eq '+' && $predicted_features{$gene_id}->{PolA}) {
            print STDOUT "$PolA_features\n";
        }
    }
}

