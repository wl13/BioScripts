#!/usr/bin/perl -w
#
#   fgenesh2gff.pl -- Convert results from fgenesh to gff3 format.
#
#   Author: Nowind
#   Created: 2011-09-18
#   Updated: 2016-05-06
#   Version: 2.0.0
#
#   Change logs:
#   Version 1.0.0 11/09/18: The initial version.
#   Version 2.0.0 15/05/06: Rewrite all codes.



=head1 NAME

process_fgenesh.pl


=head1 SYNOPSIS

  process_fgenesh.pl --help/?

=head1 DESCRIPTION

Convert results from fgenesh to other format.

=cut


use strict;

use Data::Dumper;
use Getopt::Long;
use File::Find::Rule;

use MyPerl::FileIO qw(:all);

################## Main ##################

my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '2.0.0';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";
my $SOURCE  = (scalar localtime()) . " Version: $VERSION";

my ($input, $output);
GetOptions(
            "input=s"    => \$input,
            "output=s"   => \$output,
           );

unless( $input ) {
    print <<EOF;

$0  -- Convert results from fgenesh to gff3 format.

      *Note: if the original results contain sequences with no reliable
             predictions, please remove those records or replace those
             lines with " no reliable predictions " with "//"

Version: $VERSION

Usage:   perl $0 [options]
         
         genesh2gff.pl -i fgenesh.txt -o fgenesh.gff
         
         sed 's/ no reliable predictions /\\/\\//' fgenesh.txt | \
            genesh2gff.pl -i - -o fgenesh.gff
         
Options:
    --input  <filename>
        input file generate from fgenesh, required
    --output <filename>
        output file, default to STDOUT

EOF

    exit(1);
}

$|++;


print STDERR "# $0 v$VERSION\n# " . (scalar localtime()) . "\n";

if ($output) {
    open (STDOUT, "> $output") || die $!;
}

#print STDOUT "##$CMDLINE\n##$SOURCE\n";
#print STDOUT "#ID\tExon_Num\tTSS_Start\tPolA_Start\tCDS_Len\tGene_Len\n";
print STDOUT "##gff-version 3\n";

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
    my ($info, @seq) = (split /\>/, $str);
    
    my @infos  = split /\n/, $info;
    
    my %genes    = ();
    for (my $i=0; $i<@infos; $i++)
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
        
        if ($infos[$i] =~ /(PolA|TSS)/) {
            my (undef, $G, $Str, $Feature, $Start, $Score) = (split /\s+/, $infos[$i]);
            ###$genes{$G}->{Strand}                    = $Str;
            ###$genes{$G}->{$Feature}->{Start}         = $Start;
            ###$genes{$G}->{$Feature}->{End}           = '.';
            ###$genes{$G}->{$Feature}->{Score}         = $Score;
            ###
            ###print "$contig_id\tFgenesh\t$Feature\t$Start\t.\t$Score\t$Str\t.\t" .
            ###      "ID=$contig_id\.fgenesh$G;Parent=$contig_id\.fgenesh$G\n";
            push @{$genes{$G}->{$Feature}}, "$contig_id\tFgenesh\t$Feature\t$Start\t.\t$Score\t$Str\t.\t" .
                                            "ID=$contig_id\.fgenesh$G;Parent=$contig_id\.fgenesh$G";
        }
        elsif ($infos[$i] =~ /CDS/) {
            my ($G, $Str, $order, $Feature, $Start, $End, $Score) = (split /\s+/, $infos[$i])[1..5,7,8];
            ###$genes{$G}->{Strand}       = $Str;
            ###$genes{$G}->{CDS}->{Start} = $Start;
            ###$genes{$G}->{CDS}->{End}   = $End;
            ###$genes{$G}->{CDS}->{Score} = $Score;
            ###
            ###print "$contig_id\tFgenesh\tCDS\t$Start\t$End\t$Score\t$Str\t.\t" .
            ###      "ID=$contig_id\.fgenesh$G:cds\_$order;Parent=$contig_id\.fgenesh$G\n";
            
            push @{$genes{$G}->{CDS}}, "$contig_id\tFgenesh\tCDS\t$Start\t$End\t$Score\t$Str\t.\t" .
                                       "ID=$contig_id\.fgenesh$G:cds\_$order;Parent=$contig_id\.fgenesh$G";
        }
    }
    
    ## parse mRNA
    my @mRNA   = grep {$_ =~ /\[mRNA\]/} @seq;
    for (my $i=0; $i<@mRNA; $i++)
    {
        $mRNA[$i] =~ m{
                       (\d+)\s+                # gene id
                       (\d+)\s+exon\s+\(s\)    # exon numbers
                       \s+(\d+)\s+\-\s+(\d+)   # positions
                       \s+(\d+)\s+bp\,\s+      # length
                       chain\s+(\+|\-)            # strand
                      }gx;
        
        my ($mRNA_id, $exon_num, $start, $end, $len, $strand) = ($1, $2, $3, $4, $5, $6);
        
        my $mRNA_feature   = "$contig_id\tFgenesh\tmRNA\t$start\t$end\t.\t$strand\t.\t" .
                           "ID=$contig_id\.fgenesh$mRNA_id;Name=$contig_id\.fgenesh$mRNA_id;Parent=$contig_id\.fgenesh$mRNA_id";
    
        my $cds_features   = join "\n", @{$genes{$mRNA_id}->{CDS}};
        my $TSS_features   = join "\n", @{$genes{$mRNA_id}->{TSS}}  if $genes{$mRNA_id}->{TSS};
        my $PolA_features  = join "\n", @{$genes{$mRNA_id}->{PolA}} if $genes{$mRNA_id}->{PolA};
        
        if ($strand eq '+' && $genes{$mRNA_id}->{TSS}) {
            print STDOUT "$TSS_features\n";
        }
        
        if ($strand eq '-' && $genes{$mRNA_id}->{PolA}) {
            print STDOUT "$PolA_features\n";
        }
        
        print STDOUT "$mRNA_feature\n$cds_features\n";
        
        if ($strand eq '-' && $genes{$mRNA_id}->{TSS}) {
            print STDOUT "$TSS_features\n";
        }
        
        if ($strand eq '+' && $genes{$mRNA_id}->{PolA}) {
            print STDOUT "$PolA_features\n";
        }
    }
}

