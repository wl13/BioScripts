#!/usr/bin/perl -w
#
#   gff2tables.pl -- Convert gff3 file to table-delimited file
#                          
#
#   Author: Nowind
#   Created: 2012-05-31
#   Updated: 2017-05-22
#   Version: 1.0.5
#
#   Change logs:
#   Version 1.0.0 13/04/14: The initial version.
#   Version 1.0.1 13/05/15: Add locus id and strand info in output results.
#   Version 1.0.2 13/06/28: Add options "--all-features", "--up-length" and "--down-length" to
#                           output all features including intergenic regions.
#   Version 1.0.3 13/06/29: Add function annotations to output; bug fixed: parent ids not assigned
#                           for each mRNA record while "--remove-alt" is omited.
#   Version 1.0.4 15/12/14: Bug fixed: skip header lines; Updated: output promoter regions.
#   Version 1.0.5 17/05/22: Bug fixed: skip features without ID present.







=head1 NAME

gff2tables.pl


=head1 SYNOPSIS

  gff2tables.pl --help/?

=head1 DESCRIPTION

Convert gff3 file to table-delimited file

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
my $VERSION = '1.0.5';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";
my $SOURCE  = (scalar localtime()) . " Version: $VERSION";


my $promoter_len = 500;
my ($gff_file, $output, $remove_alt, $out_all_features,
    $no_functions, @flank_lens,  $show_help);
GetOptions(
            "i|gff=s"            => \$gff_file,

            "output=s"           => \$output,

            "remove-alt"         => \$remove_alt,
            
            "all-features"       => \$out_all_features,
            
            "flank-length=i{,}"  => \@flank_lens,
            
            "promoter=i"         => \$promoter_len,
            
            "no-funcs"           => \$no_functions,
            
            "help|?"             => \$show_help,
           );

unless( !$show_help && $gff_file ) {
    print <<EOF;

$0  -- Convert gff3 file to table-delimited file.

Version: $VERSION

Usage:   perl $0 [options]

Options:
    -i, --gff     <filename>
        annotation file in gff3 format, required

    -o, --output  <dirname>
        output filename, default to STDOUT
    
    -r, --remove-alt
        remove alternative splicings, only remain the longest one as the
        representative transcript
    
    -a, --all-features
        output all features include intergenic regions, upsteam and downstream
        regions of genes
    
    -p, --promoter <number>
        regions within specified regions upstream of transcription start sites
        will be annotated as promoter regions, default: 500bp
        
    -f, --flank-length <numbers>
        length(s) of flanking regions, default: 2000bp
        this options valid only while -a option is specified
    
    -n, --no-funcs
        do not output function annotations
    
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

unless( @flank_lens > 0 ) {
    @flank_lens = (5000);
}


my %Locus_ids            = ();
my %representative_trans = ();
print STDERR ">> Start checking alternative splicings ... ";
get_represent_trans($gff_file);
print STDERR "done!\n";


print STDERR ">> Start converting $gff_file ... ";
convert_gff2tables($gff_file);
print STDERR "done!\n";

print STDERR "# " . (scalar localtime()) . "\n";

######################### Sub #########################



=head2 get_represent_trans

    About   : Get the longest transcript as the representative transcript while 
              alternative splicing present.
    Usage   : get_represent_trans($gff_file);
    Args    : Annotation file in gff format.
    Returns : Null

=cut
sub get_represent_trans
{
    my ($in) = @_;
    
    my %mRNAs_all = ();
    
    my $fh = getInputFilehandle($in);    
    while (<$fh>)
    {
        next if (/^\#/ || /^\s+$/);
        
        chomp;
        my ($chrom, $source, $feature, $start, $end,
            $score, $strand, $frame, $attribute) = (split /\t/);
        
        
        if ( $feature =~ /(mRNA|pseudogenic\_transcript)/ ) { ## mRNA
            m{
                ^(.*?)\s+.*?                   # Chromosome ID
                \s+(\d+)                       # Start Position
                \s+(\d+)\s+.*                  # End Position
                \s+(\-|\+)\s+\.                # Strand
                \s+ID\=(.*?)\;.*               # ID
                (Parent|Locus_id)\=(.*?)(;|$)  # Parent
            }x;
            
            my ($ID, $Parent) = ($5, $7);
            
            $Locus_ids{$ID} = $Parent;
            
            unless( $ID ) {
                #print STDERR "Error: No ID info found for line $.: $_";
                next;
            }
            
            ###unless( $Parent ) {
            ###    print STDERR "Error: No Parent ID info found for line $.: $_"; exit;
            ###}
            
            my $length = $end - $start + 1;
            
            if ($mRNAs_all{$Parent}->{represent}) {
                if ($mRNAs_all{$Parent}->{length} < $length) {
                    $mRNAs_all{$Parent}->{represent} = $ID;
                    $mRNAs_all{$Parent}->{length}    = $length;
                }
            }
            else {
                $mRNAs_all{$Parent}->{represent} = $ID;
                $mRNAs_all{$Parent}->{length}    = $length;
            }
        }
    }
    
    for my $mRNA (keys %mRNAs_all)
    {
        my $represent_id = $mRNAs_all{$mRNA}->{represent};
        
        $representative_trans{$represent_id} = $mRNA;
    }
}



=head2 convert_gff2tables

    About   : Convert annotation file from  gff format to table-delimited format.
    Usage   : convert_gff2tables($gff_file, \%representative_trans);
    Args    : Annotation file in gff3 format.
    Returns : Null
    Note    : This function can be easily replaced by using awk, e.g.
              awk -F'\t' 'BEGIN{OFS="\t"} {print $1,$3,$4,$5,$9;}' $gff_file | \
                  sed 's/Parent=//' | \
                  sed -r 's/ID=(.*);Name.*Locus_id=(\w+);.*/\1/g'

=cut
sub convert_gff2tables
{
    my ($in) = @_;

    print "$HEADER##" . (scalar localtime()) . "\n";
    print "#chrom\tfeature\tstart\tend\tID\tstrand\tLocus\n";

    my $fh = getInputFilehandle($in);    
    while (<$fh>)
    {
        next if (/^\#/ || /^\s+$/);
        
        chomp;
        my ($chrom, $source, $feature, $start, $end,
            $score, $strand, $frame, $attribute) = (split /\t/);
        
        next if ($feature eq 'gene');
        
        my $ID = '';
        
        if ($feature eq 'mRNA' && $attribute =~ /^ID\=(.*?)\;/) {       ## mRNA
            $ID = $1;
        }
        elsif ($attribute =~ /Parent\=(.*?)(;|$)/) {
            $ID = $1;
        }
        
        unless( $ID && $Locus_ids{$ID} ) {
            #print STDERR "Error: No ID info found for line $.: $_\n";
            next;
        }
        
        next if ($remove_alt && !$representative_trans{$ID});
        
        if ($feature eq 'mRNA') {
            if ($out_all_features) {
                for my $f_len (@flank_lens)
                {
                    my $promoter_start = $end + 1;
                    my $promoter_end   = $start - 1;
                    
                    my $aft_start = $end + 1;
                    my $bef_end   = $start - 1;
                    
                    
                    if ($strand eq '+') {
                        $promoter_start = $promoter_end - $promoter_len + 1;
                        $bef_end        = $promoter_start - 1;
                    }
                    else {
                        $promoter_end   = $promoter_start + $promoter_len - 1;
                        $aft_start = $promoter_end + 1;
                    }
                    
                    my $bef_start = $bef_end - $f_len + 1;
                    my $aft_end   = $aft_start + $f_len - 1;
                        
                        $bef_start = 1 if ($bef_start <= 0);
                    
                    
                    my $bef_flag  = ($strand eq '+') ? 'upstream'   : 'downstream';
                    my $aft_flag  = ($strand eq '+') ? 'downstream' : 'upstream';
                    
                    $bef_flag .= '_' . $f_len/1000 . 'k';
                    $aft_flag .= '_' . $f_len/1000 . 'k';
                    
                    print "$chrom\t$bef_flag\t$bef_start\t$bef_end\t$ID\t$strand\t$Locus_ids{$ID}\n";
                    print "$chrom\tpromoter_$promoter_len\t$promoter_start\t$promoter_end\t$ID\t$strand\t$Locus_ids{$ID}\n";
                    print "$chrom\t$aft_flag\t$aft_start\t$aft_end\t$ID\t$strand\t$Locus_ids{$ID}\n";
                }
            }
            
            my @attrs = (split /;/, $attribute);
            my @funcs = ();
            for my $attr (@attrs)
            {
                if ($attr =~ /(Note=|GO=|InterPro=|CGSNL)/){
                    push @funcs, $attr;
                }
            }
            
            if (@funcs > 0 && !$no_functions) {
                my $funcs = join ";", @funcs;
                print "$chrom\t$feature\t$start\t$end\t$ID\t$strand\t$Locus_ids{$ID}\t$funcs\n";
                next;
            }
        }
        
        print "$chrom\t$feature\t$start\t$end\t$ID\t$strand\t$Locus_ids{$ID}\n";
    }   
}



