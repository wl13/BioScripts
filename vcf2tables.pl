#!/usr/bin/perl -w
#
#   vcf2tables.pl -- convert vcf file to tabular file
#
#   Author: Nowind
#   Created: 2012-02-21
#   Updated: 2016-09-28
#   Version: 1.2.1
#
#   Change logs:
#   Version 1.0.0 16/07/30: The initial version.
#   Version 1.1.0 16/08/12: Updated: add option "--compact" to specify output format.
#   Version 1.2.0 16/09/21: Updated: add option "--depth" to generate overall allele depths.
#   Version 1.2.1 16/09/28: Updated: add square brackets around alleles.




use strict;

use Data::Dumper;
use Getopt::Long;

use MyPerl::FileIO qw(:all);

##################### Main ####################


my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.2.1';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";
my $SOURCE  = (scalar localtime()) . " Version: $VERSION";

my $min_out_depth = 1;
my ($vcf_file, $output, @query_info_fields, @query_sample_fields,
    $show_sample_alleles, $count_allele_depth, $compact_out);
GetOptions(
            "vcf=s"               => \$vcf_file,
            "output=s"            => \$output,
            
            "info=s{,}"           => \@query_info_fields,
            "sample-info=s{,}"    => \@query_sample_fields,
            
            "use-nucl"            => \$show_sample_alleles,
            
            "depth-count"         => \$count_allele_depth,
            "min-depth=i"         => \$min_out_depth,
            
            "compact"             => \$compact_out,
           );

unless( $vcf_file ) {
    print <<EOF;

$0  -- convert vcf file to tabular file

Version: $VERSION

Usage:   perl $0 [options]

Options:
    -v, --vcf     <filename>
        input vcf file, required
    -o, --output  <filename>
        output filename, default to STDOUT

    -i, --info    <strings>
        overall INFO fields to output, can have multiple values
    -s, --sample-info
        sample infos to output, can have multiple values [default: GT]
    
    -u, --use-nucl
        use nucleotides instead of GT code (e.g. 0,1,...) for each sample,
        the REF and ALT fields will be omitted
    
    -c, --compact
        make results more compact by output all samples in a single row
    
    -d, --depth-count
        output allele depth, require "AD" field
    -m, --min-depth <int>
        ignore alleles with total depth below this threshold, [default: 1]

EOF

    exit(1);
}

$|++;

if ($output) {
    open (STDOUT, "> $output") || die $!;
}

unless(@query_sample_fields > 0) {
    push @query_sample_fields, 'GT';
}
 
print STDERR "# $0 v$VERSION\n# " . (scalar localtime()) . "\n";

print STDERR ">> Start processing $vcf_file ... ";
print STDOUT "##source=$SOURCE $CMDLINE\n";
processVCF($vcf_file);
print STDERR "done!\n";

print STDERR "# " . (scalar localtime()) . "\n";


######################### Sub #########################

sub processVCF
{
    my ($in) = @_;
    
    my @Samples_ids = ();
    my $fh = getInputFilehandle($in);
    while (<$fh>)
    {
        if (/#CHROM/) {
            my @lines = (split /\s+/);
            for (my $i=9; $i <@lines; $i++)
            {
                push @Samples_ids, $lines[$i];
            }
            
            my $out_sample_ids = join "\t", @Samples_ids;
            my $query_sample_header = join "\;", @query_sample_fields;
            
            
            if ($compact_out) {
                $out_sample_ids = "SAMPLES($query_sample_header)";
            }
            
            my $out_header = ($show_sample_alleles) ? "#CHROM\tPOS" : "#CHROM\tPOS\tREF\tALT";
            
            
            if (@query_info_fields > 0) {
                my $out_infos = join "\t", @query_info_fields;
                
                $out_header .= "\t$out_infos";
            }
            
            if ($count_allele_depth) {
                $out_header .= "\tOverall_Allele_Depth";
            }
            
            print STDOUT "$out_header\t$out_sample_ids\n";
        }
        
        next if (/\#/ || /^\s+$/);
        
        my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, 
            $FILTER, $INFO, $FORMAT, @SAMPLES) = (split /\s+/);
        
        
        ##
        ## get query overall infos
        ##
        my @infos = (split /\;/, $INFO);
        my %infos = ();
        for (my $i=0; $i<@infos; $i++)
        {
            if($infos[$i] =~ /(\w+)\=(.*)/) {
                $infos{$1} = $2;
            }
        }
        
        my @out_infos = ();
        
        for my $info_id (@query_info_fields)
        {
            if ($infos{$info_id}) {
                push @out_infos, $infos{$info_id};
            }
            else {
                push @out_infos, '-';
            }
        }
        
        
        ##
        ## get query sample infos
        ##
        my @vars = ($REF, (split /\,/, $ALT));
        
        my @tags = (split /\:/, $FORMAT);
        my %tags = ();
        for (my $i=0; $i<@tags; $i++) { $tags{$tags[$i]} = $i; }
        
        my %total_allele_depths = ();
        
        my @out_samples = ();
        for (my $i=0; $i<@SAMPLES; $i++)
        {
            my @sample_infos = (split /\:/, $SAMPLES[$i]);
            
            my @out_sample_info = ();
            
            for (my $j=0; $j<@query_sample_fields; $j++)
            {
                if ($show_sample_alleles && $query_sample_fields[$j] eq 'GT') {
                    ## use nucleotides instead of numbers
                    if ($sample_infos[$tags{$query_sample_fields[$j]}] =~ /(\d)(\/|\|)(\d)/) {
                        $sample_infos[$tags{$query_sample_fields[$j]}] = $vars[$1] . $2 . $vars[$3];
                    }
                }
                
                if (exists($tags{$query_sample_fields[$j]}) && $sample_infos[$tags{$query_sample_fields[$j]}]) {
                    push @out_sample_info, $sample_infos[$tags{$query_sample_fields[$j]}];
                }
                else {
                    push @out_sample_info, '-';
                }
            }
            
            
            ## count allele depths
            if ($count_allele_depth) {
                if ($tags{AD} && ($sample_infos[$tags{AD}]) && ($sample_infos[$tags{AD}] ne '.')) {
                    my @depths = split /\,/, $sample_infos[$tags{AD}];
                    
                    for (my $k=0; $k<@depths; $k++)
                    {
                        $total_allele_depths{$vars[$k]} += $depths[$k];
                    }
                }
            }
            
            unless(@out_sample_info > 0) {
                print STDERR Dumper(@out_sample_info);exit;
            }
            
            my $out_sample_info = join ':', @out_sample_info;
            
            if ($compact_out) {
                push @out_samples, "$Samples_ids[$i]($out_sample_info)";
            }
            else {
                push @out_samples, $out_sample_info;
            }
        }
        
        my $out_samples = join "\t", @out_samples;
        
        if ($compact_out) {
            $out_samples = join ";", @out_samples;
        }
        
        
        
        ##
        ## count overall allele depth
        ##
        my @allele_depth_info = ();
        
        if ($count_allele_depth) {
            for my $nt (@vars)
            {
                if ($total_allele_depths{$nt} && $total_allele_depths{$nt} >= $min_out_depth) {
                    push @allele_depth_info, "[$nt]:$total_allele_depths{$nt}";
                }
            }
        }
        
        my $allele_depth_info = (@allele_depth_info > 0) ? (join ',', @allele_depth_info) : '-';
        
        my $out_line = $show_sample_alleles ? "$CHROM\t$POS" : "$CHROM\t$POS\t$REF\t$ALT";
        
        if (@query_info_fields > 0) {
            my $out_infos = join "\t", @out_infos;
            
            $out_line .= "\t$out_infos";
        }
        
        if ($count_allele_depth) {
            $out_line .= "\t$allele_depth_info";
        }
        
        print STDOUT "$out_line\t$out_samples\n";
    }
}
