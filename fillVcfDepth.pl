#!/usr/bin/perl -w
#
#   fillVcfDepth.pl -- fill the depth info of each sample each locus in the given vcf file
#                      through the results from VarScan readcounts
#
#   Author: Nowind
#   Created: 2012-02-21
#   Updated: 2019-03-20
#   Version: 1.2.0
#
#   Change logs:
#   Version 1.0.0 12/12/21: The initial version.
#   Version 1.0.1 13/01/04: Change the way to import functions from MyPerl::FileIO.
#   Version 1.0.2 13/08/01: Fix GT tag for missing alleles.
#   Version 1.0.3 13/12/04: Add support for counting reads from more than one file;
#                           remove records with no DPA found in results; rearrange
#                           output results.
#   Version 1.0.4 14/05/15: Add support for samples contain different libraries.
#   Version 1.0.5 14/06/04: Add some descriptions of LAD tag.
#   Version 1.0.6 14/11/25: Correct NUMBER descriptions of newly add FORMAT fields
#                           in output.
#   Version 1.0.7 14/12/03: Now the LN and LAD field will not be written if no
#                           multiple libraries is found.
#   Version 1.1.0 15/08/16: Bug fixed in processing mixed loci.
#   Version 1.2.0 19/03/30: Bug fixed: failed to fill ungenotyped sites.


use strict;

use Data::Dumper;
use Getopt::Long;

use MyPerl::FileIO qw(:all);
use MyPerl::Vcf qw(:all);
##################### Main ####################


my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.2.0';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";

my $SOURCE  = (scalar localtime()) . " Version: $VERSION";

my ($vcf_file, $depth_file_list, $output, @filters,
    $skip_unphased, $user_type, $update_AD, $minimu_vcf);
GetOptions(
            "vcf=s"            => \$vcf_file,
            "output=s"         => \$output,
            
            "list=s"           => \$depth_file_list,
            
            "filter=s{,}"      => \@filters,
            "phased"           => \$skip_unphased,
            
            "update-AD"        => \$update_AD,
            
            "type=s"           => \$user_type,
            
            "minimum-vcf"      => \$minimu_vcf,
           );

my $show_help = ($CMDLINE =~ /\-help/) ? 0 : 1;

unless( $vcf_file && $depth_file_list && $show_help ) {
    print <<EOF;

$0  -- convert between vcf file and fimg file

Version: $VERSION

Usage:   perl $0 [options]

Options:
    -v, --vcf    <filename>
        input vcf file, required
    -o, --output <filename>
        output filename, default to STDOUT

    -l, --list   <filename>
        file contains list of depth files and related sample names in the
        format:
        
        sample_id library_id filename 
        
        delimited by tab(s) or space(s), required
        
    -f, --filter <strings>
        skip filter loci, can have multiple values, separate by blanks, e.g.
        "LowQual SNPFilter" ...    [default: no filtering]
    -p, --phased
        skip unphased sites
    
    -t, --type   <string>
        set "snp" to process snp sites only, or set "indel" to process indels
        only
        
    -u, --update-AD
        update AD field according to the results of readcounts
    
    -m, --minimum-vcf
        remove original INFO and FORMAT fields
    
EOF

    exit(1);
}




$|++;

if ($output) {
    open (STDOUT, "> $output") || die $!;
}

 
print STDERR "# $0 v$VERSION\n# " . (scalar localtime()) . "\n";

my $filter_str = join '|', @filters;


##
## read files contain number of reads in each strand
##
print STDERR ">> Start parsing depth files in $depth_file_list ... ";
my %read_counts  = ();
my %sample_libs  = ();
my $list_fh = getInputFilehandle($depth_file_list);
while (<$list_fh>)
{
    next if (/\#/ || /^\s+$/);
    
    my ($sample_id, $library_id, $file) = (split /\s+/);
    
    $sample_libs{$sample_id}->{$library_id} = 1;
    
    parse_depth_file(\%read_counts, $sample_id, $library_id, $file);
}
print STDERR "done!\n";



##
## check whether multiple libraries is present
##
my $max_lib_cnt  = 1;
for my $sample (sort keys %sample_libs)
{
    my $lib_cnt = scalar (keys %{$sample_libs{$sample}});
    
    $max_lib_cnt = $lib_cnt if $lib_cnt > $max_lib_cnt;
}

###print STDERR Dumper($max_lib_cnt);exit;


##
## start filling
##
print STDERR ">> Start filling input vcf file $vcf_file ... ";
fillVcfDepth(\%read_counts, $vcf_file);
print STDERR "done!\n";

print STDERR "# " . (scalar localtime()) . "\n";


######################### Sub #########################


=head2 parse_depth_file

    About   : Get reads depth for each base
    Usage   : parse_depth(\%reads_depth, $depth_file);
    Args    : A hash to save all counts
              A file contains reads depth informations
    Returns : Null

=cut
sub parse_depth_file
{
    my ($rh_depth, $sample_id, $library_id, $in) = @_;
    
    my $fh   = getInputFilehandle($in);
    my $head = <$fh>;
    while (<$fh>)
    {
        my ($chrom, $pos, $ref, $depth, $qual_depth, @counts) = (split /\s+/);
        
        $rh_depth->{total}->{"$chrom:$pos"}++;
        
        for (my $i=0; $i<@counts; $i++)
        {
            my ($base, $reads, $strands, $avg_qual, $map_qual,
                $plus_reads, $minus_reads) = (split /\:/, $counts[$i]);
            
            $rh_depth->{$sample_id}->{"$chrom:$pos"}->{$library_id}->{$base}->{plus}  += $plus_reads;
            $rh_depth->{$sample_id}->{"$chrom:$pos"}->{$library_id}->{$base}->{minus} += $minus_reads;
        }
    }
}

=head2 fillVcfDepth

    About   : fill the vcf file with depth infos
    Usage   : fillVcfDepth(\%reads_depth, $vcf_file);
    Args    : Hash of all saved depth info
              Vcf file needs to be processed
    Returns : Null

=cut
sub fillVcfDepth
{
    my ($rh_depth, $vcf_file) = @_;
    
    
    my $out_header  = '';
    my @names       = ();
    
    my $vcf_fh = getInputFilehandle($vcf_file);
    while (<$vcf_fh>)
    {
        if (/#CHROM/) {
            (undef, undef, undef, undef, undef, undef, undef,
             undef, undef, @names) = (split /\s+/);
            
            if ($max_lib_cnt > 1) {
                $out_header .= <<EOF;
##FORMAT=<ID=LN,Number=1,Type=Integer,Description="number of pooled libraries">
##FORMAT=<ID=LAD,Number=.,Type=Integer,Description="allele depth of each library in the order: allele1-lib1,allele1-lib2,allele2-lib1,allele2-lib2,...">
EOF
            }
            
            $out_header .= <<EOF;
##FORMAT=<ID=RC,Number=.,Type=Integer,Description="read counts in the order: ref-forward,ref-reverse,alt1-forward,alt1-reverse,...">
##source=$SOURCE $CMDLINE
EOF
            print "$out_header";
            print "$_";
            
            next;
        }
        elsif (/\#\#/ || /^\s+$/) {
            $out_header .= $_; next;
        }
        
        next if ($skip_unphased && /\//);         ## only considering phased genotypes
        next if ($filter_str && /$filter_str/i);  ## keep variants in repeat regions
        
        my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, 
            $FILTER, $INFO, $FORMAT, @SAMPLES) = (split /\s+/);
        
        
        next unless( $rh_depth->{total}->{"$CHROM:$POS"} ); ## exclude records with no DPA info available
        
        ##
        ## filtering by variant type
        ##
        my $var_type = get_var_type($REF, $ALT);
        
        next if ($user_type && $var_type !~ /$user_type/);
        
        my @vars = ($REF, (split /\,/, $ALT));
        my @tags = (split /\:/, $FORMAT);
        my %tags = ();
        for (my $i=0; $i<@tags; $i++) { $tags{$tags[$i]} = $i; }
        
        ###my @genotyped_samples = grep { $_ ne '.' && $_ ne './.' } @SAMPLES;
        
        if( $update_AD && !$tags{AD} ) {
            $FORMAT .= ":AD";
        }
        
        my @filled_samples = ();
        for (my $i=0; $i<@SAMPLES; $i++)
        {
            ###if ($SAMPLES[$i] eq '.' || $SAMPLES[$i] eq './.') {
            ###    push @filled_samples, './.'; next;
            ###}
            
            my $sample_id  = $names[$i];
            my %detail_dps = ();
            my @libraries  = sort keys %{$rh_depth->{$sample_id}->{"$CHROM:$POS"}};
            my $lib_cnt    = scalar @libraries;
            for (my $k=0; $k<@libraries; $k++)
            {
                my $lib_id = $libraries[$k];
                
                for (my $j=0; $j<@vars; $j++)
                {
                    my $base = $vars[$j];
                    
                    my $count_plus  = $rh_depth->{$sample_id}->{"$CHROM:$POS"}->{$lib_id}->{$base}->{plus}  ?
                                      $rh_depth->{$sample_id}->{"$CHROM:$POS"}->{$lib_id}->{$base}->{plus}  : 0;
                    my $count_minus = $rh_depth->{$sample_id}->{"$CHROM:$POS"}->{$lib_id}->{$base}->{minus} ?
                                      $rh_depth->{$sample_id}->{"$CHROM:$POS"}->{$lib_id}->{$base}->{minus} : 0;
                    
                    $detail_dps{all}->{$j}->{plus}  += $count_plus;
                    $detail_dps{all}->{$j}->{minus} += $count_minus;
                    
                    $detail_dps{lib}->{$j}->{$k} = $count_plus + $count_minus;
                }
            }
            
            ##
            ## count all reads in plus and minus strand for each allele
            ##
            my @detail_dps = ();
            my @allele_dps = ();
            for my $allele (sort {$a <=> $b} keys %{$detail_dps{all}})
            {
                push @detail_dps, $detail_dps{all}->{$allele}->{plus};
                push @detail_dps, $detail_dps{all}->{$allele}->{minus};
                
                push @allele_dps, ($detail_dps{all}->{$allele}->{plus} + $detail_dps{all}->{$allele}->{minus});
            }
            
            my $read_counts   = join ',', @detail_dps;
               $read_counts ||= '.';
            
            
            ##
            ## count allele depth in each library
            ##
            my @lib_allele_dps = ();
            for my $allele (sort {$a <=> $b} keys %{$detail_dps{lib}})
            {
                for my $lib (sort {$a <=> $b} keys %{$detail_dps{lib}->{$allele}})
                {
                    push @lib_allele_dps, $detail_dps{lib}->{$allele}->{$lib};
                }
            }
            my $lib_allele_dps   = join ',', @lib_allele_dps;
               $lib_allele_dps ||= '.';
            
            ##
            ## update sample infos
            ##
            my @fields = (split /\:/, $SAMPLES[$i]);
            
            if ($minimu_vcf) {
                @fields = ($fields[$tags{GT}]);
            }
            
            if ($update_AD) {
                my $AD_new   = join ',', @allele_dps;
                   $AD_new ||= '.';
                
                if ($tags{AD} && !$minimu_vcf) {
                    $fields[$tags{AD}] = $AD_new;
                }
                else {
                    push @fields, $AD_new;
                }
            }
            
            $SAMPLES[$i] = join ':', @fields;
            
            if ($max_lib_cnt > 1) {
                push @filled_samples, "$SAMPLES[$i]:$read_counts:$lib_cnt:$lib_allele_dps";
            }
            else {
                push @filled_samples, "$SAMPLES[$i]:$read_counts";
            }
        }
        
        my $filled_samples = join "\t", @filled_samples;
        
        my $FORMAT_new = "$FORMAT:RC:LN:LAD";
        
        if ($minimu_vcf) {
            if ($update_AD) {
                $FORMAT_new = ($max_lib_cnt > 1) ? "GT:AD:RC:LN:LAD" : "GT:AD:RC";
            }
            else {
                $FORMAT_new = ($max_lib_cnt > 1) ? "GT:RC:LN:LAD" : "GT:RC";
            }
        }
        
        print "$CHROM\t$POS\t$ID\t$REF\t$ALT\t$QUAL\t$FILTER\t.\t$FORMAT_new\t$filled_samples\n";
    }
}
