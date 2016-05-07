#!/usr/bin/perl -w
#
#   count_vcf_pairs.pl --  Count informative sites between each pairs.
#
#
#   Author: Nowind
#   Created: 2011-01-08
#   Updated: 2016-02-10
#
#   Change logs:
#   Version 1.0.0 16/02/10: The initial version.





=head1 NAME

count_vcf_pairs.pl


=head1 SYNOPSIS

count_vcf_pairs.pl --help/?

  
=head1 DESCRIPTION

Count informative sites between each pairs.

=cut


use strict;

use File::Find::Rule;
use Getopt::Long;
use Data::Dumper;
use Data::Random qw(:all);
use Parallel::ForkManager;

use MyPerl::FileIO qw(:all);
use MyPerl::Vcf qw(:all);
use Statistics::Descriptive;


############################# Main ###########################

my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.0.0';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";


my %options = ();
my $max_threads     = 1;
my $window_size     = 250;
my ($vcf_file, @filters, $interval_file, $step_size, $show_help, $user_type,
    $skip_unphased, $output,
    $ref_for_missing, $hete_as_missing, $hete_as_alt, $ref_sample);
GetOptions(
            "vcf=s"              => \$vcf_file,
            "intervals=s"        => \$interval_file,
            
            "length-file=s"      => \$options{length_file},
            
            "L|window-size=i"    => \$window_size,
            "S|step-size=i"      => \$step_size,
            
            "filter=s{,}"        => \@filters,
            "phased"             => \$skip_unphased,
            
            "R|ref-for-missing"  => \$ref_for_missing,
            "H|hete-as-missing"  => \$hete_as_missing,
            "A|hete-as-alt"      => \$hete_as_alt,
            
            "output=s"           => \$output,
            
            
            "cmp-to-sample=s"    => \$ref_sample,
            
            "threads=i"          => \$max_threads,
            
            "help|?"             => \$show_help,
            
            "V|var-type=s"       => \$user_type,
            
            "exclude=s{,}"       => \@{$options{exclude_ids}},
           );

unless( !$show_help && $vcf_file ) {
    print <<EOF;

$0  -- Calculate diversities within or between different groups.

Version: $VERSION

Usage:   perl $0 [options]

Options:
    -v, --vcf       <filename>
        multiple-Samples vcf file, compressed by bgzip and index by tabix,
        required
        
    --intervals     <filename>
        file contans one or more genomic intervals over which to operate,
        calculate diversities only in those specified intervals instead
        of slide-windows along the whole chromosome, each interval contains
        4 values "id chrom start end", delimited by space or tabs e.g.
        interval01 chr01 100 200
        interval01 chr01 500 700
        ...
    
    --informative   <float>
        calculating pairwise diversity by dividing informative sites rather
        than interval length, note this option requires the vcf file contain
        all informative sites, only windows with effective sites over the
        specified percentage would be considered as informative
    
    -l, --length-file  <filename>
        a file contains chromosome length should be specified here in the
        format:

        #CHROM LENGTH
        chr01 43270923
        chr02 35937250
        chr03 36413819
        chr04 35502694
        ...
    
        required while chromosome info is absent in vcf header
        
    -e, --exclude <strings>
        exclude unwanted chromosomes or scaffolds while calculating, all
        chromosomes with ids match strings specified here would be ignored 
        
    -o, --output   <filename>
        output filename, default to STDOUT
    

    -L, --window-size  <int>
        window size(bp) [default:250]
    -S, --step-size    <int>
        step size(bp), default same as window size, indicate non-overlapping
        windows



    -c, --cmp-to-sample  <string>
        compare all other samples to one certain sample specified

    -f, --filter  <strings>
        skip filter loci, can have multiple values, separate by blanks, e.g.
        "LowQual SNPFilter" ... [default: no filtering]
    -p, --phased
        skip unphased sites

    -R, --ref-for-missing
        use the REF allele instead of the default missing genotype. Because it
        is not obvious what ploidy should be used, a user-defined string is used
        instead (e.g. 0/0)
    -H, --hete-as-missing
        assume no heterozygous sites should be present here and treat those
        heterozygous sites as missing alleles
    -A, --hete-as-alt
        assume heterozygous calls as homozygous alternative calls
    
    -t, --threads  <int>
        how many data threads should be allocated to running this analysis
        [default: 1]

    -V, --var-type    <string>
        set "snp" to process snp sites only, or set "indel" to process indels only

    -?, --help
        show this help message

EOF

    exit(1);
}

$|++;


print STDERR "# $0 v$VERSION\n# " . (scalar localtime()) . "\n";

$step_size = $window_size if (!$step_size);

if ($output) {
    open (STDOUT, "> $output") || die $!;
}

my %filters    = ();
my $filter_str = join '|', @filters;


unless( -f "$vcf_file.tbi" ) {
    print STDERR "[Error]: no index found for vcf file!\n"; exit(2);
}



##
## get sample names
##
print STDERR ">> Start parsing vcf file headers ... ";
my %CHROM_LENGTHS  = ();
my @CHROM_IDS      = ();
my @SAMPLE_NAMES   = ();
my $CONTIG_INFO    = '';
my %pair_ids       = ();
my @id_orders      = ();
parse_vcf_header($vcf_file);
if ($options{length_file}) {
    get_genome_length(\@CHROM_IDS, \%CHROM_LENGTHS, $options{length_file}, \@{$options{exclude_ids}});
}
print STDERR "done!\n";

if (!$interval_file && @CHROM_IDS == 0 && !$options{length_file}) {
    print STDERR "Error: Chromosome info not found, please check vcf header or specify a length file!\n";
    exit(2);
}



##
## parse intervals
##
print STDERR ">> Generating intervals needed to calculate diversities ... ";
my $ra_intervals = parse_intervals($interval_file);
print STDERR "done!\n";



## get all results from each childs
my %pair_infos_all = ();
my $pm = new Parallel::ForkManager($max_threads);
if ($max_threads > 1) {
    $pm->run_on_finish(
        sub{
            my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data) = @_;
            
            my ($interval, $rh_pair_infos) = @{$data};
            
            $pair_infos_all{$interval} = $rh_pair_infos;
        }
    );
}


my $total_jobs_cnt = scalar @{$ra_intervals};
my $curr_job_num   = 0;

for (my $i=0; $i<@{$ra_intervals}; $i++)
{
    my $interval = $ra_intervals->[$i];
    
    $curr_job_num++;
    
    print STDERR "\r>> Start calculating diversities ... $curr_job_num/$total_jobs_cnt";
    
    my $pid  = $pm->start and next if ($max_threads > 1);
    
    ##
    ## calculate diversities
    ##
    ##print STDERR ">> Start calculating pairwise distances ... ";
    my $rh_pair_infos = count_pairwise_informative($vcf_file, $interval);
    ##print STDERR "done!\n";
    
    if ($max_threads > 1) {
        $pm->finish(0, [$interval, $rh_pair_infos]);
    }
    else {
        $pair_infos_all{$interval} = $rh_pair_infos;
    }
}
$pm->wait_all_children;
print STDERR "\tdone!\n";


##
## generate final results
##
print STDERR ">> Start writing results ... ";
my @out_orders = ();
if ($ref_sample) {
    for my $id (@id_orders)
    {
        next if ($id eq $ref_sample);
        my $pair_id = join "\-", (sort ($id, $ref_sample));
        push @out_orders, $pair_id;
    }
}
else {
    for (my $i=0; $i<@id_orders; $i++)
    {
        for (my $j=$i+1; $j<@id_orders; $j++) {
            my $pair_id = join "\-", ($id_orders[$i], $id_orders[$j]);
            push @out_orders, $pair_id;
        } 
    }
}

###print STDERR Dumper(@id_orders);
###print STDERR Dumper(@out_orders);exit;

my $out_pair_ids = join "\t", @out_orders;

print STDOUT "$HEADER##" . (scalar localtime()) . "\n";
print STDOUT "#BIN_ID\tCHROM\tBIN_START\tBIN_END\t$out_pair_ids\n";
for my $interval (@{$ra_intervals})
{
    generate_results($pair_infos_all{$interval}, $interval);
}

print STDERR "done!\n";

print STDERR "# " . (scalar localtime()) . "\n";

########################## Sub #########################


=head2 parse_vcf_header

    About   : Parse vcf header informations.
    Usage   : parse_vcf_header($interval_file);
    Args    : Vcf file.
    Returns : Null

=cut
sub parse_vcf_header
{
    my ($in) = @_;

    ## exclude unwanted chromosomes or scaffolds while simulating, all
    ## chromosomes with ids match strings specified here would be ignored 
    my $exclude_str = '';
    if ($options{exclude_ids} && @{$options{exclude_ids}} > 0) {
        $exclude_str = join '|', @{$options{exclude_ids}};
    }
    
    my $fh = getInputFilehandle($in);
    while (<$fh>)
    {
        if (/\#\#contig=<ID=(.*?),length=(\d+)>/) {
            my ($chrom_id, $chrom_len) = ($1, $2);
            
            if ($options{exclude_ids} && @{$options{exclude_ids}} > 0) {
                next if ($chrom_id =~ /($exclude_str)/);
            }
            
            $CHROM_LENGTHS{$chrom_id} = $chrom_len;

            push @CHROM_IDS, $chrom_id;
            
            $CONTIG_INFO .= $_;
        }
        elsif (/#CHROM/) {
            my @line = (split /\s+/);
            for (my $i=9; $i<@line; $i++)
            {
                $line[$i] =~ s/-/_/g; ## replace sample id if contains "-"
                
                $SAMPLE_NAMES[$i-9] = $line[$i];
            }
            
            @id_orders = @line[9..$#line];
            
            return 0;
        }
    }
}



=head2 parse_intervals

    About   : Parse interval list or generate windows along the whole genome.
    Usage   : parse_intervals($interval_file);
    Args    : File contains all intervals (optional).
    Returns : Array reference to all intervals or windows.

=cut
sub parse_intervals
{
    my ($interval_file) = @_;
    
    my @intervals = ();
    
    ## read intervals from a file
    if ($interval_file) {
        my $fh = getInputFilehandle($interval_file);
        
        while(<$fh>)
        {
            next if (/\#/ || /^\s+$/);
            
            my ($id, $chrom, $start, $end) = (split /\s+/);
            
            push @intervals, "$id\t$chrom\t$start\t$end";
        }
    }
    else {
        ## use sliding-windows
        for my $chrom (@CHROM_IDS)
        {
            my $bin_id    = 1;
            my $bin_start = 1;
            
            while ($bin_start <= $CHROM_LENGTHS{$chrom})
            {
                my $bin_end = $bin_start + $window_size - 1;
                   $bin_end = $CHROM_LENGTHS{$chrom} if ($bin_end > $CHROM_LENGTHS{$chrom});  ## the last window
        
                push @intervals, "$bin_id\t$chrom\t$bin_start\t$bin_end";
                
                $bin_start += $step_size;
                $bin_id    ++;
            }
        }
    }
    
    return \@intervals;
}


=head2 count_pairwise_informative

    About   : Calculate pairwise distances.
    Usage   : my $rh_pairwise_divers = count_pairwise_informative(\%pairwise_differs, $interval);
    Args    : Vcf file;
              Interval string.
    Returns : Reference to a hase contains pairwise diversities
    
=cut
sub count_pairwise_informative
{
    my ($in, $interval) = @_;
    
    my ($id, $chrom, $start, $end) = (split /\t/, $interval);
    
    my $region = "$chrom:$start-$end";
    
    my %info_sites = ();
    
    open (my $fh, "tabix $in $region |") or die $!;
    while (<$fh>)
    {
        next if (/\#/ || /^\s+$/);
        
        next if ($skip_unphased && /\//);         ## only considering phased genotypes
        next if ($filter_str && /$filter_str/i);  ## keep variants in repeat regions
        
        my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO,
            $FORMAT, @Samples) = (split /\s+/);
        
        my @vars = ($REF, (split /\,/, $ALT));

        my $var_type = get_var_type($REF, $ALT);
        
        next if ($user_type && $var_type !~ /$user_type/); ## process snp or indel only
        
        my @tags = (split /\:/, $FORMAT);
        my %tags = ();
        for (my $i=0; $i<@tags; $i++) { $tags{$tags[$i]} = $i; }
        
        for (my $i=0; $i<@Samples; $i++)
        {
            my $sample1 = $SAMPLE_NAMES[$i];
            
            my $GT1 = (split /\:/, $Samples[$i])[$tags{GT}];
            
            next if ($GT1 eq '.' || $GT1 eq './.');  ## missing calls
            
            ## pairwise comparison
            for (my $j=$i+1; $j<@Samples; $j++)
            {
                my $sample2 = $SAMPLE_NAMES[$j];
                
                next if ($ref_sample && ($sample1 ne $ref_sample)       ## all samples compare to a reference sample
                                     && ($sample2 ne $ref_sample));
                
                my $GT2 = (split /\:/, $Samples[$j])[$tags{GT}];
                
                next if ($GT2 eq '.' || $GT2 eq './.');
                
                $info_sites{$chrom}->{$id}->{"$sample1\-$sample2"}++;
            }
        }
    }

    return \%info_sites;
}




=head2 generate_results

    About   : Generate results.
    Usage   : my $rh_group_divers = generate_results(\%group_divers, $interval);
    Args    : Hash of group diversities;
              Interval string.
    Returns : Null
    
=cut
sub generate_results
{
    my ($rh_infos, $interval) = @_;

    my ($id, $chrom, $start, $end) = (split /\t/, $interval);
    
    my @pair_info_sites = ();
    
    for my $pair_id (@out_orders)
    {
        my $pair_infos = $rh_infos->{$chrom}->{$id}->{$pair_id} ?
                         $rh_infos->{$chrom}->{$id}->{$pair_id} : 0;
        push @pair_info_sites, $pair_infos;
    }
    
    my $pair_info_sites = join "\t", @pair_info_sites;
    
    print STDOUT "$id\t$chrom\t$start\t$end\t$pair_info_sites\n";
}



