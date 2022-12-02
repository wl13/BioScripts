#!/usr/bin/perl -w
#
#   calc_vcf_diversity.pl --  Calculate diversities within or among different groups.
#
#
#   Author: Nowind
#   Created: 2011-01-08
#   Updated: 2022-12-02
#
#   Change logs:
#   Version 1.0.0 12/12/03: The initial version.
#   Version 1.0.1 12/12/04: Bug fixed: return an invalid value while paring interval file in function
#                           parse_intervals.
#   Version 1.1.0 12/12/18: Change input group file to a sample panel file; add option "--non-pairwise"
#                           to calculate diversities compare to reference sample.
#   Version 1.1.1 12/12/21: Bug fixed while parsing non-informatic GT fields.
#   Version 1.1.2 13/01/25: Change the way to import functions from MyPerl::FileIO.
#   Version 1.1.3 13/02/04: Bug fixed: last window is not included while the window size could not reach
#                           the default one.
#   Version 1.2.0 13/03/21: Add options to handling heterozygous calls and missing calls; add option
#                           "--no-group" to calculate diversities between all samples; change to a new
#                           algorithm to reduce memory usage and increase calculating efficiency; add
#                           support for multi-threading.
#   Version 1.2.1 13/03/27: Remove missing alleles from comparation.
#   Version 1.2.2 13/04/02: Add option "--cmp-to-sample" to specify a sample as compare reference.
#   Version 1.2.3 13/04/03: Add option "--within" to only calculate within diversities.
#   Version 1.2.4 13/05/09: Add option "--var-type" to count in only snps or indels.
#   Version 1.2.5 13/10/06: Bug fixed while sample ids contain character "-".
#   Version 1.2.6 13/10/10: Rearrange orders of output pairs according to the input sample panel file
#                           while compare to an reference sample.
#   Version 1.2.7 13/10/11: Bug fixed while output samples in order.
#   Version 1.2.8 14/01/16: Bug fixed in parsing GT field.
#   Version 1.2.9 14/03/05: Fix some spelling errors.
#   Version 1.3.0 14/06/22: Do not fork child processes while not run in parallel.
#   Version 1.3.1 14/06/26: Bug fixed while --model option receive parameters as numbers not strings. 
#   Version 1.3.2 14/12/18: Bug fixed while calculating diversity compared to reference by useing
#                           "--non-pairwise" option.
#   Version 1.3.3 14/12/19: Add option "--length-file" to read chromosome infos from a file.
#   Version 1.3.4 14/12/21: Bug fixed: no results returned while using "--cmp-to-sample" without "--group"
#                           option.
#   Version 1.3.5 15/10/31: Bug fixed in processing mixed loci.
#   Version 1.4.0 15/11/17: Update: use Statistics::Descriptive instead of MyPerl::Statistics; add option
#                           "--exclude"; add some detailed explanations.
#   Version 2.0.0 16/02/10: Update: change script name to "calc_vcf_diversity.pl"; add option "--informative"
#                           to use informative sites in calculating pairwise diversity; add option
#                           "--weight-for-het".
#   Version 2.0.1 17/05/28: Bug fixed: duplicate results while vcf header and length file both present.
#   Version 2.1.0 22/12/02: Update: (1) remove multi-allelic sites as they are not supported;
#                                   (2) change diff(R/A, R/A) from 0 to 0.5 when using "--weight-for-het" to better estimating haploid diversity;
#                                   (3) add related notes in descriptions.





=head1 NAME

calc_vcf_diversity.pl


=head1 SYNOPSIS

calc_vcf_diversity.pl --help/?

=head2 example1: calculate pairwise diversity in non-overlapping windows

calc_vcf_diversity.pl --vcf test.vcf.gz --window-size 500000 \
    --no-group --no-within --output test.pairs_divs.csv


=head2 example2: calculate pairwise diversity in specified regions

calc_vcf_diversity.pl --vcf test.vcf.gz --intervals test.regions \
    --threads 2 --no-within --no-group --output test.pairs_divs.csv


=head2 example3: calculate within and between group diversity in non-overlapping windows

calc_vcf_diversity.pl --vcf test.vcf.gz --window-size 500000 \
    --group test.sample_panel --all-pairs --output test.group_divs.csv


=head2 example4: calculate within and between group diversity in specified regions

calc_vcf_diversity.pl --vcf test.vcf.gz --intervals test.regions \
    --group test.sample_panel --all-pairs --output test.group_divs.csv

  
=head1 DESCRIPTION

Calculate pairwise diversities within or among different groups.

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
my $VERSION = '2.0.1';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";


my %options = ();
my $max_threads     = 1;
my $window_size     = 250;
my $distance_model  = 'p-distance';
my ($vcf_file, @filters, $interval_file, $step_size, $no_group, $show_help, $user_type,
    $skip_unphased, $sample_panel, $output, $all_pair, $non_pairwise, $infomative_perc,
    $ref_for_missing, $hete_as_missing, $hete_as_alt, $no_within, $only_within, $ref_sample,
    $weight_for_het);
GetOptions(
            "vcf=s"              => \$vcf_file,
            "intervals=s"        => \$interval_file,
            
            "group=s"            => \$sample_panel,
            
            "length-file=s"      => \$options{length_file},
            
            "L|window-size=i"    => \$window_size,
            "S|step-size=i"      => \$step_size,
            
            "filter=s{,}"        => \@filters,
            "phased"             => \$skip_unphased,
            
            "R|ref-for-missing"  => \$ref_for_missing,
            "H|het-as-missing"   => \$hete_as_missing,
            "A|het-as-alt"       => \$hete_as_alt,
            
            "model=s"            => \$distance_model,
            "output=s"           => \$output,
            
            "N|no-group"         => \$no_group,
            "B|no-within"        => \$no_within,
            "W|within"           => \$only_within,
            
            "all-pairs"          => \$all_pair,
            "D|non-pairwise"     => \$non_pairwise,
            
            "cmp-to-sample=s"    => \$ref_sample,
            
            "threads=i"          => \$max_threads,
            
            "help|?"             => \$show_help,
            
            "V|var-type=s"       => \$user_type,
            
            "informative=f"      => \$infomative_perc,
            
            "exclude=s{,}"       => \@{$options{exclude_ids}},
            
            "weight-for-het"     => \$weight_for_het,
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
    
    *Note: this script only use bi-allelic sites, and multi-allelic sites will
        be ignored. Although variants other than SNPs are accepted, they will
        be treated as same as SNPs (i.e., variant size is not considered)
    
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
        
    -g, --group    <filename>
        file contains group constructions, one line per member, with member
        id and group id, delimited by space or tabs, e.g.
        M1  G1
        M2  G2
        M3  G3
        ...
        if this option is not specified, all Samples will be treated as one
        single group, and calculate the mean distances between each pair of
        samples
        
    -o, --output   <filename>
        output filename, default to STDOUT
    

    -L, --window-size  <int>
        window size(bp) [default:250]
    -S, --step-size    <int>
        step size(bp), default same as window size, indicate non-overlapping
        windows

    -m, --model   <string>
        model for estimating distances, "p-distance" or "Jukes-Cantor"
        p-distance:   proportion of nucleotide sites at which two sequences
        being compared are different, obtained by dividing the number of
        nucleotide differences by the total number of nucleotides compared;
        Jukes-Cantor: Jukes and Cantor (1969) model, assumes an equality of
        substitution rates among sites 
        [default: p-distance]

    -N, --no-group
        treat each sample individually instead of groups
    -W, --within
        only calculate within-group diversities
    -B, --no-within
        only calculate between-group diversities

    -a, --all-pairs
        calculating average values by dividing the total number of the compared
        pairs, otherwise those pairs with no difference will be ignored
        default: ignore pairs with no difference
        *Note: this default behavior is intended to exclude pairs with "false"
        zero differences (e.g no difference seen due to low coverage, etc.),
        in most time this option should be specified
    
    -D, --non-pairwise
        calculating all diversities against reference sample instead of pairwise
        comparison

    -c, --cmp-to-sample  <string>
        compare all other samples to one certain sample specified

    -f, --filter  <strings>
        skip filter loci, can have multiple values, separate by blanks, e.g.
        "LowQual SNPFilter" ... [default: no filtering]
    -p, --phased
        skip unphased sites

    -R, --ref-for-missing
        use the REF allele instead of the default missing genotype. Because it
        is not obvious what ploidy should be used, a user-defined string is 
        used instead (e.g. 0/0)
    -H, --het-as-missing
        assume no heterozygous sites should be present here and treat those
        heterozygous sites as missing alleles
    -A, --het-as-alt
        assume heterozygous calls as homozygous alternative calls
    
    --weight-for-het
        by default this script assumes a haploid or homozygous diploid model,
        by which heterozygous genotype (R/A) is treated as different allele
        from both reference (R/R) and alternative allele (A/A), i.e.
        
            diff(R/A, R/R) = (R/A, A/A) = 1, and diff(R/A, R/A) = 0
            
        specifying this option will consider a more realistic diploid model,
        by which difference between homozygous and heterozygous, as well as
        between heterozygous and heterozygous will be counted as 0.5 instead
        of 1, i.e.
        
            diff(R/A, R/A) = diff(R/A, R/R) = diff(R/A, A/A) = 0.5
        
        Higher ploidy is not supported currently
    
    -t, --threads  <int>
        how many data threads should be allocated to running this analysis
        [default: 1]

    -V, --var-type    <string>
        set "snp" to process snp sites only, or set "indel" to process indels
        only

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
## get group structures
##
my %group_members = ();
my %group_count   = ();
my @id_orders     = ();
if ($sample_panel) {
    my $g_fh = getInputFilehandle($sample_panel);
    while (<$g_fh>)
    {
        next if (/\#/ || /^\s+$/);
        my ($member_id, $group_id) = (split /\s+/)[0,1];
        
        $member_id =~ s/-/_/g; ## replace sample id if contains "-"
        $group_id  =~ s/-/_/g;
        
        $group_id = 'all' if $no_group;
        
        $group_members{$member_id} = $group_id;
        $group_count{$group_id} ++ ;
        
        push @id_orders, $group_id;
    }
}

###print STDERR Dumper(@member_orders);exit;



##
## get sample names
##
print STDERR ">> Start parsing vcf file headers ... ";
my %CHROM_LENGTHS  = ();
my @CHROM_IDS      = ();
my @SAMPLE_NAMES   = ();
my $CONTIG_INFO    = '';
parse_vcf_header($vcf_file);
if ($options{length_file} && (scalar @CHROM_IDS == 0)) {
    get_genome_length(\@CHROM_IDS, \%CHROM_LENGTHS, $options{length_file}, \@{$options{exclude_ids}});
}
print STDERR "done!\n";

if (!$interval_file && (scalar @CHROM_IDS == 0) && !$options{length_file}) {
    print STDERR "Error: Chromosome info not found, please check vcf header or specify a length file!\n";
    exit(2);
}

###print STDERR Dumper(%group_members);exit;


##
## count compared members
##
my %group_pairs = ();
my @group_ids   = sort (keys %group_count);
for (my $i=0; $i<@group_ids; $i++)
{
    next if ($ref_sample && ($group_ids[$i] ne $ref_sample));  ## choose this sample as a reference
    
    ###print STDERR "$ref_sample\n";
    
    if ($non_pairwise) {
        my $group_pair = join "\-", (sort ($group_ids[$i], $group_ids[$i]));
        $group_pairs{id}->{$group_pair} ++;
        $group_pairs{num}->{$group_pair} = 1;
        next;
    }
    
    for (my $j=0; $j<@group_ids; $j++)
    {
        next if ($ref_sample && ($group_ids[$j] eq $ref_sample));  ## choose this sample as a reference
        
        my $group_pair = join "\-", (sort ($group_ids[$i], $group_ids[$j]));
        
        next if ($only_within && ($group_ids[$i] ne $group_ids[$j]));
        next if ($no_within   && ($group_ids[$i] eq $group_ids[$j]));
        
        $group_pairs{id}->{$group_pair}++;
        
        if ($all_pair) {
            if ($group_ids[$i] eq $group_ids[$j]) {
                $group_pairs{num}->{$group_pair} = $group_count{$group_ids[$i]} *
                                                  ($group_count{$group_ids[$i]} - 1) / 2;
            }
            else {
                $group_pairs{num}->{$group_pair} = $group_count{$group_ids[$i]} *
                                                   $group_count{$group_ids[$j]};
            }
        }
    }
}

###print STDERR Dumper($ref_sample, %group_pairs);exit;


##
## parse intervals
##
print STDERR ">> Generating intervals needed to calculate diversities ... ";
my $ra_intervals = parse_intervals($interval_file);
print STDERR "done!\n";



## get all results from each childs
my %group_divers_all = ();
my $pm = new Parallel::ForkManager($max_threads);
if ($max_threads > 1) {
    $pm->run_on_finish(
        sub{
            my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data) = @_;
            
            my ($interval, $rh_group_divers) = @{$data};
            
            $group_divers_all{$interval} = $rh_group_divers;
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
    my $rh_pairwise_divers = calc_pairwise_diversity($vcf_file, $interval);
    ##print STDERR "done!\n";
    
    ##
    ## calculate average diversities
    ##
    ##print STDERR ">> Start calculating group diversities ... ";
    my $rh_group_divers = calc_group_divers($rh_pairwise_divers, $interval);
    ##print STDERR "done!\n";
    
    if ($max_threads > 1) {
        $pm->finish(0, [$interval, $rh_group_divers]);
    }
    else {
        $group_divers_all{$interval} = $rh_group_divers;
    }
}
$pm->wait_all_children;
print STDERR "\tdone!\n";


##
## generate final results
##
print STDERR ">> Start writing results ... ";



my @out_orders = sort keys %{$group_pairs{id}};
if ($ref_sample) {
    @out_orders = ();
    for my $id (@id_orders)
    {
        next if ($id eq $ref_sample);
        my $pair_id = join "\-", (sort ($id, $ref_sample));
        push @out_orders, $pair_id;
    }
}

###print STDERR Dumper(@id_orders);
###print STDERR Dumper(@out_orders);exit;

my $group_pair_ids = join "\t", @out_orders;

print STDOUT "$HEADER##" . (scalar localtime()) . "\n";
print STDOUT "#BIN_ID\tCHROM\tBIN_START\tBIN_END\t$group_pair_ids\n";
for my $interval (@{$ra_intervals})
{
    generate_results($group_divers_all{$interval}, $interval);
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
            
            ## assume all samples in a single group if no group list specified
            unless($sample_panel) {
                for my $member_id (@SAMPLE_NAMES)
                {
                    my $group_id = $no_group ? $member_id : 'all';
                    $group_members{$member_id} = $group_id;
                    $group_count{$group_id} ++;
                }
                
                @id_orders = @line[9..$#line];
            }
            
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


=head2 parse_genotype

    About   : Parsing GT field.
    Usage   : parse_genotype($GT, @variants);
    Args    : INFO field for each sample;
              Array of all alleles.
    Returns : Variant and genotype string.

=cut
sub parse_genotype
{
    my ($GT, @VARS) = @_;
    
    $GT =~ /((\d)(\/|\|)(\d))/;
    
    my ($geno, $ref_geno, $alt_geno) = ($1, $2, $4);
    
    ###print STDERR "$GT\t$geno\t$ref_geno\t$alt_geno\n";exit;
    
    if ($GT eq '.' || $GT eq './.') {
        my $ret_val = $ref_for_missing ? "ref:$VARS[0]" : "?";
        return $ret_val;
    }
    elsif ($GT =~ /0(\/|\|)0/ ) {
        return "ref:$VARS[$alt_geno]";
    }
    elsif ($ref_geno == $alt_geno) {
        return "alt:$VARS[$alt_geno]";
    }
    else {
        my $ret_val = "het:$VARS[$alt_geno]";
        if ($hete_as_missing) {
            $ret_val = "?";
        }
        elsif ($hete_as_alt) {
            $ret_val = "alt:$VARS[$alt_geno]";
        }
        
        return $ret_val;
    }
}



=head2 calc_pairwise_diversity

    About   : Calculate pairwise distances.
    Usage   : my $rh_pairwise_divers = calc_pairwise_diversity(\%pairwise_differs, $interval);
    Args    : Vcf file;
              Interval string.
    Returns : Reference to a hase contains pairwise diversities
    
=cut
sub calc_pairwise_diversity
{
    my ($in, $interval) = @_;
    
    my %pair_divers = ();
    
    my ($id, $chrom, $start, $end) = (split /\t/, $interval);
    
    my $region = "$chrom:$start-$end";
    
    my %differs    = ();
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
        
        next unless (scalar @vars == 2);          ## remove multi-allelic sites as they are not supported
        
        my $var_type = get_var_type($REF, $ALT);
        
        next if ($user_type && $var_type !~ /$user_type/); ## process snp or indel only
        
        my @tags = (split /\:/, $FORMAT);
        my %tags = ();
        for (my $i=0; $i<@tags; $i++) { $tags{$tags[$i]} = $i; }
        
        for (my $i=0; $i<@Samples; $i++)
        {
            my $sample1 = $SAMPLE_NAMES[$i];
            
            next unless( $group_members{$sample1} );             ## unused sample
            
            my $GT1 = (split /\:/, $Samples[$i])[$tags{GT}];
            
            my $genotype1 = parse_genotype($GT1, @vars);
            
            next if ($genotype1 eq "?");  ## missing calls
            
            if ($non_pairwise) {          ## compare all samples with reference sample
                my $genotype2 = "ref:$REF";

                if ($infomative_perc) {
                    $info_sites{"$SAMPLE_NAMES[$i]\-$SAMPLE_NAMES[$i]"} ++;
                }
                
                if ($weight_for_het && ($genotype1 =~ /het/ || $genotype2 =~ /het/)) {
                    $differs{"$SAMPLE_NAMES[$i]\-$SAMPLE_NAMES[$i]"} += 0.5;
                } elsif ($genotype1 ne $genotype2) {
                    $differs{"$SAMPLE_NAMES[$i]\-$SAMPLE_NAMES[$i]"} ++;
                }
            }
            else {                       ## pairwise comparison
                for (my $j=$i+1; $j<@Samples; $j++)
                {
                    my $sample2 = $SAMPLE_NAMES[$j];
                    
                    next unless( $group_members{$sample2} );                ## unused sample
                    next if ($only_within && ($group_members{$sample1} ne   ## skip samples belong to the different groups
                                              $group_members{$sample2}));
                    next if ($no_within   && ($group_members{$sample1} eq   ## skip samples belong to the same group
                                              $group_members{$sample2}));
                    next if ($ref_sample && ($sample1 ne $ref_sample)       ## all samples compare to a reference sample
                                         && ($sample2 ne $ref_sample));
                    
                    
                    my $GT2 = (split /\:/, $Samples[$j])[$tags{GT}];
                    my $genotype2 = parse_genotype($GT2, @vars);
                    
                    next if ($genotype2 eq "?");  ## missing calls
                    
                    if ($infomative_perc) {
                        $info_sites{"$sample1\-$sample2"} ++;
                    }
                    
                    ###my $pair = join "\t", (sort ($Names[$i], $Names[$j]));
                    
                    if ($weight_for_het && ($genotype1 =~ /het/ || $genotype2 =~ /het/)) {
                        $differs{"$sample1\-$sample2"} += 0.5;
                    } elsif ($genotype1 ne $genotype2) {
                        $differs{"$sample1\-$sample2"} ++;
                    }
                }
            }
        }
    }
    
    ###print STDERR Dumper(%info_sites);
    
    my $len = $end - $start + 1;
    
    for my $pair (keys %differs)
    {
        my $distance = $differs{$pair} / $len;
        
        if ($infomative_perc) {
            $distance = (100 * ($info_sites{$pair} / $len) >= $infomative_perc) ?
                        $differs{$pair} / $info_sites{$pair} : -1;
                        
            ###print STDERR "$chrom\t$start\t$end\t$differs{$pair}\t$info_sites{$pair}\n";
        }
        
        ###print STDERR "Choose model: $distance_model\n";
        
        if ($distance_model eq 'Jukes-Cantor') {
            ###print STDERR "Use Jukes-Cantor, before: $distance\n";
            $distance = ($distance >= 3/4) ? -1 : (0 - (3/4) * log(1 - (4/3*$distance)));
            ###print STDERR "Use Jukes-Cantor, after: $distance\n";
        }
        
        $pair_divers{$chrom}->{$id}->{$pair} = $distance;
    }


    return \%pair_divers;
}



=head2 calc_group_divers

    About   : Calculate differences within or between groups
    Usage   : my $rh_group_divers = calc_group_divers(\%pairwise_divers, $interval);
    Args    : Hash of pairwise distances;
              Interval string.
    Returns : Reference to a hase contains mean value of group diversities
    
=cut
sub calc_group_divers 
{
    my ($rh_pairwise_divers, $interval) = @_;
    
    my %total_divers = ();

    my ($id, $chrom, $start, $end) = (split /\t/, $interval);
    
    my %group_divers = ();
    
    for my $pair (keys %{$rh_pairwise_divers->{$chrom}->{$id}})
    {
        my ($sample1, $sample2) = (split /\-/, $pair);
        
        ###print STDERR "$pair\t$sample1\t$sample2\t$group_members{$sample1}\t$group_members{$sample2}\n";
        
        my $group_pair = join "\-", (sort ($group_members{$sample1}, $group_members{$sample2}));
        
        push @{$group_divers{$group_pair}}, $rh_pairwise_divers->{$chrom}->{$id}->{$pair};
    }
    
    for my $group_pair (keys %group_divers)
    {
        my $stats = Statistics::Descriptive::Full->new();
           $stats->add_data(@{$group_divers{$group_pair}});
           
        my $mean  = $all_pair ? ($stats->sum / $group_pairs{num}->{$group_pair}) : ($stats->mean);
        
        $total_divers{$chrom}->{$id}->{$group_pair} = $mean;
    }

    return \%total_divers;
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
    my ($rh_divers, $interval) = @_;

    my ($id, $chrom, $start, $end) = (split /\t/, $interval);
    
    my @group_divers = ();
    
    for my $group_pair (@out_orders)
    {
        my $group_divers = $rh_divers->{$chrom}->{$id}->{$group_pair} ?
                           $rh_divers->{$chrom}->{$id}->{$group_pair} : 0;
        push @group_divers, $group_divers;
    }
    
    my $divers = join "\t", @group_divers;
    
    print STDOUT "$id\t$chrom\t$start\t$end\t$divers\n";
}



