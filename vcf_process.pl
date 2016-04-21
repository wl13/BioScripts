#!/usr/bin/perl -w
#
#   vcf_process.pl -- Process of vcf file.
#
#
#   Author: Nowind
#   Created: 2012-05-30
#   Updated: 2016-02-10
#   Version: 2.4.4
#
#   Change logs:
#   Version 1.0.0 12/09/28: The initial version.
#   Version 1.0.1 12/10/24: Add option "--remove" to remove useless sample rows;
#                           add option "--missing" and "--no-hete" to filter sites
#                           with too much missing alleles.
#   Version 1.0.2 12/12/12: Add option "--unique" and "--quality".
#   Version 1.0.3 12/12/17: Add option "--min-alleles" and "--max-alleles" to set
#                           range of different alleles; add option "--downsample"
#                           to downsampling file randomly.
#   Version 1.0.4 12/12/20: Bug fixed while no missing filtering should be applied.
#   Version 1.0.5 13/01/03: Change the way to import functions from MyPerl::FileIO;
#                           add option "--no-snp-within" to remove snps around indels.
#   Version 1.0.6 13/01/08: Change options: "--missing" to "--max-missing", "--no-hete"
#                           to "--hete-as-missing", "--remove" to "--remove-samples";
#                           add options "--min-non-ref-ac", "--max-non-ref-ac" and
#                           "--var-type".
#   Version 1.0.7 13/01/10: Bug fixed in parsing missing alleles.
#   Version 1.0.8 13/01/11: Bug fixed while no missing alleles is permitted for each site.
#   Version 1.0.9 13/01/17: Add several comments in vcf header.
#   Version 1.1.0 14/04/02: Bug fixed while using option "--var-type" to filter SNPs or INDELs.
#   Version 1.1.1 14/04/08: Add option "--hete-samples".
#   Version 1.1.2 14/04/09: Add options "--min-ref-ac" and "--max-ref-ac".
#   Version 1.1.3 14/04/10: Add options "--min-het-ac" and "--max-het-ac".
#   Version 1.1.4 14/04/16: Add some comments.
#   Version 1.2.0 14/05/10: Add some code for statistics.
#   Version 1.2.1 14/05/11: Bug fixed in counting alleles after regenotyping.
#   Version 1.2.2 14/05/12: Bug fixed in updating GT fields; add process of rare variants;
#                           add option "--minimum-vcf".
#   Version 1.2.3 14/05/14: Bug fixed while no allele available when filtering rare variants.
#   Version 1.2.4 14/05/15: Add option "--check-miss-AD".
#   Version 1.2.5 14/05/16: Bug fixed while no AD tags available when "--check-miss-AD" is specified.
#   Version 1.3.0 14/05/19: Bug fixed in considering true and false heterozygous variants as missing;
#                           add function overlap_vcf and find_segregate_loci.
#   Version 1.4.0 14/05/20: Correct a format error in segregating results; add function markers2blocks.
#   Version 1.4.1 14/05/21: Bug fixed while collecting results from child processes.
#   Version 1.4.2 14/05/23: Bug fixed in parsing sample names while --remove-samples is specified.
#   Version 1.4.3 14/05/27: Add option "--no-ref-rare"; correct a typo.
#   Version 1.5.0 14/05/30: Add function to count base changes; add check of vcf header to support
#                           multiple vcf records from a pipe.
#   Version 1.5.1 14/06/03: Add option "--ref-depth" to use reference allele depth instead minor allele
#                           depth in stats output; add several options in regenotyping processes; change
#                           optiont "--hete-as-missing" to "--max-pseudo-het".
#   Version 1.5.2 14/06/04: Change the threshold used to regentyping heterozygous samples to a range.
#   Version 1.5.3 14/06/09: Output non-reference depth ratio instead of reference depth ratio while
#                           "--ref-depth" is used; set reversed genotyped alleles to missing.
#   Version 1.5.4 14/06/10: Add printing of process status while merged blocks.
#   Version 1.5.5 14/06/11: Add total markers number to merged blocks results; bug fixed in check sample
#                           depth.
#   Version 1.5.6 14/06/17: Revise output vcf header; add options "--default-sample-type" and
#                           "--pseudo-as-missing".
#   Version 1.5.7 14/06/19: Add option "--remain-samples"; rearrange several codes.
#   Version 1.5.8 14/07/02: Add options of allele frequency based filtering in different group of sample types.
#   Version 1.5.9 14/07/06: Set heterozygous alleles in homozygous samples to missing while "--pseudo-as-missing"
#                           even when no regenotyping option is specified.
#   Version 1.6.0 14/07/08: Add function collect_metrics; add variant type filtering in other functions.
#   Version 1.6.1 14/07/10: Add option "--fill-trans" in clustering process.
#   Version 1.6.2 14/07/20: Add support for single sample compare in function find_segregate_loci.
#   Version 1.6.3 14/07/27: Add some comments.
#   Version 1.7.0 14/11/08: Add new function count_diagnose to summary results from GATK DiagnoseTargets;
#                           bug fixed while no GT field present; bug fixed in filtering by FILTER field.
#   Version 1.7.1 14/11/17: Add option "--pass-only" to filter all except "PASS" records; bug fixed:
#                           "--pseudo-as-missing" not effect while AD tag is missing.
#   Version 1.7.2 14/11/18: Add option "--out-genotype-stats" and "--out-locus-stats" to generate various stats
#                           of variants.
#   Version 2.0.0 14/11/18: Move all sub functions to MyPerl::Vcf, re-structure codes in a OO favor.
#   Version 2.0.1 14/11/19: Add option "--check-gt-depth" to toggle regenotyping operations at certain loci;
#                           add some explantions about the regenotyping operations.
#   Version 2.1.0 14/11/24: Add several options to inspect sequence context of variants; add some comments.
#   Version 2.2.0 14/12/12: Major updates in clustering functions: much faster, less memory consumption and more
#                           accurate borders.
#   Version 2.2.1 14/12/19: Add filtering options for alt count.
#   Version 2.2.2 15/01/04: Add option "--exclude-samples" for filtering rare variants.
#   Version 2.2.3 15/01/19: Add options to filter indels by length.
#   Version 2.3.0 15/08/13: Add functions to generate statistics of distance between adjacent markers; add more detailed
#                           results of marker infos in clustered blocks.
#   Version 2.3.1 15/08/16: Bug fixed in processing mixed loci.
#   Version 2.4.0 15/09/26: Add options to combine vcf files; update Vcf.pm to version 1.6.0.
#   Version 2.4.1 15/09/29: Add option "--gt-diff-as-missing" to handle samples with GT fail regenotyping thresholds;
#                           update Vcf.pm to version 1.6.2.
#   Version 2.4.2 15/10/01: Add options to set combine tags while combine two vcf files; update Vcf.pm to version 1.6.3.
#   Version 2.4.3 15/11/04: Add options to set combine rows while combine two vcf files; update Vcf.pm to version 1.6.4.
#   Version 2.4.4 16/02/10: Updated: correct some explanation of options.


=head1 NAME

vcf_process.pl


=head1 SYNOPSIS

  vcf_process.pl --help/?

=head1 DESCRIPTION

Vcf format file related processes.

=cut


use strict;

use Data::Dumper;
use Getopt::Long;
use File::Find::Rule;
use File::Basename;
use Parallel::ForkManager;
use Statistics::Descriptive;
use Statistics::PointEstimation;

use MyPerl::FileIO qw(:all);
use MyPerl::Vcf qw(:all);

################### Main #################

my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '2.4.4';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";
my $SOURCE  = (scalar localtime()) . " Version: $VERSION";


my %options = ();
   $options{source_line} = "##source=$SOURCE $CMDLINE";
   $options{sample_type} = 'het';
   $options{threads} = 1;
   $options{snp_extend}   = 5;
   $options{indel_extend} = 10;
GetOptions(
            "vcf=s"                  => \$options{vcf},
            "depth-file=s"           => \$options{sample_depth_file},
            
            "secondary-vcf=s"        => \$options{secondary_vcf},
            "primary-tag=s"          => \$options{primary_tag},
            "secondary-tag=s"        => \$options{secondary_tag},
            "intersect-tag=s"        => \$options{intersect_tag},
            "combine-rows=i{,}"      => \@{$options{combine_rows}},
            
            "output=s"               => \$options{output},
            
            "out-metrics"            => \$options{out_metrics},
            "metrics=s{,}"           => \@{$options{metrics}},
            
            "contig=s"               => \$options{length_file},
            
            "minimum-vcf"            => \$options{remove_info},
            
            "pass-only"              => \$options{pass_only},
            "filter=s{,}"            => \@{$options{filters}},
            "phased"                 => \$options{skip_unphased},
            
            "var-type=s"             => \$options{var_type},
            
            "remain-samples=s{,}"    => \@{$options{samples_remain}},
            "remove-samples=s{,}"    => \@{$options{samples_remove}},
            "max-missing=f"          => \$options{max_missing},
            "max-hom-missing=i"      => \$options{max_hom_missing},
            "max-het-missing=i"      => \$options{max_het_missing},
            "pseudo-as-missing"      => \$options{pseduo_het_as_missing},
            "missing-as-ref"         => \$options{missing_as_ref},
            
            
            "max-pseudo-het=i"       => \$options{max_pseduo_het},
            "min-pseudo-het=i"       => \$options{min_pseduo_het},
            "hete-samples=s{,}"      => \@{$options{het_samples}},
            "homo-samples=s{,}"      => \@{$options{hom_samples}},
            "default-sample-type=s"  => \$options{sample_type},
            
            "unique"                 => \$options{unique},
            "quality=f"              => \$options{min_quality},
            
            "min-sample-depth=i"     => \$options{min_sample_dp},
            "max-sample-depth=i"     => \$options{max_sample_dp},
            
            "min-alleles=i"          => \$options{min_allele_types},
            "max-alleles=i"          => \$options{max_allele_types},
            
            "min-ref=f"              => \$options{min_ref_count},
            "max-ref=f"              => \$options{max_ref_count},
            
            "min-alt=f"              => \$options{min_alt_count},
            "max-alt=f"              => \$options{max_alt_count},
            
            "min-hom-ref=f"          => \$options{min_hom_ref_count},
            "max-hom-ref=f"          => \$options{max_hom_ref_count},
            "min-het-ref=f"          => \$options{min_het_ref_count},
            "max-het-ref=f"          => \$options{max_het_ref_count},
            
            "min-non-ref=f"          => \$options{min_non_ref_count},
            "max-non-ref=f"          => \$options{max_non_ref_count},
            
            "min-het=f"              => \$options{min_het_count},
            "max-het=f"              => \$options{max_het_count},
            
            "rare-only=f"            => \$options{max_rare_count},
            "no-rare-ref"            => \$options{no_rare_ref},
            "exclude-samples=s{,}"   => \@{$options{samples_exclude}},
            
            "downsample=f"           => \$options{random_size},
            
            "no-snp-within=i"        => \$options{min_indel_range},
            
            "stats-outfile=s"        => \$options{stats_outfile},
            "out-genotype-stats"     => \$options{genotype_stats},
            "out-locus-stats"        => \$options{locus_stats},
            "stats-only"             => \$options{stats_only},
            "ref-depth"              => \$options{use_ref_depth},
            "regenotype-hom=f"       => \$options{hom_threshold},
            "regenotype-het=s"       => \$options{het_threshold},
            "gt-diff-as-missing"     => \$options{gt_diff_as_missing},
            "check-gt-depth"         => \$options{check_gt_dp},
            
            "stat-var-dist"          => \$options{stat_var_dist},
            
            "check-miss-AD"          => \$options{check_miss_AD},
            
            "overlap-gts"            => \$options{overlap_genotypes},
            
            "segregating=s"          => \$options{segregating},
            "out-symbols=s"          => \$options{out_symbols},
            
            "out-blocks"             => \$options{cluster_markers},
            "source-tag=s"           => \$options{source_tag},
            "type2char=s"            => \$options{type2char_str},
            "char2type=s"            => \$options{char2type_str},
            "fill-gaps"              => \$options{fill_gaps},
            "min-frag-length=i"      => \$options{min_frag_length},
            "min-frag-markers=i"     => \$options{min_frag_markers},
            "min-frag-density=f"     => \$options{min_frag_density},
            "min-major-perc=f"       => \$options{min_major_perc},
            "min-block-density=f"    => \$options{min_block_density},
            "fill-trans=s"           => \$options{trans_fill_char},
            
            
            "threads=i"              => \$options{threads},
            
            "sum-diagnose"           => \$options{summary_diagnose},
            
            "base-changes"           => \$options{base_changes},
            "GT-types=s{,}"          => \@{$options{query_GTs}},
            
            
            "check-context"          => \$options{check_context},
            "fasta=s"                => \$options{fasta},
            "extend-snp=i"           => \$options{snp_extend},
            "extend-indel=i"         => \$options{indel_extend},
            
            "min-indel-len=i"        => \$options{min_indel_len},
            "max-indel-len=i"        => \$options{max_indel_len},
           );

unless( $options{vcf} ) {
    print <<EOF;

$0  -- filter vcf format files

Version: $VERSION

Usage:   perl $0 [--vcf FILE | STDIN] [Filtering Options] [Output Options]

Input Options:

    --vcf   <filename>
        input vcf file, support compressed file with suffix *.vcf.gz or
        *.vcf.bz2, required

    --secondary-vcf <filename>
        set a secondary vcf file, which will be combined into the primary
        vcf file specified use "--vcf" option, only new records in secondary
        file will be added to the final vcf file
    
        
    --fasta <filename>
        sequence file in fasta format, required while --context-check option
        is specified
    
    --hete-samples <strings>
        specify heterozygous samples
    --homo-samples <strings>
        specify homozygous samples 
    --default-sample-type <string>
        all sample will be treated as type specifed here if neither of
        preivous two options is specified, default: het
        
        
    --contig  <filename>
        a file contains chromosome names and lengths in the format:
        
        #CHROM LENGTH
        chr01 43270923
        chr02 35937250
        chr03 36413819
        chr04 35502694
        ...
        
        required if no chromosome info found in the markers file header
        
    --metrics <strings>
        metrics to be collected, can have multiple values
    
    
Output Options:

    --output    <filename>
        output file, default to STDOUT

    --minimum-vcf
        remove original info fields in the output vcf file

    --stats-out <filename>
        do some statistics and write the results to this file

    --out-genotype-stats
        output overall stats of genotypes for each sample
    --out-locus-stats
        output count of variants in each locus
        
    --stats-only
        only output statistics results, while this option is specified, the
        output will be redirected to STDOUT if "--stats-out" is not specified
    --ref-depth
        use non-reference allele depth instead minor allele depth in stats
        output

    --out-symbols <strings>
        only output samples not in the segregating groups, assume those samples
        are all from the two segregating groups and assign a symbol to each
        sample stand for the source of this sample, the symbols used could
        be specified here using strings like "G1;G2"
    
    --out-blocks
        try to cluster markers into large blocks

    --base-changes
        output stats of base changes, only bi-allelic snp loci will be counted
    --GT-types    <strings>
        only count base changes of sample with the specified GT types

    --out-metrics
        output variant metrics for SNP or INDELs
    
    --stat-var-dist
        output statistics of adjacent variants distances for each sample
    
Filtering Options:

    --pass-only
        only retain "PASS" loci
    --filter  <strings>
        skip filter loci, can have multiple values, separate by blanks, 
        e.g. "LowQual SNPFilter" ... [default: no filtering]
    --phased
        skip unphased sites
    --quality <float>
        filter sites with quality below this threshold.
    
    --remain-samples <strings>
        only remain specified samples
    --remove-samples <strings>
        remove infos of given samples

    --var-type
        filter by variant types, e.g. "snp" or "indel", mixed locus (with
        both snp and indel) will be considered both as "snp" and "indel",
        MNP will be considered as "snp" if with equal size or "indel" if
        with unequal size
        
    --unique
        retain only unique records, remove those sites with same chr:pos

    --no-snp-within  <int>
        remove snps around indels which have a high possibility of false
        positive, set to x means no snp in the upstream x bp and downstream
        x bp of the indel, set to 0 means only remove snps overlapped with
        indels.

        
    --max-missing <int>
        maxmimum allowed missing numbers, sites with missing alleles more
        than this will be filtered
    --max-hom-missing <int>
        maxmimum allowed missing numbers in homozygous samples, sites with
        missing alleles more than this will be filtered
    --max-het-missing <int>
        maxmimum allowed missing numbers in heterozygous samples, sites with
        missing alleles more than this will be filtered
    --check-miss-AD
        treat allele as missing while neither AD or NR and NV fields are
        available
    --missing-as-ref
        treat missing alleles as reference alleles
    --min-pseudo-het <int>
    --max-pseudo-het <int>
        allowed pseduo-heterozygous alleles in homozygous samples within this
        range
    
    --min-alleles <int>
    --max-alleles <int>
        include only sites with a number of different alleles within the
        specified range

    --min-ref <int>
    --max-ref <int>
        include only sites with all Reference Sample Counts within the
        specified range

    --min-hom-ref <int>
    --max-hom-ref <int>
        include only sites with all Reference Sample Counts within the
        specified range in homozygous samples
    --min-het-ref <int>
    --max-het-ref <int>
        include only sites with all Reference Sample Counts within the
        specified range in heterozygous samples
        
        
        
    --min-non-ref <int>
    --max-non-ref <int>
        include only sites with all Non-Reference Sample Counts within the
        specified range

    --min-het <int>
    --max-het <int>
        include only sites with all Heterozygous Sample Counts within the
        specified range

    --rare-only       <int>
        only output loci contain rare alleles, the rare allele was defined
        as non-reference allele with a frequency no more than the value
        specified here
    --exclude-samples <strings>
        loci shared with those samples will be skipped while scanning for
        rare variants
    
    
    --downsample <float>
        downsampling file, randomly extract n loci(while n > 1) or n
        percentage of total loci(while n <= 1)

    --min-sample-depth <int>
    --max-sample-depth <int>
        samples with depth not within than this range will be treated as
        missing, DP or AD tag required
    --depth-file <filename>
        read sample depth filtering criteria from a file, the file
        should contain following rows
        
        "sample_id" "min_depth" "max_depth"
        
        set to a value less than 0 indicates no filtering
        
    --min-indel-len <int>
        filtering indels with size smaller than this value
    --max-indel-len <int>
        filtering indels with size larger than this value
        
    *Note: for multi-allelic indel loci, the indel length refers to the
        largest allele length
        
        
Group Comparison Options:

    --segregating <strings>
        screen out loci which segregating alleles between two groups, groups
        should be specied here in the format:
        
        "Group1_1,Group1_2,Group1_3,...;Group2_1,Group2_2,Group2_3,..."
        
        can use an un-defined sample name in second group, thus the other
        samples would been determined only through comparison with the
        first group

Genotype manipulation Options:

    --overlap-gts
        screen out loci contain concordant GT fields in all samples
    
    --check-gt-depth
        make sure the genotype do have supporting reads, otherwise set them
        to missing
    --regenotype-hom <float>
        regnotype variants in homozygous samples according to the allele 
        depth calculated above, the homozygousity would be called only when 
        the minor allele count was less than the specified fraction of total 
        reliable mapped reads
    --regenotype-het <string>
        regnotype variants in heterozygous samples according to the calculated
        monior allele depth frequency, set a range here like 0.1,0.3, which
        indicates only minor allele depth frequency below the first threshold
        will be called as homozgyous, while frequency above second threshold
        will be called as heterozygous
    --gt-diff-as-missing
        do not actually regenotype those samples inconsistent with previous
        thresholds, but set those samples as missing instead
        
        
    --pseudo-as-missing
        set heterozygous alleles in homozygous samples to missing

    *Notes for Genotype counting:
        most of the genotype related operations require the present of
        covered depth for each allele, use AD field from GATK UnifiedGenotyper
        or HaplotyperCaller is recommand; while using NR and NV fields from
        Platypus, only bi-allelic loci will be processed as the NR field
        contains depths for each allele which could overlap each other

        
Clustering Options:

    --source-tag <string>
        tag of source type for each sample used to form blocks, can be set
        as "GT", "SC" or others, otherwise will use all sample field as
        the source type

    --fill-gaps
        trying to fill gaps between fragments, the original blocks with
        continuous markers were treated as fragments here, regions with no
        markers or failed to pass minimum threshold required for a reliable
        fragment would be treated as gaps, those gaps will be filled as
        describled below:
        
        merged into flanking fragments while they are in same type or
        merged by extending flanking fragments while they are in different
        types, e.g.
            AAAAAAAAAAAAXXXXAAAAAAAAAAA => AAAAAAAAAAAAaaaaAAAAAAAAAAA
            AAAAAAAAAAAAXXXXBBBBBBBBBBB => AAAAAAAAAAAAaabbBBBBBBBBBBB
        final results also extend to regions with no markers, which use a
        half-half rule to split into flanking blocks
        
        
    --min-frag-length   <int>
        minimum length of a fragment
    --min-frag-markers  <int>
        minimum number of markers in a fragment
    --min-frag-density  <float>
        minimum fraction of markers in a fragment
        
        Those three options are used to determine whether a fragment will be
        treated as "seeds" which used to extend, or "gaps" which need to be
        filled
        
    --min-major-perc    <float>
        blocks with (major markers / total markers) * 100 less than this value
        will be annotated as "Complex"
    --min-block-density <float>
        blocks with (total markers / block length) smaller than this value will
        be annotated as "LowDensity"
        
        Those two options only add tags to output
    
    *********************** deprecated options ***********************
    --type2char <string>
        specify symbols to overide default set, in the format:
        "Source_A:A;Source_B:B;...", required
    --char2type <string>
        reverse convert chars to types, in the format:
        "A:Source_A;B:Source_B;...",
        default will use the reverse of the values "--type2char" set
    
    --fill-trans <string>
        fill gaps between different types of fragments (denoted as transfrom
        regions here) using certain specified character instead of previous
        "half-half" rule, character set here must have a correspond type
        defind in "--type2char" option
    
    *Notes for deprecated options:
        The last three options were only remained for backward compatiblity,
        which use an old string match clustering algorithm, that is
        magnitudes slower than current used one, and could also cause the
        memory to run out if too much samples is given, thus is deprecated
        currently;
        However, if you do want to use the old algorithm, also be cautious
        characters from "X" to "Z" is reserved for special use, do not set
        those characters unless you known what you're doing.
    ******************************************************************


Combine Options:
    
    --primary-tag   <string>
        Tag used for records only from master vcf file in combined vcf file,
        default: Primary
    --secondary-tag <string>
        Tag used for records only from secondary vcf file in combined vcf
        file, default: Secondary
    --intersect-tag <string>
        Tag used for intersection records in combined vcf file,
        default: Intersection

    --combine-rows  <numbers>
        Specify rows (0-based) used as compare keys, these rows determine
        whether two records are same (intersection) or not (primary or
        secondary), the CHROM and POS rows should always be specified,
        default only use CHROM and POS fields (e.g. 0 1)
        
Other Options:

    --threads  <int>
        how many data threads should be allocated to running this analysis,
        currently only valid while doing clustering [default: 1]
        
    *Note:
        1) As forking new thread is relatively expensive, using multiple
        threads on simple jobs could cost more time, thus this option is only
        suggested while doing processes such as fill gaps, which contains
        steps much more time consuming;
        2) Currently clustering process if much faster, threading is thus
        unnecessary.

    --sum-diagnose
        indicate input vcf file is from GATK DiagnoseTargets, output will be
        summary of those diagnose results
        
    
    --check-context
        check sequence context within or around variants loci
    --extend-snp    <int>
        extend regions to inspect flanking sequence context for SNPs
        [default: 5bp]
    --extend-indel  <int>
        extend regions to inspect flanking sequence context for INDELs
        [default: 10bp]

    *Notes for context checking:
        1) Only bi-allelic loci is supported while analysis sequence context,
        one needs to break multi-alleles into multi bi-allele records before
        using this function;
        2) Extension here is different for SNPs and INDELs, e.g. upstream
        5bp and downstream 5bp for SNPs, while only downstream 10bp for
        INDELs, thus the INDELs are assumed to be already left aligned

EOF

    exit(0);
}

$|++;

print STDERR "# $0 v$VERSION\n# " . (scalar localtime()) . "\n";


if ($options{output}) {
    open (STDOUT, "> $options{output}") || die $!;
}


my $vcf_process = MyPerl::Vcf->new(%options);


if ($options{out_metrics}) {
    if (@{$options{metrics}} < 1) {
        print STDERR "Error: at least 1 metric should be specified!\n";
        exit(2);
    }
    
    ## collect metrics
    print STDERR ">> Start collecting metrics in $options{vcf} ... ";
    $vcf_process->collect_metrics();
    print STDERR "\tdone!\n";
}
elsif ($options{overlap_genotypes}) {
    ## checking overlapped genotypes
    print STDERR ">> Start check overlapping in $options{vcf} ... ";
    $vcf_process->get_overlapped_loci();
    print STDERR "\tdone!\n";
}
elsif ($options{segregating}) {
    ## filter loci segregated in two groups
    print STDERR ">> Start filtering $options{vcf} ... ";
    $vcf_process->find_segregate_loci();
    print STDERR "\tdone!\n";
}
elsif ($options{cluster_markers}) {
    ## cluster markers
    if ($options{type2char_str}) {
        $vcf_process->markers2blocks_str();
    }
    else {
        $vcf_process->markers2blocks();
    }
}
elsif ($options{base_changes}) {
    ## count change of all bases
    print STDERR ">> Start counting $options{vcf} ... ";
    $vcf_process->count_base_changes();
    print STDERR "\tdone!\n";
}
elsif ($options{summary_diagnose}) {
    ## process diagnose results from GATK DiagnoseTargets
    print STDERR ">> Start parsing $options{vcf} ... ";
    $vcf_process->count_diagnose();
    print STDERR "\tdone!\n";
}
elsif ($options{check_context}) {
    ## check sequence context
    print STDERR ">> Start checking $options{vcf} ... ";
    $vcf_process->check_variants_context();
    print STDERR "\tdone!\n";
}
elsif ($options{stat_var_dist}) {
    ## check variants distance
    print STDERR ">> Start counting $options{vcf} ... ";
    $vcf_process->stat_vcf_dist();
    print STDERR "\tdone!\n";
}
elsif ($options{secondary_vcf}) {
    ## combine vcf files
    print STDERR ">> Start combine vcf files ... ";
    $vcf_process->combine_vcfs();
    print STDERR "\tdone!\n";
}
else {
    ## general filtering and processing processes
    print STDERR ">> Start filtering $options{vcf} ... ";
    $vcf_process->filter_vcf();
    print STDERR "\tdone!\n";
}


print STDERR "# " . (scalar localtime()) . "\n";

######################### Sub #########################



