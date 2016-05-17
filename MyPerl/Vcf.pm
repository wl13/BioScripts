#!/usr/bin/perl -w
#
#  Vcf.pm -- Process of vcf file
#
#  Author: Nowind
#  Created: 2014-11-08
#  Updated: 2016-05-16
#  Version: 1.6.8
#
#  Change logs:
#  Version 1.0.0 14/11/08: The initial version.
#  Version 1.0.1 14/11/19: Bug fixed: no results written while "--stats-out" is specified;
#                          add support for "NR" and "NV" tags used in Platypus. 
#  Version 1.1.0 14/11/24: Add function check_variants_context.
#  Version 1.2.0 14/12/12: Complete rewrite of function markers2blocks use new algorithm which
#                          is magnitude faster than previous method with less memory consumption,
#                          also peformed better in determine borders of blocks.
#  Version 1.2.1 14/12/19: Add filtering for alt count.
#  Version 1.2.2 15/01/04: Add filtering for rare variants.
#  Version 1.3.0 15/01/19: Bug fixed in parsing NR fields; Add two functions to process AD infos.
#  Version 1.4.0 15/08/13: Add functions stat_vars_dist and stat_vcf_dist to generate statistics of 
#                          distance between adjacent markers; add more detailed output results of 
#                          marker infos in function markers2blocks.
#  Version 1.5.0 15/08/16: Add funtion to detect variant types; add funtions to export.
#  Version 1.6.0 15/09/26: Add funtion combine_vcfs to combine two vcf files.
#  Version 1.6.1 15/09/28: Bug fixed: contig info lost while writing results.
#  Version 1.6.2 15/09/29: Now one can set inconsistent genotypes to missing without regenotype it.
#  Version 1.6.3 15/10/01: Update function combine_vcfs to support user defined combine tags.
#  Version 1.6.4 15/11/04: Update function combine_vcfs to support user defined combine rows.
#  Version 1.6.5 15/11/13: Updated: move function get_genome_length to MyPerl::FileIO.
#  Version 1.6.6 15/11/29: Bug fixed: error abort due to no available infos present for some contigs.
#  Version 1.6.7 16/05/04: Updated: add support for merge records with same master fields but different
#                          sub fields in function combine_vcfs.
#  Version 1.6.8 16/05/16: Bug fixed in function combine_vcfs:
#                               1) uninitialized value while combine duplicate records;
#                               2) results lost if contig infos missed or incomplete in vcf header



=head1 NAME

MyPerl::Vcf - Local perl module for vcf file processes


=head1 SYNOPSIS

  use MyPerl::Vcf qw(:all);

=head1 DESCRIPTION

Local perl module used for reading, parsing and analysis vcf files

=cut

package MyPerl::Vcf;

use strict;
use Data::Dumper;
use Statistics::Descriptive;

use MyPerl::FileIO qw(:all);

require Exporter;

##
## Global Constants and Variables
##
use vars qw(
  @ISA
  %EXPORT_TAGS
  @EXPORT_OK
  @EXPORT
);

@ISA = qw(Exporter);

%EXPORT_TAGS = (
    'all' => [
        qw(
            filter_vcf
            get_overlapped_loci
            collect_metrics
            find_segregate_loci
            markers2blocks
            count_base_changes
            count_diagnose
            fix_AD_fields
            NRNV2AD
            get_var_type
            check_variants_context
            stat_vcf_dist
            stat_vars_dist
            combine_vcfs
        )
    ]
);

@EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );
@EXPORT    = qw();


$MyPerl::Vcf::VERSION = '1.6.8';


=head1 METHODS

=head2 new

    About   : Create a new object.
    Usage   : my $aln = MyPerl::Align->new(%settings);
    Args    : Using defined settings about:
                input:    input file
                output:   output file
    Returns : An object.

=cut
sub new
{
    my ($class, %options)  = @_;
    
    return bless \%options, $class;
}


=head2 check_var_type

    About   : Return variants types (snp or indels).
    Usage   : check_var_type($ref, $alts);
    Args    : Reference allele;
              Alternative alleles.
    Returns : Variants types.

=cut
sub get_var_type
{
    my ($ref, $alts) = @_;
    
    my @var_types = ();
    for my $alt (split(/,/,$alts))
    {
        my $len = length($ref) - length($alt);
        if ( $len == 0 ) { push @var_types, 'snp'; }
        else { push @var_types, 'indel'; } 
    }
    
    my $var_types = join ',', @var_types;
    
    return $var_types;
}


=head2 filter_vcf

    About   : Filter vcf file accroding to certain criteria.
    Usage   : filter_vcf($opts->{vcf});
    Args    : Vcf file.
    Returns : Null

=cut
sub filter_vcf
{
    my ($opts) = @_;
    

    ##
    ## parse options
    ##
    my %sample_types     = ();
       $sample_types{$_} = 'het' for @{$opts->{het_samples}};
       $sample_types{$_} = 'hom' for @{$opts->{hom_samples}};

    my %remove_samples     = ();
       $remove_samples{$_} = 1 for @{$opts->{samples_remove}};
       $remove_samples{$_} = 0 for @{$opts->{samples_remain}};
        
       $opts->{filter_str} = join '|', @{$opts->{filters}};
    
    
    my %exclude_samples    = ();
       $exclude_samples{$_} = 1 for @{$opts->{samples_exclude}};
    ##
    ## read depth criteria for each sample from  a file
    ##
    my %samples_depth = ();
    if ($opts->{sample_depth_file}) {
        print STDERR ">> Start reading depth range from $opts->{sample_depth_file} ... ";
        my $fh = getInputFilehandle($opts->{sample_depth_file});
        while (<$fh>)
        {
            next if (/^\#/ || /^\s+$/);
        
            my ($sample, $min_depth, $max_depth) = (split /\s+/);
            
            $samples_depth{$sample}->{min} = $min_depth;
            $samples_depth{$sample}->{max} = $max_depth;
        }
        print STDERR "done!\n";
        
        ###print STDERR Dumper(%samples_depth);exit;
    }

    ##
    ## write stats results to a file
    ##
    if ($opts->{stats_outfile}) {
        open (STATS, "> $opts->{stats_outfile}") || die $!;
        print STATS "$opts->{source_line}\n";
    }
    
    if ($opts->{stats_only} && !$opts->{stats_outfile}) { ## redirect the output to STDOUT for pipe
        print STDOUT "$opts->{source_line}\n";
    }
    
    ## file header for loci stats
    if ($opts->{locus_stats}) {
        if ($opts->{stats_outfile}) {
            print STATS "#chrom\tpos\tref_count\talt_count\thet_count\tmissing_count\tlow_depth_count\tover_covered_count\n";
        }
        else {
            print STDOUT "#chrom\tpos\tref_count\talt_count\thet_count\tmissing_count\tlow_depth_count\tover_covered_count\n";
        }
    }
   
    ##
    ## record all loci in vcf file
    ##
    my %records_count = ();
    my $fh = getInputFilehandle($opts->{vcf});
    if ($opts->{unique} || $opts->{random_size} || defined $opts->{min_indel_range}) {
        while (<$fh>)
        {
            next if (/^\#/ || /^\s+$/);
    
            next if ($opts->{skip_unphased} && /\//);         ## only considering phased genotypes
            next if ($opts->{filter_str} && /$opts->{filter_str}/i);
            
            my ($chrom, $pos, $id, $ref, $alts) = (split /\s+/, $_);
            
            my $var_type = get_var_type($ref, $alts);
    
            $records_count{all}->{"$chrom\t$pos"}++;
            
            if ((defined $opts->{min_indel_range}) && ($var_type =~ /indel/)) {
                my $start = $pos + 1 - $opts->{min_indel_range};
                my $end   = $pos + (length $ref) + $opts->{min_indel_range};
            
                $records_count{indel}->{"$chrom\t$_"}++ for ($start..$end);
            }
            
        }
        
        $fh = getInputFilehandle($opts->{vcf});        
    }
    
    
    ##
    ## generate random sets
    ##
    my %random_set = ();
    if ($opts->{random_size}) {
        my $total_loci_num = scalar (keys %{$records_count{all}});
        $opts->{random_size} = ($opts->{random_size} > 1) ? $opts->{random_size} :
                       int($total_loci_num * $opts->{random_size});
                       
        my $random_count = 0;
        while(1)
        {
            my $rand_pos = int(rand($total_loci_num));
            
            $random_set{$rand_pos}++;
            
            last if ((keys %random_set) >= $opts->{random_size});
        }
    }

    
    ##
    ## main process part
    ##
    my %duplicates    = ();
    my @remain_rows   = ();
    my @Samples_ids   = ();
    my %Samples_stats = ();
    while (<$fh>)
    {
        ##
        ## parse vcf header
        ##
        if (/#CHROM/) {
            next if (@Samples_ids > 0);  ## support for multi vcf files to pipe in
            
            my @line   = (split /\s+/);
            
            if (@{$opts->{samples_remain}} > 0) {
                for my $sample_id (@line[9..$#line])
                {
                    unless(defined $remove_samples{$sample_id}) {
                        $remove_samples{$sample_id} = 1;
                    }
                }
            }
            
            @remain_rows = grep { !$remove_samples{$line[$_]} } (0..$#line);
            @Samples_ids = grep { !$remove_samples{$_} } @line[9..$#line];
            
            for my $sample_id (@Samples_ids)
            {
                unless ($sample_types{$sample_id}) {
                    if (@{$opts->{hom_samples}}> 0) {
                        $sample_types{$sample_id} = 'het';
                    }
                    elsif (@{$opts->{het_samples}} > 0) {
                        $sample_types{$sample_id} = 'hom';
                    }
                    else {
                        $sample_types{$sample_id} = $opts->{sample_type};
                    }
                }
            }
            
            my $line = join "\t", @line[@remain_rows];
            
            ###print STDERR Dumper(@remain_rows, $line);exit;
            
            unless ($opts->{stats_only}) {
                if ($opts->{max_rare_count}) {
                print STDOUT <<EOF;
##INFO=<ID=RAREALLELE,Number=1,Type=String,Description="Non-reference allele with a minimum frequency">
##INFO=<ID=RAREFQ,Number=1,Type=Integer,Description="Frequency of rare allele">
##INFO=<ID=RARESAMPLES,Number=1,Type=String,Description="Samples contain the rare allele">                    
EOF
                }
                print STDOUT "$opts->{source_line}\n";
                print STDOUT "$line\n";
            }
            
            ###print STDERR Dumper(%sample_types);exit
            
            next;
        }
        if (/^\#/ || /^\s+$/) {
            next if (@Samples_ids > 0);  ## support for multi vcf files to pipe in
            
            unless($opts->{stats_only}) {
                print STDOUT;
            }
            next;
        }
        
        ###last if ($. > 10000);
        
        next if ($opts->{skip_unphased} && /\//);         ## only considering phased genotypes
        
        my @line   = (split /\s+/);
        my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT) = @line[0..8];
        
        
        next if ($opts->{pass_only} && $FILTER ne 'PASS');
        next if ($opts->{filter_str} && $FILTER =~ /$opts->{filter_str}/i);
        
        
        ##
        ## remove low-quality records
        ##
        next if ($opts->{min_quality} && (($QUAL eq '.') || ($QUAL < $opts->{min_quality})));
        
        ##
        ## filtering by variant type
        ##
        my $var_type = get_var_type($REF, $ALT);
        
        next if ($opts->{var_type} && ($var_type !~ /$opts->{var_type}/));
        
        
        ##
        ## remove records with the same chr:pos to avoid collisions
        ##
        next if ($opts->{unique} && $records_count{all}->{"$CHROM\t$POS"} > 1);
        
        ##
        ## remove snps overlapped with or around indels
        ##
        next if (defined $opts->{min_indel_range} && ($var_type =~ /snp/) && $records_count{indel}->{"$CHROM\t$POS"});
        
        ##
        ## remove duplicate records
        ##
        next if ($duplicates{"$CHROM\t$POS\t$REF\t$ALT"});
                 $duplicates{"$CHROM\t$POS\t$REF\t$ALT"}++;
        
        ##
        ## random sampling
        ##
        next if ($opts->{random_size} && !$random_set{"$CHROM\t$POS"});
        
        ##
        ## filtering by numbers of different alleles
        ##
        my @alleles = ($REF, (split /\,/, $ALT));
        
        next if (($opts->{min_allele_types} && @alleles < $opts->{min_allele_types}) ||
                 ($opts->{max_allele_types} && @alleles > $opts->{max_allele_types}));
        
        
        ##
        ## filtering by indel length
        ##
        my $indel_len = 0;
        for my $allele (@alleles)
        {
            my $len = abs(length($allele) - length($REF));
            
            $indel_len = $len if $len > $indel_len;
        }
        
        next if (defined $opts->{min_indel_len} && $indel_len < $opts->{min_indel_len});
        next if (defined $opts->{max_indel_len} && $indel_len > $opts->{max_indel_len});
        
        
        my @tags = (split /\:/, $FORMAT);
        my %tags = ();
        for (my $i=0; $i<@tags; $i++) { $tags{$tags[$i]} = $i; }
            
        @line = @line[@remain_rows];
        
        
        ##
        ## genotype related processes
        ##
        my $ref_cnt      = 0;
        my $alt_cnt      = 0;
        my $het_cnt      = 0;
        my $pseudo_het   = 0;
        my $hom_ref_cnt  = 0;
        my $het_ref_cnt  = 0;
        my $hom_alt_cnt  = 0;
        my $het_alt_cnt  = 0;
        my $low_dp_cnt   = 0;
        my $high_dp_cnt  = 0;
        my %missing_cnts = ();
        my %allele_freqs = ();
        
        if (defined($tags{GT})) {
            for (my $i=9; $i<@line; $i++)
            {
                my $GT = (split /:/, $line[$i])[$tags{GT}];
                
                my $sample_id = $Samples_ids[$i-9];
                
                if ($GT !~ /((\d)(\/|\|)(\d))/) {
                    $missing_cnts{$sample_types{$sample_id}} ++;
                    $Samples_stats{missing}->{$sample_id} ++;
                    next;
                }
                
                ## parse GT tag
                $GT =~ /((\d)(\/|\|)(\d))/;
                
                my ($geno, $allele1, $separator, $allele2) = ($1, $2, $3, $4);
                
                my $geno_tag = '';
                if ($allele1 == 0 && $allele2 == 0) {
                    $geno_tag = 'ref';
                }
                elsif ($allele1 != $allele2) {
                    $geno_tag = 'het';
                }
                else {
                    $geno_tag = 'alt';
                }
                
                ##
                ## check depth by each sample
                ##
                if (defined $opts->{min_sample_dp} || defined $opts->{max_sample_dp} || $opts->{sample_depth_file}) {
                    my $DP = 0;
                    if ($tags{DP}) {
                        $DP = (split /:/, $line[$i])[$tags{DP}];
                    }
                    elsif ($tags{NR}) {
                        ## NR and NV is two tags used in Platypus
                        ## FORMAT=<ID=NR,Number=.,Type=Integer,Description="Number of reads covering variant location in this sample">
                        my $NRs = (split /:/, $line[$i])[$tags{NR}];
                        my @NRs = (split /\,/, $NRs);
                        
                        $DP = (sort {$a <=> $b} @NRs)[-1]; ## choose the maximum depth in multi-allele locus
                    }
                    elsif ($tags{AD}) {
                        my $AD = (split /:/, $line[$i])[$tags{AD}];
                        
                        my @dps = (split /,/, $AD);
                        
                        $DP += $_ for @dps;
                        
                        $DP = 0 if $AD eq '.';
                    }
                    
                    $DP = 0 if $DP eq '.';
                    
                    my $fail_depth_check = 0;
                    
                    my $min_depth = $opts->{min_sample_dp} ? $opts->{min_sample_dp} : -1;
                    my $max_depth = $opts->{max_sample_dp} ? $opts->{max_sample_dp} : -1;
                    
                    ## if a file is specified, use the depth range from the file in priority
                    if ($opts->{sample_depth_file} && $samples_depth{$sample_id}) {
                        $min_depth = $samples_depth{$sample_id}->{min};
                        $max_depth = $samples_depth{$sample_id}->{max};
                    }
                    
                    if ($min_depth >= 0 && $DP < $min_depth) {
                        $fail_depth_check = 1;
                        $low_dp_cnt ++;
                    }
                    
                    if ($max_depth >= 0 && $DP > $max_depth) {
                        $fail_depth_check = 1;
                        $high_dp_cnt ++;
                    }
                    
                    if ($fail_depth_check) {
                        my $GT_new = '.' . $separator . '.';
                        
                        $line[$i] =~ s/$GT/$GT_new/;
                        
                        $missing_cnts{$sample_types{$sample_id}} ++;
                        $Samples_stats{missing}->{$sample_id} ++;
                        ###print STDERR "$CHROM\t$POS\t$DP\t$line[$i]\n";exit;
                        next;
                    }
                }
                
                ##
                ## treat sample without AD info as missing alleles
                ##
                if ($opts->{check_miss_AD}) {
                    if ($tags{AD}) {
                        my $AD = (split /:/, $line[$i])[$tags{AD}];
                        
                        if (!$AD || ($AD eq '.')) {
                            my $GT_new = '.' . $separator . '.';
                            
                            $line[$i] =~ s/$GT/$GT_new/;
                        
                            $missing_cnts{$sample_types{$sample_id}} ++;
                            $Samples_stats{missing}->{$sample_id} ++;
                            next;
                        }
                    }
                    elsif ($tags{NR} && $tags{NV}) {
                        ## NR and NV is two tags used in Platypus
                        ## FORMAT=<ID=NR,Number=.,Type=Integer,Description="Number of reads covering variant location in this sample">
                        ## FORMAT=<ID=NV,Number=.,Type=Integer,Description="Number of reads containing variant in this sample">
                        ## do nothing if these two tags is present
                    }
                    else {
                        my $GT_new = '.' . $separator . '.';
                        
                        $line[$i] =~ s/$GT/$GT_new/;
                        
                        $missing_cnts{$sample_types{$sample_id}} ++;
                        $Samples_stats{missing}->{$sample_id} ++;
                        next;
                    }
                }
                
                ##
                ## genotyping processes
                ##
                if ($opts->{genotype_stats} || defined $opts->{hom_threshold} || defined $opts->{het_threshold}) {
                    
                    ##
                    ## get read depth for each allele
                    ##
                    my @dps = ();
                    
                    if ($tags{AD}) {
                        ## use AD tag if present
                        my $AD = (split /:/, $line[$i])[$tags{AD}];
                        
                        if ($AD && ($AD ne '.')) {
                            @dps = (split /\,/, $AD);
                        }
                    }
                    elsif ($tags{NR} && $tags{NV} && (@alleles == 2)) {
                        ## use NR and NV tags if present
                        ## did not support multi-allele loci
                        my $NR  = (split /:/, $line[$i])[$tags{NR}];
                        my $NV  = (split /:/, $line[$i])[$tags{NV}];
                        
                        @dps = ($NR - $NV, $NV);
                        
                        ###my @NVs = (split /\,/, $NV);
                        ###$dp_ref -= $_ for @NVs;
                        ###@dps = ($dp_ref, @NVs);
                        
                        ###print STDERR "$CHROM\t$POS\t$sample_id\t@dps\n"; exit;
                    }
                    
                    if (@dps > 1) {
                        my $dp_sum = 0;
                           $dp_sum += $_ for @dps;
                        
                        my $ref_dp       = $dps[0];
                        my @allele_sort  = sort {$dps[$b] <=> $dps[$a]} (0..$#dps);
                        my $major_allele = $allele_sort[0];
                        my $minor_allele = $allele_sort[1];
                        my $minor_dp     = $dps[$minor_allele];
                        
                        my $ref_perc   = ($dp_sum > 0) ? (100 * $ref_dp / $dp_sum) : -1;
                        my $minor_perc = ($dp_sum > 0) ? (100 * $minor_dp / $dp_sum) : -1;
                        
                        my $bin  = $opts->{use_ref_depth} ? int(100 - $ref_perc) : int($minor_perc);
                        
                        $Samples_stats{GT}->{$geno_tag}->{$bin}->{$sample_id} ++;
                        $Samples_stats{GT}->{All}->{$bin}->{$sample_id} ++;
                        
                        ##
                        ## reverse genotyped alleles such as "GT:AD  0/0:0,2"
                        ##
                        if ($opts->{check_gt_dp}) {
                            if (int($ref_perc) == 100 && ($geno_tag eq 'alt')) {
                                my $GT_new = '.' . $separator . '.';
                                
                                $line[$i] =~ s/$GT/$GT_new/;
                                
                                $missing_cnts{$sample_types{$sample_id}} ++;
                                $Samples_stats{missing}->{$sample_id} ++;
                                next;
                            }
                            
                            if (int($ref_perc) == 0 && ($geno_tag eq 'ref')) {
                                my $GT_new = '.' . $separator . '.';
                                
                                $line[$i] =~ s/$GT/$GT_new/;
                                
                                $missing_cnts{$sample_types{$sample_id}} ++;
                                $Samples_stats{missing}->{$sample_id} ++;
                                next;
                            }
                        }

                    
                        ##
                        ## regenotype variants in homozygous samples
                        ##
                        if ((defined $opts->{hom_threshold}) && ($sample_types{$sample_id} eq 'hom')) {
                            my $regenotyping = 0;
                            
                            if (($geno_tag eq 'het') && ($minor_perc >= 0 && $minor_perc <= 100 * $opts->{hom_threshold})) {
                                ###print STDERR "prev:$CHROM\t$POS\t$line[$i]\n";
                                ## heterozygous => homozygous
                                my $GT_new = $major_allele . $separator . $major_allele;
                                
                                if ($GT_new =~ /(0(\/|\|)0)/) {
                                    $geno_tag = 'ref';
                                }
                                else {
                                    $geno_tag = 'alt';
                                }
                                
                                $line[$i] =~ s/$GT/$GT_new/;
                                
                                $GT = $GT_new;  ## update GT fields
                                ###print STDERR "curr:$CHROM\t$POS\t$line[$i]\n";
                                
                                $regenotyping = 1;
                            }
                            elsif (($geno_tag ne 'het') && ($minor_perc > 100 * $opts->{hom_threshold})) {
                                ###print STDERR "prev:$CHROM\t$POS\t$line[$i]\n";
                                ## homozygous => heterozygous
                                my $GT_new = ($minor_allele > $major_allele) ?
                                             ($major_allele . $separator . $minor_allele) :
                                             ($minor_allele . $separator . $major_allele);
                                
                                $geno_tag = 'het';
                                
                                $line[$i] =~ s/$GT/$GT_new/;
                                
                                $GT = $GT_new;  ## update GT fields
                                ###print STDERR "curr:$CHROM\t$POS\t$line[$i]\n";
                                
                                $regenotyping = 1;
                            }
                            
                            ## set as missing if $minor_perc < 0 or --gt-diff-as-missing option is specified
                            if ($minor_perc < 0 || ($regenotyping && $opts->{gt_diff_as_missing})) {
                                my $GT_new = '.' . $separator . '.';
                                
                                $line[$i] =~ s/$GT/$GT_new/;
                                
                                $missing_cnts{$sample_types{$sample_id}} ++;
                                $Samples_stats{missing}->{$sample_id} ++;
                                next;
                            }
                            
                            $Samples_stats{GT_new}->{$geno_tag}->{$bin}->{$sample_id} ++;
                            $Samples_stats{GT_new}->{All}->{$bin}->{$sample_id} ++;
                        }
                        
                        
                        ##
                        ## regenotype variants in heterozygous samples
                        ##
                        if ((defined $opts->{het_threshold}) && ($sample_types{$sample_id} eq 'het')) {
                            my ($het_to_hom, $hom_to_het) = (split /\,/, $opts->{het_threshold});
                            my $regenotyping = 0;
                            
                            if (($geno_tag eq 'het') && ($minor_perc >= 0 && $minor_perc <= 100 * $het_to_hom)) {
                                ## heterozygous => homozygous
                                my $GT_new = $major_allele . $separator . $major_allele;
                                
                                if ($GT_new =~ /(0(\/|\|)0)/) {
                                    $geno_tag = 'ref';
                                }
                                else {
                                    $geno_tag = 'alt';
                                }
                                
                                $line[$i] =~ s/$GT/$GT_new/;
                                
                                $GT = $GT_new;  ## update GT fields
                                
                                $regenotyping = 1;
                            }
                            elsif (($geno_tag ne 'het') && ($minor_perc >= 100 * $hom_to_het)) {
                                ## homozygous => heterozygous
                                my $GT_new = ($minor_allele > $major_allele) ?
                                             ($major_allele . $separator . $minor_allele) :
                                             ($minor_allele . $separator . $major_allele);
                                
                                $geno_tag = 'het';
                                
                                $line[$i] =~ s/$GT/$GT_new/;
                                
                                $GT = $GT_new;  ## update GT fields
                                
                                $regenotyping = 1;
                            }
                            
                            ## set as missing if $minor_perc < 0 or --gt-diff-as-missing option is specified
                            if (($minor_perc < 0) || ($minor_perc > 100 * $het_to_hom && $minor_perc < 100 * $hom_to_het)
                                                  || ($regenotyping && $opts->{gt_diff_as_missing})) {
                                my $GT_new = '.' . $separator . '.';
                                
                                $line[$i] =~ s/$GT/$GT_new/;
                                
                                $missing_cnts{$sample_types{$sample_id}} ++;
                                $Samples_stats{missing}->{$sample_id} ++;
                                next;
                            }
                            
                            $Samples_stats{GT_new}->{$geno_tag}->{$bin}->{$sample_id} ++;
                            $Samples_stats{GT_new}->{All}->{$bin}->{$sample_id} ++;
                        }
                    }
                }
                
                ## set pseudo-hete allele to missing
                if ($opts->{pseduo_het_as_missing} && ($sample_types{$sample_id} eq 'hom') && ($geno_tag eq 'het')) {
                    my $GT_new = '.' . $separator . '.';
                    
                    $line[$i] =~ s/$GT/$GT_new/;
                    
                    $missing_cnts{$sample_types{$sample_id}} ++;
                    $Samples_stats{missing}->{$sample_id} ++;
                    next;
                }
                
                
                $GT =~ /((\d)(\/|\|)(\d))/;
                
                ($geno, $allele1, $separator, $allele2) = ($1, $2, $3, $4);
                
                if ($allele1 == 0 && $allele2 == 0) {
                    $ref_cnt ++;
                    if ($sample_types{$sample_id} eq 'hom') {
                        $hom_ref_cnt ++;
                    }
                    else {
                        $het_ref_cnt ++;
                    }
                }
                elsif ($allele1 != $allele2) {
                    $het_cnt ++;
                    if ($sample_types{$sample_id} eq 'hom') {
                        $pseudo_het ++;
                    }
                }
                else {
                    $alt_cnt ++;
                }
                
                ## count allele frequency
                push @{$allele_freqs{$allele1}}, $sample_id;
                if ($allele1 ne $allele2) {
                    push @{$allele_freqs{$allele2}}, $sample_id;
                }
            }
            
        }
        
        ##
        ## filtering by non reference allele counts
        ##
        my $non_ref_cnt = ($het_cnt + $alt_cnt);
        my $hom_missing = $missing_cnts{hom} ? $missing_cnts{hom} : 0;
        my $het_missing = $missing_cnts{het} ? $missing_cnts{het} : 0;
        my $missing_cnt = $hom_missing + $het_missing;
        
        if ($opts->{missing_as_ref}) {
            $hom_ref_cnt += $hom_missing;
            $het_ref_cnt += $het_missing;
            $ref_cnt     += $missing_cnt;
        }
        
        if ($opts->{locus_stats}) {
            if ($opts->{stats_outfile}) {
                print STATS "$CHROM\t$POS\t$ref_cnt\t$alt_cnt\t$het_cnt\t$missing_cnt\t$low_dp_cnt\t$high_dp_cnt\n";
            }
            else {
                print STDOUT "$CHROM\t$POS\t$ref_cnt\t$alt_cnt\t$het_cnt\t$missing_cnt\t$low_dp_cnt\t$high_dp_cnt\n";
            }
            
        }
            
        if ($opts->{stats_only}) {
            next;
        }
        
        
        ###print STDERR;
        ###print STDERR "ref:$ref_cnt\thom_ref:$hom_ref_cnt\tnon-ref:$non_ref_cnt\n";exit;
        
        ##
        ## filtering by various frequency
        ##
        next if (defined $opts->{min_ref_count} && ($ref_cnt < $opts->{min_ref_count}));
        next if (defined $opts->{max_ref_count} && ($ref_cnt > $opts->{max_ref_count}));
        next if (defined $opts->{max_het_count} && ($het_cnt > $opts->{max_het_count}));
        next if (defined $opts->{min_het_count} && ($het_cnt < $opts->{min_het_count}));
        next if (defined $opts->{min_alt_count} && ($alt_cnt < $opts->{min_alt_count}));
        next if (defined $opts->{max_alt_count} && ($alt_cnt > $opts->{max_alt_count}));
        
        next if (defined $opts->{min_hom_ref_count} && ($hom_ref_cnt < $opts->{min_hom_ref_count}));
        next if (defined $opts->{max_hom_ref_count} && ($hom_ref_cnt > $opts->{max_hom_ref_count}));
        next if (defined $opts->{min_het_ref_count} && ($het_ref_cnt < $opts->{min_het_ref_count}));
        next if (defined $opts->{max_het_ref_count} && ($het_ref_cnt > $opts->{max_het_ref_count}));
        
        next if (defined $opts->{min_non_ref_count} && ($non_ref_cnt < $opts->{min_non_ref_count}));
        next if (defined $opts->{max_non_ref_count} && ($non_ref_cnt > $opts->{max_non_ref_count}));
        
        next if (defined $opts->{max_hom_missing} && ($hom_missing > $opts->{max_hom_missing}));
        next if (defined $opts->{max_het_missing} && ($het_missing > $opts->{max_het_missing}));
        next if (defined $opts->{max_missing}     && ($missing_cnt > $opts->{max_missing}));
        
        next if (defined $opts->{min_pseduo_het} && ($pseudo_het < $opts->{min_pseduo_het}));
        next if (defined $opts->{max_pseduo_het} && ($pseudo_het > $opts->{max_pseduo_het}));
        
        
        if ($opts->{remove_info}) {
            $line[7] = '.';
        }
        
        ##
        ## screen out rare alleles
        ##
        if (defined $opts->{max_rare_count}) {
            my @alleles_order = sort { (scalar @{$allele_freqs{$a}}) <=>
                                       (scalar @{$allele_freqs{$b}}) } (keys %allele_freqs);
            my $min_allele = $alleles_order[0];
            
            ###unless (defined $min_allele) {
            ###    print STDERR "$CHROM\t$POS\tmiss:$missing_cnt\talt:$alt_cnt\thet:$het_cnt\tref:$ref_cnt\t$min_allele\n";
            ###    print STDERR Dumper(%allele_freqs); exit;
            ###}
            
            next unless((defined $min_allele) && $allele_freqs{$min_allele});
            
            next if ($opts->{no_rare_ref} && ($min_allele == 0));
            
            my $rare_count = scalar @{$allele_freqs{$min_allele}};
            
            next if ($rare_count > $opts->{max_rare_count});
            
            
            ##
            ## remove loci shared with non-desired samples
            ##
            my $is_exclude = 0;
            
            for my $sample (@{$allele_freqs{$min_allele}})
            {
                if ($exclude_samples{$sample}) {
                    $is_exclude = 1; last;
                }
            }
            
            next if ($is_exclude);
            
            my $rare_info = "RAREALLELE=$alleles[$min_allele];RAREFQ=$rare_count;RARESAMPLES="
                           . join ',', (@{$allele_freqs{$min_allele}});
            
            if ($opts->{remove_info}) {
                $line[7] = "$rare_info";
            }
            else {
                $line[7] .= ";$rare_info";
            }
        }
        
        my $line = join "\t", @line;
        print STDOUT "$line\n";  
    }
    
    
    if ($opts->{genotype_stats}) {
        my $out_samples = join "\t", @Samples_ids;
        
        if ($opts->{stats_outfile}) {
            print STATS "##DP_Minor(%): Percentage of minor allele depth, 0 means [0%, 1%), etc. " .
                        "while -1 means no reliable reads found for both reference and non-reference alleles, e.g. AD=0,0\n";
            print STATS "#Geno\tDP_Minor(%)\t$out_samples\n";
        }
        else {
            print STDOUT "##DP_Minor(%): Percentage of minor allele depth, 0 means [0%, 1%), etc. " .
                         "while -1 means no reliable reads found for both reference and non-reference alleles, e.g. AD=0,0\n";
            print STDOUT "#Geno\tDP_Minor(%)\t$out_samples\n";
        }

        
        
        ###print STATS Dumper($Samples_stats{GT});
        
        ## calculate minor allele depth ratio using each sample's AD values,
        ## note only reads with mapping quality upon a certain value was
        ## counted in AD fields
        my $max_range = $opts->{use_ref_depth} ? 100 : 50;
        for my $tag (qw(GT GT_new))
        {
            next unless($Samples_stats{$tag});
            
            for my $geno (qw(ref alt het))
            {
                for my $perc (-1..$max_range)
                {
                    my @freqs = ();
                    for my $sample (@Samples_ids)
                    {
                        my $freq = $Samples_stats{$tag}->{$geno}->{$perc}->{$sample} ?
                                   $Samples_stats{$tag}->{$geno}->{$perc}->{$sample} : 0;
                        push @freqs, $freq;
                    }
                    
                    my $freqs = join "\t", @freqs;
                    
                    if ($opts->{stats_outfile}) {
                        print STATS "$tag:$geno\t$perc\t$freqs\n";
                    }
                    else {
                        print STDOUT "$tag:$geno\t$perc\t$freqs\n";
                    }
                    
                }
            }
        }
        
        if ($opts->{stats_outfile}) {
            print STATS "#Tag\tsample\tmissing\n";
        }
        else {
            print STDOUT "#Tag\tsample\tmissing\n";
        }
        
        for my $sample (@Samples_ids)
        {
            my $sample_missing = $Samples_stats{missing}->{$sample} ?
                                 $Samples_stats{missing}->{$sample} : 0;
                                 
            if ($opts->{stats_outfile}) {
                print STATS "MISSING\t$sample\t$sample_missing\n";
            }
            else {
                print STDOUT "MISSING\t$sample\t$sample_missing\n";
            }
        }
        
        close STATS;
    }
}



=head2 get_overlapped_loci

    About   : Screen out loci with samples which all could be overlapped in
              all libraries.
    Usage   : get_overlapped_loci(\%options);
    Args    : Input vcf file;
              Source line to be written to output vcf header (optional).
    Returns : Null

=cut
sub get_overlapped_loci
{
    my ($opts) = @_;
    
    my $fh = getInputFilehandle($opts->{vcf});
    my @Samples_ids = ();
    my @all_ids     = ();
    while (<$fh>)
    {
        if (/^\#CHROM/) {
            next if (@Samples_ids > 0);  ## support for multi vcf files to pipe in
            
            my @line   = (split /\s+/);
            
            @all_ids = @line[9..$#line];
            
            for (my $i=9; $i<@line; $i+=3)
            {
                my ($id, $tag) = split /\./, $line[$i];
                push @Samples_ids, $id;
            }
            
            my $line = join "\t", (@line[0..8], @Samples_ids);
            
            print STDOUT "$opts->{source_line}\n";
            print STDOUT "$line\n";
            
            next;
        }
        if (/^\#/ || /^\s+$/) {
            next if (@Samples_ids > 0);  ## support for multi vcf files to pipe in
            
            print STDOUT;
            next;
        }
        
        my @line   = (split /\s+/);
        my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT) = @line[0..8];
        

        my @tags = (split /\:/, $FORMAT);
        my %tags = ();
        for (my $i=0; $i<@tags; $i++) { $tags{$tags[$i]} = $i; }
            
        my %sample_gts = ();
        my @out_rows   = ();
        for (my $i=9; $i<@line; $i++)
        {
            my $GT = (split /:/, $line[$i])[$tags{GT}];
            
            my ($id, $tag) = split /\./, $all_ids[$i-9];
            
            $sample_gts{$id}->{$GT} ++;
            
            if ($tag eq 'sample') {
                push @out_rows, $line[$i];
            }
        }
        
        my $discordant = 0;
        for my $id (keys %sample_gts)
        {
            if (scalar (keys %{$sample_gts{$id}}) != 1) {
                $discordant = 1;
            }
        }
        
        next if ($discordant);
        
        my $line = join "\t", (@line[0..8], @out_rows);
        print STDOUT "$line\n";  
    }  
}



=head2 collect_metrics

    About   : Collect multiple metics for variants analysis
    Usage   : collect_metrics(\%options);
    Args    : Vcf file needs to be processed;
              Query metrics;
              Other options.
    Returns : Null

=cut
sub collect_metrics
{
    my ($opts) = @_;
    
    my @query_metrics = @{$opts->{metrics}};
    
    my $fh = getInputFilehandle($opts->{vcf});
    my %variants_metrics  = ();
    while (<$fh>)
    {
        if (/^\#/ || /^\s+$/) {
            next;
        }
        
        my @line   = (split /\s+/);
        my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, @SAMPLES) = @line;
        
        ##
        ## filtering by variant type
        ##
        my $var_type = get_var_type($REF, $ALT);
        
        next if ($opts->{var_type} && ($var_type !~ /$opts->{var_type}/));
        
        for my $metric (@query_metrics)
        {
            if ($INFO =~ m/$metric=(.*?)(\;|$)/x) {
                ###print STDERR "$metric\t$1\n";
                
                push @{$variants_metrics{$metric}}, $1;
            }
        }
        
        ###last if $.>10000;
    }
    
    my %metrics_stats = ();
    for my $metric (@query_metrics)
    {
        my @values = @{$variants_metrics{$metric}};
        
        if (@values > 0) {
            my $stat = Statistics::Descriptive::Full->new();
               $stat->add_data(@values);
            
            my $count  = $stat->count();
            my $mean   = $stat->mean();
            my $min    = $stat->min();
            my $max    = $stat->max();
            my $stdev  = $stat->standard_deviation();
            
            push @{$metrics_stats{count}}, $count;
            push @{$metrics_stats{mean}}, $mean;
            push @{$metrics_stats{min}}, $min;
            push @{$metrics_stats{max}}, $max;
            push @{$metrics_stats{stdev}}, $stdev;
            
            my @counts = ();
            for (my $i=1; $i<=100; $i++)
            {
                my $perc_cnt = $stat->percentile($i);
                push @{$metrics_stats{percentile}->{$i}}, $perc_cnt;
            }
        }
        else {
            push @{$metrics_stats{count}}, 0;
            push @{$metrics_stats{mean}}, 0;
            push @{$metrics_stats{min}}, 0;
            push @{$metrics_stats{max}}, 0;
            push @{$metrics_stats{stdev}}, 0;
            
            my @counts = ();
            for (my $i=0; $i<=100; $i++)
            {
                push @{$metrics_stats{percentile}->{$i}}, 0;
            }
        }
    }
    
    my $out_metrics = join "\t", @query_metrics;
    
    
    print STDOUT "$opts->{source_line}\n";
    print STDOUT "#Stats\t$out_metrics\n";
    for my $key (qw(count mean min max stdev))
    {
        my $metrics_vals = join "\t", @{$metrics_stats{$key}};
        print "$key\t$metrics_vals\n";
    }
    
    print STDOUT "#Percentile\t$out_metrics\n";
    for my $i (sort {$a <=> $b} keys %{$metrics_stats{percentile}})
    {
        my $perc_vals = join "\t", @{$metrics_stats{percentile}->{$i}};
        print STDOUT "$i\t$perc_vals\n";
    }
    
}





=head2 find_segregate_loci

    About   : Screen out loci 
    Usage   : find_segregate_loci($opts->{vcf});
    Args    : Vcf file.
    Returns : Null

=cut
sub find_segregate_loci
{
    my ($opts) = @_;
    
    my ($group1, $group2) = (split /\;/, $opts->{segregating});
    
    my @group1_members = split /\,/, $group1;
    my @group2_members = split /\,/, $group2;
    
    my %groups = ();
       $groups{$_} = 1 for @group1_members;
       $groups{$_} = 2 for @group2_members;
    
    my @group1_rows = ();
    my @group2_rows = ();
    my @other_rows  = ();
    
    my @symbols = split /\;/, $opts->{out_symbols} if ($opts->{out_symbols});
    
    my $fh = getInputFilehandle($opts->{vcf});
    my @Samples_ids = ();
    while (<$fh>)
    {
        if (/#CHROM/) {
            next if (@Samples_ids > 0);  ## support for multi vcf files to pipe in
            
            my @line   = (split /\s+/);
            
            @Samples_ids = @line[9..$#line];
            
            @group1_rows = grep { $groups{$line[$_]} && $groups{$line[$_]} == 1 } (9..$#line);
            @group2_rows = grep { $groups{$line[$_]} && $groups{$line[$_]} == 2 } (9..$#line);
            @other_rows  = grep { !$groups{$line[$_]} } (9..$#line);
            
            unless(@group1_rows > 0) {
                print STDERR "Error: no specified samples found in the first group, exit!\n"; exit(0);
            }
            
            my $line = join "\t", @line;
            
            if ($opts->{out_symbols}) {
                $line = join "\t", @line[0..8, @other_rows];
                
            print STDOUT <<EOF;
##INFO=<ID=$symbols[0],Number=1,Type=String,Description="Corresponding variant of the symbol">
##INFO=<ID=$symbols[1],Number=1,Type=String,Description="Corresponding variant of the symbol">
##FORMAT=<ID=SC,Number=1,Type=String,Description="Allele source">
EOF
            }
            
            print STDOUT "$opts->{source_line}\n";
            print STDOUT "$line\n";
            next;
        }
        if (/^\#/ || /^\s+$/) {
            next if (@Samples_ids > 0);  ## support for multi vcf files to pipe in
            
            print STDOUT;
            next;
        }
        
        my @line   = (split /\s+/);
        my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT) = @line[0..8];
        

        my @tags = (split /\:/, $FORMAT);
        my %tags = ();
        for (my $i=0; $i<@tags; $i++) { $tags{$tags[$i]} = $i; }
        
        my @alleles = ($REF, (split /\,/, $ALT));
        
        my @group1_gts = @line[@group1_rows];
        
        my %group1_gts = ();
        for (@group1_gts)
        {
            my $GT = (split /:/, $_)[$tags{GT}];
            $group1_gts{$GT} ++ 
        }
        my $group1_gt = (keys %group1_gts)[0];
        my $group2_gt = '';
        
        if (@group2_rows > 0) {
            my @group2_gts = @line[@group2_rows];
            my %group2_gts = ();
            for (@group2_gts)
            {
                my $GT = (split /:/, $_)[$tags{GT}];
                $group2_gts{$GT} ++ 
            }
            
            next if ((scalar (keys %group1_gts) > 1) || (scalar (keys %group2_gts) > 1)); ## not fixed in two groups
            
            $group2_gt = (keys %group2_gts)[0];
            
            next if ($group1_gt eq $group2_gt); ## not segregated in two groups
        }
        
        
        my $line = join "\t", @line;
        
        if ($opts->{out_symbols}) {
            my $group1_allele = (split /\/|\|/, $group1_gt)[1];
            my $group2_allele = (@group2_rows > 0) ? (split /\/|\|/, $group2_gt)[1] : (1-$group1_allele);
            
            for (my $i=0; $i<@other_rows; $i++)
            {
                my $GT = (split /:/, $line[$other_rows[$i]])[$tags{GT}];
                
                if ($GT =~ /((\d)(\/|\|)(\d))/)
                {
                    my ($geno, $allele1, $separator, $allele2) = ($1, $2, $3, $4);
                    
                    my $allele1_source = 0;
                    if ($allele1 == $group2_allele) {
                        $allele1_source = 1;
                    }
                    
                    my $allele2_source = 0;
                    if ($allele2 == $group2_allele) {
                        $allele2_source = 1;
                    }
                    
                    if ($separator eq '/') { ## sort allele if non-phased
                        $line[$other_rows[$i]] = join "$separator", (sort ($symbols[$allele1_source], $symbols[$allele2_source]));
                    }
                    else {
                        $line[$other_rows[$i]] = join "$separator", ($symbols[$allele1_source], $symbols[$allele2_source]);
                    }
                }
                else {
                    $line[$other_rows[$i]] = './.';
                }
            }
            
            $line[7] = "$symbols[0]=$alleles[$group1_allele];$symbols[1]=$alleles[$group2_allele]";
            $line[8] = 'SC';
            $line = (join "\t", @line[0..8, @other_rows]);
        }
        
        print STDOUT "$line\n";  
    }  
}



=head2 stat_vars_dist

    About   : Stats of distance between two adjacent variant sites.
    Usage   : $dist_stats = stat_vars_dist(\@vars_pos);
    Args    : Array reference of variants positions within certain range.
    Returns : String of stat results.

=cut
sub stat_vars_dist
{
    my ($ra_marker_pos) = @_;

    if (scalar @{$ra_marker_pos} < 2) {
        return "-,-,-,-,-,-";
    }
    
    my @marker_dist = ();
    for (my $i=1; $i<@{$ra_marker_pos}; $i++)
    {
        push @marker_dist, $ra_marker_pos->[$i] - $ra_marker_pos->[$i-1] - 1;
        ###print "$chrom\t$sample\t$marker_pos[$i-1]\t$marker_pos[$i]\n";
    }
    ###print Dumper(@marker_dist);
    
    my $stat = Statistics::Descriptive::Full->new();
       $stat->add_data(@marker_dist);
    
    my $count  = $stat->count();
    my $mean   = sprintf("%.2f", $stat->mean());
    my $median = $stat->median();
    my $min    = $stat->min();
    my $max    = $stat->max();
    my $stdev  = sprintf("%.2f", $stat->standard_deviation());
    my $stderr = sprintf("%.2f", $stdev / sqrt($count));
    
    return "$min,$median,$mean,$max,$stdev,$stderr";
    ###print "$sample\t$chrom\t($mean,$median,$stdev)[$min,$max]\n";
}



=head2 stat_vcf_dist

    About   : Stats of distance between two adjacent variant sites.
    Usage   : $dist_stats = stat_vars_dist(\@vars_pos);
    Args    : Array reference of variants positions within certain range.
    Returns : String of stat results.

=cut
sub stat_vcf_dist
{
    my ($opts) = @_;
    
    my %sample_vars_all = ();
    
    my $fh = getInputFilehandle($opts->{vcf});
    my @Samples_ids = ();
    while (<$fh>)
    {
        next if (/\#\#/ || /^\s+$/);
        
        if (/#CHROM/) {
            next if (@Samples_ids > 0);  ## support for multi vcf files to pipe in
            
            my @line   = (split /\s+/);
            
            @Samples_ids = @line[9..$#line];
            
            next;
        }
        
        ###print STDERR; exit;
        my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, @SAMPLES) = (split /\s+/);
        
        my @tags = (split /\:/, $FORMAT);
        my %tags = ();
        for (my $i=0; $i<@tags; $i++) { $tags{$tags[$i]} = $i; }
        
        for (my $i=0; $i<@SAMPLES; $i++)
        {
            my $type = $opts->{source_tag} ? ((split /:/, $SAMPLES[$i])[$tags{$opts->{source_tag}}]) : $SAMPLES[$i];
            
            next if($type eq './.' || $type eq '.');
            
            push @{$sample_vars_all{$Samples_ids[$i]}->{$CHROM}}, $POS;
        }
    }
    
    print STDOUT "$opts->{source_line}\n";
    print STDOUT "#sample\tchrom\tmin\tmedian\tmean\tmax\tstdev\tstderr\n";
    for my $sample (sort keys %sample_vars_all)
    {
        for my $chrom (sort keys %{$sample_vars_all{$sample}})
        {
            if (exists $sample_vars_all{$sample}->{$chrom}) {
                my $dist_stat = stat_vars_dist($sample_vars_all{$sample}->{$chrom});
                   $dist_stat =~ s/\,/\t/g;
                print STDOUT "$sample\t$chrom\t$dist_stat\n";
            }
            else {
                print STDOUT "$sample\t$chrom\t-\t-\t-\t-\t-\t-\n";
            }
        }
    }
}



=head2 markers2blocks

    About   : Link markers into large blocks
    Usage   : markers2blocks($opts->{vcf});
    Args    : Vcf file.
    Returns : Null

=cut
sub markers2blocks
{
    my ($opts) = @_;
    
    ##
    ## get ids and lengths of each chromosome from a file
    ##
    my @CHROM_IDs     = ();
    my %CHROM_LENGTHs = ();
    if ($opts->{length_file}) {
        print STDERR ">> Start parsing $opts->{length_file} ... ";
        get_genome_length(\@CHROM_IDs, \%CHROM_LENGTHs, $opts->{length_file});
        print STDERR "done!\n";
    }

    print STDERR ">> Start parsing $opts->{vcf} ... ";
    my @Samples_ids = ();
    my %MARKERS     = ();
    my $fh = getInputFilehandle($opts->{vcf});
    while (<$fh>)
    {
        if (/#CHROM/) {
            next if (@Samples_ids > 0);  ## support for multi vcf files to pipe in
            
            my @lines = (split /\s+/);
            
            @Samples_ids = @lines[9..$#lines];
        }
        elsif (/\#\#contig=<ID=(.*?),length=(\d+)>/ && !$opts->{length_file}) {
            next if (@Samples_ids > 0);  ## support for multi vcf files to pipe in
            
            push @CHROM_IDs, $1;
            $CHROM_LENGTHs{$1} = $2; 
        }
        
        next if (/^\#/ || /^\s+$/);
        
        next if ($opts->{filter_str} && /$opts->{filter_str}/i); 
        
        my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, 
            $FILTER, $INFO, $FORMAT, @SAMPLES) = (split /\s+/);

        ##
        ## filtering by variant type
        ##
        my $var_type = get_var_type($REF, $ALT);
        
        next if ($opts->{var_type} && ($var_type !~ /$opts->{var_type}/));
            
        my @tags = (split /\:/, $FORMAT);
        my %tags = ();
        for (my $i=0; $i<@tags; $i++) { $tags{$tags[$i]} = $i; }
        
        for (my $i=0; $i<@SAMPLES; $i++)
        {
            my $type = $opts->{source_tag} ? ((split /:/, $SAMPLES[$i])[$tags{$opts->{source_tag}}]) : $SAMPLES[$i];
            
            next if($type eq './.' || $type eq '.');
            
            my $prev_type = ($MARKERS{$Samples_ids[$i]}->{$CHROM}->{prev_type}) ?
                             $MARKERS{$Samples_ids[$i]}->{$CHROM}->{prev_type} : '';
            
            push @{$MARKERS{$Samples_ids[$i]}->{$CHROM}->{pos_all}}, $POS;
            
            if ($type eq $prev_type) {
                $MARKERS{$Samples_ids[$i]}->{$CHROM}->{linked}->[-1]->[-3] = $POS;
                
                $MARKERS{$Samples_ids[$i]}->{$CHROM}->{prev_count} ++;
                $MARKERS{$Samples_ids[$i]}->{$CHROM}->{linked}->[-1]->[-1] =
                $MARKERS{$Samples_ids[$i]}->{$CHROM}->{prev_count};
            }
            else {
                push @{$MARKERS{$Samples_ids[$i]}->{$CHROM}->{linked}}, [$POS, $POS, $type, 1];
                
                $MARKERS{$Samples_ids[$i]}->{$CHROM}->{prev_type}  = $type;
                $MARKERS{$Samples_ids[$i]}->{$CHROM}->{prev_count} = 1;
            }
        }
    }
    print STDERR "done!\n";
    
    ###print Dumper(%MARKERS);exit;
    
    #print STDERR ">> Start process blocks ... \n";
    
    ##
    ## get all results from each childs
    ##
    my %final_blocks = ();
    my $pm = new Parallel::ForkManager($opts->{threads});
    if ($opts->{threads} > 1) {
        $pm->run_on_finish(
            sub{
                my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data) = @_;
                
                my ($sample, $chrom, $ra_blocks) = @{$data};
                
                $final_blocks{$sample}->{$chrom} = $ra_blocks;
            }
        );
    }
    
    my $total_jobs  = 0;
    my $curr_job_id = 0;
    for my $sample (@Samples_ids)
    {
        for my $chrom (sort keys %{$MARKERS{$sample}})
        {
            $total_jobs ++;
        }
    }
    
    print STDOUT "$opts->{source_line}\n";
    if ($opts->{fill_gaps}) {
        print STDOUT "##Format of Marker_Details: Total_Markers;Marker_Density;Marker_Stats;Marker_Consists;Major_Perc\n";
        print STDOUT "##Total_Markers: total number of markers\n";
        print STDOUT "##Marker_Density(x1000): (number of markers / linked block length), (number of markers / extened block length)\n";
        print STDOUT "##Marker_Stats: stats of distance between adjacent markers in the order min,median,mean,max,stdev,stderr\n";
        print STDOUT "##Marker_Consists: (Marker type : Number of markers in this type)\n";
        print STDOUT "##Major_Perc(%): (Number of markers support the block type / Total markers)\n";
        print STDOUT "#Sample\tChrom\tType\tStart\tEnd\tLength\tMarker_Start\tMarker_End\tLength\tMarker_Details\n";
    }
    else {
        print STDOUT "#Sample\tChrom\tType\tStart\tEnd\tLength\tMarkers_Count\n";
    }
    for my $sample (@Samples_ids)
    {
        for my $chrom (sort keys %{$MARKERS{$sample}})
        {
            $curr_job_id ++;
            
            my $pid = $pm->start and next if ($opts->{threads} > 1);
            
            print STDERR "\r>> Start initial processing $sample:$chrom ... $curr_job_id/$total_jobs";
            
            #########################################################################
            ## child processes
            my @block_results = ();
            
            ##
            ## gap filling processes
            ##
            if ($opts->{fill_gaps}) {
                ##
                ## classify linked markers into reliable fragments (seeds) or unreliable gaps
                ##
                for (my $i=0; $i<@{$MARKERS{$sample}->{$chrom}->{linked}}; $i++)
                {
                    my ($bg_start, $bg_end, $bg_type, $markers_cnt) = @{$MARKERS{$sample}->{$chrom}->{linked}->[$i]};
                    
                    my $bg_len = $bg_end - $bg_start + 1;
                  
                    my $is_gap = 0;
                    
                    ## gaps failing minimum length requirement
                    $is_gap = 1 if ($opts->{min_frag_length} && $bg_len < $opts->{min_frag_length});
                    ## gaps failing minimum markers requirement
                    $is_gap = 1 if ($opts->{min_frag_markers} && $markers_cnt < $opts->{min_frag_markers});
                    ## gaps failing minimum density requirement
                    $is_gap = 1 if ($opts->{min_frag_density} && ($markers_cnt/$bg_len) < $opts->{min_frag_density});
                    
                    
                    if ($is_gap) {
                        push @{$MARKERS{$sample}->{$chrom}->{gaps}}, $i;
                    }
                    else {
                        if ($MARKERS{$sample}->{$chrom}->{blocks}) {
                            my ($prev_start, $prev_end, $prev_type, $rh_prev_markers) = @{$MARKERS{$sample}->{$chrom}->{blocks}->[-1]};
                            
                            if ($bg_type eq $prev_type) {
                                ##
                                ## merge adjacent fragments with same type
                                ##
                                $MARKERS{$sample}->{$chrom}->{blocks}->[-1]->[1] = $bg_end;
                                $MARKERS{$sample}->{$chrom}->{blocks}->[-1]->[3]->{$bg_type} += $markers_cnt;
                                
                                for my $gap_index (@{$MARKERS{$sample}->{$chrom}->{gaps}})
                                {
                                    my ($type, $cnt) = @{$MARKERS{$sample}->{$chrom}->{linked}->[$gap_index]}[2,3];
                                    
                                    $MARKERS{$sample}->{$chrom}->{blocks}->[-1]->[3]->{$type} += $cnt;
                                }
                                
                                ## empty previous merged gaps
                                @{$MARKERS{$sample}->{$chrom}->{gaps}} = ();
                            }
                            else {
                                ##
                                ## extend previous fragment and turn to a new fragment with different type
                                ##
                                
                                ###print STDERR "markers:@{$MARKERS{$sample}->{$chrom}->{gaps}}\n";
                                
                                ## find the maximum extendable position of previous fragments
                                my $prev_ext_index     = -1;
                                my %gaps_markers       = ();
                                   $gaps_markers{add}  = 0;
                                   $gaps_markers{prev} = 0;
                                   $gaps_markers{perc} = 0;
                                for my $gap_index (@{$MARKERS{$sample}->{$chrom}->{gaps}})
                                {
                                    my ($gap_start, $gap_end, $type, $cnt) = @{$MARKERS{$sample}->{$chrom}->{linked}->[$gap_index]};
                                    $gaps_markers{add} += $cnt;
                                    
                                    if ($type eq $prev_type) {
                                        $gaps_markers{prev} += $cnt;
                                    }
                                    
                                    my $prev_perc = 100 * $gaps_markers{prev} / $gaps_markers{add};
                                    
                                    ###print STDERR "$prev_ext_index\t$prev_perc\t$gaps_markers{perc}\n";
                                    
                                    if ($prev_perc >= 50 && $prev_perc > $gaps_markers{perc}) {
                                        $prev_ext_index = $gap_index;
                                        
                                        $gaps_markers{perc} = $prev_perc;
                                    }
                                }
                                
                                ##
                                ## extend previous fragments
                                ##
                                my @prev_ext_gaps = grep { $_ <= $prev_ext_index } (@{$MARKERS{$sample}->{$chrom}->{gaps}});
                                
                                
                                
                                my $prev_ext_end  = $prev_end;
                                for my $gap_index (@prev_ext_gaps)
                                {
                                    my ($gap_start, $gap_end, $type, $cnt) = @{$MARKERS{$sample}->{$chrom}->{linked}->[$gap_index]};
                                    
                                    $MARKERS{$sample}->{$chrom}->{blocks}->[-1]->[3]->{$type} += $cnt;
                                    
                                    $prev_ext_end = $gap_end;
                                }
                                
                                $MARKERS{$sample}->{$chrom}->{blocks}->[-1]->[1] = $prev_ext_end;
                                
                                
                                ##
                                ## extend next fragments
                                ##
                                push @{$MARKERS{$sample}->{$chrom}->{blocks}}, [$bg_start, $bg_end, $bg_type, {$bg_type => $markers_cnt}];
                                
                                my @next_ext_gaps  = grep { $_ >  $prev_ext_index } (@{$MARKERS{$sample}->{$chrom}->{gaps}});
                                
                                ###print STDERR "final:$prev_ext_index\tprev:@prev_ext_gaps\tnext:@next_ext_gaps\n";
                                
                                my $next_ext_start = $bg_start;
                                for my $gap_index (@next_ext_gaps)
                                {
                                    my ($gap_start, $gap_end, $type, $cnt) = @{$MARKERS{$sample}->{$chrom}->{linked}->[$gap_index]};
                                    
                                    $MARKERS{$sample}->{$chrom}->{blocks}->[-1]->[3]->{$type} += $cnt;
                                    
                                    $next_ext_start = $gap_start if ($gap_start < $next_ext_start);
                                }
                                
                                $MARKERS{$sample}->{$chrom}->{blocks}->[-1]->[0] = $next_ext_start;
                                
                                ## empty previous merged gaps
                                @{$MARKERS{$sample}->{$chrom}->{gaps}} = ();
                            }
                        }
                        else {
                            push @{$MARKERS{$sample}->{$chrom}->{blocks}}, [$bg_start, $bg_end, $bg_type, {$bg_type => $markers_cnt}];
                            
                            ##
                            ## check if some gaps occurs at the begining of the chromosome, simply merge it to the first fragment
                            ##
                            if ($MARKERS{$sample}->{$chrom}->{gaps} && scalar @{$MARKERS{$sample}->{$chrom}->{gaps}} > 0) {
                                my $ext_start = $bg_start;
                                for my $gap_index (@{$MARKERS{$sample}->{$chrom}->{gaps}})
                                {
                                    my ($gap_start, $gap_end, $type, $cnt) = @{$MARKERS{$sample}->{$chrom}->{linked}->[$gap_index]};
                                    
                                    $MARKERS{$sample}->{$chrom}->{blocks}->[-1]->[3]->{$type} += $cnt;
                                    
                                    $ext_start = $gap_start if ($gap_start < $ext_start);
                                }
                                
                                $MARKERS{$sample}->{$chrom}->{blocks}->[-1]->[0] = $ext_start;
                                
                                ## empty previous merged gaps
                                @{$MARKERS{$sample}->{$chrom}->{gaps}} = ();
                            }
                        }
                    }
                }
                
                ##
                ## check if any gaps remains at the end of the chromosome, simply merge it to the last fragment
                ##
                if ($MARKERS{$sample}->{$chrom}->{gaps} && scalar @{$MARKERS{$sample}->{$chrom}->{gaps}} > 0) {
                    next unless(defined $MARKERS{$sample}->{$chrom}->{blocks});
                    
                    my $ext_end = $MARKERS{$sample}->{$chrom}->{blocks}->[-1]->[1];
                    
                    for my $gap_index (@{$MARKERS{$sample}->{$chrom}->{gaps}})
                    {
                        my ($gap_start, $gap_end, $type, $cnt) = @{$MARKERS{$sample}->{$chrom}->{linked}->[$gap_index]};
                        
                        $MARKERS{$sample}->{$chrom}->{blocks}->[-1]->[3]->{$type} += $cnt;
                        
                        $ext_end = $gap_end if ($gap_end > $ext_end);
                    }
                    
                    $MARKERS{$sample}->{$chrom}->{blocks}->[-1]->[1] = $ext_end;
                    
                    ## empty previous merged gaps
                    @{$MARKERS{$sample}->{$chrom}->{gaps}} = ();
                }
                
                ###print STDERR Dumper(scalar @{$MARKERS{$sample}->{$chrom}->{blocks}});exit;
                
                for (my $i=0; $i < scalar @{$MARKERS{$sample}->{$chrom}->{blocks}}; $i++)
                {
                    my ($bg_start, $bg_end, $bg_type, $rh_markers) = @{$MARKERS{$sample}->{$chrom}->{blocks}->[$i]};
                    
                    my $bg_len = $bg_end - $bg_start + 1;
                    
                    ##
                    ## further extend blocks to cover whole chromosome
                    ##
                    my $ext_start = $bg_start;
                    my $ext_end   = $bg_end;
                    if ($i == 0) {
                        $ext_start = 1;
                    }
                    if ($i+1 == scalar @{$MARKERS{$sample}->{$chrom}->{blocks}}) {
                        $ext_end = $CHROM_LENGTHs{$chrom};
                    }
                    
                    if ($i > 0) {
                        my ($prev_start, $prev_end) = @{$MARKERS{$sample}->{$chrom}->{blocks}->[$i-1]}[0,1];
                        
                        $ext_start = int(($prev_end + $bg_start) / 2) + 1;
                    }
                    elsif ($i+1 < scalar @{$MARKERS{$sample}->{$chrom}->{blocks}}) {
                        my ($next_start, $next_end) = @{$MARKERS{$sample}->{$chrom}->{blocks}->[$i+1]}[0,1];
                        
                        $ext_end = int(($bg_end + $next_start) / 2);
                    }
                    
                    my $ext_len = $ext_end - $ext_start + 1;
                    
                    ##
                    ## check marker consists
                    ##
                    my @marker_consists = ();
                    my $markers_total   = 0;
                    for my $type (sort keys %{$rh_markers})
                    {
                        push @marker_consists, "$type:$rh_markers->{$type}";
                        $markers_total += $rh_markers->{$type};
                    }
                    
                    my $marker_density     = $markers_total ? sprintf("%.3f", 1000 * $markers_total / $bg_len) : '-';
                    my $marker_density_ext = $markers_total ? sprintf("%.3f", 1000 * $markers_total / $ext_len) : '-';
                    my $marker_consists    = (@marker_consists == 0) ? "-" : (join ",", @marker_consists);
                    
                    ## count of major types, to be rewrite for better representation
                    my $major_cnt  = $rh_markers->{$bg_type} ? $rh_markers->{$bg_type} : 0;
                    my $major_perc = $markers_total ? sprintf("%.2f", 100 * $major_cnt/$markers_total) : '-';
                    
                    
                    ##
                    ## calculate stats of distance between two adjancent markers
                    ##
                    my @marker_pos  = grep { $_ >= $bg_start && $_ <= $bg_end } @{$MARKERS{$sample}->{$chrom}->{pos_all}};
                    
                    my $marker_stats = stat_vars_dist(\@marker_pos);
                    
                    my $markers_details = "$markers_total;$marker_density,$marker_density_ext;$marker_stats;$marker_consists;$major_perc";
                    
                    
                    ##
                    ## soft filter (mask only) of clustered blocks
                    ##
                    my %out_filters = ();
                    if ($opts->{min_major_perc} && $major_perc < $opts->{min_major_perc}) {
                        $out_filters{Complex} = 1;
                    }
                    if ($opts->{min_block_density} && ($markers_total/$ext_len) < $opts->{min_block_density}) {
                        $out_filters{LowDensity} = 1;
                    }
                    
                    my @out_filters = sort keys %out_filters;
                    my $out_filters = join ",", @out_filters;
                    
                    if (@out_filters > 0) {
                        $bg_type .= "($out_filters)";
                    }
                    
                    my $out_line = "$sample\t$chrom\t$bg_type\t$ext_start\t$ext_end\t$ext_len\t$bg_start\t$bg_end\t$bg_len\t$markers_details";
                    
                    if ($opts->{threads} > 1) {
                        push @block_results, $out_line;
                    }
                    else {
                        print STDOUT "$out_line\n";
                    }
                }
            }
            else {
                for my $ra_block (@{$MARKERS{$sample}->{$chrom}->{linked}})
                {
                    my ($bg_start, $bg_end, $bg_type, $markers_cnt)   = @{$ra_block};
                    
                    my $bg_len = $bg_end - $bg_start + 1;
                    
                    my $out_line = sprintf("%s\t%s\t%s" . "\t%d" x 4, $sample, $chrom, $bg_type, $bg_start, $bg_end, $bg_len, $markers_cnt);
                    if ($opts->{threads} > 1) {
                        push @block_results, $out_line;
                    }
                    else {
                        print STDOUT "$out_line\n";
                    }
                } 
            }
            #########################################################################
            
            if ($opts->{threads} > 1) {
                $pm->finish(0, [$sample, $chrom, \@block_results]);
            }
        }
    }
    $pm->wait_all_children;
    
    print STDERR "\tdone!\n";
    
    ###print Dumper(%MARKERS);exit;
    
    return 0 unless($opts->{threads} > 1);
    
    print STDERR ">> Start collecting results ... ";
    for my $sample (sort keys %final_blocks)
    {
        for my $chrom (sort keys %{$final_blocks{$sample}})
        {
            my @blocks_results = @{$final_blocks{$sample}->{$chrom}};
            
            for my $block_detail (@blocks_results)
            {
                print STDOUT "$block_detail\n";
            }
        }
    }
    print STDERR "\tdone!\n";
    
    return 0;
}



sub markers2blocks_str
{
    my ($opts) = @_;
    
    ## prepare conversion symbols
    my %type2char = ();
    for my $pair (split /\;/, $opts->{type2char_str})
    {
        my ($type, $char) = (split /\:/, $pair);
        
        $type2char{$type} = $char;
    }
    
    my %char2type     = reverse %type2char;
    my %char2type_new = %char2type;
    if ($opts->{char2type_str}) {
        for my $pair (split /\;/, $opts->{char2type_str})
        {
            my ($char, $type) = (split /\:/, $pair);
            
            $char2type_new{$char} = $type;
        }  
    }
    
    ## generate match patterns
    my @matches = ();
    for my $char (keys %char2type)
    {
        push @matches, "$char\+";
    }
    my $match_str = join "\|", @matches;
    

    ##
    ## get ids and lengths of each chromosome from a file
    ##
    my @CHROM_IDs     = ();
    my %CHROM_LENGTHs = ();
    if ($opts->{length_file}) {
        print STDERR ">> Start parsing $opts->{length_file} ... ";
        get_genome_length(\@CHROM_IDs, \%CHROM_LENGTHs, $opts->{length_file});
        print STDERR "done!\n";
    }

    print STDERR ">> Start parsing $opts->{vcf} ... ";
    my @Samples_ids = ();
    my %MARKERS     = ();
    my $fh = getInputFilehandle($opts->{vcf});
    while (<$fh>)
    {
        if (/#CHROM/) {
            next if (@Samples_ids > 0);  ## support for multi vcf files to pipe in
            
            my @lines = (split /\s+/);
            
            @Samples_ids = @lines[9..$#lines];
        }
        elsif (/\#\#contig=<ID=(.*?),length=(\d+)>/ && !$opts->{length_file}) {
            next if (@Samples_ids > 0);  ## support for multi vcf files to pipe in
            
            push @CHROM_IDs, $1;
            $CHROM_LENGTHs{$1} = $2; 
        }
        
        next if (/^\#/ || /^\s+$/);
        
        next if ($opts->{filter_str} && /$opts->{filter_str}/i); 
        
        my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, 
            $FILTER, $INFO, $FORMAT, @SAMPLES) = (split /\s+/);

        ##
        ## filtering by variant type
        ##
        my $var_type = get_var_type($REF, $ALT);
        
        next if ($opts->{var_type} && ($var_type !~ /$opts->{var_type}/));
            
        my @tags = (split /\:/, $FORMAT);
        my %tags = ();
        for (my $i=0; $i<@tags; $i++) { $tags{$tags[$i]} = $i; }

        for (my $i=0; $i<@SAMPLES; $i++)
        {
            my $type = $opts->{source_tag} ? ((split /:/, $SAMPLES[$i])[$tags{$opts->{source_tag}}]) : $SAMPLES[$i];
            
            next unless($type2char{$type});
                
            push @{$MARKERS{$Samples_ids[$i]}->{$CHROM}->{char}}, $type2char{$type};
            push @{$MARKERS{$Samples_ids[$i]}->{$CHROM}->{pos}}, $POS;
            
            $MARKERS{$Samples_ids[$i]}->{$CHROM}->{$POS} = $type;
        }
    }
    print STDERR "done!\n";
    
    print STDERR ">> Start process blocks ... \n";
    
    ##
    ## get all results from each childs
    ##
    my %final_blocks = ();
    my $pm = new Parallel::ForkManager($opts->{threads});
    if ($opts->{threads} > 1) {
        $pm->run_on_finish(
            sub{
                my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data) = @_;
                
                my ($sample, $chrom, $ra_blocks) = @{$data};
                
                $final_blocks{$sample}->{$chrom} = $ra_blocks;
            }
        );
    }

    my $total_jobs  = 0;
    my $curr_job_id = 0;
    for my $sample (@Samples_ids)
    {
        for my $chrom (sort keys %{$MARKERS{$sample}})
        {
            $total_jobs ++;
        }
    }
    
    print STDOUT "$opts->{source_line}\n";
    print STDOUT "#Sample\tChrom\tType\tStart\tEnd\tLength\tMarker_Start\tMarker_End\tMarkers_Count\n";
    for my $sample (@Samples_ids)
    {
        for my $chrom (sort keys %{$MARKERS{$sample}})
        {
            $curr_job_id ++;
            
            my $pid = $pm->start and next if ($opts->{threads} > 1);
            
            print STDERR ">> Start initial processing $sample:$chrom ... $curr_job_id/$total_jobs\n";
            
            #########################################################################
            ## child processes
            my $marker_seq = join '', @{$MARKERS{$sample}->{$chrom}->{char}};
            my @genome_pos = @{$MARKERS{$sample}->{$chrom}->{pos}};
            
            while ($marker_seq =~ m/(Z+)/gx)
            {
                my $unknown = $1;
                
                my $end     = (pos $marker_seq) - 1;
                my $len     = length $unknown;
                my $start   = $end - $len + 1;
                
                my $prev_type = substr($marker_seq, $start-1, 1);
                my $next_type = substr($marker_seq, $end+1, 1);
                
                if ($start < 1) {
                    substr($marker_seq, $start, $len, $next_type x $len);
                }
                elsif (!$next_type) {
                    substr($marker_seq, $start, $len, $prev_type x $len);
                }
                elsif ($prev_type eq $next_type) {
                    substr($marker_seq, $start, $len, $prev_type x $len);
                }
                else {
                    if ($opts->{trans_fill_char}) {
                        substr($marker_seq, $start, $len, $opts->{trans_fill_char} x $len);
                    }
                    else {
                        substr($marker_seq, $start, int($len/2), $prev_type x int($len/2));
                        substr($marker_seq, $start+int($len/2), $len-int($len/2), $next_type x ($len-int($len/2)));
                    }
                }
            }
            ###print STDERR ">> Finished step1 \n";
            
            my $chrom_seq  = 'Z' x $CHROM_LENGTHs{$chrom} if ($opts->{fill_gaps});
            
            ###print STDERR ">> Start step2 ... \n";
            my @block_results = ();
            while ($marker_seq =~ m/($match_str)/g)
            {
                my $block = $1;
                
                my $char    = substr($block, 0, 1);
                my $end     = (pos $marker_seq) - 1;
                my $len     = length $block;
                my $start   = $end - $len + 1;
                
                my $rel_start = $genome_pos[$start];
                my $rel_end   = $genome_pos[$end];
                my $rel_len   = $rel_end - $rel_start + 1;
                   
                my $type = $char2type_new{$char};
                
                
                ## find gaps between fragments and trying to filling them
                if ($opts->{fill_gaps}) {
                    if ($rel_start < $opts->{min_frag_length}) {
                        $rel_start = 1;
                    }
                    if ($CHROM_LENGTHs{$chrom} - $rel_end < $opts->{min_frag_length}) {
                        $rel_end = $CHROM_LENGTHs{$chrom};
                    }
            
                    my $rel_len  = $rel_end - $rel_start + 1;
                    
                    my $is_gap = 0;
                    
                    $is_gap = 1 if ($opts->{min_frag_length} && $rel_len < $opts->{min_frag_length});          ## gaps failing minimum length requirement
                    $is_gap = 1 if ($opts->{min_frag_markers} && $len < $opts->{min_frag_markers});            ## gaps failing minimum markers requirement
                    $is_gap = 1 if ($opts->{min_frag_density} && ($len/$rel_len) < $opts->{min_frag_density}); ## gaps failing minimum density requirement
                    
                    $char = 'X' if $is_gap;
                    
                    substr($chrom_seq, $rel_start-1, $rel_len, $char x $rel_len);
                }
                else {
                    my $line = sprintf("%s\t%s\t%s" . "\t%d" x 6,
                                       $sample, $chrom, $type, $rel_start, $rel_end, 
                                       $rel_len, ($start+1), ($end+1), $len);
                    if ($opts->{threads} > 1) {
                        push @block_results, $line;
                    }
                    else {
                        print STDOUT "$line\n";
                    }
                }
            }
            ###print STDERR ">> Finished step2 ... \n";
            
            if ($opts->{fill_gaps}) {
                ###
                ### Merge gaps into adjacent fragments
                ###
                while ($chrom_seq =~ m/((X|Z)+)/gx)
                {
                    my $frag = $1;
    
                    my $frag_end     = (pos $chrom_seq);
                    my $frag_len     = length $frag;
                    my $frag_start   = $frag_end - $frag_len + 1;
    
                    my $prev_char    = ($frag_start > 2) ? substr($chrom_seq, $frag_start-2, 1) : "BEGIN";
                    my $next_char    = ($frag_end < $CHROM_LENGTHs{$chrom}) ? substr($chrom_seq, $frag_end, 1) : "END";
                    
                    if ($prev_char eq $next_char) {
                        substr($chrom_seq, $frag_start-1, $frag_len, $prev_char x $frag_len);
                    }
                    elsif ($prev_char eq "BEGIN") {
                        substr($chrom_seq, $frag_start-1, $frag_len, $next_char x $frag_len);
                    }
                    elsif ($next_char eq "END") {
                        substr($chrom_seq, $frag_start-1, $frag_len, $prev_char x $frag_len);
                    }
                    else {
                        if ($opts->{trans_fill_char}) {
                            substr($chrom_seq, $frag_start-1, $frag_len, $opts->{trans_fill_char} x $frag_len);
                        }
                        else {
                            substr($chrom_seq, $frag_start-1, int($frag_len/2), $prev_char x int($frag_len/2));
                            substr($chrom_seq, $frag_start+int($frag_len/2)-1,
                                   $frag_len - int($frag_len/2),
                                   $next_char x ($frag_len - int($frag_len/2)));
                        }
                    }
                }
                
                
                print STDERR ">> Start block refinement $sample:$chrom ... $curr_job_id/$total_jobs\n";
                while ($chrom_seq =~ m/($match_str)/g)
                {
                    my $bg_block = $1;
                    
                    my $bg_char    = substr($bg_block, 0, 1);
                    my $bg_end     = (pos $chrom_seq);
                    my $bg_len     = length $bg_block;
                    my $bg_start   = $bg_end - $bg_len + 1;
                    
                    my $bg_type  = $char2type{$bg_char} ? $char2type{$bg_char} : "Unknown";
                    
                    my $out_line = '';
                    
                    ##
                    ## check marker consists
                    ##
                    my @markers_all = ();
                    
                    ## original marker consists
                    my @original_markers = grep { $MARKERS{$sample}->{$chrom}->{$_} } ($bg_start..$bg_end);
                        
                    if (@original_markers == 0) {
                        $out_line = "$sample\t$chrom\t$bg_type\t$bg_start\t$bg_end\t$bg_len\t-\t-\t0";
                    }
                    else {
                        my %counts  = ();
                        for my $pos (@original_markers)
                        {
                            my $origin_type = $MARKERS{$sample}->{$chrom}->{$pos};
                            $counts{$origin_type}++;
                            $counts{total}++;
                        }
                        
                        my @marker_consists = ();
                        for my $origin_type (sort keys %counts)
                        {
                            next if ($origin_type eq 'total');
                            next if ($counts{$origin_type} <= 0);
                            
                            push @marker_consists, "$origin_type:$counts{$origin_type}";
                        }
                        
                        my $marker_consists = (@marker_consists == 0) ? "-" : (join ",", @marker_consists);
                        
                        ## count of major types, to be rewrite for better representation
                        my $major_cnt  = $counts{$bg_type} ? $counts{$bg_type} : 0;
                        my $major_perc = $counts{total} ? sprintf("%.2f", 100 * $major_cnt/$counts{total}) : '-';
                        
                        push @markers_all, "$counts{total};$marker_consists;$major_perc";
                        
                        my $markers_all = join "\t", @markers_all;
                        
                        $out_line = "$sample\t$chrom\t$bg_type\t$bg_start\t$bg_end\t$bg_len\t$original_markers[0]\t$original_markers[-1]\t$markers_all";
                    }

                    
                    if ($opts->{threads} > 1) {
                        push @block_results, $out_line;
                    }
                    else {
                        print STDOUT "$out_line\n";
                    }
                }
                
            }
            #########################################################################
            
            if ($opts->{threads} > 1) {
                $pm->finish(0, [$sample, $chrom, \@block_results]);
            }
        }
    }
    $pm->wait_all_children;
    
    return 0 unless($opts->{threads} > 1);
    
    print STDERR ">> Start collecting results ... ";
    for my $sample (sort keys %final_blocks)
    {
        for my $chrom (sort keys %{$final_blocks{$sample}})
        {
            my @blocks_results = @{$final_blocks{$sample}->{$chrom}};
            
            for my $block_detail (@blocks_results)
            {
                print STDOUT "$block_detail\n";
            }
        }
    }
    print STDERR "\tdone!\n";
    
    return 0;
}




=head2 count_base_changes

    About   : Count base changes
    Usage   : count_base_changes($opts->{vcf});
    Args    : Vcf file needs to be processed
    Returns : Null

=cut
sub count_base_changes
{
    my ($opts) = @_;
    
    my %query_GTs = ();
       $query_GTs{$_} = 1 for @{$opts->{query_GTs}};
    
    my $fh = getInputFilehandle($opts->{vcf});
    my @Samples_ids = ();
    my %counts_all  = ();
    while (<$fh>)
    {
        if (/#CHROM/) {
            next if (@Samples_ids > 0);  ## support for multi vcf files to pipe in
            
            my @line   = (split /\s+/);
            
            @Samples_ids = @line[9..$#line];
            
            my $Samples_ids = join "\t", @Samples_ids;
            
            #print STDOUT "##source=$SOURCE $CMDLINE\n";
            print STDOUT "#changes\t$Samples_ids\n";
            
            next;
        }
        if (/^\#/ || /^\s+$/) {
            next;
        }
        
        my @line   = (split /\s+/);
        my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, @SAMPLES) = @line;
        
        my @alleles  = ($REF, (split /\,/, $ALT));
        my $var_type = get_var_type($REF, $ALT);
        next if ($var_type ne 'snp');
        
        next if (($opts->{min_allele_types} && @alleles < $opts->{min_allele_types}) ||
                 ($opts->{max_allele_types} && @alleles > $opts->{max_allele_types}));
        
        my @tags = (split /\:/, $FORMAT);
        my %tags = ();
        for (my $i=0; $i<@tags; $i++) { $tags{$tags[$i]} = $i; }

        for (my $i=0; $i<@SAMPLES; $i++)
        {
            my $GT = (split /:/, $SAMPLES[$i])[$tags{GT}];
            
            next unless($GT =~ /((\d)(\/|\|)(\d))/);
            
            next if(@{$opts->{query_GTs}} > 0 && !$query_GTs{$GT});
            
            my ($geno, $allele1, $separator, $allele2) = ($1, $2, $3, $4);
            
            my $non_ref = $allele1 > $allele2 ? $allele1 : $allele2;
            
            next if ($non_ref == 0);
            
            $counts_all{"$REF>$alleles[$non_ref]"}->{$Samples_ids[$i]} ++;
        }
    }
    
    for my $change (sort keys %counts_all)
    {
        my @samples_counts = ();
        for my $sample (@Samples_ids)
        {
            my $count = $counts_all{$change}->{$sample} ? $counts_all{$change}->{$sample} : 0;
            
            push @samples_counts, $count;
        }
        
        my $samples_counts = join "\t", @samples_counts;
        
        print "$change\t$samples_counts\n";
    }
}



=head2 count_diagnose

    About   : Summary results from GATK DiagnoseTargets 
    Usage   : count_diagnose($opts->{vcf});
    Args    : Vcf file needs to be processed
    Returns : Null

=cut
sub count_diagnose
{
    my ($opts) = @_;
    
    my $fh = getInputFilehandle($opts->{vcf});
    my @Samples_ids     = ();
    my %Samples_summary = ();
    my $interval_total  = 0;
    while (<$fh>)
    {
        if (/#CHROM/) {
            next if (@Samples_ids > 0);  ## support for multi vcf files to pipe in
            
            my @line   = (split /\s+/);
            
            @Samples_ids = @line[9..$#line];
            
            print STDOUT "$opts->{source_line}\n";
            print STDOUT "#SUMMARY_TYPE\tID\tTOTAL\tPASS\tFILTER\n";
            next;
        }
        if (/^\#/ || /^\s+$/) {
            next if (@Samples_ids > 0);  ## support for multi vcf files to pipe in
            
            print STDOUT if (/\#\#FILTER/);
            next;
        }
        
        my @line   = (split /\s+/);
        my ($CHROM, $START, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT) = @line[0..8];
        
        my @tags = (split /\:/, $FORMAT);
        my %tags = ();
        for (my $i=0; $i<@tags; $i++) { $tags{$tags[$i]} = $i; }
        
        my ($END) = ($INFO =~ /END=(\d+)/);
        
        my %interval_summary = ();
        for (my $i=9; $i<@line; $i++)
        {
            my $sample_flt_field = defined($tags{FT}) ? (split /:/, $line[$i])[$tags{FT}] : "FT_MISSING";
               $sample_flt_field =~ s/\;/,/g;
            
            $Samples_summary{$Samples_ids[$i-9]}->{$sample_flt_field} ++;
            $interval_summary{$sample_flt_field} ++;
        }
        
        ##
        ## summary diagnose results in each interval
        ##
        my @sample_details = ();
        for my $flt (sort keys %interval_summary)
        {
            next if ($flt eq 'PASS');
            push @sample_details, "$flt:$interval_summary{$flt}";
        }
        my $sample_details = (@sample_details >= 1) ? (join ";", @sample_details) : '-';
        
        my $sample_total = scalar @Samples_ids;
        my $sample_pass  = $interval_summary{PASS} ? $interval_summary{PASS} : 0;
        
        print STDOUT "INTERVAL\t$CHROM:$START\-$END\t$sample_total\t$sample_pass\t$sample_details\n";
    
        $interval_total ++;
    }

    ##
    ## summary diagnose results in each sample
    ##
    for my $sample_id (sort keys %Samples_summary)
    {
        my @interval_details = ();
        for my $flt (sort keys %{$Samples_summary{$sample_id}})
        {
            next if ($flt eq 'PASS');
            
            push @interval_details, "$flt:$Samples_summary{$sample_id}->{$flt}";
        }
        
        my $interval_pass  = $Samples_summary{$sample_id}->{PASS} ? $Samples_summary{$sample_id}->{PASS} : 0;
        
        my $interval_details = (@interval_details >= 1) ? (join ";", @interval_details) : '-';
        
        print STDOUT "SAMPLE\t$sample_id\t$interval_total\t$interval_pass\t$interval_details\n";
    }
    
}


=head2 check_variants_context

    About   : Sequence context analysis within and around variants loci.
    Usage   : check_variants_context($file);
    Args    : File contains variants info.
    Returns : Null

=cut
sub check_variants_context
{
    my ($opts) = @_;

    if ($opts->{output}) {
        open (STDOUT, "> $opts->{output}") || die $!;
    }
  
    my %SEQs = ();
    parse_fasta_SEQs(\%SEQs, $opts->{fasta});
    
    my %counts_all = ();
    
    my $fh = getInputFilehandle($opts->{vcf});
    while (<$fh>)
    {
        if (/^\#\#/) {
            print STDOUT; next;
        }
        elsif (/^\#/) {
            print STDOUT <<EOF;
##INFO=<ID=Context,Number=1,Type=String,Description="Local sequence context within and around variants loci">
EOF
            print STDOUT "$opts->{source_line}\n";
            print STDOUT; next;
        }
        elsif (/^\s+$/) {
            next;
        }

        my @line   = (split /\s+/);
        my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT) = @line[0..8];
        
        ###print STDERR "$CHROM\t$POS\n";exit;
        
        my $ref_len = length $REF;
        my $alt_len = length $ALT;
        
        my $prev_base = substr($REF, 0, 1);
        
        my $seq_context = '';
        
        if ($ref_len != $alt_len) {
            ##
            ## check INDEL context
            ##
            
            ##
            ## get inserted or deleted bases
            ##
            my $indel = '';
            
            substr($REF, 0, 1, '');
            substr($ALT, 0, 1, '');
            
            if ($ref_len > $alt_len) {
                $indel = $REF;
                if ($alt_len >= 2) {
                    $indel =~ s/$ALT//;
                }
            }
            elsif ($ref_len < $alt_len) {
                $indel = $ALT;
                if ($ref_len >= 2) {
                    $indel =~ s/$REF//;
                }
            }
            
            ##
            ## check context of indel itself
            ##
            my @nts = split '', $indel;
            my %nts = ();
               $nts{$_}++ for (@nts);
            
            my @bases = (keys %nts);
            
            my $cmp_seq = (@bases == 1) ? $bases[0] : $indel;
            my $ref_seq = uc(substr($SEQs{$CHROM}, $POS-1, $opts->{indel_extend}));
            
            ##
            ## check flanking sequence context
            ##
            my $ref_cnt   = 0;
            if ($ref_seq =~ m/^($prev_base)(($cmp_seq)+)/gx)
            {
                $ref_cnt = (length $2) / (length $cmp_seq);
            }
            
            my $indel_cnt = 0;
            if ($indel =~ m/(($cmp_seq)+)/gx)
            {
                $indel_cnt = (length $1) / (length $cmp_seq);
            }
            
            $seq_context = "[$cmp_seq]$indel_cnt,$ref_cnt";
        }
        else {
            ##
            ## check SNP context
            ##
            my $ref_seq = uc(substr($SEQs{$CHROM}, $POS-$opts->{snp_extend}-1, $opts->{snp_extend}*2));
            my $ref_cnt = 0;
            if ($ref_seq =~ m/[^$REF]($REF+)[^$REF]/gx)
            {
                $ref_cnt = (length $1);
            }
            
            my $alt_seq = $ref_seq;
            
            substr($alt_seq, $opts->{snp_extend}-1, 1, $ALT);
            
            my $alt_cnt   = 0;
            if ($alt_seq =~ m/[^$ALT]($ALT+)[^$ALT]/gx)
            {
                $alt_cnt = (length $1);
            }
            
            $seq_context = "[$REF]$ref_cnt,[$ALT]$alt_cnt";
        }
        
        if ($line[7] eq '.') {
            $line[7] = "Context=$seq_context";
        }
        else {
            $line[7] .= ";Context=$seq_context";
        }
        
        my $out_line = join "\t", @line;
        print STDOUT "$out_line\n";
    }
}




=head2 fix_AD_fields

    About   : Fix unmatch AD fields after merge vcf files.
    Usage   : fix_AD_fields($AD, $GT, $allele_number);
    Args    : Previous AD info;
              Sample's GT info;
              Total alleles in this locus.
    Returns : Fixed AD info.

=cut
sub fix_AD_fields
{
    my ($AD, $GT, $allele_number) = @_;
    
    my @prev_dps = split /\,/, $AD;
    
    ## fix AD tags if some fileds is missing
    my @GTs = split /\/|\|/, $GT;
    my %GT_dps = ();
       
    for (my $j=0; $j<@GTs; $j++)
    {
        $GT_dps{$GTs[$j]} = $prev_dps[$j];
    }
    
    my @fixed_dps = ();
    
    for (my $j=0; $j<$allele_number; $j++)
    {
        if (defined $GT_dps{$j}) {
            push @fixed_dps, $GT_dps{$j};
        }
        else {
            push @fixed_dps, 0;
        }
    }
    
    my $AD_fixed = join ',', @fixed_dps;
    
    return $AD_fixed;
}



=head2 NRNV2AD

    About   : Convert NR NV infos to AD info
    Usage   : NRNV2AD($sample, $NR_order, $NV_order);
    Args    : Infos for each sample;
              Order of NR field in each sample infos;
              Order of NV field in each sample infos;
    Returns : New AD info

=cut
sub NRNV2AD
{
    my ($info, $NR_order, $NV_order) = @_;
    
    my $NRs  = (split /:/, $info)[$NR_order];
    my $NVs  = (split /:/, $info)[$NV_order];
    
    ## NR and NV is two tags used in Platypus
    ## FORMAT=<ID=NR,Number=.,Type=Integer,Description="Number of reads covering variant location in this sample">
    ## FORMAT=<ID=NV,Number=.,Type=Integer,Description="Number of reads containing variant in this sample">
    ## *Note: the NV fields possible counts all reads differ from the reference sequence!
    my @NRs = (split /\,/, $NRs);
    my @NVs = (split /\,/, $NVs);
    
    my $NR_max = (sort {$a <=> $b} @NRs)[-1];
    my $NV_max = (sort {$a <=> $b} @NVs)[-1];
    
    
    ## use "max(NR) - max(NV)" as the reference allele depth, may be inaccurate
    my $ref_dp = $NR_max - $NV_max;
    
    my $AD = join ',', ($ref_dp, @NVs);
    
    return $AD;
}





=head2 combine_vcfs

    About   : Combine two vcf files
    Usage   : combine_vcfs(\%options);
    Args    : Primary vcf file used as master records;
              Secondary vcf file.
    Returns : Null

=cut
sub combine_vcfs
{
    my ($opts) = @_;
    
    my $primary_tag   = $opts->{primary_tag}   ? $opts->{primary_tag}   : 'Primary';
    my $secondary_tag = $opts->{secondary_tag} ? $opts->{secondary_tag} : 'Secondary';
    my $intersect_tag = $opts->{intersect_tag} ? $opts->{intersect_tag} : 'Intersection';
    
    my @combine_rows  = @{$opts->{combine_rows}} > 0 ? @{$opts->{combine_rows}} : (0, 1);
    my @compare_rows  = @{$opts->{compare_rows}} > 0 ? @{$opts->{compare_rows}} : ();
    
    ###print STDERR Dumper(@combine_rows);exit;
    
    my %headers = ();
       $headers{combine} = 0;
       $headers{compare} = 0;
    
    ## step1: read through primary file
    my $primary_fh   = getInputFilehandle($opts->{vcf});
    my %vcf_records  = ();
    my @vcf_contigs  = ();
    my $vcf_header   = '';
    my @vcf_header   = ();
    while (<$primary_fh>)
    {
        if (/^\#\#contig=<ID=(.*?),length=\d+>/) {
            push @vcf_contigs, $1; print; next;
        }
        elsif (/^\#CHROM/) {
            chomp($vcf_header = $_);
            @vcf_header = split /\s+/, $vcf_header;
            next;
        }
        elsif (/^\#\#/ || /^\s+$/) {
            if (/ID=Combine/) {
                $headers{combine} = 1;
            }
            if (/ID=SDIFF/) {
                $headers{compare} = 1;
            }
            print; next;
        }
        
        chomp(my $record = $_);
        
        my @fields = (split /\s+/, $record);
        
        ## rows for combine, if those rows differ, each will give a unique record in new vcf file
        my ($CHROM, $POS, @append_keys) = @fields[@combine_rows];
        
        my $append_keys = (@append_keys > 0) ? (join "\t", @append_keys) : "FIELDS";
        
        ## use array to ensure no duplicate records would get removed
        $vcf_records{$CHROM}->{$POS}->{$append_keys}->{primary} = 1; 
        push @{$vcf_records{$CHROM}->{$POS}->{$append_keys}->{record}}, $record;
        
        ## rows only for compare, records with same combine rows but different compare rows will
        ## be merge into one single record
        if (@compare_rows > 0) {
            push @{$vcf_records{$CHROM}->{$POS}->{$append_keys}->{compare}}, join "\t", @fields[@compare_rows];
        }
    }
    
    my $secondary_fh   = getInputFilehandle($opts->{secondary_vcf});
    my %merged_records = ();
    while (<$secondary_fh>)
    {
        if (/^\#/ || /^\s+$/) {
            next;
        }
        
        chomp(my $record = $_);
        
        my @fields = (split /\s+/, $record);
        
        ## rows for combine, if those rows differ, each will give a unique record in new vcf file
        my ($CHROM, $POS, @append_keys) = @fields[@combine_rows];
        
        my $append_keys = (@append_keys > 0) ? (join "\t", @append_keys) : "FIELDS";
        
        if ($vcf_records{$CHROM}->{$POS}->{$append_keys}->{primary}) {
            $merged_records{$CHROM}->{$POS}->{$append_keys} ++;
            
            for (my $i=0; $i<@{$vcf_records{$CHROM}->{$POS}->{$append_keys}->{record}}; $i++)
            {
                push @{$vcf_records{$CHROM}->{$POS}->{$append_keys}->{tag}}, $intersect_tag;
                
                ## rows only for compare, records with same combine rows but different compare rows will
                ## be merge into one single record
                if ($vcf_records{$CHROM}->{$POS}->{$append_keys}->{compare}->[$i]) {
                    my @primary_cmp    = split "\t", $vcf_records{$CHROM}->{$POS}->{$append_keys}->{compare}->[$i];
                    my @secondary_cmp  = @fields[@compare_rows];
                    
                    my @secondary_diffs = ();
                    for (my $j=0; $j<@primary_cmp; $j++)
                    {
                        if ($primary_cmp[$j] ne $secondary_cmp[$j]) {
                            push @secondary_diffs, "$vcf_header[$compare_rows[$j]]:$secondary_cmp[$j]";
                        }
                    }
                    
                    if (@secondary_diffs > 0) {
                        $vcf_records{$CHROM}->{$POS}->{$append_keys}->{diff}->{$i} = join "|", @secondary_diffs;
                    }
                }
            }
            
            if ($merged_records{$CHROM}->{$POS}->{$append_keys} > scalar @{$vcf_records{$CHROM}->{$POS}->{$append_keys}->{record}}) {
                push @{$vcf_records{$CHROM}->{$POS}->{$append_keys}->{record}}, $record;
            }
        }
        else {
            push @{$vcf_records{$CHROM}->{$POS}->{$append_keys}->{record}}, $record;
            push @{$vcf_records{$CHROM}->{$POS}->{$append_keys}->{tag}},    $secondary_tag;
        }
    }
    
    ## vcf contigs info miss or imcomplete
    my @present_contigs = sort (keys %vcf_records);
    if (scalar @vcf_contigs < scalar @present_contigs) {
        @vcf_contigs = @present_contigs;
    }
    
    if (!$headers{combine}) {
        print STDOUT "##INFO=<ID=Combine,Number=1,Type=String,Description=\"Source VCF for the merged record\">\n";
    }
    if (!$headers{compare}) {
        print STDOUT "##INFO=<ID=SDIFF,Number=1,Type=String,Description=\"Different fields in secondary file\">\n";
    }
    print STDOUT "$opts->{source_line}\n";
    print STDOUT "$vcf_header\n";
    for my $CHROM (@vcf_contigs)
    {
        if ($vcf_records{$CHROM})
        {
            for my $POS (sort {$a <=> $b} keys %{$vcf_records{$CHROM}})
            {
                for my $append_keys (sort keys %{$vcf_records{$CHROM}->{$POS}})
                {
                    for (my $i = 0; $i < @{$vcf_records{$CHROM}->{$POS}->{$append_keys}->{record}}; $i++)
                    {
                        my $tag = $vcf_records{$CHROM}->{$POS}->{$append_keys}->{tag} ?
                                  $vcf_records{$CHROM}->{$POS}->{$append_keys}->{tag}->[$i] : $primary_tag;
                        
                        my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, @Samples)
                         = (split /\s+/, $vcf_records{$CHROM}->{$POS}->{$append_keys}->{record}->[$i]);
                        
                        if ($vcf_records{$CHROM}->{$POS}->{$append_keys}->{diff}->{$i}) {
                            my $diff_str = $vcf_records{$CHROM}->{$POS}->{$append_keys}->{diff}->{$i};
                            $tag .= "(SDIFF=$diff_str)";
                        }
                        
                        $INFO .= ";Combine=$tag";
                        
                        my $Samples = join "\t", @Samples;
                        print STDOUT "$CHROM\t$POS\t$ID\t$REF\t$ALT\t$QUAL\t$FILTER\t$INFO\t$FORMAT\t$Samples\n";
                    }
                }
            }
        }
    }
}



1;

=head1 VERSION

1.6.8

=head1 AUTHOR

Nowind, noitulove9891@gmail.com

=head1 COPYRIGHT

Copyright (c) Nowind's Area. All rights reserved. This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself. 


=cut


