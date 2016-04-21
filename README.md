# BioScripts
Scripts for bioinformatics processing and analysis. Supplement to other more useful tools like bcftools, vcftools and etc.

### INSTALL
Add "MyPerl" folder to PERL5LIB, something like "export PERL5LIB=$HOME/folder_contain_MyPerl", or copy it to a pre-exist folder like "perl/site/lib", or just copy it to the same folder where you run the script


#### calc_vcf_diversity.pl
Calculate within- or between-groups diversities from a vcf file, an alternative choice of "vcftools --window-pi".

    Version: 2.0.0

    Usage:   perl calc_vcf_diversity.pl [options]

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
        difference between homozygous and heterozygous will be counted as 0.5
        instead of 1
    
    -t, --threads  <int>
        how many data threads should be allocated to running this analysis
        [default: 1]

    -V, --var-type    <string>
        set "snp" to process snp sites only, or set "indel" to process indels
        only

    -?, --help
        show this help message



#### fasta_process.pl
Query, extract and processing fasta sequences.

#### gff2fasta.pl
Extract sequences from gff3 file to fasta file.

#### paintGenomeBlocks.pl
Plot blocks.

#### reference_align.pl
Align sequences to a reference sequence.

#### vcf_process.pl
Vcf format file related processing.
