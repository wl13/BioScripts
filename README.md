# BioScripts

### INTRODUCTION
Scripts for bioinformatics processing and analysis. Supplement to other more useful tools like bcftools (http://samtools.github.io/bcftools/), vcftools (https://vcftools.github.io/index.html) and etc.

### INSTALL
Add "MyPerl" folder to PERL5LIB, something like "export PERL5LIB=$HOME/folder_contain_MyPerl", or copy it to a pre-exist folder like "perl/site/lib", or just copy it to the same folder where you run the script

### USAGE
Simply type "perl certain_script.pl" for details of each option or read through the descriptions below ..


#### calc_vcf_diversity.pl
> Calculate within- or between-groups diversities from a vcf file, an alternative choice of "vcftools --window-pi".

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
> Query, extract and processing fasta sequences.

    Version: 1.1.2
    
    Usage:   perl fasta_process.pl [options]
    
    Input/Output Options:
    
        --fasta  <filename>
            sequences file(s) in fasta format, required
        --output <filename>
            output file, default to STDOUT
            
        --out-format <string>
            specify output format, currently support:
            
            "fasta":    standard fasta format;
            "tabular":  tabulated results with id and sequence in each line,
            "seq":      output sequences only
            
            [default: fasta]
            
        --out-order <string>
            specify orders of output results,
            
            "input":    same order as input file;
            "alpha":    sort output sequences by ids using alphabetic sorting;
            "num":      sort output sequences by ids using numeric sorting;
            
            or set a file name, then the output ids will be sorted accroding
            to the order in the refer list
            
            [default: input]
    
        --no-found-order <string>
            while sort sequences by refer list, put those sequences not found in
            refer list to the "top" or "bottom" of all others, [default: top] 
            
        --wordwrap  <int>
            line feed for print, only valid while output in fasta format
    
    
    Manipulation Options:
    
        --sort-by-list  <filename>
            sort sequences by the order of ids in another file, the file contains
            sequence ids same as the fasta file, one id per line, otherwise will
            use numeric sort by default
            
        --reverse
            reverse sequences before output        
        --complement
            complement sequences before output
    
            
    Extracting Options:
     
        --query  <filename>
            query file with sequence ids need to be extracted
            
        --rows   <int>
            specify the input fields, can have multiple values (0-based),
            the first will be used as the row of the query IDs, others
            will be treated as append infomations [default: 0]
        --subset    <numbers>
            extract only part of the query sequences, set 2 values to
            indicate start and end row (0-based)
    
        --delimiter <numbers>
            specify a delimiter while reading input files, such as ",",
            "\t", multiple delimiters can be set such as ",|\t"
            [default: "\s+"]
        --match     <string>
            only considering lines matching a pattern, support perl
            regular expression
    
        --lower     <int>
            only output sequences with length above this value
        --upper     <int>
            only output sequences with length below this value
            
        --replace <string>
            replace reference nucleotide with the specified nucleotide
            at specified position for each sequence in the format:
            
            "row of position,row of reference nucleotides,row of substitution
            nucleotides"
    
    Other options:

        --count-nucl
            count "dinucleotide" or "triplet" nucleotide context
            


#### gff2fasta.pl
> Extract sequences from gff3 file to fasta file.

    Version: 1.1.0
    
    Usage:   perl gff2fasta.pl [options]
    
    Options:
        -i, --gff     <filename>
            annotation file in gff3 format, required
        -s, --seqs    <filename>
            reference sequences in fasta format, required
            
        -o, --output  <dirname>
            output filename, default to STDOUT
        
        -F, --format  <string>
            output file format, can be set to fasta or tabular,
            default: fasta
        
        -R, --features
            output features, gene or cds, default: gene
    
        -w, --word-wrap  <int>
            line feed for print
        
        -d, --details
            output verbose info in header line
        
        -?, --help
            show this help message


#### paintGenomeBlocks.pl
> Plot blocks.

    Version: 1.2.8
    
    Usage:   perl paintGenomeBlocks.pl [options]
    
    Options:
        -i, --input   <filename>
            input file(s), required
        -o, --output  <string>
            output file prefix, default named as graph.* to current directory
        -l, --length  <filename>
            a file contains contig name and length in the format:
            
            #CHROM LENGTH
            chr01 43270923
            chr02 35937250
            chr03 36413819
            chr04 35502694
            ...
            
            required while no length info can be found in the
            input file
    
        -r, --rows   <int>
            specify the input fields, can have multiple values (0-based), in
            the order "block id", "chromosome", "type", "start position",
            "end position", ["RGB values"], if the last row is specified,
            user can set local defined colors directly from RGB values in the
            format "R,G,B" [default: 0 1 2 3 4]
        
        -f, --format    <string>
            image type, currently support png and svg [default: svg]
            *currently bugs could not be fixed while output png format,
            please use svg format instead
        
        -C, --colors <strings>
            specify colors for certain block types, in the format:
            
            "type1:color1;type2:color2;..."
            
            otherwise all types will be set to white by default, current support
            colors: white, black, red, blue, yellow, grey, light_grey, green,
            cyan, light_cyan
        
        -s, --sort-blocks  <strings>
            speicfy output blocks' orders, default use alphabetic sort
        
        -w, --width     <int>
            width of the image [default: 1000]
        -h, --height    <int>
            height of the image [default: 1000]
        -t, --thickness <int>
            height of each rectangle [default: 200]
        -D, --chrom-distance <int>
            vertical distance between each chromosomes [default: 500]
        -b, --block-distance  <int>
            vertical distance between each blocks [default: 200]
    
        -a, --add-border
            add border to image
    
        --split-pair
            draw chromosome pairs, this need the type row contain two
            types delimited by "/"



#### reference_align.pl
> Mapping with clustalw2(http://www.clustal.org/clustal2/) or muscle(http://www.drive5.com/muscle/).

Align sequences to a reference sequence, this was done by 2 steps:

**Step1:** align the first sequence to the reference sequence, and get
a expanded reference sequence with gaps inserted, then align the
second sequence to the new reference sequence, iterate this process
to generate a reference sequence expand all query sequences;

**Step2:** build a consensus sequence with nucleotide replaced in reference
sequence, and re-align all sequences to the consensus sequence.

    Version: 1.1.1
    
    Usage:   perl reference_align.pl [options]
    
    Options:
        -i, --input   <filename>
            input file contains at least two sequences in fasta
            format, the first sequence appeared in this file will
            be used as the reference sequence, required
    
        -o, --output  <filename>
            output filename, output extracted sequences in fasta
            format, default to STDOUT
        
        -c, --consensus
            output consensus sequence
            
        -a, --aligner <string>
            choose aligner, 'clustalw2' or 'muscle', [default: clustalw2]
        
        -p, --params  <string>
            change parameters for specified aligner
            default: '-gapopen=15 -gapext=6.66' [clustalw2]
                     '-quiet' [muscle]
    
        -m, --maxiters <int>
            maximum number of iterations for muscle [default: 3]
    

#### vcf_process.pl
> Vcf format file related processing.

    Version: 2.4.4
    
    Usage:   perl vcf_process.pl [--vcf FILE | STDIN] [Filtering Options] [Output Options]
    
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



More to add ...
