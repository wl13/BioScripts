# BioScripts

## INTRODUCTION
Scripts for bioinformatics processing and analysis. Supplement to other more useful tools like bcftools (http://samtools.github.io/bcftools/), vcftools (https://vcftools.github.io/index.html) and etc. Most scripts were designed to support "pipe-in" and "pipe-out".


## INSTALL
Add "MyPerl" folder to PERL5LIB, something like "export PERL5LIB=$HOME/folder_contain_MyPerl", or copy it to a pre-exist folder like "perl/site/lib", or just copy it to the same folder where you run the script


## USAGE
Simply type "perl certain_script.pl" or "perl certain_script.pl -h" for details of each option 


### calc_vcf_diversity.pl
> Calculate within- or between-groups diversities from a vcf file, an alternative choice of "vcftools --window-pi".


### fasta_process.pl
> Query, extract and processing fasta sequences.


### gff2fasta.pl
> Extract sequences (mRNA or CDS) from gff3 file to fasta file.


### paintGenomeBlocks.pl
> Plot blocks along each chromosome.


### reference_align.pl
> Mapping with clustalw2(http://www.clustal.org/clustal2/) or muscle(http://www.drive5.com/muscle/).

Align sequences to a reference sequence, this was done by 2 steps:

**Step1:** align the first sequence to the reference sequence, and get
a expanded reference sequence with gaps inserted, then align the
second sequence to the new reference sequence, iterate this process
to generate a reference sequence expand all query sequences;

**Step2:** build a consensus sequence with nucleotide replaced in reference
sequence, and re-align all sequences to the consensus sequence.



### vcf_process.pl
> VCF format file related processing.

This script does quite a lot things, including filtering, combining, clustering and etc., seems I put too many functions here ...
However, since the VCF format generated from different caller varies, this script was manily tailored for vcf file generated from GATK (UnifiedGenotyper or HaplotypeCaller, http://www.broadinstitute.org/gatk/), some functions require the AD (allele depth) field, so it may not perform very well for VCF files generated from other caller.


#### Clustering

The clustering function is used to identify genome blocks through certain type of markers. This was done by fisrt search for the reliable seeds (segments with consecutive markers of the same type and pass the criteria, the "seeding" stage), then merge adjacent seeds with same type to form blocks (the "extension" stage), the boundary between blocks of different type was determined according to the markers present between two blocks or use the middle point while no more markers present. 
The "seeding-and-extension" algorithm was borrowed from "Wijnker, E. et al. The genomic landscape of meiotic crossovers and gene conversions in Arabidopsis thaliana. eLife 2, e01426 (2013)", which used for identify recombinat blocks.
  
      The source type could be genotypes (use GT field) or other user defined types (e.g. 
      ancestral status) an example of markers.vcf.gz could be:
      #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  sample
      chr01   161     .       C       A       54.12   .       .       GT:SC   0/0:A/A
      chr01   431     .       C       T       44.2    .       .       GT:SC   0/0:A/A
      chr01   1641    .       G       A       64.15   .       .       GT:SC   1/1:B/B
      chr01   4165    .       C       A       34.31   .       .       GT:SC   1/1:B/B
      ...
  
*No clustering, just output blocks with consecutive markers with same source type, 
    
    vcf_process.pl --vcf markers.vcf.gz --out-blocks --source-tag "SC" > markers.blocks.csv

*Clustering
    
    vcf_process.pl --vcf markers.vcf.gz --out-blocks --source-tag "SC" --fill-gaps \
        --min-frag-length 10000 --min-frag-markers 25 > markers.blocks.l10km25.csv
    
*Use paintGenomeBlocks.pl to visually compare two results
    
    awk 'BEGIN{OFS="\t";} !/\#/ {$1 = $1"-original"; print;}' markers.blocks.csv | \
        cat markers.blocks.l10km25.csv - | \
        paintGenomeBlocks.pl --input - \
        --width 1600 --height 3000 --thickness 10 --chrom-distance 20 --block-distance 2 \
        --output markers.blocks.cmp --length reference_genome.fasta.fai \
        --colors "type1:strong_red2;B:strong_blue2" --sort-blocks sample-original sample --format png
    



### fgenesh2gff.pl
> Convert results from fgenesh to gff3 format


More to add ...
