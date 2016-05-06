# BioScripts

### INTRODUCTION
Scripts for bioinformatics processing and analysis. Supplement to other more useful tools like bcftools (http://samtools.github.io/bcftools/), vcftools (https://vcftools.github.io/index.html) and etc. Most scripts were designed to support "pipe-in" and "pipe-out".

### INSTALL
Add "MyPerl" folder to PERL5LIB, something like "export PERL5LIB=$HOME/folder_contain_MyPerl", or copy it to a pre-exist folder like "perl/site/lib", or just copy it to the same folder where you run the script

### USAGE
Simply type "perl certain_script.pl" or "perl certain_script.pl -h" for details of each option 


#### calc_vcf_diversity.pl
> Calculate within- or between-groups diversities from a vcf file, an alternative choice of "vcftools --window-pi".



#### fasta_process.pl
> Query, extract and processing fasta sequences.


#### gff2fasta.pl
> Extract sequences from gff3 file to fasta file.


#### paintGenomeBlocks.pl
> Plot blocks.


#### reference_align.pl
> Mapping with clustalw2(http://www.clustal.org/clustal2/) or muscle(http://www.drive5.com/muscle/).

Align sequences to a reference sequence, this was done by 2 steps:

**Step1:** align the first sequence to the reference sequence, and get
a expanded reference sequence with gaps inserted, then align the
second sequence to the new reference sequence, iterate this process
to generate a reference sequence expand all query sequences;

**Step2:** build a consensus sequence with nucleotide replaced in reference
sequence, and re-align all sequences to the consensus sequence.



#### vcf_process.pl
> Vcf format file related processing.


#### fgenesh2gff.pl
> Convert results from fgenesh to gff3 format


More to add ...
