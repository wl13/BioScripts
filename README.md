# BioScripts

## INTRODUCTION
Scripts for bioinformatics processing and analysis. Supplement to other more useful tools like bcftools (http://samtools.github.io/bcftools/), vcftools (https://vcftools.github.io/index.html) and etc. Most scripts were designed to support "pipe-in" and "pipe-out".


## INSTALL
Add "MyPerl" folder to PERL5LIB, something like "export PERL5LIB=$HOME/folder_contain_MyPerl", or copy it to a pre-exist folder like "perl/site/lib", or just copy it to the same folder where you run the script


## USAGE
Simply type "perl certain_script.pl" or "perl certain_script.pl -h" for details of each option 


## File Process

### fasta_process.pl
> Query, extract and processing fasta sequences.

* Query a single gene

		echo "${gene_id}" | fasta_process.pl --fasta all.seq --query - --rows 0 > gene.fa

* Query multiple genes

		fasta_process.pl --rows 0 --query query_genes.list --fasta all_genes.fasta > query_genes.fas

* Extract a single region

		echo "${query_region}" | fasta_process.pl --rows 0 1 2 --subset 1 2 --query - \
        		--fasta genome.fasta > query.fa

* Extract multiple sequences from the genome

		fasta_process.pl --rows 0 1 2 --query regions.list --fasta genome.fasta --subset 1 2 > query_regions.fas
    
    
### convert_fastq_quality.pl
> Convert fastq encodings 

* Convert Sanger encoding to Illumina 1.5+

		zcat reads.fq.gz | convert_fastq_quality.pl -i - -c sanger2illumina | gzip -c > reads.converted.fq.gz


### vcf_process.pl
> VCF format file related processing.

This script does quite a lot things, including filtering, combining, clustering and etc., seems I put too many functions here ...
However, since the VCF format generated from different caller varies, this script was manily tailored for vcf file generated from GATK (UnifiedGenotyper or HaplotypeCaller, http://www.broadinstitute.org/gatk/), some functions require the AD (allele depth) field, so it may not perform very well for VCF files generated from other caller.


#### Use vcf_process.pl to clustering markers (genetically linked regions)

The clustering function is used to identify genome blocks through certain type of markers. This was done by fisrt search for the reliable seeds (segments with consecutive markers of the same type and pass the criteria, the "seeding" stage), then merge adjacent seeds with same type to form blocks (the "extension" stage), the boundary between blocks of different type was determined according to the markers present between two blocks or use the middle point while no more markers present. 
The "seeding-and-extension" algorithm was borrowed from "Wijnker, E. et al. The genomic landscape of meiotic crossovers and gene conversions in Arabidopsis thaliana. eLife 2, e01426 (2013)", which used for identify recombinat blocks.

* Example about the input vcf file. The source type could be genotypes (use GT field) or other user defined types (e.g. ancestral status) an example of markers.vcf.gz could be:

		#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  sample
		chr01   161     .       C       A       54.12   .       .       GT:SC   0/0:A/A
		chr01   431     .       C       T       44.2    .       .       GT:SC   0/0:A/A
		chr01   1641    .       G       A       64.15   .       .       GT:SC   1/1:B/B
		chr01   4165    .       C       A       34.31   .       .       GT:SC   1/1:B/B
		...
  
* No clustering, just output blocks with consecutive markers with same source type, 

		vcf_process.pl --vcf markers.vcf.gz --out-blocks --source-tag "SC" > markers.blocks.csv

* Clustering

		vcf_process.pl --vcf markers.vcf.gz --out-blocks --source-tag "SC" --fill-gaps \
		    --min-frag-length 10000 --min-frag-markers 25 > markers.blocks.l10km25.csv

* Use paintGenomeBlocks.pl to visually compare two results

		awk 'BEGIN{OFS="\t";} !/\#/ {$1 = $1"-original"; print;}' markers.blocks.csv | \
		    cat markers.blocks.l10km25.csv - | \
		    paintGenomeBlocks.pl --input - \
		    --width 1600 --height 3000 --thickness 10 --chrom-distance 20 --block-distance 2 \
		    --output markers.blocks.cmp --length reference_genome.fasta.fai \
		    --colors "type1:strong_red2;B:strong_blue2" --sort-blocks sample-original sample --format png




### fgenesh2gff.pl
> Convert results from fgenesh to gff3 format

* All sequences have predictions

		fgenesh2gff.pl -i fgenesh.txt -o fgenesh.gff

* Contain sequences did not have reliable predictions

		sed 's/ no reliable predictions /\/\//' fgenesh.txt | fgenesh2gff.pl -i - -o fgenesh.gff


### gff2fasta.pl
> Extract sequences (mRNA or CDS) from gff3 file to fasta file.

* Extract mRNA sequences

		gff2fasta.pl --gff example.gff --seqs example.fasta --features mRNA > example.mRNA.fasta

* Extract CDS sequences

		gff2fasta.pl --gff example.gff --seqs example.fasta --features cds > example.cds.fasta


### gff2tables.pl
> Convert gff to tabular format


### fasta2tabular.pl
> Convert fasta to tabular format


### intervals2bed
> Convert intervals to bed format, e.g. chr01:1-1000 -> chr01	0	1000


### maskedSEQ2bed.pl
> Extract masked (lowercase and Ns) regions in bed format, quite slow ...



### bam2fasta.pl
> Convert SAM/BAM to FASTA/Q format

* Extract unmapped reads in BAM file and send it to blast+

		bam2fasta.pl -b example.bam -s "-f 4" | \
            		blastall -p blastn -d genome.fasta -e 1e-70 -a 1 -m 8 -F F -n T -o unmapped.blast.csv


### sam2fastq.pl
> Convert SAM/BAM to FASTQ format

* Extract paired reads

		sam2fastq.pl -i example.sam -o prefix


### extract_bam_pairs.pl
> Extract reads in pairs from SAM/BAM file

* Extract reads mapped to region "chr03 4000 5000"

		echo "chr03 4000 5000" | \
		    extract_bam_pairs.pl --input - --extend 2000 --samtools "-f2 -F3852 -q20" \
		    --min-insert 450 --max-insert 550 --max-clipping 2 --format fastq \
		    --bam sample1.bam sample2.bam --output extracted


**Note:** Difference between **bam2fasta.pl**, **sam2fastq.pl** and **extract_bam_pairs.pl**

		bam2fasta.pl ignore paired infos;
		sam2fastq.pl convert in pairs;
		extract_bam_pairs.pl extract all reads overlap certain regions and search for its pairs.

They can be actually integrated, but why 3 scripts? Because I forget the previous one when I started write a new one, and finally I got three ...

## Calculation

### calc_vcf_diversity.pl
> Calculate within- or between-groups diversities from a vcf file, an alternative choice of "vcftools --window-pi".

* Calculate pairwise diversity in non-overlapping windows

		calc_vcf_diversity.pl --vcf test.vcf.gz --window-size 500000 --no-group --no-within > pairs_divs.csv

* Calculate pairwise diversity in specified regions

		calc_vcf_diversity.pl --vcf test.vcf.gz --intervals test.regions --no-within --no-group > pairs_divs.csv


* Calculate within and between group diversity in non-overlapping windows

		calc_vcf_diversity.pl --vcf test.vcf.gz --window-size 500000 \
			--group test.sample_panel --all-pairs  --threads 2 > group_divs.csv

* Calculate within and between group diversity in specified regions

		calc_vcf_diversity.pl --vcf test.vcf.gz --intervals test.regions \
    			--group test.sample_panel --all-pairs --threads 2 > group_divs.csv


**Note:** values calculated directly from VCF file tends to underestimate the true diversity as it divide all sites rather than informative sites, one workaround is to use a three step calculation which also get the number of informative sites

* First count all informative sites (this requires the VCF file also contain the non-variant allele infos)

		count_vcf_pairs.pl --vcf sites.vcf.gz --intervals example.regions > pair_infos.csv

* Second, calculate differences between each pairs

		calc_vcf_diversity.pl --vcf snp.vcf.gz --no-group --no-within --intervals example.regions > pair_diffs.csv

* Finally, obtain the adjusted diversity by diffs / informative_sites

		calc_pop_divs.pl --query pair_diffs.csv --subject pair_infos.csv --min-info-perc 50 > pairs_divs.mean.csv


### count_vcf_pairs.pl
> A company script of calc_vcf_diversity.pl to count informative sites


### calc_pop_divs.pl
> Post-calculation script of results from calc_vcf_diversity.pl and count_vcf_pairs.pl


### slide_windows_count.pl
> Count by slide windows

* Count in 100kb windows, excluding scaffolds

		slide_windows_count.pl -i example.csv -w 100000 -l genome.fa.fai -e scaffold > example.100k.csv

* Count with different types

		slide_windows_count.pl -i example.csv -w 100000 -l genome.fa.fai -t > example.100k.csv

**Note:** This script use two different algorithms while counting non-overlapping and overlapping (when a step size is specified) windows. The first algorithm is relatively fast and does not require too much memory as it count on-the-fly while reading file; however, count non-overlapping windows needs first read the entire file, and the subsequent counting rely on the grep function is not
quite efficient, thus could be very slow and memory-intensive while process large files.


## Visualization

### paintGenomeBlocks.pl
> Plot blocks along each chromosome.

* Plot recombinat blocks with centromeres

		cat centromere.csv blocks.csv | \
		    paintGenomeBlocks.pl --input - --rows 0 1 2 3 4 --split-pairs \
		    --width 1600 --height 12000 --thickness 10 --chrom-distance 20 --block-distance 2 \
		    --output blocks --length genome.fa.fai --format svg \
		    --colors "A:strong_red2;B:strong_blue2;Centromere:grey" --sort-blocks Centromere sample1 sample2

* Plot variants (just extend it to a "block")

		awk 'BEGIN{OFS="\t"} !/\#/ {print "Example",$1,"SNP",$2-10,$2+10;}' example.vcf | \
		    paintGenomeBlocks.pl --input - \
		    --width 800 --height 600 --thickness 20 --chrom-distance 20 --block-distance 2 \
		    --output example.snp --length genome.fa.fai --colors "SNP:red"


### ggplot-windows-count.R
> Plot distributions with ggplot2

* Prepare examples

		bgzip -dc example1.vcf.gz | slide_windows_count.pl --input - --window-size 200000 > example1.w200k.csv
		bgzip -dc example2.vcf.gz | slide_windows_count.pl --input - --window-size 200000 > example2.w200k.csv

* Reformat

		cat example1.w200k.csv | \
		    awk 'BEGIN{OFS="\t"} !/##/ {if(/#/){print "Sample\tCHROM\tINTERVAL\tBIN_ID\tCOUNT";}
		    	else{print "Example1",$1,$2,$3,$4;}}' \
		    > example.w200k.txt
		
		cat example2.w200k.csv | \
		    awk 'BEGIN{OFS="\t"} !/#/ {print "Example2",$1,$2,$3,$4;}' \
		    >> example.w200k.txt
	
* Plot

		Rscript ggplot-windows-count.R example.w200k.txt example.w200k.tiff
    

### ggplot-windows-multi.R
> Alternative to ggplot-windows-count.R, designed for different input format

* Prepare examples

		echo -e "chr01\t122\tSNP\nchr01\t142\tSNP\nchr01\t254\tINDEL\n" | \
		    slide_windows_count.pl -i - -w 500 -t | cut -f 4 --complement | \
		    perl -ne 'next if (/^\#\#/); s/\#//g; print;' \
		    > counts.w500.txt

* Plot

		Rscript ggplot-windows-multi.R counts.w500.txt counts.w500.tiff


### ggplot-liner-fit.R
> Liner fit with ggplot2


### plot-scatter_with_lines.R
> Scatter plot


## Miscellaneous

### reference_align.pl
> Mapping with clustalw2(http://www.clustal.org/clustal2/) or muscle(http://www.drive5.com/muscle/).

Align sequences to a reference sequence, this was done by 2 steps:

**Step1:** align the first sequence to the reference sequence, and get
a expanded reference sequence with gaps inserted, then align the
second sequence to the new reference sequence, iterate this process
to generate a reference sequence expand all query sequences;

**Step2:** build a consensus sequence with nucleotide replaced in reference
sequence, and re-align all sequences to the consensus sequence.

* Align all subsequent sequences to the first sequence use clustalw2

		reference_align.pl -i example.fa > example.aln.fas

* Align with muscle

		reference_align.pl -i example.fa -a muscle > example.aln.fas



### parallel_baseml.pl
> Wrapper to run PAML (baseml) in parallel

* Run all models in a shot

		parallel_baseml.pl --input aln.fas --model 0 1 2 3 4 5 6 7 8 --align-length 100 > aln.csv


### parallel_codeml.pl
> Wrapper to run PAML (codeml) in parallel

* Run with multiple threads

		parallel_codeml.pl --ref ref.tbl --samples samples.tbl --query query.tsv --threads 6 > query.codeml.csv


### ascii2num
> Convert ASCII characters to Numbers

* Mainly used to inspect Fastq qualities

		ascii2num "AAAA"


These scripts are free softwares; you can redistribute it and/or modify it under the same terms as Perl itself.
