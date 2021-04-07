#!/usr/bin/perl -w
#
#   parallel_codeml.pl -- A wrapper to run codeml in parallel.
#                          
#
#   Author: Nowind
#   Created: 2012-05-31
#   Updated: 2021-04-07
#   Version: 1.0.1
#
#   Change logs:
#   Version 1.0.0 18/03/30: The initial version.
#   Version 1.0.1 21/04/07: Correct some descriptions.




=head1 NAME

parallel_codeml.pl


=head1 SYNOPSIS

  parallel_codeml.pl --help/?

=head1 DESCRIPTION

A wrapper to run codeml in parallel.

=cut




use strict;

use Data::Dumper;
use Getopt::Long;
use File::Find::Rule;
use File::Basename;
use File::Temp;
use Parallel::ForkManager;

use MyPerl::FileIO qw(:all);
use MyPerl::Convert;
use MyPerl::Compare;

######################### Main #########################

my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.0.1';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";
my $SOURCE  = (scalar localtime()) . " Version: $VERSION";


my $max_threads  = 1;
my @models       = ();
my $min_aln_len  = 30;
my ($aln_file, $output, $out_directly, $show_help);
GetOptions(
            "input=s"            => \$aln_file,
            
            "output=s"           => \$output,
            
            "threads=i"          => \$max_threads,
            
            "models=i{,}"        => \@models,
            
            "align-length=i"     => \$min_aln_len,
            
            "help|?"             => \$show_help,
           );

unless( !$show_help && $aln_file ) {
    print <<EOF;

$0  -- Run codeml in parallel.

Version: $VERSION

Usage:   perl $0 [options]

Options:
    -i, --input    <filename>
        input alignment files in fasta format, required

    -o, --output  <dirname>
        output filename, default to STDOUT
    
    -m, --model    <numbers>
        currently not used
    
    -a, --align-length <int>
        minimum alignment length to process, [default: 30 (bp)]
    
    -t, --threads  <int>
        how many data threads should be allocated to running this analysis
        [default: 1]
    
    -?, --help
        show this help message
    
EOF

    exit(1);
}

$|++;


print STDERR "# $0 v$VERSION\n# " . (scalar localtime()) . "\n";


if ($output) {
    open (STDOUT, "> $output") || die $!;
}

#unless(@models > 0) {
    @models = (0);
#}



##
## read reference sequences
##
my $pm = new Parallel::ForkManager($max_threads) if $max_threads > 1;

print STDOUT "$HEADER##" . (scalar localtime()) . "\n";
print "#Sample1\tSample2\taln_len\tref_len\tcmp_len\tt\tS\tN\tdNdS\tdN\tdS\n";


my $curr_job_num   = 0;
my $aln_fh = getInputFilehandle($aln_file);
while (<$aln_fh>)
{
    chomp(my $ref_id  = $_);
    chomp(my $ref_seq = <$aln_fh>);
    chomp(my $cmp_id  = <$aln_fh>);
    chomp(my $cmp_seq = <$aln_fh>);
    
    $ref_id =~ s/^\>//;
    $cmp_id =~ s/^\>//;
    
    my $aln_len = ($ref_seq =~ tr/ATGCatgc-/ATGCatgc-/);
    my $ref_len = $aln_len - ($ref_seq =~ tr/-/-/);
    my $cmp_len = $aln_len - ($cmp_seq =~ tr/-/-/);
    
    next if ($min_aln_len && $aln_len < $min_aln_len);
    
    $curr_job_num++;
    
    my $pid  = $pm->start and next if ($max_threads > 1);
    
    print STDERR "\r>>Start calculating ... $curr_job_num";
    
    my $tmp_dir   = File::Temp->newdir(CLEANUP => 1);
    my $tmp_in_fh = File::Temp->new(DIR => $tmp_dir, UNLINK => 1);
    my $tmp_in    = $tmp_in_fh->filename;
    
    open (my $tmp_fh, "> $tmp_in") || die $!;
    
    print {$tmp_fh} "2\t$aln_len\n";
    print {$tmp_fh} "$ref_id\n      $ref_seq\n";
    print {$tmp_fh} "$cmp_id\n      $cmp_seq\n";
    
    my @results = ();
    for my $model (@models)
    {
        write_codeml_ctl("$tmp_in.ctl", $tmp_in, $model);
        
        my $run_codeml = "codeml $tmp_in.ctl";
        
        my $return = system "$run_codeml >/dev/null 2>&1";
        
        my $result_fh = getInputFilehandle("$tmp_in.out");
        while (<$result_fh>)
        {
            ##print STDOUT;
            
            next unless(/t=\s+(.*?)\s+S=\s+(.*?)\s+N=\s+(.*?)\s+dN\/dS=\s+(.*?)\s+dN\s+=\s+(.*?)\s+dS\s+=\s+(.*?)\s+/);
            
            my ($t, $S, $N, $dNdS, $dN, $dS) = ($1, $2, $3, $4, $5, $6);
            
            push @results, $1;
            
            print STDOUT "$ref_id\t$cmp_id\t$aln_len\t$ref_len\t$cmp_len\t$t\t$S\t$N\t$dNdS\t$dN\t$dS\n";
            
            last;
        }
    }
    
    if ($max_threads > 1) {
        $pm->finish;
    }
}
print STDERR "\tdone!\n";


######################### Sub #########################


sub write_codeml_ctl
{
    my ($config_file, $seqfile, $select_model) = @_;
    
    my $treefile = $seqfile . ".tree";
    my $outfile  = $seqfile . ".out";
    
    open (OUT, "> $config_file") || die $!;
    print OUT <<EOF;
      seqfile = $seqfile  * sequence data filename
     treefile = $treefile * tree structure file name
      outfile = $outfile  * main result file name

        noisy = 0  * 0,1,2,3,9: how much rubbish on the screen
      verbose = 0  * 0: concise; 1: detailed, 2: too much
      runmode = -2  * 0: user tree;  1: semi-automatic;  2: automatic
                   * 3: StepwiseAddition; (4,5):PerturbationNNI; -2: pairwise

      seqtype = 1  * 1:codons; 2:AAs; 3:codons-->AAs
    CodonFreq = 2  * 0:1/61 each, 1:F1X4, 2:F3X4, 3:codon table

*        ndata = 10
        clock = 0  * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis
       aaDist = 0  * 0:equal, +:geometric; -:linear, 1-6:G1974,Miyata,c,p,v,a
   aaRatefile = dat/jones.dat  * only used for aa seqs with model=empirical(_F)
                   * dayhoff.dat, jones.dat, wag.dat, mtmam.dat, or your own

        model = 2
                   * models for codons:
                       * 0:one, 1:b, 2:2 or more dN/dS ratios for branches
                   * models for AAs or codon-translated AAs:
                       * 0:poisson, 1:proportional, 2:Empirical, 3:Empirical+F
                       * 6:FromCodon, 7:AAClasses, 8:REVaa_0, 9:REVaa(nr=189)

      NSsites = 0  * 0:one w;1:neutral;2:selection; 3:discrete;4:freqs;
                   * 5:gamma;6:2gamma;7:beta;8:beta&w;9:beta&gamma;
                   * 10:beta&gamma+1; 11:beta&normal>1; 12:0&2normal>1;
                   * 13:3normal>0

        icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below
        Mgene = 0
                   * codon: 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff
                   * AA: 0:rates, 1:separate

    fix_kappa = 0  * 1: kappa fixed, 0: kappa to be estimated
        kappa = 2  * initial or fixed kappa
    fix_omega = 0  * 1: omega or omega_1 fixed, 0: estimate 
        omega = .4 * initial or fixed omega, for codons or codon-based AAs

    fix_alpha = 1  * 0: estimate gamma shape parameter; 1: fix it at alpha
        alpha = 0. * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0  * different alphas for genes
        ncatG = 8  * # of categories in dG of NSsites models

        getSE = 0  * 0: don't want them, 1: want S.E.s of estimates
 RateAncestor = 1  * (0,1,2): rates (alpha>0) or ancestral states (1 or 2)

   Small_Diff = .5e-6
    cleandata = 1  * remove sites with ambiguity data (1:yes, 0:no)?
*  fix_blength = -1  * 0: ignore, -1: random, 1: initial, 2: fixed
       method = 0  * Optimization method 0: simultaneous; 1: one branch a time

* Genetic codes: 0:universal, 1:mammalian mt., 2:yeast mt., 3:mold mt.,
* 4: invertebrate mt., 5: ciliate nuclear, 6: echinoderm mt., 
* 7: euplotid mt., 8: alternative yeast nu. 9: ascidian mt., 
* 10: blepharisma nu.
* These codes correspond to transl_table 1 to 11 of GENEBANK.
EOF
    close OUT;
}



