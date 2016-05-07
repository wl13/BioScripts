#!/usr/bin/perl -w
#
#   parallel_codeml.pl -- A wrapper to run codeml in parallel.
#                          
#
#   Author: Nowind
#   Created: 2012-05-31
#   Updated: 2015-05-01
#   Version: 1.0.0
#
#   Change logs:
#   Version 1.0.0 15/05/01: The initial version.




=head1 NAME

parallel_baseml.pl


=head1 SYNOPSIS

  parallel_baseml.pl --help/?

=head1 DESCRIPTION

A wrapper A wrapper to run baseml in parallel.

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
my $VERSION = '1.0.0';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";
my $SOURCE  = (scalar localtime()) . " Version: $VERSION";


my $max_threads  = 1;
my @models       = ();
my $min_aln_len  = 100;
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

$0  -- Run baseml in parallel.

Version: $VERSION

Usage:   perl $0 [options]

Options:
    -i, --input    <filename>
        input alignment files in fasta format, required

    -o, --output  <dirname>
        output filename, default to STDOUT
    
    -m, --model    <numbers>
        specifies one or more models of nucleotide substitution, 0:JC69,
        1:K80, 2:F81, 3:F84, 4:HKY85, 5:T92, 6:TN93, 7:REV, 8:UNREST
        [default: 0]
    
    -a, --align-length <int>
        minimum alignment length to process, [default: 100 (bp)]
    
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

unless(@models > 0) {
    @models = (0);
}

my %substitution_models = (
                            0 => "JC69",
                            1 => "K80",
                            2 => "F81",
                            3 => "F84",
                            4 => "HKY85",
                            5 => "T92",
                            6 => "TN93",
                            7 => "REV",
                            8 => "UNREST",
                           );

my @model_names = ();
for my $model (@models)
{
    push @model_names, $substitution_models{$model};
}
my $model_names = join "\t", @model_names;


##
## read reference sequences
##
my $pm = new Parallel::ForkManager($max_threads) if $max_threads > 1;

print STDOUT "$HEADER##" . (scalar localtime()) . "\n";
print "#ref_id\tcmp_id\tref_len\tcmp_len\ttaln_length\t$model_names\n";

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
    
    my $aln_len = length($ref_seq);
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
        write_baseml_ctl("$tmp_in.ctl", $tmp_in, $model);
        
        my $run_baseml = "baseml $tmp_in.ctl";
        
        my $return = system "$run_baseml >/dev/null 2>&1";
        
        my $result_fh = getInputFilehandle("$tmp_in.out");
        while (<$result_fh>)
        {
            next unless(m/$cmp_id\s+(\d+\.\d+)(\(|$)/);
            
            push @results, $1;
            
            last;
        }
    }
    
    my $results = join "\t", @results;
    
    print "$ref_id\t$cmp_id\t$ref_len\t$cmp_len\t$aln_len\t$results\n";

    
    if ($max_threads > 1) {
        $pm->finish;
    }
}
print STDERR "\tdone!\n";


######################### Sub #########################


sub write_baseml_ctl
{
    my ($config_file, $seqfile, $select_model) = @_;
    
    my $treefile = $seqfile . ".tree";
    my $outfile  = $seqfile . ".out";
    
    open (OUT, "> $config_file") || die $!;
    print OUT <<EOF;
      seqfile = $seqfile
*    treefile = $treefile

      outfile = $outfile       * main result file
        noisy = 0   * 0,1,2,3: how much rubbish on the screen
      verbose = 0   * 1: detailed output, 0: concise output
      runmode = 2   * 0: user tree;  1: semi-automatic;  2: automatic
                    * 3: StepwiseAddition; (4,5):PerturbationNNI 

        model = $select_model   * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85
                                * 5:T92, 6:TN93, 7:REV, 8:UNREST, 9:REVu; 10:UNRESTu

        Mgene = 0   * 0:rates, 1:separate; 2:diff pi, 3:diff kapa, 4:all diff

*        ndata = 100
        clock = 0   * 0:no clock, 1:clock; 2:local clock; 3:CombinedAnalysis
    fix_kappa = 0   * 0: estimate kappa; 1: fix kappa at value below
        kappa = 5  * initial or fixed kappa

    fix_alpha = 0   * 0: estimate alpha; 1: fix alpha at value below
        alpha = 0.5   * initial or fixed alpha, 0:infinity (constant rate)
       Malpha = 0   * 1: different alpha's for genes, 0: one alpha
        ncatG = 5   * # of categories in the dG, AdG, or nparK models of rates
        nparK = 0   * rate-class models. 1:rK, 2:rK&fK, 3:rK&MK(1/K), 4:rK&MK 

        nhomo = 0   * 0 & 1: homogeneous, 2: kappa for branches, 3: N1, 4: N2
        getSE = 0   * 0: don't want them, 1: want S.E.s of estimates
 RateAncestor = 0   * (0,1,2): rates (alpha>0) or ancestral states

   Small_Diff = 7e-6
    cleandata = 1  * remove sites with ambiguity data (1:yes, 0:no)?
*        icode = 0  * (with RateAncestor=1. try "GC" in data,model=4,Mgene=4)
*  fix_blength = -1  * 0: ignore, -1: random, 1: initial, 2: fixed
       method = 0  * Optimization method 0: simultaneous; 1: one branch a time
EOF
    close OUT;
}



