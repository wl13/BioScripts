#!/usr/bin/perl -w
#
#   reference_align.pl -- Align sequences to a reference sequence.
#
#
#   Author: Nowind
#   Created: 2012-02-21
#   Updated: 2016-04-19
#   Version: 1.1.1
#
#   Change logs:
#   Version 1.0.0 12/12/31: The initial version.
#   Version 1.0.1 13/07/02: Update usage infos.
#   Version 1.0.2 13/07/04: Correct id name for consensus sequence.
#   Version 1.0.3 13/09/07: Add option "--consensus" to output consensus sequence.
#   Version 1.1.0 16/03/17: Updated: 1) add support for choose muscle as an alternative aligner;
#                                    2) add options to set parameters.
#   Version 1.1.1 16/04/19: Updated: 1) add more comments; 2) remove some useless codes.


use strict;

use Data::Dumper;
use Getopt::Long;

use MyPerl::FileIO qw(:all);
use MyPerl::Align;

##################### Main ####################


my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.1.1';
my $HEADER  = "# $CMDLINE\n# Version: $VERSION\n";

my $SOURCE  = (scalar localtime()) . " Version: $VERSION";

my $aligner  = 'clustalw2';
my $maxiters = 3;
my ($input, $output, $params, $out_cns);
GetOptions(
            "input=s"          => \$input,
            "output=s"         => \$output,
            
            "aligner=s"        => \$aligner,
            "params=s"         => \$params,
            "maxiters=i"       => \$maxiters,
            
            "consensus"        => \$out_cns,
           );

my $show_help = ($CMDLINE =~ /\-help/) ? 0 : 1;

unless( $input && $show_help ) {
    print <<EOF;

$0  -- Align sequences to a reference sequence, this was done by 2 steps,

    Step1: align the first sequence to the reference sequence, and get
    a expanded reference sequence with gaps inserted, then align the
    second sequence to the new reference sequence, iterate this process
    to generate a reference sequence expand all query sequences;
    
    Step2: build a consensus sequence with nucleotide replaced in reference
    sequence, and re-align all sequences to the consensus sequence.

Version: $VERSION

Usage:   perl $0 [options]

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
    
EOF

    exit(1);
}




$|++;

if ($output) {
    open (STDOUT, "> $output") || die $!; 
}

unless($params) {
    if ($aligner eq 'clustalw2') {
        $params = '-gapopen=15 -gapext=6.66';
    }
    if ($aligner eq 'muscle') {
        $params = '-quiet';
    }
}


print STDERR "# $0 v$VERSION\n# " . (scalar localtime()) . "\n";


print STDERR ">> Start align sequences in $input ... ";
pairwise_align($input);
print STDERR "done!\n";

print STDERR "# " . (scalar localtime()) . "\n";


######################### Sub #########################


=head2 pairwise_align

    About   : Align each sequence to a reference sequence
    Usage   : pairwise_align($fasta_file);
    Args    : Sequence file in fasta format
    Returns : Null

=cut
sub pairwise_align
{
    my ($in) = @_;
    
    my @SEQs  = ();
    
    my @ids = parse_fasta_SEQs(\@SEQs, $in);
    
    ###print format_fasta_SEQs($ids[0], \$SEQs[0]);
    
    ## step1: generate a longest gapped reference sequence
    my $tmp_ref = $SEQs[0];
    
    for (my $i=1; $i<@ids; $i++)
    {
        my $aln = MyPerl::Align->new(prog     => $aligner,
                                     type     => 'DNA',
                                     params   => $params,
                                     maxiters => $maxiters);
        
        my $rh_aln_seqs = $aln->align_seqs($tmp_ref, $SEQs[$i]);
        
        $tmp_ref = $rh_aln_seqs->[0];
        
        ###print format_fasta_SEQs($ids[0],  \$rh_aln_seqs->[0]);
        ###print format_fasta_SEQs($ids[$i], \$rh_aln_seqs->[1]);
    }
    
    print format_fasta_SEQs($ids[0], \$tmp_ref);
    
    ## step2: re-align all sequences to the longest reference sequence
    my %aln_nts = ();
    
    for (my $i=1; $i<@ids; $i++)
    {
        my $aln = MyPerl::Align->new(prog     => $aligner,
                                     type     => 'DNA',
                                     params   => $params,
                                     maxiters => $maxiters);
        
        my $rh_aln_seqs = $aln->align_seqs($tmp_ref, $SEQs[$i]);
        
        my @nts = split //, $rh_aln_seqs->[1];
        
        ## count the frequncy of each nucleotide in position j
        for (my $j=0; $j<@nts; $j++)
        {
            next if ($nts[$j] eq '-');
            $aln_nts{$j}->{$nts[$j]}++;
        }
        
        ###print format_fasta_SEQs($ids[0],  \$rh_aln_seqs->[0]);
        ###print format_fasta_SEQs($ids[$i], \$rh_aln_seqs->[1]);
        ###print "###\n";
    }
    
    ## step3: choose the nucleotide with highest frequncy in each position to build
    ## a consensus sequence, if absent, use reference instead
    my @ref_nts = split //, $tmp_ref;
    my @cns_nts = ();
    for (my $j=0; $j<(length $tmp_ref); $j++)
    {
        my $major_nt = $ref_nts[$j];
        if( $aln_nts{$j} ) {
            $major_nt   = (sort {$aln_nts{$j}->{$a} <=> $aln_nts{$j}->{$b}} (keys %{$aln_nts{$j}}))[-1];
        }
        
        push @cns_nts, $major_nt;
    }
    
    my $cns_seq = join '', @cns_nts;
    
    
    if ($out_cns) {
        print format_fasta_SEQs("$ids[0]-cns", \$cns_seq);
    }
    
    
    ## step4: re-align each sequence to the consensus sequence
    for (my $i=1; $i<@ids; $i++)
    {
        my $aln = MyPerl::Align->new(prog     => $aligner,
                                     type     => 'DNA',
                                     params   => $params,
                                     maxiters => $maxiters);
        
        my $rh_aln_seqs = $aln->align_seqs($cns_seq, $SEQs[$i]);
        
        ###print format_fasta_SEQs($ids[0],  \$rh_aln_seqs->[0]);
        print format_fasta_SEQs($ids[$i], \$rh_aln_seqs->[1]);
        ###print "###\n";
    }
}


