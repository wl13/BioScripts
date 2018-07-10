#!/usr/bin/perl -w
#
#  Align.pm -- Sequences alignment operations.
#
#  Author: Nowind
#  Created: 2010-10-09
#  Updated: 2016-03-17
#  Version: 1.1.3
#
#  Change logs:
#  Version 1.0.0 10/10/09: The initial version.
#  Version 1.1.0 12/12/31: Rearrange full constructions. 
#  Version 1.1.1 12/01/02: Add option "parameters" to set additional program specified
#                          parameters.
#  Version 1.1.2 12/03/30: Use File::Temp to handle temporary files safely in order to
#                          aviod collisions while two processes access same temporary file.
#  Version 1.1.3 16/03/17: Bug fixed: could not write result while muscle is specified without
#                                     "--quiet" option.
#                          Updated: Add module "Data::Dumper".





=head1 NAME

MyPerl::Align - Local perl module for alignment operations


=head1 SYNOPSIS

  use MyPerl::Align;
  
  ## creat a new align object
  my $aln = MyPerl::Align->new(
                                input  => "input.fasta",
                                output => "output.fasta",
                                prog   => "clustalw2",
                                type   => "DNA",
                               );
  
  ## start alignment
  $aln->start;  
  
  ## directly align sequences
  my $rh_aln_seqs = $aln->align_seqs($seq1, $seq2, ...)

=head1 DESCRIPTION

Local perl module invoke clustalw2 and muscle to align sequences in fasta format

=cut

package MyPerl::Align;

use strict;

require Exporter;

use File::Temp;
use Bio::SeqIO;
use Data::Dumper;

use MyPerl::FileIO qw(:all);


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
            Align_DNA
            Align_By_Codon
        )
    ]
);

@EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );
@EXPORT    = qw();


$MyPerl::FileIO::VERSION = '1.1.3';


=head1 METHODS

=head2 new

    About   : Create an object
    Usage   : my $aln = MyPerl::Align->new(%settings);
    Args    : Using defined settings about:
                input:    input file
                output:   output file
                prog:     align program, clustalw2 or muscle
                type:     sequence type, DNA or PROTEIN
                dir:      output directory, default to current directory
                maxiters: max iterations(only valid while using muscle)
                log:      log file
                params:   additional parameters
    Returns : An object to run the alignment

=cut
sub new
{
    my ($class, %user_set)  = @_;
    
    my $default_set = {
        'input'        => undef,
        'output'       => undef,
        'prog'         => 'clustalw2', 
        'type'         => 'DNA',
        'dir'          => '.',
        'maxiters'     => 3,
        'log'          => '',
    };
    
    for (keys %{$default_set})
    {
        unless (defined $user_set{$_}) {
            $user_set{$_} = $default_set->{$_};
        }
    }    

    unless( -d $user_set{dir} ) { mkdir $user_set{dir} || die $!; }

    return bless \%user_set, $class;
}


=head2 start

    About   : Run alignment
    Usage   : my $obj = $aln->start;
    Args    : Null
    Returns : Null

=cut
sub start
{
    my ($settings) = @_;
    
    my $aln_in   = $settings->{input};
    my $aln_out  = $settings->{output};
    my $rh_ids   = $settings->{id};
    my $rh_seqs  = $settings->{seq};
    my $prog     = $settings->{prog};
    my $type     = $settings->{type};
    my $dir      = $settings->{dir};
    my $maxiters = $settings->{maxiters};
    my $log      = $settings->{log};
    my $params   = $settings->{params};
    
    my $run_clustaw2 = "clustalw2 -infile=$aln_in -type=$type "
                     . "-align  -outfile=$aln_out -output=fasta "
                     . "-outorder=input";
    
    my $run_muscle   = "muscle -in $aln_in -out $aln_out "
                     . "-maxiters $maxiters";        

    my $cmdline = ($prog eq 'clustalw2') ? $run_clustaw2 : $run_muscle;
       $cmdline = "$cmdline $params " if $params;
    
    if ($log) {
        $cmdline .= " >>$log 2>&1"
    }
    else {
        $cmdline .= " >/dev/null 2>&1"
    }
    
    ###print STDERR "$cmdline\n";exit;
    
    system "$cmdline";
    
    return 0;
}


=head2 align_seqs

    About   : Directly align sequences
    Usage   : my $ra_aln_seqs = $aln->align_seqs($seq1, $seq2, ...);
    Args    : Sequences need to be aligned, require at least two sequences
    Returns : Reference to an array contains aligned sequences

=cut
sub align_seqs
{
    my ($obj, @seqs) = @_;
    
    my $tmp_dir   = File::Temp->newdir();
    my $tmp_in_fh = File::Temp->new( DIR => $tmp_dir->dirname );
    my $tmp_in    = $tmp_in_fh->filename;
    my $tmp_ot    = $tmp_in . "_aln";
    
    for (my $i=0; $i<@seqs; $i++)
    {
        print $tmp_in_fh ">$i\n$seqs[$i]\n";
    }
    
    $obj->{input}  = $tmp_in;
    $obj->{output} = $tmp_ot;
    
    my $flag = $obj->start($obj);
    
    return -1 if ($flag);
    
    my @aln_seqs = ();
    parse_fasta_SEQs(\@aln_seqs, $tmp_ot);
    ###print STDERR Dumper($tmp_ot);
    return \@aln_seqs;
}



##
## Align Codons
##
sub Align_By_Codon
{
    my ($seqobj1, $seqobj2, $aln_out) = @_;
    
    my $gene1  = $seqobj1->id;
    my $gene2  = $seqobj2->id;
    
    my @nt1 = split //, $seqobj1->seq;
    my @nt2 = split //, $seqobj2->seq;
    
    my $tmp_dir   = File::Temp->newdir();
    my $tmp_in_fh = File::Temp->new( DIR => $tmp_dir->dirname );
    my $tmp_in    = $tmp_in_fh->filename;
    
    my $seqout    = Bio::SeqIO->new( -file   => "> $tmp_in", -format => 'fasta' );

    my $trans1 = $seqobj1->translate();
    my $trans2 = $seqobj2->translate();
    $seqout->write_seq($trans1);
    $seqout->write_seq($trans2);        
    
    my $aln_file = $tmp_in . "_aln";
    my $command  = "clustalw2 -INFILE=$tmp_in -TYPE=PROTEIN -ALIGN -OUTFILE=$aln_file -OUTPUT=FASTA -QUIET";
    system "$command";
    
    
    my %aln_content = ();
    my $aln_in  = Bio::SeqIO->new( -file => "$aln_file", -format => 'fasta' );
    my $aln_seq1 = $aln_in->next_seq();
    my $aln_seq2 = $aln_in->next_seq();
    $aln_content{$gene1} = [split //, $aln_seq1->seq];
    $aln_content{$gene2} = [split //, $aln_seq2->seq];
    
    
    my $rseq1 = Recover_NT_Sequence(\@nt1, $aln_content{$gene1});
    my $rseq2 = Recover_NT_Sequence(\@nt2, $aln_content{$gene2});

    if ($aln_out) {
        open (FA, ">> $aln_out") or die $!;
        print FA ">$gene1\n";
        print FA $$rseq1 . "\n";
        print FA ">$gene2\n";
        print FA $$rseq2 . "\n";
        close FA;
    }
    
    return ($rseq1, $rseq2);
}

#Recover Nucleotide Sequence
sub Recover_NT_Sequence
{
    my ($nt, $aa) = @_;
    
    my $seq = '';
    
    for (my $i=0; $i<=$#{$aa};$i++)
    {
        if ($aa->[$i] eq '-') {
            $seq .= '-' x 3;
        } else {
            $seq .= shift @{$nt};
            $seq .= shift @{$nt};
            $seq .= shift @{$nt};
        }
    }
    
    return \$seq;
}

1;


=head1 VERSION

1.1.3

=head1 AUTHOR

Nowind, noitulove9891@gmail.com

=head1 COPYRIGHT

Copyright (c) Nowind's Area. All rights reserved. This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself. 


=cut
