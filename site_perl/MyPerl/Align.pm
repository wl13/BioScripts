#!/usr/bin/perl -w
#
#  Align.pm -- Sequences alignment operations.
#
#  Author: Nowind
#  Created: 2010-10-09
#  Updated: 2016-09-14
#  Version: 1.2.2
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
#  Version 1.2.0 16/09/04: Updated: remove dependency on Bioperl.
#  Version 1.2.1 16/09/13: Updated: append the stop codon while recovering the original DNA
#                          sequences through aligned proteins.
#  Version 1.2.2 16/09/14: Bug fixed: failed to generate tmp results while align with muscle.


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
use Data::Dumper;

use MyPerl::FileIO qw(:all);
use MyPerl::Convert qw(:all);

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
        )
    ]
);

@EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );
@EXPORT    = qw();


$MyPerl::FileIO::VERSION = '1.2.2';


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
                     . "-maxiters $maxiters -quiet";        

    my $cmdline = ($prog eq 'clustalw2') ? $run_clustaw2 : $run_muscle;
       $cmdline = "$cmdline $params " if $params;
    
    if ($log) {
        $cmdline .= " >>$log 2>&1";
    }
    elsif ($prog eq 'clustalw2') {
        $cmdline .= " >/dev/null 2>&1";
    }
    else {
        $cmdline .= " 2>/dev/null";
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
    my $tmp_in_fh = File::Temp->new( DIR => $tmp_dir );
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



=head2 align_codons

    About   : Codon-guided sequences alignments
    Usage   : my $ra_aln_seqs = $aln->align_seqs($seq1, $seq2, ...);
    Args    : Sequences need to be aligned, require at least two sequences
    Returns : Reference to an array contains aligned sequences

=cut
sub align_codons
{
    my ($obj, @seqs) = @_;
    
    my $tmp_dir   = File::Temp->newdir();
    my $tmp_in_fh = File::Temp->new( DIR => $tmp_dir );
    my $tmp_in    = $tmp_in_fh->filename;
    my $tmp_ot    = $tmp_in . "_aln";

    my %original_nts = ();
    for (my $i=0; $i<@seqs; $i++)
    {
        @{$original_nts{$i}} = split //, $seqs[$i];
        
        my $prot = Translate($seqs[$i]);
        print $tmp_in_fh ">$i\n$prot\n";
    }
    
    $obj->{input}  = $tmp_in;
    $obj->{output} = $tmp_ot;
    $obj->{type}   = 'PROTEIN';
    
    my $flag = $obj->start($obj);
    
    return -1 if ($flag);
    
    my @aln_prots = ();
    parse_fasta_SEQs(\@aln_prots, $tmp_ot);
    
    my @aln_seqs  = ();
    
    my %aligned_aas = ();
    for (my $i=0; $i<@aln_prots; $i++)
    {
        @{$aligned_aas{$i}} = split //, $aln_prots[$i];
        
        my $rs_aln_seq = recover_codon2nt(\@{$original_nts{$i}}, $aligned_aas{$i});
        
        push @aln_seqs, $$rs_aln_seq;
    }
    
    return \@aln_seqs;
}


=head2 recover_codon2nt

    About   : Recover Nucleotide Sequence
    Usage   : my $recovered_seqs = recover_codon2nt($nt, $aa);
    Args    : Array reference to original nucleotides;
              Array reference to proteins.
    Returns : Reference to recovered sequences.

=cut
sub recover_codon2nt
{
    my ($ra_nt, $ra_aa) = @_;
    
    my $seq = '';
    
    for (my $i=0; $i<=$#{$ra_aa};$i++)
    {
        if ($ra_aa->[$i] eq '-') {
            $seq .= '-' x 3;
        } else {
            $seq .= shift @{$ra_nt};
            $seq .= shift @{$ra_nt};
            $seq .= shift @{$ra_nt};
        }
    }
    
    if (@{$ra_nt} > 0) {
        $seq .= join '', @{$ra_nt};
    }
    
    return \$seq;
}

1;


=head1 VERSION

1.2.2

=head1 AUTHOR

Nowind, noitulove9891@gmail.com

=head1 COPYRIGHT

Copyright (c) Nowind's Area. All rights reserved. This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself. 


=cut