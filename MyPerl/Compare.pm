#!/usr/bin/perl -w
#
#  Compare.pm -- Sequences comparation
#
#  Author: Nowind
#  Created: 2010-10-09
#  Updated: 2012-11-04
#  Version: 1.1.0
#
#  Change logs:
#  Version 1.1.0 12/11/04: Add function check_codon_TsTv; remove dependency on bioperl.





package MyPerl::Compare;

use MyPerl::Convert;

@ISA    = qw(Exporter);

@EXPORT = qw(Clip_Seq_By_Ref Compare_Codon Count_Indel
             Format_Aln_Sequence Count_Diff_Num check_codon_TsTv);


$VERSION = '1.1.0';


##
## Compare codons in two aligned sequences
##
sub Compare_Codon
{
    my ($ref, $cmp) = @_;
    
    my @cmp = ();
    my $pos = 0;
    
    while ((length $ref) >3)
    {
        ++$pos;
        
        my $code1 = substr($ref, 0, 3, '');
        my $code2 = substr($cmp, 0, 3, '');
        
        if ($code1 ne $code2) {
            my $p1 = Codon2AA($code1);
            my $p2 = Codon2AA($code2);
            
            if ($p1 ne -1 && $p2 ne -1) {
                if ($p1 eq $p2) {
                    push @cmp, "$pos\t$code1->$code2\tSyn";
                }
                else {
                    push @cmp, "$pos\t$code1->$code2\tNon";
                }
            }
        }
    }

    return \@cmp;
}

##
## Find indel in sequences
##
sub Count_Indel
{
    my ($seq) = shift;
    
    my @indel = ();
    
    while ($seq =~ m/\-+/g)
    {
        my $pos = pos $seq;
        my $len = length $&;
        
        my $start = $pos - $len + 1;
        push @indel, "$start\t$len";
    }
    
    return \@indel;
}

##
## Use one of the sequences as the reference, and get rid of all '-' in the reference
##
sub Clip_Seq_By_Ref
{
    my ($rseq1, $rseq2, $limit) = @_;
    
    while ($$rseq1 =~ m/\-+/g)
    {
        my $pos = pos($$rseq1);
        my $len = length $&;
        
        next if ($limit && $len <= $limit);
        
        substr($$rseq1, $pos-$len, $len, '');
        substr($$rseq2, $pos-$len, $len, '');
    }
}

##
## Count the difference number of bases
##
sub Count_Diff_Num
{
    my ($rseq1, $rseq2) = @_;
    
    my @nt1 = split //, $$rseq1;
    my @nt2 = split //, $$rseq2;
    
    my $aln_length    = 0;
    my $diff_num      = 0;

    for (my $i=0; $i<=$#nt1; $i++)
    {
        my $cmp1 = $nt1[$i];
        my $cmp2 = $nt2[$i];
        
        if ($cmp1 ne '-' && $cmp2 ne '-') {
            $aln_length++;
            
            if($cmp1 ne $cmp2) {
                $diff_num++;
            }
        }
        
    }
    
    return ($diff_num, $aln_length);
}

## Get rid of '-' at the head after alignment
sub Format_Aln_Sequence
{
    my ($rseq1, $rseq2) = @_;
    
    my $pos1 = pos($$rseq1) if $$rseq1 =~ m/\w/g;
    my $pos2 = pos($$rseq2) if $$rseq2 =~ m/\w/g;
    
    my $pos  = $pos1 >= $pos2 ? $pos1 : $pos2;
    
    return 1 if ($pos == 1); # Neither of sequence has '-' in head
    
    substr($$rseq1, 0, $pos-1, '');
    substr($$rseq2, 0, $pos-1, '');
    
    Format_Aln_Sequence($rseq1, $rseq2);
}

## Return the shorter length
sub cmp_seq_length
{
    my ($gene1, $gene2, $cds_seq) = @_;
    
    my $length1 = length $cds_seq->{$gene1}->seq;
    my $length2 = length $cds_seq->{$gene2}->seq;
    
    my $short = $length1 >= $length2 ? $length2 : $length1;
    
    return $short;
}


sub check_codon_TsTv     # Transition or transverion
{
    my ($cmp1, $cmp2) = @_;
    
    my %Ts = (
                 'A'  => 'G',
                 'G'  => 'A',
                 'T'  => 'C',
                 'C'  => 'T',
              );
    
    my @cmp1 = split //, $cmp1;
    my @cmp2 = split //, $cmp2;
    
    my $diff = 0;
    my $type = '';
    
    for my $n (0,1,2)
    {
        if ($cmp1[$n] ne $cmp2[$n]) { 
            $diff++;
            
            return 'N' unless exists($Ts{$cmp1[$n]});
            
            if ($Ts{$cmp1[$n]} ne $cmp2[$n]) { # Transversion
                $type = 'Tv';
            }
            else {                             # Transition
                $type = 'Ts';
            }
        }
    }
    
    return 'N' unless $diff == 1; # Direct Neighbour
    
    return $type;    
}


1;