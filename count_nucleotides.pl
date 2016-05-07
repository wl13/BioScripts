#!/usr/bin/perl -w
#
#   count_nucleotides.pl -- Count of nucleotides contents.
#
#
#   Author: Nowind
#   Created: 2012-05-30
#   Updated: 2015-10-20
#   Version: 2.0.0
#
#   Change logs:
#   Version 1.0.0 15/02/06: The initial version.
#   Version 1.0.1 15/04/24: Remove unused options.
#   Version 1.0.2 15/05/09: Reduce verbose; add test for too short (<3 bp) sequences.
#   Version 1.0.3 15/08/18: Add match strings to skip scaffolds.
#   Version 2.0.0 15/10/20: Add support for stats of dinucleotide contents.



=head1 NAME

count_nucleotides.pl


=head1 SYNOPSIS

  count_nucleotides.pl --help/?

=head1 DESCRIPTION

Analysis of codon changes.

=cut


use strict;

use Data::Dumper;
use Getopt::Long;
use File::Find::Rule;
use File::Basename;
use Parallel::ForkManager;
use Statistics::Descriptive;

use MyPerl::FileIO qw(:all);
use MyPerl::Convert;

################### Main #################

my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '2.0.0';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";
my $SOURCE  = (scalar localtime()) . " Version: $VERSION";


my $count_type  = 'triplet';
my ($fasta_file, $output);
GetOptions(
            "fasta=s"                => \$fasta_file,
            
            "O|output=s"             => \$output,
            
            "context=s"              => \$count_type,
           );

unless( $fasta_file ) {
    print <<EOF;

$0  -- Analysis of codon changes

Version: $VERSION

Usage:   perl $0 [--fasta FILE | STDIN] [Options]

Options:

    -f, --fasta    <filename>
        sequence file in fasta format, required
    -O, --output   <filename>
        output file, default to STDOUT
        
    -c, --context  <string>
        nucleotide context to count, can be set as "triplet" or
        "dinucleotide", default: triplet
        
EOF

    exit(0);
}

$|++;

print STDERR "# $0 v$VERSION\n# " . (scalar localtime()) . "\n";





if ($output) {
    open (STDOUT, "> $output") || die $!;
}


print STDERR ">> Start reading $fasta_file ... ";
my @SEQs = ();
my @IDs  = parse_fasta_SEQs(\@SEQs, $fasta_file);
print STDERR "done!\n";

if ($count_type eq 'dinucleotide') {
    count_dinucleotide($fasta_file);
}
elsif ($count_type eq 'triplet') {
    count_triplets($fasta_file);
}


print STDERR "# " . (scalar localtime()) . "\n";

######################### Sub #########################


=head2 count_dinucleotide

    About   : Count dinucleotide content in fasta file.
    Usage   : count_dinucleotide($fasta_file);
    Args    : Sequences file in fasta format.
    Returns : Null

=cut
sub count_dinucleotide
{
    my ($in) = @_;
    
    ###my $seq = 'AACAATACGA';
    ###my @dints = unpack("(A2)*", substr($seq, 2));
    ###
    ###my $comp_seq = $seq;
    ###   $comp_seq =~ tr/ATGC/TACG/;
    ###
    ###my @dints2 = unpack("(A2)*", $comp_seq);
    ###
    ###print STDERR Dumper(@dints, @dints2);exit;
    
    my %dints_all = ();
    for (my $i=0; $i<@IDs; $i++)
    {
        next if ($IDs[$i] =~ /(mitochondria|chloroplast|scaffold)/);
        
        print STDERR "\r>> Start counting ... $IDs[$i]";
        
        next if (length($SEQs[$i]) < 2);
        
        my $seq = uc ($SEQs[$i]);
        
        my @dints_forward1 = unpack("(A2)*", $seq);
        my @dints_forward2 = unpack("(A2)*", substr($seq, 1));
        
        $dints_all{$_}->{forward1} ++ for @dints_forward1;
        $dints_all{$_}->{forward2} ++ for @dints_forward2;
        
        my $comp_seq = $seq;
           $comp_seq =~ tr/ATGC/TACG/;
           $comp_seq = reverse $comp_seq;
           
        my @dints_reverse1 = unpack("(A2)*", $comp_seq);
        my @dints_reverse2 = unpack("(A2)*", substr($comp_seq, 1));
        
        $dints_all{$_}->{reverse1} ++ for @dints_reverse1;
        $dints_all{$_}->{reverse2} ++ for @dints_reverse2;
    }
    print STDERR "\tdone!\n";

    
    
    print STDERR ">> Start generating results ... ";
    print STDOUT "##$CMDLINE\n##$SOURCE\n";
    print STDOUT "#dinucleotides\tforward1\tforward2\treverse1\treverse2\n";
    for my $dint (sort keys %dints_all)
    {
        next unless (length($dint) == 2);
        
        my @counts = ();
        for my $tag (qw(forward1 forward2 reverse1 reverse2))
        {
            my $dint_num = $dints_all{$dint}->{$tag} ? $dints_all{$dint}->{$tag} : 0;
            
            push @counts, $dint_num;
        }
        
        my $counts = join "\t", @counts;
        print STDOUT "$dint\t$counts\n";
    }
    print STDERR "done!\n";
}



=head2 count_triplets

    About   : Count triplets content in fasta file.
    Usage   : count_triplets($fasta_file);
    Args    : Sequences file in fasta format.
    Returns : Null

=cut
sub count_triplets
{
    my ($in) = @_;
    
    ###my $seq = 'AACAATACGA';
    ###my @codons = unpack("(A3)*", substr($seq, 2));
    ###
    ###my $comp_seq = $seq;
    ###   $comp_seq =~ tr/ATGC/TACG/;
    ###
    ###my @codons2 = unpack("(A3)*", $comp_seq);
    ###
    ###print STDERR Dumper(@codons, @codons2);exit;
    
    my %triples_all = ();
    for (my $i=0; $i<@IDs; $i++)
    {
        next if ($IDs[$i] =~ /(mitochondria|chloroplast|scaffold)/);
        
        print STDERR "\r>> Start counting ... $IDs[$i]";
        
        next if (length($SEQs[$i]) < 3);
        
        my $seq = uc ($SEQs[$i]);
        
        my @codons_forward1 = unpack("(A3)*", $seq);
        my @codons_forward2 = unpack("(A3)*", substr($seq, 1));
        my @codons_forward3 = unpack("(A3)*", substr($seq, 2));
        
        $triples_all{$_}->{forward1} ++ for @codons_forward1;
        $triples_all{$_}->{forward2} ++ for @codons_forward2;
        $triples_all{$_}->{forward3} ++ for @codons_forward3;
        
        my $comp_seq = $seq;
           $comp_seq =~ tr/ATGC/TACG/;
           $comp_seq = reverse $comp_seq;
           
        my @codons_reverse1 = unpack("(A3)*", $comp_seq);
        my @codons_reverse2 = unpack("(A3)*", substr($comp_seq, 1));
        my @codons_reverse3 = unpack("(A3)*", substr($comp_seq, 2));
        
        $triples_all{$_}->{reverse1} ++ for @codons_reverse1;
        $triples_all{$_}->{reverse2} ++ for @codons_reverse2;
        $triples_all{$_}->{reverse3} ++ for @codons_reverse3;
    }
    print STDERR "\tdone!\n";

    
    
    print STDERR ">> Start generating results ... ";
    print STDOUT "##$CMDLINE\n##$SOURCE\n";
    print STDOUT "#triplets\tforward1\tforward2\tforward3\treverse1\treverse2\treverse3\n";
    for my $codon (sort keys %triples_all)
    {
        next unless (length($codon) == 3);
        
        my @counts = ();
        for my $tag (qw(forward1 forward2 forward3 reverse1 reverse2 reverse3))
        {
            my $codon_num = $triples_all{$codon}->{$tag} ? $triples_all{$codon}->{$tag} : 0;
            
            push @counts, $codon_num;
        }
        
        my $counts = join "\t", @counts;
        print STDOUT "$codon\t$counts\n";
    }
    print STDERR "done!\n";
}

