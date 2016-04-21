#!/usr/bin/perl -w
#
#   fasta_process.pl -- query fasta sequences
#                          
#
#   Author: Nowind
#   Created: 2012-05-31
#   Updated: 2016-04-06
#   Version: 1.1.2
#
#   Change logs:
#   Version 1.0.0 14/11/10: The initial version.
#   Version 1.1.0 15/12/24: Add functions to count sequence context.
#   Version 1.1.1 16/03/10: Update: add support for multiple fasta files and fasta records from pipeline.
#   Version 1.1.2 16/04/06: Bug fixed: direct exit if no entrance were extracted.



=head1 NAME

fasta_process.pl


=head1 SYNOPSIS

  fasta_process.pl --help/?

=head1 DESCRIPTION

Fasta format file related processes.

=cut


use strict;

use Data::Dumper;
use Getopt::Long;
use File::Find::Rule;
use File::Basename;

use MyPerl::FileIO qw(:all);

################################# Main ###############################


my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.1.2';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";
my $SOURCE  = (scalar localtime()) . " Version: $VERSION";

my $delimiter  = '\s+';
my $out_format = 'fasta';
my $out_order  = 'input';
my $no_found_order = 'top';
my (@fasta_files, $output, $word_wrap,
    $query_file,  @query_rows, $match_str, @sub_set, $substitution,
    $min_length, $max_length, $reverse_seq, $complement_seq,
    $uppercase_seq, $lowercase_seq, $count_nucl);
GetOptions(
            "fasta=s{,}"       => \@fasta_files,
            "output=s"         => \$output,
            
            "query=s"          => \$query_file,
            "rows=i{,}"        => \@query_rows,
            
            "out-format=s"     => \$out_format,
            "out-order=s"      => \$out_order,
            
            "no-found-order"   => \$no_found_order,
            
            "wordwrap=i"       => \$word_wrap,         ## Line feed for print
            "delimiter=s"      => \$delimiter,
            "match=s"          => \$match_str,
            "subset=i{2}"      => \@sub_set,
            "lower=i"          => \$min_length,
            "upper=i"          => \$max_length,
            
            "reverse"          => \$reverse_seq,
            "complement"       => \$complement_seq,
            "uppercase"        => \$uppercase_seq,
            "lowercase"        => \$lowercase_seq,
            
            "count-nucl=s"     => \$count_nucl,
            
            "replace=s"        => \$substitution,
           );

unless( @fasta_files > 0 ) {
    print <<EOF;

$0  -- query fasta sequences

Version: $VERSION

Usage:   perl $0 [options]

Input/Output Options:

    --fasta  <filename>
        sequences file(s) in fasta format, required
    --output <filename>
        output file, default to STDOUT
        
    --out-format <string>
        specify output format, currently support:
        
        "fasta":    standard fasta format;
        "tabular":  tabulated results with id and sequence in each line,
        "seq":      output sequences only
        
        [default: fasta]
        
    --out-order <string>
        specify orders of output results,
        
        "input":    same order as input file;
        "alpha":    sort output sequences by ids using alphabetic sorting;
        "num":      sort output sequences by ids using numeric sorting;
        
        or set a file name, then the output ids will be sorted accroding
        to the order in the refer list
        
        [default: input]

    --no-found-order <string>
        while sort sequences by refer list, put those sequences not found in
        refer list to the "top" or "bottom" of all others, [default: top] 
        
    --wordwrap  <int>
        line feed for print, only valid while output in fasta format


Manipulation Options:

    --sort-by-list  <filename>
        sort sequences by the order of ids in another file, the file contains
        sequence ids same as the fasta file, one id per line, otherwise will
        use numeric sort by default
        
    --reverse
        reverse sequences before output        
    --complement
        complement sequences before output

        
Extracting Options:
 
    --query  <filename>
        query file with sequence ids need to be extracted
        
    --rows   <int>
        specify the input fields, can have multiple values (0-based),
        the first will be used as the row of the query IDs, others
        will be treated as append infomations [default: 0]
    --subset    <numbers>
        extract only part of the query sequences, set 2 values to
        indicate start and end row (0-based)

    --delimiter <numbers>
        specify a delimiter while reading input files, such as ",",
        "\\t", multiple delimiters can be set such as ",|\\t"
        [default: "\\s+"]
    --match     <string>
        only considering lines matching a pattern, support perl
        regular expression

    --lower     <int>
        only output sequences with length above this value
    --upper     <int>
        only output sequences with length below this value
        
    --replace <string>
        replace reference nucleotide with the specified nucleotide
        at specified position for each sequence in the format:
        
        "row of position,row of reference nucleotides,row of substitution
        nucleotides"

Other options:

    --count-nucl
        count "dinucleotide" or "triplet" nucleotide context
        
EOF

    exit(1);
}

$|++;


print STDERR "# $0 v$VERSION\n# " . (scalar localtime()) . "\n";

if ($output) {
    open (STDOUT, "> $output") || die $!;
}

unless(@query_rows){ @query_rows = (0) };


## read into all sequences
my @input_seqs = ();
my @input_ids  = ();
for my $in (@fasta_files)
{
    print STDERR ">> Start reading sequences from $in ... ";
    my @ids  = parse_fasta_SEQs(\@input_seqs, $in);
    push @input_ids, @ids;
    print STDERR "done!\n";
}

## convert array to hash
my %input_seqs = ();
for (my $i=0; $i<@input_ids; $i++)
{
    $input_seqs{$input_ids[$i]} = $input_seqs[$i];
}


my @output_ids  = @input_ids;
my @output_seqs = @input_seqs;

if (defined $query_file) {
    print STDERR ">> Start querying $query_file ... ";
    my %query_records = ();
    extract_seqs(\%query_records, $query_file);
    print STDERR "done!\n";
    
    if($query_records{id}) {
        @output_ids  = @{$query_records{id}};
        @output_seqs = @{$query_records{seq}};
    }
    else {
        print STDERR "No sequence extracted!\n"; exit(0);
    }

}

##
## sort output results
##
my $sort_by_list = 0;
my %output_seqs  = ();
if ($out_order ne 'input') {
    my ($ra_sorted, $sort_flag) = sort_fasta(\@output_ids);
    
    $sort_by_list = $sort_flag;
    
    if ($sort_flag) {
        $output_seqs{$output_ids[$_]} = $output_seqs[$_] for (0..$#output_ids);
        @output_ids   = @{$ra_sorted};
    }
    else {
        @output_ids   = @output_ids[@{$ra_sorted}];
        @output_seqs  = @output_seqs[@{$ra_sorted}];
    }
}


##
## generating results
##
if ($count_nucl) {
    print STDOUT "##$CMDLINE\n##$SOURCE\n";
    if ($count_nucl eq 'triplet') {
        print STDOUT "#seq_id\ttriplets\tforward1\tforward2\tforward3\treverse1\treverse2\treverse3\n";
    }
    elsif ($count_nucl eq 'dinucleotide') {
        print STDOUT "#seq_id\tdinucleotides\tforward1\tforward2\treverse1\treverse2\n";
    }
}


print STDERR ">> Start writing results ... ";
for (my $i=0; $i < @output_ids; $i++)
{
    my $out_id = $output_ids[$i];
    
    if ($sort_by_list) {
        unless($input_seqs{$out_id}) {
            print STDERR "Error: $out_id not found!\n";
            exit(2);
        }
    }
    
    my $seq    = $sort_by_list ? $input_seqs{$out_id} : $output_seqs[$i];
    
    next if ($min_length && (length $seq) < $min_length);
    next if ($max_length && (length $seq) > $max_length);
    
    if ($reverse_seq) {
        $seq = reverse $seq;
    }
    if ($complement_seq) {
        $seq =~ tr/ATGCatgc/TACGtacg/;
    }
    
    if ($count_nucl) {
        if ($count_nucl eq 'triplet') {
            count_triplets($out_id, \$seq);
        }
        elsif ($count_nucl eq 'dinucleotide') {
            count_dinucleotide($out_id, \$seq);
        }
    }
    elsif ($out_format eq 'fasta') {
        print STDOUT format_fasta_SEQs($out_id, \$seq, $word_wrap);
    }
    elsif ($out_format eq 'tabular') {
        print STDOUT "$out_id\t$seq\n";
    }
    elsif ($out_format eq 'seq') {
        print STDOUT "$seq\n";
    }
}
print STDERR "done!\n";

print STDERR "# " . (scalar localtime()) . "\n";




######################### Sub #########################


=head2 extract_seqs

    About   : Query sequence data.
    Usage   : extract_seqs($rh_query_records, $fasta_file);
    Args    : Hash reference to query results;
              Source fasta file.
    Returns : Null

=cut
sub extract_seqs
{
    my ($rh_query_records, $in) = @_;
    
    
    my $fh = getInputFilehandle($in);
    while (<$fh>)
    {
        next if (/\#/ || /^\s+$/);
    
        next if ($match_str && !/$match_str/);
        
        my @rows = (split /$delimiter/, $_)[@query_rows];
        
        my ($query_id, @desc)  = @rows;
        
        if ($input_seqs{$query_id}) {
            my $seq  =  $input_seqs{$query_id};
               $seq  =~ s/\*$//;
    
            my $del_len = 0;
            if ($substitution) {
                my ($sub_pos_row, $ref_row, $sub_nt_row) = (split /\,/, $substitution);
                
                my $sub_pos = $rows[$sub_pos_row];
                my $ref_nt  = $rows[$ref_row];
                my $sub_nt  = $rows[$sub_nt_row];
                
                $del_len = length($ref_nt) - length($sub_nt);
                
                $del_len = 0 if ($del_len < 0);
            }
            
               
            if (@sub_set > 0) {
                my ($start, $end) = (split /$delimiter/, $_)[@sub_set];
                
                if ($start > $end) {             ## reverse complement
                    $start = $start + $del_len;  ## adjust end positions to ensure flanking length while deletion found
                    
                    $seq = substr($seq, $end-1, $start-$end+1);
                    $seq =~ tr/ATGCatgc/TACGtacg/;
                    $seq = reverse $seq;
                    
                    push @desc, "reverse_complement";
                }
                else {
                    $end = $end + $del_len;      ## adjust end positions to ensure flanking length while deletion found
                    
                    $seq = substr($seq, $start-1, $end-$start+1);
                }
            }
            
            if ($substitution) {
                my ($sub_pos_row, $ref_row, $sub_nt_row) = (split /\,/, $substitution);
                
                my $sub_pos = $rows[$sub_pos_row];
                my $ref_nt  = $rows[$ref_row];
                my $sub_nt  = $rows[$sub_nt_row];
                
                substr($seq, $sub_pos-1, length($ref_nt), $sub_nt);
            }
            
            my $info = join "\_", @desc;
    
            my $out_id  = $query_id;
               $out_id .= "\_$info" if ($info);
    
            push @{$rh_query_records->{id}},  $out_id;
            push @{$rh_query_records->{seq}}, $seq;
        }
        else {
            print STDERR "no record found: $query_id\n";
        }
    }
}


=head2 sort_fasta

    About   : Sort sequences.
    Usage   : my ($ra_sorted_ids, $ra_sorted_index) = sort_fasta(\@input_ids);
    Args    : Array reference to input sequence ids.
    Returns : Array reference to sorted ids or index;
              Flag of whether sorting by a reference list or not.

=cut
sub sort_fasta
{
    my ($ra_input_ids) = @_;
    
    if ($out_order eq 'alpha') {
        my @sorted_index = sort {$ra_input_ids->[$a] cmp $ra_input_ids->[$b]} (0..$#{$ra_input_ids});
        
        return ([@sorted_index], 0);
    }
    elsif ($out_order eq 'num') {
        my @sorted_index = sort {$ra_input_ids->[$a] <=> $ra_input_ids->[$b]} (0..$#{$ra_input_ids});
        
        return ([@sorted_index], 0);
    }
    elsif (-f $out_order) {
        my @ref_ids = ();
        
        print STDERR ">> Start parsing $out_order ... ";
        my $fh = getInputFilehandle($out_order);
        while (<$fh>)
        {
            next if (/^\#/ || /^\s+$/);
            
            my ($id) = (split /\s+/);
            
            push @ref_ids, $id;
        }
        print STDERR "done!\n";
        
        my %counts   = ();
        my @no_order = ();
        if (@{$ra_input_ids} > @ref_ids) {
            $counts{$_}++ for @{$ra_input_ids};
            $counts{$_}++ for @ref_ids;
            
            for my $id (sort keys %counts)
            {
                if ($counts{$id} == 1) {
                    push @no_order, $id;
                }
            }
        }
        
        if ($no_found_order eq 'top') {
            return ([(@no_order, @ref_ids)], 1);
        }
        else {
            return ([(@ref_ids, @no_order)], 1);
        }
    }
}



=head2 count_dinucleotide

    About   : Count dinucleotide content in fasta file.
    Usage   : count_dinucleotide($seq_id, $rs_seq);
    Args    : Sequence id;
              Scalar reference to sequence.
    Returns : Null

=cut
sub count_dinucleotide
{
    my ($id, $rs_seq) = @_;
    
    return if (length($$rs_seq) < 2);
    
    my $seq = uc ($$rs_seq);
    
    my @dints_forward1 = unpack("(A2)*", $seq);
    my @dints_forward2 = unpack("(A2)*", substr($seq, 1));
    
    my %dints = ();
    
    $dints{$_}->{forward1} ++ for @dints_forward1;
    $dints{$_}->{forward2} ++ for @dints_forward2;
    
    my $comp_seq = $seq;
       $comp_seq =~ tr/ATGC/TACG/;
       $comp_seq = reverse $comp_seq;
       
    my @dints_reverse1 = unpack("(A2)*", $comp_seq);
    my @dints_reverse2 = unpack("(A2)*", substr($comp_seq, 1));
    
    $dints{$_}->{reverse1} ++ for @dints_reverse1;
    $dints{$_}->{reverse2} ++ for @dints_reverse2;
    
    for my $dint (sort keys %dints)
    {
        my $rel_len = ($dint =~ tr/ATGC/ATGC/);
        next unless ($rel_len == 2);
        
        my @counts = ();
        for my $tag (qw(forward1 forward2 reverse1 reverse2))
        {
            my $dint_num = $dints{$dint}->{$tag} ? $dints{$dint}->{$tag} : 0;
            
            push @counts, $dint_num;
        }
        
        my $counts = join "\t", @counts;
        print STDOUT "$id\t$dint\t$counts\n";
    }
}


=head2 count_triplets

    About   : Count triplets content in fasta file.
    Usage   : count_triplets($seq_id, $rs_seq);
    Args    : Sequence id;
              Scalar reference to sequence.
    Returns : Null

=cut
sub count_triplets
{
    my ($id, $rs_seq) = @_;
    
    return if (length($$rs_seq) < 3);
    
    my $seq = uc ($$rs_seq);
    
    my @codons_forward1 = unpack("(A3)*", $seq);
    my @codons_forward2 = unpack("(A3)*", substr($seq, 1));
    my @codons_forward3 = unpack("(A3)*", substr($seq, 2));
    
    my %triples = ();
    
    $triples{$_}->{forward1} ++ for @codons_forward1;
    $triples{$_}->{forward2} ++ for @codons_forward2;
    $triples{$_}->{forward3} ++ for @codons_forward3;
    
    my $comp_seq = $seq;
       $comp_seq =~ tr/ATGC/TACG/;
       $comp_seq = reverse $comp_seq;
       
    my @codons_reverse1 = unpack("(A3)*", $comp_seq);
    my @codons_reverse2 = unpack("(A3)*", substr($comp_seq, 1));
    my @codons_reverse3 = unpack("(A3)*", substr($comp_seq, 2));
    
    $triples{$_}->{reverse1} ++ for @codons_reverse1;
    $triples{$_}->{reverse2} ++ for @codons_reverse2;
    $triples{$_}->{reverse3} ++ for @codons_reverse3;


    for my $codon (sort keys %triples)
    {
        my $rel_len = ($codon =~ tr/ATGC/ATGC/);
        next unless ($rel_len == 3);
        
        my @counts = ();
        for my $tag (qw(forward1 forward2 forward3 reverse1 reverse2 reverse3))
        {
            my $codon_num = $triples{$codon}->{$tag} ? $triples{$codon}->{$tag} : 0;
            
            push @counts, $codon_num;
        }
        
        my $counts = join "\t", @counts;
        
        print STDOUT "$id\t$codon\t$counts\n";
    }
}
