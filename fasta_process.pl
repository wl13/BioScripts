#!/usr/bin/perl -w
#
#   fasta_process.pl -- Fasta-format file related operations.
#                          
#
#   Author: Nowind
#   Created: 2012-05-31
#   Updated: 2023-05-25
#   Version: 2.2.0
#
#   Change logs:
#   Version 1.0.0 14/11/10: The initial version.
#   Version 1.1.0 15/12/24: Add functions to count sequence context.
#   Version 1.1.1 16/03/10: Updated: add support for multiple fasta files and fasta records from pipeline.
#   Version 1.1.2 16/04/06: Bug fixed: direct exit if no entrance were extracted.
#   Version 2.0.0 16/05/08: Updated: add function to translate sequences; add function to split fasta file;
#                           support more output format; revise some descriptions.
#   Version 2.0.1 16/05/09: Bug fixed: remove unused option "--sort-by-list".
#   Version 2.0.2 17/03/21: Updated: add "--split" option in usage.
#   Version 2.1.0 19/11/20: Updated: add support for counting homopolymer runs; add support for removing
#                           unwanted sequences.
#   Version 2.2.0 23/05/25: Updated: add support for counting informative length; add option "--quite" to
#                           suppress progress infos.




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
use MyPerl::Convert qw(:all);

################################# Main ###############################


my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '2.2.0';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";
my $SOURCE  = (scalar localtime()) . " Version: $VERSION";

my $delimiter    = '\s+';
my $out_format   = 'fasta';
my $out_order    = 'input';
my $data_type    = 'DNA';
my $indel_symbol = '-';
my $no_found_order = 'top';
my (@fasta_files, $output, $word_wrap, $out_dir, $query_file,  @query_rows, @exclude_ids,
    $match_str, @sub_set, $substitution, $min_length, $max_length, $reverse_seq, $complement_seq,
    $uppercase_seq, $lowercase_seq, $count_nucl, $count_overall, $translate, $seperate, $numeric,
    $split_fasta, $new_id, $id_as_folder, $omit_id, $prefix, $suffix, $quite_mode);
GetOptions(
            "fasta=s{,}"       => \@fasta_files,
            "output=s"         => \$output,

            "query=s"          => \$query_file,
            "rows=i{,}"        => \@query_rows,
            "exclude=s{,}"     => \@exclude_ids,
            
            "out-format=s"     => \$out_format,
            "out-order=s"      => \$out_order,
            
            "numeric"          => \$numeric,
            "seperate"         => \$seperate,
            
            "data-type=s"      => \$data_type,
            "symbol=s"         => \$indel_symbol,
            
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
            "count-all"        => \$count_overall,
            
            "replace=s"        => \$substitution,
            
            "translate"        => \$translate,
            
            "split"            => \$split_fasta,
            
            "outdir=s"         => \$out_dir,
            "prefix=s"         => \$prefix,
            "suffix=s"         => \$suffix,
            
            "new-id=s"         => \$new_id,
            
            "id-as-folder"     => \$id_as_folder,
            "omit-id"          => \$omit_id,
            
            "quite"            => \$quite_mode,
           );

unless( @fasta_files > 0 ) {
    print <<EOF;

$0  -- Process fasta sequences

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
        "mega":     mega format
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

    --reverse
        reverse sequences before output        
    --complement
        complement sequences before output

    --translate
        translate nucleotides to proteins
    
    --numeric
        convert nucleotide bases to numbers according to the following
        conversions:
        A:1, T:2, G:3, C:4, -:5, N*:6
        *N stands for all unkown bases, used for tabular output

    --seperate
        seperate characters using comma, used for tabular output


    --outdir     <filename>
        output directory while splitting fasta file, default to current
        working directory
    
    
    --split
        split each fasta record to a seperate file
    --prefix     <string>
        prefix of output filename after splitting
    --suffix     <string>
        suffix of output filename after splitting
    --new-id     <string>
        use a new id to replace original id in all output files after
        splitting
    
    --id-as-folder
        use original sequence id as folder name and put each splitted file
        in related folder
    --omit-id
        omit ids in output filename, if this option is specified, at least
        one of the "--prefix" or "--suffix" options should be setted


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

    --replace <string>
        replace reference nucleotide with the specified nucleotide
        at specified position for each sequence in the format:
        
        "row of position,row of reference nucleotides,row of substitution
        nucleotides"

    --exclude <strings>
        exclude sequences whose ids match the specified strings
    
Filtering Options:
    
    --match     <string>
        only considering lines matching a pattern, support perl
        regular expression
    
    --lower     <int>
        only output sequences with length above this value
    --upper     <int>
        only output sequences with length below this value


 
Other options:

    --count-nucl <string>
        count "informative" length, or count "dinucleotide", "triplet",
        "homopolymer" nucleotide context
    --count-all
        sum up all counts across given sequences, note this option will first
        join all sequences into a single sequence by padding "N" between each
        sequence, so the counting results for di- or tri- nt context for each
        frame will be slightly different from directly sum up of all nts from
        each sequence though the sum of all frames will be identical
    
    --data-type  <string>
        specify data type of input fasta file if output in mega format,
        can be set to DNA or Protein, default: DNA
    --symbol     <string>
        symbol for indels used for mega output, default: "-"
        
    
    --quite
        no verbosity
    
   *Note: some options have orders, for example you can first extract then
    sort the extracted sequences, while the opposite will not work; simply
    break it into two steps by first sort, then extract. Some options could
    be combined, some could not.

EOF

    exit(1);
}

$|++;


print STDERR "# $0 v$VERSION\n# " . (scalar localtime()) . "\n" unless($quite_mode);

if ($split_fasta) {
    unless ( $out_dir ) { $out_dir = '.'; }
    unless( -e $out_dir ) { mkdir $out_dir; }
    
    $out_dir =~ s/\/$//g;
}
elsif ($output) {
    open (STDOUT, "> $output") || die $!;
}

unless(@query_rows){ @query_rows = (0) };


## read into all sequences
my @input_seqs = ();
my @input_ids  = ();
for my $in (@fasta_files)
{
    print STDERR ">> Start reading sequences from $in ... " unless($quite_mode);
    my @ids  = parse_fasta_SEQs(\@input_seqs, $in);
    push @input_ids, @ids;
    print STDERR "done!\n" unless($quite_mode);
}

## convert array to hash
my %input_seqs  = ();
for (my $i=0; $i<@input_ids; $i++)
{
    $input_seqs{$input_ids[$i]}  = $input_seqs[$i];
}


##
## copy input to output
##
my @output_ids  = ();
my @output_seqs = ();

if (@exclude_ids > 0) {
    my $exclude_str = join '|', @exclude_ids;
    
    for (my $i=0; $i<@input_ids; $i++)
    {
        next if ($input_ids[$i] =~ /$exclude_str/);
        
        push @output_ids,  $input_ids[$i];
        push @output_seqs, $input_seqs[$i];
    }
}
else {
    @output_ids  = @input_ids;
    @output_seqs = @input_seqs;
}




##
## query sequences
##
if (defined $query_file) {
    print STDERR ">> Start querying $query_file ... " unless($quite_mode);
    my %query_records = ();
    extract_seqs(\%query_records, $query_file);
    print STDERR "done!\n" unless($quite_mode);
    
    if($query_records{id}) {
        ## update the output arrays
        @output_ids  = @{$query_records{id}};
        @output_seqs = @{$query_records{seq}};
    }
    else {
        print STDERR "No sequence extracted!\n"; exit(0);
    }

}


##
## translate sequences
##
if ($translate) {
    print STDERR "\r>> Start translating sequences ... " unless($quite_mode);
    my @translated_seqs = ();
    for (my $i=0; $i<@output_ids; $i++)
    {
        push @translated_seqs, Translate($output_seqs[$i]);
    }
    
    @output_seqs = @translated_seqs;
    print STDERR "done!\n" unless($quite_mode);
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
        ## store the sorted ids and sequences in a new hash
        $output_seqs{$output_ids[$_]} = $output_seqs[$_] for (0..$#output_ids);
        @output_ids   = @{$ra_sorted};
    }
    else {
        ## only update output arrays
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
    elsif ($count_nucl eq 'homopolymer') {
        print STDOUT "#seq_id\thomopolymers\tlength\tcount\n";
    }
    elsif ($count_nucl eq 'informative') {
        print STDOUT "#seq_id\tinformative_length\n";
    }
}

if ($out_format eq 'mega') {
    print STDOUT "#mega\n";
    print STDOUT "!Title ;\n";
    print STDOUT "!Format DataType=$data_type indel=$indel_symbol;\n\n";
}

##
## sum up all sequences when counting
##
if ($count_nucl && $count_overall) {
    my $combined_seq = join "N", @output_seqs;
    @output_seqs = ($combined_seq);
    @output_ids  = qw(Overall);
}

print STDERR ">> Start writing results ... " unless($quite_mode);
for (my $i=0; $i < @output_ids; $i++)
{
    my $out_id = $output_ids[$i];
    
    if ($sort_by_list) {
        unless($output_seqs{$out_id}) {
            print STDERR "Error: $out_id not found!\n";
            exit(2);
        }
    }
    
    my $seq = $sort_by_list ? $output_seqs{$out_id} : $output_seqs[$i];
    
    next if ($min_length && (length $seq) < $min_length);
    next if ($max_length && (length $seq) > $max_length);
    
    if ($reverse_seq) {
        $seq = reverse $seq;
    }
    if ($complement_seq) {
        $seq =~ tr/ATGCatgc/TACGtacg/;
    }
    
    if ($split_fasta) {
        my $out_name  = $prefix ? $prefix . $out_id : $out_id;
           $out_name .= $suffix if ($suffix);
           
        if ($omit_id) {
            if ($prefix && $suffix) {
                $out_name = $prefix . $suffix;
            }
            elsif ($prefix) {
                $out_name = $prefix;
            }
            elsif ($suffix) {
                $out_name = $suffix;
            }
            else {
                print STDERR "Error: no valid output filename specified!\n";
                exit(2);
            }
        }
        
        my $out_id  = $new_id ? $new_id : $out_id;
        
        if ($id_as_folder) {
            $out_dir = "$out_dir\/$out_id";
            
            unless(-e $out_dir){ mkdir $out_dir };
        }
        
        open (my $ot, "> $out_dir\/$out_name") || die $!;
        print $ot format_fasta_SEQs($out_id, \$seq, $word_wrap);
        
    }
    elsif ($count_nucl) {
        if ($count_nucl eq 'triplet') {
            count_triplets($out_id, \$seq);
        }
        elsif ($count_nucl eq 'dinucleotide') {
            count_dinucleotide($out_id, \$seq);
        }
        elsif ($count_nucl eq 'homopolymer') {
            count_polymers($out_id, \$seq);
        }
        elsif ($count_nucl eq 'informative') {
            my $info_len = ($seq =~ tr/ATGCatgc/ATGCatgc/);
            print STDOUT "$out_id\t$info_len\n";
        }
    }
    elsif ($out_format eq 'fasta') {
        print STDOUT format_fasta_SEQs($out_id, \$seq, $word_wrap);
    }
    elsif ($out_format eq 'tabular') { ## write in tabular format
        my $out_str = $seq;
        
        if ($numeric) {                ## use numbers instead of nucleotides
            $out_str =~ tr/ATGC\-N?*ZM/123456/;
        }
        
        $out_str = $seperate ? (join ',', (split //, $out_str)) : $out_str;
        
        print STDOUT "$out_id\t$out_str\n";
    }
    elsif ($out_format eq 'mega') {
        my $formated_seq = format_fasta_SEQs($out_id, \$seq, $word_wrap);
           $formated_seq =~ s/\>/#/g;
        
        print STDOUT "$formated_seq\n";
    }
    elsif ($out_format eq 'seq') {
        print STDOUT "$seq\n";
    }
}
print STDERR "done!\n" unless($quite_mode);

print STDERR "# " . (scalar localtime()) . "\n" unless($quite_mode);




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
        
        print STDERR ">> Start parsing $out_order ... " unless($quite_mode);
        my $fh = getInputFilehandle($out_order);
        while (<$fh>)
        {
            next if (/^\#/ || /^\s+$/);
            
            my ($id) = (split /\s+/);
            
            push @ref_ids, $id;
        }
        print STDERR "done!\n" unless($quite_mode);
        
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


=head2 count_polymers

    About   : Count homopolymer runs in fasta file.
    Usage   : count_polymers($seq_id, $rs_seq);
    Args    : Sequence id;
              Scalar reference to sequence.
    Returns : Null

=cut
sub count_polymers
{
    my ($id, $rs_seq) = @_;
    
    return if (length($$rs_seq) < 3);
    
    my $seq = uc ($$rs_seq);
    

    my %polymer_runs = ();
    
    for my $nt (qw(A T G C))
    {
        my $rep_seq = $seq;
           $rep_seq =~ s/[^$nt]/N/g;
        
        my @nts = (split /N+/, $rep_seq);
        
        $polymer_runs{$_} ++ for @nts;
    }
    
    for my $polymer (sort keys %polymer_runs)
    {
        my $rel_len = ($polymer =~ tr/ATGC/ATGC/);
        
        next if ($rel_len < 3);
        
        my $count = join "\t", $polymer_runs{$polymer};
        
        print STDOUT "$id\t$polymer\t$rel_len\t$count\n";
    }
}




