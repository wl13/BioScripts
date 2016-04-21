#!/usr/bin/perl -w
#
#  FileIO.pm -- Process of files
#
#  Author: Nowind
#  Created: 2010-10-09
#  Updated: 2015-11-13
#  Version: 2.1.0
#
#  Change logs:
#  Version 1.0.0 11/11/12: Add function ParseSNPfile, ParseVCFSNPs, ParseShoreSNPs, ParseOtherSNPs.
#  Version 1.1.0 11/11/13: Add function ParseMergedSNPs.
#  Version 1.2.0 11/11/13: Remove function ParseOtherSNPs.
#  Version 1.3.0 11/11/13: Remove function ParseSNPfile, ParseVCFSNPs, ParseShoreSNPs, ParseMergedSNPs.
#  Version 1.3.1 11/12/04: Change function Parse_GFF to ParseGFF3. 
#  Version 1.3.2 11/12/14: Change the way in reading files in function Parse_Fasta.
#  Version 1.4.0 12/06/15: Add function getInputFilehandle.
#  Version 1.4.1 12/06/19: Add gzipped and bzipped file support in reading fasta files.
#  Version 1.4.2 12/07/03: Modify match patterns in function ParseGFF3; change some
#                          data structures to reduce memory use.
#  Version 1.4.3 12/07/04: Modify data structures in function ParseGFF3.
#  Version 1.4.4 12/07/28: Add @ids as return values in function Parse_Fasta.
#  Version 1.5.0 12/07/30: Remove dependency of "Bio::SeqIO"; remove function Get_Seq_Hash.
#  Version 1.5.1 12/08/09: Add existence checking in function ParseGFF3.
#  Version 1.5.2 12/08/14: Bug fixed in ParseGFF3 in match patterns.
#  Version 1.5.3 12/10/21: Modify match patterns in ParseGFF3.
#  Version 1.5.4 12/11/04: Add checking for gene type in ParseGFF3.
#  Version 1.5.5 12/11/12: Add support for reading from stdin in getInputFilehandle.
#  Version 1.5.6 12/12/20: Bug fixed in reading stdin.
#  Version 2.0.0 12/12/28: Change to a safer way for exporting functions; change serveral functions'
#                          names: Format_Fasta -> format_fasta_SEQs, Parse_Fasta -> parse_fasta_SEQs,
#                          ParseGFF3 -> parse_gff3_file; remove unused functions Get_Seq_Hash,
#                          Get_cds_info and Parse_Variant; complete the documents.
#  Version 2.0.1 13/01/21: Add annotations of intergenic regions in parse_gff3_file(uncompleted*).
#  Version 2.0.2 15/01/06: Bug fixed: "Use of uninitialized value $id" while id descriptions contian
#                          character ">".
#  Version 2.1.0 15/11/13: Updated: add function get_genome_length.



=head1 NAME

MyPerl::FileIO - Local perl module for file-related operations


=head1 SYNOPSIS

  use MyPerl::FileIO qw(:all);
  
  my $fh = getInputFilehandle( $filename );

  my $genome_size = get_genome_length(\%chrom_ids, \%chrom_lengths, $length_file, \@exclude_chroms);
  
  my @SEQ_IDs = parse_fasta_SEQs(\%SEQs, $filename);
  my @SEQ_IDs = parse_fasta_SEQs(\@SEQs, $filename, 'full');
  
  my $formated_str = format_fasta_SEQs($id, \$seq, $word_wrap);
  
  parse_gff3_file(\%GFF_info, $filename);

=head1 DESCRIPTION

Local perl module used for reading, parsing and format files used in bioinformatics

=cut

package MyPerl::FileIO;

use strict;

require Exporter;

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
            getInputFilehandle
            get_genome_length
            format_fasta_SEQs
            parse_fasta_SEQs
            parse_gff3_file
        )
    ]
);

@EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );
@EXPORT    = qw();


$MyPerl::FileIO::VERSION = '2.1.0';


=head1 METHODS

=head2 getInputFilehandle

    About   : Open and return filehandles
    Usage   : my $fh = getInputFilehandle($filename);
    Args    : Filename
    Returns : Filehandle to the opened file

=cut
sub getInputFilehandle
{
    my ($in) = shift;
    
    my $expr = "-";  
    
    if (!defined $in || $in eq "-") { ## read from STDIN
        $expr = "-";   
    }
    elsif ($in =~ /\.tar\.gz$/) {     ## read from a tar gzip ball
        $expr = "tar -zxf $in -O |";
    }
    elsif ($in =~ /\.tar\.bz2$/) {    ## read from a tar bzip2 ball
        $expr = "tar -jxf $in -O |";
    }
    elsif ($in =~ /\.gz$/) {          ## read from a gzip file
        $expr = "gzip -dc $in |";
    }
    elsif ($in =~ /\.bz2$/) {         ## read from a bzip2 file
        $expr = "bzip2 -dc $in |";
    }
    elsif ($in =~ /\.zip$/) {         ## read from a zip file
        $expr = "unzip -p $in |";
    }
    else {
        $expr = "< $in";
    }
    
    open (my $fh, $expr) || die $!;
    
    return $fh;
}


=head2 parse_fasta_SEQs

    About   : Reading sequences from a file in fasta format
    Usage   : my @SEQ_IDs = parse_fasta_SEQs(\%SEQs, $filename, $description);
    Args    : Reference to a hash or array to hold all the sequences;
              Fasta filename;
              Descriptions used for parsing sequence IDs.
    Returns : Array of sequence IDs.

=cut
sub parse_fasta_SEQs
{
    my ($rf_seq, $in, $desc) = @_;
    
    my $match_str = '^(.*?)\s+';
    
    if ($desc) {
        if ($desc eq 'full') {
            $match_str = '^(.*?)\n' 
        }
        elsif ($desc eq 'gbk') {
            $match_str = 'gb\|(.*?)\.';
        }        
    }
    
    my $fh  = getInputFilehandle($in);
    my $fas = do { local $/; <$fh> };
    
    my @fas = split /\n\>/, $fas;
    
    $fas[0] =~ s/^\>//;
    #shift @fas;
    
    my @ids = ();
    
    for my $str (@fas)
    {
        $str =~ /$match_str/;
        my $id = $1;
        
        $str =~ s/.*\n//;
        #$str =~ s/\>//;
        $str =~ s/\s+//g;
        
        if (ref($rf_seq) eq 'HASH') {     ## read into a hash
            $rf_seq->{$id} = $str;
        }
        elsif (ref($rf_seq) eq 'ARRAY') { ## read into an array
            push @{$rf_seq}, $str;
        }
        
        push @ids, $id;
    }
    
    return (@ids);
}


=head2 format_fasta_SEQs

    About   : Format sequence for output in fasta format
    Usage   : my $formated_str = format_fasta_SEQs($id, \$seq, $word_wrap);
    Args    : ID of the sequence;
              Reference to sequence;
              Maximum columns before line wraps.
    Returns : Formated string.

=cut
sub format_fasta_SEQs
{
    my ($id, $rs_seq, $wrap) = @_;
    
    my $seq = $$rs_seq;
    my $str = ">$id\n";
    
    if ($wrap) {
        while (length $seq > $wrap)
        {
            $str .= substr($seq, 0, $wrap, '');
            $str .= "\n";
        }        
    }
    $str .= "$seq\n";
    
    return $str;
}



=head2 parse_gff3_file

    About   : Parsing annotation informations in gff3 format file
    Usage   : parse_gff3_file(\%GFF_info, $filename);
    Args    : Reference to a hash to hold all the annotation informations;
              Filename of the gff3 format file.
    Returns : Null

=cut
sub parse_gff3_file
{
    my ($rh_GFF, $in) = @_;
    
    my $fh = getInputFilehandle($in);

    while (<$fh>)
    {
        next if (/#/ || /^\s+$/);
        
        chomp;
        
        my ($chr, $desc, $type) = (split /\s+/)[0..2];
        
        if ( $type =~ /chromosome/ ) { ## chromosome 
            m{
                ^(.*?)\s+.*?        # Chromosome ID
                chromosome\s+(\d+)  # Start Position
                \s+(\d+)            # End Position
            }x;
            
            next unless $1;
            
            my ($chr, $start, $end) = ($1, $2, $3);
            
            $rh_GFF->{$chr}->{pos} = "$start\-$end";
        }
        elsif ( $type =~ /(.*?gene)/ ) { ## gene
            my $gene_type = $1;
            
            m{
                ^(.*?)\s+.*?        # Chromosome ID
                gene\s+(\d+)        # Start Position
                \s+(\d+)\s+\.       # End Position
                \s+(\-|\+)\s+\.     # Strand
                \s+ID\=(.*?)\;.*    # ID
                Note\=(.*?)(;|$)    # Note
            }x;
            
            next unless $1;
            
            my ($chr, $start, $end, $Strand, $ID, $Note) =
               ($1, $2, $3, $4, $5, $6);
            
            $rh_GFF->{$chr}->{gene}->{$ID}->{type}   = $gene_type;
            $rh_GFF->{$chr}->{gene}->{$ID}->{strand} = $Strand;
            $rh_GFF->{$chr}->{gene}->{$ID}->{pos}    = "$start\-$end";
            $rh_GFF->{$chr}->{gene}->{$ID}->{note}   = $Note;
        }
        elsif ( $type =~ /(mRNA|pseudogenic\_transcript)/ ) { ## mRNA
            m{
                ^(.*?)\s+.*?                   # Chromosome ID
                \s+(\d+)                       # Start Position
                \s+(\d+)\s+.*                  # End Position
                \s+(\-|\+)\s+\.                # Strand
                \s+ID\=(.*?)\;.*               # ID
                (Parent|Locus_id)\=(.*?)(;|$)  # Parent
            }x;
            
            next unless $1;
            
            my ($chr, $start, $end, $Strand, $ID, $Parent) =
               ($1, $2, $3, $4, $5, $7);
            
            ###print "$chr\t$start\t$end\t$Strand\t$ID\t$Parent\n";exit;
            
            $rh_GFF->{$chr}->{mRNA}->{$ID}->{strand} = $Strand;
            $rh_GFF->{$chr}->{mRNA}->{$ID}->{pos}    = "$start\-$end";
            
            push @{$rh_GFF->{$chr}->{gene}->{$Parent}->{mRNA}}, $ID;
        }
        elsif ( $type =~ /(five|three)\_prime\_UTR/ ) {      ## five prime UTR 
            my $terminal = $1;
            
            m{
                ^(.*?)\s+.*?        # Chromosome ID
                UTR\s+(\d+)         # Start Position
                \s+(\d+).*?         # End Position
                Parent\=(.*)(;|$)   # Parent
            }x;
            
            next unless $1;
            
            my ($chr, $start, $end, $Parent) = ($1, $2, $3, $4);
            push @{$rh_GFF->{$chr}->{mRNA}->{$Parent}->{$terminal}}, "$start\-$end";
        }
        elsif ( $type =~ /(exon|pseudogenic\_exon)/ ) {    ## exon
            m{
                ^(.*?)\s+.*?        # Chromosome ID
                exon\s+(\d+)        # Start Position
                \s+(\d+).*?         # End Position
                Parent\=(.*)(;|$)   # Parent
            }x;
            
            next unless $1;
            
            my ($chr, $start, $end, $Parent) = ($1, $2, $3, $4);
            
            push @{$rh_GFF->{$chr}->{mRNA}->{$Parent}->{exon}}, "$start\-$end";
        }
        elsif ( $type =~ /CDS/ ) {     ## CDS 
            m{
                ^(.*?)\s+.*?        # Chromosome ID
                CDS\s+(\d+)         # Start Position
                \s+(\d+).*?         # End Position
                Parent=(.*?)(,|$)   # Parent
            }x;
            
            next unless $1;
            
            my ($chr, $start, $end, $Parent) = ($1, $2, $3, $4);
            push @{$rh_GFF->{$chr}->{mRNA}->{$Parent}->{CDS}}, "$start\-$end";
        }
        
        
        ###last if $.>40;
    }
    
    #for my $chrom (sort keys %{$rh_GFF})
    #{
    #    
    #    my @genes = sort {(split /\-/, $rh_GFF->{$chrom}->{gene}->{$a}->{pos})[0] <=>
    #                      (split /\-/, $rh_GFF->{$chrom}->{gene}->{$b}->{pos})[0]} (keys %{$rh_GFF->{$chrom}->{gene}});
    #    
    #    # intergenic region before the first gene
    #    my ($gene_id, $gene_start, $gene_end) = (split /\-/, $genes[0]); 
    #    my $intgenic_id    = $gene_id;
    #    my $intgenic_start = 1;
    #    my $intgenic_end   = $gene_start - 1;
    #    
    #    print STDOUT "$chrom\t.\tintergenic\t$intgenic_start\t$intgenic_end\t."
    #               . "\t+\t.\tID=$intgenic_id;Name=$intgenic_id;Note=N/A\n";
    #    
    #    for (my $i=0; $i<$#genes; $i++)
    #    {
    #        my ($bef_id, $bef_start, $bef_end) = (split /\-/, $genes[$i]);
    #        my ($aft_id, $aft_start, $aft_end) = (split /\-/, $genes[$i+1]);
    #        
    #        if ($bef_end < $aft_start) {
    #            $intgenic_id    = "$bef_id-$aft_id";
    #            $intgenic_start = $bef_end + 1;
    #            $intgenic_end   = $aft_start - 1;
    #            
    #            print STDOUT "$chrom\t.\tintergenic\t$intgenic_start\t$intgenic_end\t."
    #                       . "\t+\t.\tID=$intgenic_id;Name=$intgenic_id;Note=N/A\n";            
    #        }
    #        elsif ($aft_end < $bef_end) {
    #            @genes = @genes[0..$i, ($i+2)..$#genes];
    #            $i--;
    #        }
    #    }
    #    
    #    # intergenic region after the last gene
    #    ($gene_id, $gene_start, $gene_end) = (split /\-/, $genes[-1]);
    #    $intgenic_id    = $gene_id;
    #    $intgenic_start = $gene_end + 1;
    #    $intgenic_end   = $chroms{length}->[$i];;
    #    
    #    print STDOUT "$chrom\t.\tintergenic\t$intgenic_start\t$intgenic_end\t."
    #               . "\t+\t.\tID=$intgenic_id;Name=$intgenic_id;Note=N/A\n";
    #}
}


=head2 get_genome_length

    About   : Get length of each chromosomes.
    Usage   : my $genome_size = get_genome_length(\@chrom_ids, \%chrom_lengths, $length_file, \@exclude_chroms);
    Args    : Array reference to save all chromosome ids;
              Hash reference to save all chromosome lengths;
              Input file contain chromosome ids and lengthes;
              Chromosome ids to be ignored (pattern matches).
    Returns : Total size of all supplied chromosomes.

=cut
sub get_genome_length
{
    my ($ra_chrom_ids, $rh_chrom_lens, $in, $ra_exclude_chroms) = @_;
    
    ## exclude unwanted chromosomes or scaffolds while simulating, all
    ## chromosomes with ids match strings specified here would be ignored 
    my $exclude_str = '';
    if ($ra_exclude_chroms && @{$ra_exclude_chroms} > 0) {
        $exclude_str = join '|', @{$ra_exclude_chroms};
    }
    
    my $genome_size = 0;
    
    my $fh = getInputFilehandle($in);
    while (<$fh>)
    {
        next if (/^\#/ || /^\s+$/);
        
        my ($CHROM, $LENGTH) = (split /\s+/);
        
        if ($ra_exclude_chroms && @{$ra_exclude_chroms} > 0) {
            next if ($CHROM =~ /($exclude_str)/);
        }
        
        push @{$ra_chrom_ids}, $CHROM;
        $rh_chrom_lens->{$CHROM} = $LENGTH;
        
        $genome_size += $LENGTH;
    }
    
    return $genome_size;
}




1;

=head1 VERSION

2.1.0

=head1 AUTHOR

Nowind, noitulove9891@gmail.com

=head1 COPYRIGHT

Copyright (c) Nowind's Area. All rights reserved. This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself. 


=cut


