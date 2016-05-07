#!/usr/bin/perl -w
#
#   detectRecombEvents.pl -- Detect candidate recombination events.
#                          
#
#   Author: Nowind
#   Created: 2012-05-31
#   Updated: 2015-08-12
#   Version: 2.3.0
#
#   Change logs:
#   Version 1.0.0 13/02/16: The initial version.
#   Version 1.0.1 13/02/19: Exclude events in block adjacent regions.
#   Version 1.0.2 13/02/20: Add option "--min-markers", "--max-markers", "--max-length" 
#                           and "--min-length" to filter events; add some statistics in
#                           final results.
#   Version 1.1.0 13/04/07: Add option "--block" and "--rows" to set background blocks
#                           file; add detection of crossover breakpoints.
#   Version 1.1.1 13/04/08: Add option "--apply-filter" to apply filters in the final
#                           results.
#   Version 1.1.2 13/04/09: Remove vectors in crossover events.
#   Version 1.2.0 13/04/12: Reconstruct all functions; fix the processing of complex gene 
#                           conversion events; remove used options.
#   Version 1.2.1 13/04/15: Bug fixed while filtering crossover events.
#   Version 1.2.2 13/04/16: Add option "--min-major-perc" to filter mixed blocks which have
#                           a high possiblity to generate false positive events.
#   Version 1.2.3 13/04/17: Remove tests for the last event; add option "--complex-as-normal",
#                           fix markers infos for complex events, choose the maximum count of 
#                           continuous markers number as the number of complex events while
#                           this option is specified.
#   Version 2.0.0 14/01/15: Major update in detect gene conversions:
#                           (1) Add option "--contig" to get chromosome info while absent;
#                           (2) Merge gene conversions clustered together as they are more
#                               likely belong to a single event;
#                           (3) Rearrange related codes and implement an mapping-based strategy
#                               in merge and parsing gene conversions;
#                           (4) Add several other options to toggle whether events should be
#                               merged or filtered;
#                           (5) Adjust processes of complex gene conversions.
#   Version 2.0.1 14/01/17: Bug fixed in distance compare.
#   Version 2.1.0 14/05/20: Bug fixed while crossover occurred at the end of the chromosome;
#                           add support for multi-samples.
#   Version 2.2.0 14/05/21: Add option "--source-tag" to control which field describe the source
#                           type; remove option "--contig" and related functions.
#   Version 2.2.1 14/05/22: Add stats results of gene conversion events.
#   Version 2.2.2 14/06/03: Skip missing markers "./." from markers file.
#   Version 2.2.3 14/06/10: Add adjacent block infos while output crossover breakpoints.
#   Version 2.2.4 14/06/11: Add some codes for debugging; revise pattern match of missing markers.
#   Version 2.2.5 14/06/20: Add option "--linked" to toggle on previous detection method.
#   Version 2.2.6 14/06/21: Add option "--min-border-dist" to filter gene conversions.
#   Version 2.2.7 14/08/12: Skip markers whose background blocks is missing.
#   Version 2.3.0 15/08/13: Add support for detailed block infos for crossover events.


use strict;

use Data::Dumper;
use Getopt::Long;
use File::Find::Rule;
use File::Basename;

use MyPerl::FileIO qw(:all);

######################## Main ########################


my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '2.3.0';
my $HEADER  = "##$CMDLINE\n##Version: $VERSION\n";


my ($markers_file, $output, $block_file, $min_markers, $max_markers, $complex_as_normal,
    $source_tag, $min_gc_length, $max_gc_length, $min_gc_distance,
    $min_bg_length, $min_co_length, $min_diff_perc, $min_border_distance,
    @query_rows, $link_convert_markers, $out_block_details);
GetOptions(
            "vcf=s"               => \$markers_file,
            "source-tag=s"        => \$source_tag,
            
            "blocks=s"            => \$block_file,
            "rows=i{,}"           => \@query_rows,
            
            "output=s"            => \$output,
            
            "M|min-markers=i"     => \$min_markers,
            "X|max-markers=i"     => \$max_markers,

            
            "G|min-bg-len=i"      => \$min_bg_length,
            "L|min-gc-len=i"      => \$min_gc_length,
            "U|max-gc-len=i"      => \$max_gc_length,
            "D|min-gc-dist=i"     => \$min_gc_distance,
            "B|min-border-dist=i" => \$min_border_distance,
            "P|min-conv-perc=f"   => \$min_diff_perc,
            
            
            "C|min-co-len=i"      => \$min_co_length,
            "S|block-details"     => \$out_block_details,
            
            "A|complex-as-normal" => \$complex_as_normal,
            "K|linked"            => \$link_convert_markers,
           );

unless( $block_file ) {
    print <<EOF;

$0  -- Detect candidate recombination events.

Version: $VERSION

Usage:   perl $0 [options]

File Options:
    -b, --blocks  <filename>
        file contains background blocks, required
    -r, --rows    <numbers>
        specify the row fields of sample id, chromosome, block type,
        start and end position(0-based) in blocks file [default: 0 1 2 3 4]

    -v, --vcf     <filename>
        input markers file in vcf format, required in detecting gene
        conversions
        
        *Note: the sample ids in this file should match those in the blocks
        file

    --source-tag <string>
        tag of source type for each sample used to form blocks, can be set
        as "GT", "SC" or others, otherwise will use all sample field as
        the source type
        
        
    -o, --output  <filename>
        output filename, default to STDOUT

        
Options for filtering gene conversion events:
    
    -M, --min-markers <int>
        minimum number of markers each gene conversion contains
    -X, --max-markers <int>
        maximum number of markers each gene conversion contains
    
    -G, --min-bg-len      <int>
        minimum background size used to infer gene conversions
    -B, --min-border-dist <int>
        minimum distance to background borders, the border of background
        is affected by many effects and thus is not very accurate, calling
        gene conversion at border region will sometimes ended in much
        false positives
        
    -L, --min-gc-len  <int>
        minimum gene conversion length
    -U, --max-gc-len  <int>
        maximum gene conversion length
    -D, --min-gc-dist <int>
        minimum distance between two gene converions, otherwise these two
        events will be merged (in same type) or discarded (in diff types)
    
    -P, --min-conv-perc <float>
        minimum percentage of markers support the gene conversion event
    -A, --complex-as-normal
        choose the maximum count of continuous markers number as the number of
        complex events and do filtering as the normal events
    
    -K, --linked
        link adjancent converted markers into a long event

Options for processing crossover events:
    
    -C, --min-co-len  <int>
        minimum length of each background block used to infer crossover events
    -S, --block-details
        output details of blocks surrounding crossover events, user must
        specify a row of detailed block infos at the last of "--rows" option
    
EOF

    exit(1);
}

$|++;

if ($output) {
    open (STDOUT, "> $output") || die $!;
}

 
print STDERR "# $0 v$VERSION\n# " . (scalar localtime()) . "\n";

print STDOUT "$HEADER##" . (scalar localtime()) . "\n";


unless(@query_rows > 0){ @query_rows = qw(0 1 2 3 4) };


##
## read background file
##
print STDERR ">> Start parsing $block_file ... ";
my %BLOCKS = ();
parse_background_blocks(\%BLOCKS, $block_file);
print STDERR "done!\n";



##
## detect recombinations
##
if ($markers_file) {
    print STDERR ">> Start detecting candidate gene conversions ... ";
    my ($rh_events_info, $rh_Markers) = detect_gene_convers($markers_file);
    print STDERR "done!\n";
    
    print STDERR ">> Start refining and filtering detected gene conversoins ... ";
    procs_gene_convers($rh_events_info, $rh_Markers);
    print STDERR "done!\n";
}
else {
    detect_crossover();
}




print STDERR "# " . (scalar localtime()) . "\n";

######################### Sub #########################


=head2 parse_background_blocks

    About   : Parse background blocks.
    Usage   : parse_background_blocks(\%blocks, $file);
    Args    : Array reference to hold background blocks infos;
              File contains blocks infos.
    Returns : Null
    
=cut
sub parse_background_blocks
{
    my ($rh_blocks, $in) = @_;
    
    my $fh = getInputFilehandle($in);    
    while (<$fh>)
    {
        next if (/\#/ || /^\s+$/);
        
        my ($sample_id, $chrom, $type, $start, $end, $desc) = (split /\s+/)[@query_rows];
        
        push @{$rh_blocks->{$sample_id}->{$chrom}->{pos}}, [$start, $end];
        push @{$rh_blocks->{$sample_id}->{$chrom}->{type}}, $type;
        
        if ($desc) {
            push @{$rh_blocks->{$sample_id}->{$chrom}->{desc}}, $desc;
        }
    }
}




=head2 detect_crossover

    About   : Detect candidate crossover events.
    Usage   : detect_crossover();
    Args    : Null
    Returns : Null
    
=cut
sub detect_crossover
{
    my %Events = ();
    my %Stats  = ();
    
    my %Categories = ();
    
    ##
    ## detect candidate crossover events
    ##

    print STDOUT "#Sample\tChrom\tChanges\tEvent_Start\tEvent_End\tEvent_Len\tPrev_Block_Info\tNext_Block_Info\n";
    for my $sample (sort keys %BLOCKS)
    {
        for my $chrom (sort keys %{$BLOCKS{$sample}})
        {
            print STDERR "\r>> Start detecting candidate crossover events in $sample:$chrom ... ";
            
            my @bg_types  = ();
            my @block_pos = ();
            my @infos     = ();
            
            ## filtering background blocks
            if ($min_co_length) {
                for (my $i=0; $i<@{$BLOCKS{$sample}->{$chrom}->{pos}}; $i++)
                {
                    my ($block_start, $block_end) = @{$BLOCKS{$sample}->{$chrom}->{pos}->[$i]};
                    
                    my $block_len = $block_end - $block_start + 1;
                    
                    next if ($block_len < $min_co_length);
                    
                    push @bg_types,  $BLOCKS{$sample}->{$chrom}->{type}->[$i];
                    push @block_pos, $BLOCKS{$sample}->{$chrom}->{pos}->[$i];
                    
                    if ($out_block_details) {
                        push @infos, $BLOCKS{$sample}->{$chrom}->{desc}->[$i];
                    }
                }
            }
            else {
                @bg_types  = @{$BLOCKS{$sample}->{$chrom}->{type}};
                @block_pos = @{$BLOCKS{$sample}->{$chrom}->{pos}};
                
                if ($out_block_details) {
                    @infos = @{$BLOCKS{$sample}->{$chrom}->{desc}};
                }
            }
            
            for (my $i=1; $i<@bg_types; $i++)
            {
                if ($bg_types[$i] ne $bg_types[$i-1]) {
                    my $event_type  = join '<->', (sort ($bg_types[$i-1], $bg_types[$i]));
                    my $event_start = $block_pos[$i-1]->[1];
                    my $event_end   = $block_pos[$i]->[0];
                    my $event_len   = $event_end - $event_start + 1;
                    
                    my $prev_block_info = $block_pos[$i-1]->[1] - $block_pos[$i-1]->[0] + 1;
                    my $next_block_info = $block_pos[$i]->[1] - $block_pos[$i]->[0] + 1;
                    
                    if ($out_block_details) {
                        $prev_block_info .= ";" . $infos[$i-1];
                        $next_block_info .= ";" . $infos[$i];
                    }
 
                    print STDOUT "$sample\t$chrom\t$event_type\t$event_start\t$event_end\t$event_len\t$prev_block_info\t$next_block_info\n";
                    
                    $Stats{$sample}->{$chrom}->{$event_type} ++;
                    $Stats{Total}->{$chrom}->{$event_type} ++;
                    $Stats{$sample}->{All}->{$event_type} ++;
                    $Stats{Total}->{All}->{$event_type} ++;
                    $Stats{$sample}->{$chrom}->{All}  ++;
                    $Stats{Total}->{$chrom}->{All}  ++;
                    $Stats{$sample}->{All}->{All} ++;
                    $Stats{Total}->{All}->{All} ++;
                    
                    $Categories{sample}->{$sample} ++;
                    $Categories{sample}->{Total} ++;
                    $Categories{chrom}->{$chrom} ++;
                    $Categories{chrom}->{All} ++;
                    $Categories{type}->{$event_type} ++;
                    $Categories{type}->{All} ++;
                }
            }
        }
    }

    print STDERR "done!\n";
    
    
    my @sample_ids = sort keys %{$Categories{sample}};
    my $sample_ids = join "\t", @sample_ids;
    
    print STDERR ">> Start generating results ... ";
    for my $event_type (sort keys %{$Categories{type}})
    {
        print STDOUT "^#$event_type:chrom\t$sample_ids\n";
        
        for my $chrom (sort keys %{$Categories{chrom}})
        {
            my @counts = ();
            for my $sample (@sample_ids)
            {
                my $count = $Stats{$sample}->{$chrom}->{$event_type} ?
                            $Stats{$sample}->{$chrom}->{$event_type} : 0;
                push @counts, $count;
            }
            
            my $counts = join "\t", @counts;
            print STDOUT "^$event_type:$chrom\t$counts\n";
        }
    }

    print STDERR "done!\n";
}



=head2 detect_gene_convers

    About   : Detect candidate crossover events.
    Usage   : ($rh_events_info) = detect_gene_convers($file);
    Args    : File contains all markers infos.
    Returns : Reference to hash contains several infos of detected gene conversoins.
    
=cut
sub detect_gene_convers
{
    my ($in) = @_;
    
    my %Markers     = ();
    my %Block_nums  = ();
    
    my %events_all  = ();
    my @Samples_ids = ();
    my $fh = getInputFilehandle($in);
    while (<$fh>)
    {
        if (/#CHROM/) {
            my @line   = (split /\s+/);
            @Samples_ids = @line[9..$#line];
            
            ##
            ## initial hash values
            ##
            for my $sample_id (sort keys %BLOCKS)
            {
                for my $CHROM (sort keys %{$BLOCKS{$sample_id}})
                {
                    $Markers{$sample_id}->{$CHROM}->{prev}  = 0;
                    $Markers{$sample_id}->{$CHROM}->{curr}  = 0;
                    $Markers{$sample_id}->{$CHROM}->{start} = 0;
                    $Markers{$sample_id}->{$CHROM}->{end}   = 0;
                    $Markers{$sample_id}->{$CHROM}->{index} = 0;
                    
                    $Markers{CHROM}->{$CHROM} = 1;
                }
            }
        }
        
        next if (/\#/ || /^\s+$/);
        
        my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, @Samples) = (split /\s+/);
        
        next unless( $Markers{CHROM}->{$CHROM} );
        
        my @tags = (split /\:/, $FORMAT);
        my %tags = ();
        for (my $i=0; $i<@tags; $i++) { $tags{$tags[$i]} = $i; }
        

        for (my $n=0; $n < @Samples; $n++)
        {
            my $sample_id = $Samples_ids[$n];
            
            my $sample_type = $source_tag ? ((split /:/, $Samples[$n])[$tags{$source_tag}]) : $Samples[$n];
            
            next if ($sample_type =~ /\.(\/|\|)\./);
            
            $Markers{$sample_id}->{$CHROM}->{index} ++;  ## set a index for each position in each chromosome
            $Markers{$sample_id}->{$CHROM}->{type}->{$POS} = $sample_type;
            
            ##
            ## map this position to its background
            ##
            my $bg_type = '-';
            
            my $i = 0;
            while ($i < @{$BLOCKS{$sample_id}->{$CHROM}->{pos}})
            {
                ## check whether this POS in this block
                if ($POS >= $BLOCKS{$sample_id}->{$CHROM}->{pos}->[$i]->[0] &&
                    $POS <= $BLOCKS{$sample_id}->{$CHROM}->{pos}->[$i]->[1]) {
                    last;  ## found the block where the markers in
                }
                
                $i++;
            }
            
            if ($i < @{$BLOCKS{$sample_id}->{$CHROM}->{pos}}) {
                $bg_type = $BLOCKS{$sample_id}->{$CHROM}->{type}->[$i];
                
                my $block_start = $BLOCKS{$sample_id}->{$CHROM}->{pos}->[$i]->[0];
                my $block_end   = $BLOCKS{$sample_id}->{$CHROM}->{pos}->[$i]->[1];
                
                $Block_nums{$sample_id}->{$CHROM}->{$POS} = "$block_start\t$block_end";
            }
            else {
                next;
            }
            
            
            next if (!$sample_type || ($sample_type eq "."));  ## skip while background block did not exist
            
            
            ##
            ## start detection
            ##
            if ($link_convert_markers) {         ## search for converted markers in a linkage
                if ($bg_type ne $sample_type) {  ## save those linked markers
                    if ($Markers{$sample_id}->{$CHROM}->{curr} - $Markers{$sample_id}->{$CHROM}->{prev} == 1) {  
                        $Markers{$sample_id}->{$CHROM}->{end}  = $POS;
                        
                        $Markers{$sample_id}->{$CHROM}->{prev} ++;
                        $Markers{$sample_id}->{$CHROM}->{curr} ++;
                        
                        $Markers{$sample_id}->{$CHROM}->{count} ++;
                        
                        push @{$Markers{$sample_id}->{$CHROM}->{change}->{"$bg_type->$sample_type"}}, $Markers{$sample_id}->{$CHROM}->{index};
                    }
                    else {                                           ## reset all values while linkage is breaked by different changes
                        $Markers{$sample_id}->{$CHROM}->{start} = $POS;
                        $Markers{$sample_id}->{$CHROM}->{end}   = $POS;
                        
                        $Markers{$sample_id}->{$CHROM}->{count} = 1;
                        
                        $Markers{$sample_id}->{$CHROM}->{prev} = $Markers{$sample_id}->{$CHROM}->{curr};
                        $Markers{$sample_id}->{$CHROM}->{curr} ++;
                        
                        $Markers{$sample_id}->{$CHROM}->{change}->{"$bg_type->$sample_type"} = [$Markers{$sample_id}->{$CHROM}->{index}];
                        
                        $Markers{$sample_id}->{$CHROM}->{background} = $bg_type;
                        
                        ###print STDERR Dumper($Markers{$sample_id}->{$CHROM});exit;
                    }
                }
                else {  ## save previous search and start a new search while this linkage is breaked by markers in background type
                    $Markers{$sample_id}->{$CHROM}->{curr} ++;
                    
                    ##
                    ## found blocks changes in the background blocks
                    ##
                    if ($Markers{$sample_id}->{$CHROM}->{change}) {
                        
                        ## changes in adjacent regions of background blocks have a high possibility
                        ## of false positive so should not be counted in the final results
                        if ($Markers{$sample_id}->{$CHROM}->{background} eq $bg_type) {
                            my $event_start = $Markers{$sample_id}->{$CHROM}->{start};
                            my $event_end   = $Markers{$sample_id}->{$CHROM}->{end};
                            my $event_len   = $event_end - $event_start + 1;
                            
                            ##
                            ## check background
                            ##
                            next unless ($Block_nums{$sample_id}->{$CHROM}->{$event_start} &&    ## make sure the event located in only one block
                                         $Block_nums{$sample_id}->{$CHROM}->{$event_end} &&
                                         $Block_nums{$sample_id}->{$CHROM}->{$event_start} eq $Block_nums{$sample_id}->{$CHROM}->{$event_end});
                            
                            my ($bg_start, $bg_end) = split /\t/, $Block_nums{$sample_id}->{$CHROM}->{$event_start};
                            
                            push @{$events_all{$sample_id}->{$CHROM}}, [$event_start, $event_end, "$bg_type\t$bg_start\t$bg_end", 'Full'];
                        }
                        ###else {
                        ###    my $event_start = $Markers{$sample_id}->{$CHROM}->{start};
                        ###    my $event_end   = $Markers{$sample_id}->{$CHROM}->{end};
                        ###    print STDERR "$CHROM\t$sample_id\t$event_start\t$event_end\t$Block_nums{$sample_id}->{$CHROM}->{$event_start}\t"
                        ###               . "\t$Markers{$sample_id}->{$CHROM}->{background}\t$bg_type\n";
                        ###}
                        
                        ##
                        ## empty all records
                        ##
                        undef $Markers{$sample_id}->{$CHROM}->{change};
                        undef $Markers{$sample_id}->{$CHROM}->{background};
                    }
                    
                    $Markers{$sample_id}->{$CHROM}->{prev} = $Markers{$sample_id}->{$CHROM}->{curr};
                }
            }
            elsif ($bg_type ne $sample_type) { ## directly detect converted markers
                my ($bg_start, $bg_end) = split /\t/, $Block_nums{$sample_id}->{$CHROM}->{$POS};
                
                push @{$events_all{$sample_id}->{$CHROM}}, [$POS, $POS, "$bg_type\t$bg_start\t$bg_end", 'Full'];
            }
        }
    }
    
    return (\%events_all, \%Markers);
}


=head2 procs_gene_convers

    About   : Merge gene conversion events which more seems to be a single event and do some filtering
    Usage   : procs_gene_convers(\%events_all);
    Args    : Reference to hash contains several infos of detected gene conversoins.
    Returns : Null
    
=cut
sub procs_gene_convers
{
    my ($rh_events_all, $rh_Markers) = @_;
    
    my %Events_All = ();
    
    print STDOUT "##Explanation of output file headers: \n";
    print STDOUT "##Event_Type(Single): Single gene conversion event with only one type of changes\n";
    print STDOUT "##Event_Type(Multi): Multiple gene conversion event with more than one type of changes\n";
    print STDOUT "##Marker_Consists: Total number of markers;Consists of markers;percentage of all supported markers\n";
    print STDOUT "##BG_Info: Info of background blocks: type:start-end:length\n";
    print STDOUT "#Sample\tChrom\tEvent_Type\tEvent_Start\tEvent_End\tEvent_Length\tChanges"
               . "\tMarker_Consists\tBG_Info\n";
    for my $sample_id (sort keys %{$rh_events_all})
    {
        for my $chrom (sort keys %{$rh_events_all->{$sample_id}})
        {
            my @events = sort { $a->[0] <=> $b->[0] } @{$rh_events_all->{$sample_id}->{$chrom}};
            
            ##
            ## merge two events if two events been too closely
            ##
            if ($min_gc_distance) {
                for (my $i=0; $i<$#events; $i++)
                {
                    if ($events[$i+1]->[0] - $events[$i]->[1] < $min_gc_distance) {
                        if ($events[$i+1]->[2] ne $events[$i]->[2]) {
                            $events[$i]->[3]   = 'Exchange';
                            $events[$i+1]->[3] = 'Exchange';
                        }
                        else {
                            $events[$i]->[1]   = $events[$i+1]->[1];
                            $events[$i+1]->[0] = $events[$i]->[0];
                            $events[$i+1]->[2] = $events[$i]->[2];
                            
                            $events[$i]->[3]   = $events[$i]->[3] eq 'Exchange' ? 'Exchange' : 'Part';
                            $events[$i+1]->[3] = $events[$i]->[3];
                        }
                    }
                }
            }
    
            
            my %duplicated = ();
            for (my $i=0; $i<@events; $i++)
            {
                my ($event_start, $event_end, $bg_info, $tag) = @{$events[$i]};
                
                next if ($duplicated{$event_start});
                next if ($tag eq 'Exchange');
                
                $duplicated{$event_start}++;
                
                my $event_len = $event_end - $event_start + 1;
                
                ## remove events with too short or too long
                next if ($min_gc_length && $event_len < $min_gc_length);
                next if ($max_gc_length && $event_len > $max_gc_length);
                
                ##
                ## filter events in small background
                ##
                my ($bg_type, $bg_start, $bg_end) = (split /\t/, $bg_info);
                my $bg_len = $bg_end - $bg_start + 1;
                next if ($min_bg_length && $bg_len < $min_bg_length);  
                
                ##
                ## filter by distance relative to borders
                ##
                my $border_dist = ($event_start-$bg_start < $bg_end-$event_end) ?
                                  ($event_start-$bg_start) : ($bg_end-$event_end);
                next if ($min_border_distance && $border_dist < $min_border_distance); 
                
                ##
                ## check marker consists
                ##
                my @markers = grep { $rh_Markers->{$sample_id}->{$chrom}->{type}->{$_} } ($event_start..$event_end);
                my %counts  = ();
                
                for my $pos (@markers)
                {
                    my $type = $rh_Markers->{$sample_id}->{$chrom}->{type}->{$pos};
                    
                    $counts{total} ++;
                    
                    if ($type ne $bg_type) {
                        $counts{convs}->{$type} ++;
                        $counts{convs_all} ++;
                    }
                }
    
                ## remove events with too less or too much markers
                next if ($min_markers && $counts{total} < $min_markers);
                next if ($max_markers && $counts{total} > $max_markers);
    
                ## remove events with insufficient proportion of markers in different
                ## type compared to background
                my $percent  = 100 * $counts{convs_all} / $counts{total};
                next if ($min_diff_perc && $percent < $min_diff_perc);
                
                my @convs_type = sort {$counts{convs}->{$b} <=>
                                       $counts{convs}->{$a}} keys %{$counts{convs}};
                
                my $event_type = (@convs_type > 1) ? "Multi" : "Single";
                
                my @consists = ();
                for my $type (sort keys %{$counts{convs}})
                {
                    push @consists, "$type:$counts{convs}->{$type}";
                }
                
                my $consists = join ",", @consists;
                
                my @changes = ();
                for my $c_type (@convs_type)
                {
                    my $change = "$bg_type->$c_type";
                    push @changes, $change;
                }
                
                my $changes = join ';', @changes;
                
                print STDOUT "$sample_id\t$chrom\t$event_type\t$event_start\t$event_end\t$event_len"
                           . "\t$changes\t$counts{total};$consists;$percent\t$bg_type:$bg_start-$bg_end:$bg_len\n";
                
                
                push @{$Events_All{$chrom}->{$event_start}}, "$sample_id\t$changes\t$event_end\t$counts{convs_all}";
            }
        }
    }
    
    
    my %Events_Stat = ();
    my %Categories  = ();
    for my $chrom (sort keys %Events_All)
    {
        for my $start (keys %{$Events_All{$chrom}})
        {
            my @events = @{$Events_All{$chrom}->{$start}};
            my $fq = scalar @events;
            
            $Events_Stat{shared}->{$fq} ++;
            
            for my $event (@events)
            {
                my ($sample, $change, $end, $marker_num) = split /\t/, $event;
                
                $change = "Multi" if $change =~ /\;/;
                
                my $len = $end - $start + 1;
                
                $Events_Stat{change}->{$fq}->{$sample}->{$chrom}->{$change} ++;
                $Events_Stat{change}->{$fq}->{Total}->{$chrom}->{$change}  ++;
                $Events_Stat{change}->{$fq}->{$sample}->{All}->{$change}  ++;
                $Events_Stat{change}->{$fq}->{$sample}->{$chrom}->{All}  ++;
                $Events_Stat{change}->{$fq}->{Total}->{All}->{$change}  ++;
                $Events_Stat{change}->{$fq}->{Total}->{$chrom}->{All}  ++;
                $Events_Stat{change}->{$fq}->{$sample}->{All}->{All}  ++;
                $Events_Stat{change}->{$fq}->{Total}->{All}->{All}  ++;
                
                
                $Events_Stat{change}->{0}->{$sample}->{$chrom}->{$change} ++;
                $Events_Stat{change}->{0}->{Total}->{$chrom}->{$change} ++;
                $Events_Stat{change}->{0}->{$sample}->{All}->{$change} ++;
                $Events_Stat{change}->{0}->{$sample}->{$chrom}->{All} ++;
                $Events_Stat{change}->{0}->{Total}->{All}->{$change} ++;
                $Events_Stat{change}->{0}->{Total}->{$chrom}->{All} ++;
                $Events_Stat{change}->{0}->{$sample}->{All}->{All} ++;
                $Events_Stat{change}->{0}->{Total}->{All}->{All} ++;
    
                $Events_Stat{marker}->{$fq}->{$marker_num} ++;
                
                my $len_bin = int ($len / 10) + 1;
                $Events_Stat{length}->{$fq}->{$len_bin} ++;
                
                $Categories{chrom}->{$chrom} ++;
                $Categories{sample}->{$sample} ++;
                $Categories{change}->{$change} ++;
                
                $Categories{chrom}->{All} ++;
                $Categories{sample}->{Total} ++;
                $Categories{change}->{All} ++;
            }
        }
    }
    
    my @sample_ids = sort keys %{$Categories{sample}};
    my $sample_ids = join "\t", @sample_ids;
    
    print STDOUT "^##Statistics of gene conversion events\n";
    print STDOUT "^#shared:frequency\tcount\n";
    for my $fq (sort {$a <=> $b} keys %{$Events_Stat{shared}})
    {
        my $cnt = $Events_Stat{shared}->{$fq};
        print STDOUT "^shared:$fq\t$cnt\n";
    }
    
    print STDOUT "^##fq: Shared frequency of events occurred in same region, a value of 0 indicates all frequencies\n";
    print STDOUT "^#fq:change:chrom\t$sample_ids\n";
    for my $fq (sort {$a <=> $b} keys %{$Events_Stat{change}})
    {
        for my $change (sort keys %{$Categories{change}})
        {
            for my $chrom (sort keys %{$Categories{chrom}})
            {
                my @counts = ();
                for my $sample (@sample_ids)
                {
                    my $count = $Events_Stat{change}->{$fq}->{$sample}->{$chrom}->{$change} ?
                                $Events_Stat{change}->{$fq}->{$sample}->{$chrom}->{$change} : 0;
                    push @counts, $count;
                }
                
                my $counts = join "\t", @counts;
                print STDOUT "^$fq:$change:$chrom\t$counts\n";
            }
        }
    }
    
    print STDOUT "^#markers:fq\tsupport_marker_num\tfrequency\n";
    for my $fq (sort {$a <=> $b} keys %{$Events_Stat{marker}})
    {
        for my $marker_num (sort {$a <=> $b} keys %{$Events_Stat{marker}->{$fq}})
        {
            my $marker_cnt = $Events_Stat{marker}->{$fq}->{$marker_num};
            
            print STDOUT "^markers:$fq\t$marker_num\t$marker_cnt\n";
        }
    }

    
    print "^#length:fq\trange\tevent_length(x10bp)\tfrequency\n";
    for my $fq (sort {$a <=> $b} keys %{$Events_Stat{length}})
    {
        my $max_bin_end = (sort {$a <=> $b} keys %{$Events_Stat{length}->{$fq}})[-1];
        
        for my $bin (1 .. $max_bin_end)
        {
            my $num = $Events_Stat{length}->{$fq}->{$bin} ? $Events_Stat{length}->{$fq}->{$bin} : 0;
            
            my $bin_start = ($bin - 1) * 10;
            my $bin_end   = $bin_start + 10 - 1;
            
            print "^length:$fq\t$bin_start-$bin_end\t$bin\t$num\n";
        }
    }
}





=head2 check_continuous

    About   : Check whether the given numeric array is continuous, and return the
              count for each continuous group.
    Usage   : check_continuous(@numbers);
    Args    : Array of numbers.
    Returns : Array of counts of numbers in each continuous group.
    
=cut
sub check_continuous
{
    my (@numbers) = @_;
    
    if (@numbers <= 1) {
        return (1);
    }
    
    @numbers = sort {$a <=> $b} @numbers;
    
    my @groups = ();
    my $count  = 1;
    for (my $i=1; $i<@numbers; $i++)
    {
        if ($numbers[$i] - $numbers[$i-1] == 1) {
            $count ++;
        }
        else {
            push @groups, $count;
            $count = 1;
        }
        
        if ($i == (@numbers-1)) {
            push @groups, $count;
        }
    }
    
    return @groups;
}


