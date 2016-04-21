#!/usr/bin/perl -w
#
#   paintGenoGraph.pl -- slide window analysis of genotypes
#                          
#
#   Author: Nowind
#   Created: 2012-05-31
#   Updated: 2014-08-18
#   Version: 1.2.8
#
#   Change logs:
#   Version 1.0.0 12/06/13: The initial version.
#   Version 1.0.1 12/06/14: Set default output directory to current directory; add options 
#                           "--length", "--markers" and "--border".
#   Version 1.0.2 12/06/15: Modify file open operations.
#   Version 1.0.3 12/06/18: Bug fixed in output "empty" chromosomes.
#   Version 1.0.4 12/06/21: Add default color set; add option "--split".
#   Version 1.0.5 12/06/21: Add option "--contig".
#   Version 1.0.6 12/06/30: Bug fixed in read input file in function getGenomeLength;
#                           change chromosome orders same as the input file.
#   Version 1.1.0 12/07/19: Add options "--name" and "--color" to for user setting colors.
#   Version 1.1.1 12/07/20: Change option "--name" to "--tags"; add option "--prefix".
#   Version 1.1.2 12/09/23: Add option "--rows" to specify query rows of input files.
#   Version 1.1.3 13/01/06: Change the way to import functions from MyPerl::FileIO; bug fixed
#                           while block id contains character "-"; omit output of paintings
#                           with no data available.
#   Version 1.2.0 13/01/16: Remove options "--markers", "--length", "--tags", "--primary" and "--split";
#                           rename options "--contig" to "--length", "--blocks" to "--sort-blocks",
#                           "--border" to "--add-border".
#   Version 1.2.1 13/01/17: Bug fixed while outputting file in png format.
#   Version 1.2.2 13/01/27: Now user can specify colors by RGB values in the input block file; bug fixed
#                           while user set colors contain spaces.
#   Version 1.2.3 13/02/04: Add color "purple".
#   Version 1.2.4 14/01/06: Add several colors.
#   Version 1.2.5 14/05/28: Add two colors; add option "--split-pairs" to draw chromosome pairs.
#   Version 1.2.6 14/06/06: Bug fixed while pair tags were delimited with phased tag "|".
#   Version 1.2.7 14/06/09: Remove option "--prefix", now the prefix would be specified directly through
#                           "--output" option.
#   Version 1.2.8 14/08/18: Add several new colors.


use strict;

use Data::Dumper;
use Getopt::Long;
use File::Find::Rule;
use File::Basename;

use MyPerl::FileIO qw(:all);

################### Main ###################

my $CMDLINE = "perl $0 @ARGV";
my $VERSION = '1.2.8';
my $HEADER  = "# $CMDLINE\n# Version: $VERSION\n";


my $image_width    = 1000;                 ## resize fold
my $image_height   = 1000;
my $thickness      = 200;
my $block_distance = 200;
my $chrom_distance = 500;
my $image_type     = 'svg';
my (@input_files, $out_prefix, $border, $length_file,
    $user_colors, @query_rows, @block_orders, $split_pairs);
GetOptions(
            "input=s{,}"           => \@input_files,
            "output=s"             => \$out_prefix,
            "length=s"             => \$length_file,
            
            "rows=i{,}"            => \@query_rows,
            
            "C|colors=s"           => \$user_colors,

            "sort-blocks=s{,}"     => \@block_orders,
            
            "width=i"              => \$image_width,
            "height=i"             => \$image_height,
            "format=s"             => \$image_type,
            "thickness=i"          => \$thickness,
            
            "block-distance=i"     => \$block_distance,
            "D|chrom-distance=i"   => \$chrom_distance,
            
            "add-border"           => \$border,
            
            "split-pairs"          => \$split_pairs,
           );

unless( @input_files > 0 ) {
    print <<EOF;

$0  -- slide window analysis of genotypes

Version: $VERSION

Usage:   perl $0 [options]

Options:
    -i, --input   <filename>
        input file(s), required
    -o, --output  <string>
        output file prefix, default named as graph.* to current directory
    -l, --length  <filename>
        a file contains contig name and length in the format:
        
        #CHROM LENGTH
        chr01 43270923
        chr02 35937250
        chr03 36413819
        chr04 35502694
        ...
        
        required while no length info can be found in the
        input file

    -r, --rows   <int>
        specify the input fields, can have multiple values (0-based), in
        the order "block id", "chromosome", "type", "start position",
        "end position", ["RGB values"], if the last row is specified,
        user can set local defined colors directly from RGB values in the
        format "R,G,B" [default: 0 1 2 3 4]
    
    -f, --format    <string>
        image type, currently support png and svg [default: svg]
        *currently bugs could not be fixed while output png format,
        please use svg format instead
    
    -C, --colors <strings>
        specify colors for certain block types, in the format:
        
        "type1:color1;type2:color2;..."
        
        otherwise all types will be set to white by default, current support
        colors: white, black, red, blue, yellow, grey, light_grey, green,
        cyan, light_cyan
    
    -s, --sort-blocks  <strings>
        speicfy output blocks' orders, default use alphabetic sort
    
    -w, --width     <int>
        width of the image [default: 1000]
    -h, --height    <int>
        height of the image [default: 1000]
    -t, --thickness <int>
        height of each rectangle [default: 200]
    -D, --chrom-distance <int>
        vertical distance between each chromosomes [default: 500]
    -b, --block-distance  <int>
        vertical distance between each blocks [default: 200]

    -a, --add-border
        add border to image

    --split-pair
        draw chromosome pairs, this need the type row contain two
        types delimited by "/"
        
EOF

    exit(1);
}

$|++;

print STDERR "# $0 v$VERSION\n# " . (scalar localtime()) . "\n";

##
## set default values
##
unless( $out_prefix ) { $out_prefix = 'graph' ;}
unless(@query_rows){ @query_rows = (0, 1, 2, 3, 4) };



##
## paring blocks files
##
my %CHROM_LENGTH   = ();
my @CHROM_IDs      = ();
my %BLOCKS         = ();
my %BLOCK_IDs      = ();
my %user_rgb_set   = ();
for my $input (@input_files)
{
    print STDERR ">> Start parsing $input ... ";
    parse_block_file(\%BLOCKS, $input);
    print STDERR "done!\n";
}

my @BLOCK_IDs = (sort {$a cmp $b} keys %BLOCK_IDs);

if (@block_orders > 0)  {
    @BLOCK_IDs = @block_orders;
    
    if ($split_pairs) {
        @BLOCK_IDs = ();
        for my $id (@block_orders)
        {
            if ($BLOCK_IDs{$id}) {
                push @BLOCK_IDs, $id;
            }
            else {
                push @BLOCK_IDs, "$id-a";
                push @BLOCK_IDs, "$id-b";
                push @BLOCK_IDs, "$id-c";
            }
        }
    }
}


## get the length of chromosomes from a file if specified
if ($length_file) {
    getGenomeLength($length_file);
}


##
## determine which package to use
##
my $package  = ($image_type eq 'svg') ? 'GD::SVG' : 'GD';

eval "use $package";

my $image_pkg = $package . '::Image';
my $font_pkg  = $package . '::Font';


##
## set default colors
##
my %rgb_set  = (
                   'white'              => [255,255,255],
                   'black'              => [0,0,0],
                   
                   'red'                => [255,0,0],
                   'dark_red'           => [139,5,9],
                   'strong_red'         => [188,12,7],
                   'strong_red2'        => [194,21,23],   #C21517
                   'light_red1'         => [242,139,146],
                   'light_red2'         => [234,58,60], #ea3a3c
                   'red1'               => [140,0,0],
                   'red2'               => [179,0,0],
                   'red3'               => [218,0,0],
                   'red4'               => [255,8,8],
                    
                   'pink'               => [255,192,203],

                   'light_magenta'      => [222,65,255],  #de41ff
                   'strong_magenta'     => [204,0,204],   #cc00cc
                   
                   'blue'               => [0,0,255],
                   'strong_blue'        => [7,92,188],
                   'strong_blue2'       => [66,73,149],   #424995
                   'light_blue1'        => [119,128,229],
                   'light_blue2'        => [103,111,188], #686fbc
                   'slate_blue'         => [106,90,205],  #6a5acd
                   'royal_blue'         => [65,105,225],
                   'light_sky_blue'     => [135,206,250],
                   'blue1'              => [0,0,162],
                   'blue2'              => [0,0,213],
                   'blue3'              => [8,8,255],
                   'blue4'              => [60,60,255],
                   
                   'yellow'             => [255,255,0],
                   'yellow2'            => [252,179,21],  #fcb315
                   
                   'grey'               => [190,190,190],
                   'light_grey'         => [211,211,211],
                   
                   'green'              => [50,205,50],
                   'moderate_green'     => [132,205,90],  #84cd5a
                   'strong_green'       => [90,179,0],    #5ab300
                   'vivid_green'        => [185,255,40],  #b9ff28
                   
                   'cyan'               => [0,255,255],
                   'cyan2'              => [21,210,252],  #15d2fc
                   'light_cyan'         => [224,255,255],
                   'dark_grayish_cyan'  => [119,148,148], #779494
                   'moderate_cyan'      => [90,189,205],  #5abdcd
                   'strong_cyan'        => [0,179,179],   #00b3b3
                   
                   'honeydew'           => [240,255,240],
                   
                   'orange'             => [255,165,0],
                   'vivid_orange'       => [255,110,40],  #ff6e28
                   
                   'purple'             => [128,0,128],
                   
                   'moderate_violet'    => [164,90,205],  #a45acd
                   'vivid_violet'       => [110,40,255],  #6e28ff
                );


##
## set start and end positions
##
my $image_start = $image_width * 0.05;
my $image_end   = $image_width * 0.95;
my $fill_size   = $image_width * 0.90;

my $longest_chr = (sort {$a <=> $b} values %CHROM_LENGTH)[-1];

my $resize_fold = $longest_chr / $fill_size;


## create a new image
my $img   = $image_pkg->new($image_width,$image_height);


##
## allocate all colors
##
my %color_set = ();
   $color_set{white} = $img->colorAllocate(@{$rgb_set{white}});   # allocate white first as the background color
for my $color (keys %rgb_set)
{
    next if ($color eq 'white');
    $color_set{$color} = $img->colorAllocate(@{$rgb_set{$color}});
}

##
## set colors for each block type
##
my %user_color_set = ();
   $user_color_set{BLANK} = $color_set{white};
if ($user_colors) {
    $user_colors =~ s/\s+//g;
    for my $pair (split /\;/, $user_colors)
    {
        my ($type, $color) = (split /\:/, $pair);
        
        $user_color_set{$type} = $color_set{$color};
    }  
}


for my $type (keys %user_rgb_set)
{
    $user_color_set{$type} = $img->colorAllocate(@{$user_rgb_set{$type}});
}


## marks white as being transparent, make background transparent
$img->transparent($color_set{white});

## make background interlaced
$img->interlaced('true');



my $start_y = $image_height * 0.05;

for my $chrom (@CHROM_IDs)
{
    next unless( $BLOCKS{$chrom} );

    my $chrom_len = $CHROM_LENGTH{$chrom};
    my $width     = $chrom_len / $resize_fold;
    
    for my $block_id (@BLOCK_IDs)
    {   
        print STDERR "\rStart painting $chrom:$block_id ... ";
        
        draw_graphs($img, $BLOCKS{$chrom}->{$block_id},
                    $width, $start_y, $start_y + $thickness);
        
        $start_y += ($thickness + $block_distance);
    }
    print STDERR "done!\n";
    
    $start_y += ($thickness + $chrom_distance);

}




##
## write graph
##
my $outname = "$out_prefix.$image_type";

print STDERR ">> Start writing graph to $outname ... ";

open (DISPLAY, "> $outname") or die $!;

binmode DISPLAY;  # make sure we are writing to a binary stream 

my $img_data = $img->$image_type;

print DISPLAY $img_data;

close DISPLAY;



print STDERR "done!\n";

print STDERR "# " . (scalar localtime()) . "\n";

######################### Sub #########################

sub getGenomeLength
{
    my ($in) = shift;
    
    @CHROM_IDs = ();
    
    my $fh = getInputFilehandle($in);
    while (<$fh>)
    {
        next if (/\#/ || /^\s+$/);
        
        my ($CHROM, $LENGTH) = (split /\s+/);
        
        $CHROM_LENGTH{$CHROM} = $LENGTH;
    
        push @CHROM_IDs, $CHROM;
    }
}

sub draw_graphs
{
    my ($img, $ra_data, $width, $start_y, $end_y) = @_;

    $img->filledRectangle($image_start, $start_y,
                          $image_start + $width, $end_y,
                          $color_set{white});
    
    for (my $i=0; $i<=$#{$ra_data}; $i++)
    {
        my ($desc, $start, $end) = (split /;/, $ra_data->[$i]);

        my $x1 = $start / $resize_fold;
        my $x2 = $end   / $resize_fold;
        
        my $color = $user_color_set{$desc} ? $user_color_set{$desc}
                                           : $color_set{white};
        
        $img->filledRectangle($x1+$image_start, $start_y,
                              $x2+$image_start, $end_y, $color);
    }
    
    if ( $border ) {
        $img->rectangle($image_start, $start_y,
                        $image_start + $width, $end_y,
                        $color_set{black});
    }

}



=head2 parse_block_file

    About   : Parsing blocks file
    Usage   : parse_block_file(\%Blocks, $filename);
    Args    : Reference to a hash to save all block infos;
              Block filename.
    Returns : Array of sequence IDs.

=cut
sub parse_block_file
{
    my ($rh_blocks, $in) = @_;

    my $fh = getInputFilehandle($in);
    
    while (<$fh>)
    {
        if (/\#\#contig=<ID=(.*?),length=(\d+)>/) {
            next if ($CHROM_LENGTH{$1});
            
            $CHROM_LENGTH{$1} = $2;
            
            push @CHROM_IDs, $1;
        }
        
        next if (/\#/ || /^\s+$/);
        
        my ($block_id, $chrom, $desc, $start, $end, $rgb_values)
        = (split /\s+/)[@query_rows];
        
        
        if ($split_pairs && $desc =~ /(\/|\|)/) {
            my ($p1_type, $p2_type) = (split /\/|\|/, $desc);

            push @{$rh_blocks->{$chrom}->{"$block_id-a"}}, "$p1_type;$start;$end";
            push @{$rh_blocks->{$chrom}->{"$block_id-b"}}, "$p2_type;$start;$end";
            push @{$rh_blocks->{$chrom}->{"$block_id-c"}}, "BLANK;$start;$end";
            
            $BLOCK_IDs{"$block_id-a"}++;
            $BLOCK_IDs{"$block_id-b"}++;
            $BLOCK_IDs{"$block_id-c"}++;
        }
        else {
            $BLOCK_IDs{$block_id}++;
            push @{$rh_blocks->{$chrom}->{$block_id}}, "$desc;$start;$end";
        }
    
        
        if ($rgb_values) {
            $user_rgb_set{$desc} = [split /\,/, $rgb_values];
        }
    }
}

