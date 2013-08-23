#############################################################################################
#
# SliceSILAC pipeline : Step 1
#
# Parse tab files from MaxQuant to screen ratios band by band.
# Look for entries for which at least 2 ratios are following these criteria:
# - at least one ratio below 'less_value' 
# - at least one ratio above 'more_value'
#
#
# Usage:
#  perl ./gel_mobility_2-0.pl <input_file_name> <less_value> <more_value> <slice_number> output_file_name>
#
#  perl					Perl interpreter
#  input_file_name		path to the list of protein groups 
# 						 (MaxQuant output as generated with SliceSILAC experimental design).
#  less_value			threshold for SILAC H/L ratios
#  more_value			threshold for SILAC H/L ratios
#  output_file_name		path to the output file
#
#
# Example of command line:
# 
#  Select protein groups with at least one ratio below 0.8 and one above 1.0 among 18 slices, 
# from file 'proteinGroups.txt' located in the same folder. The filtered list will be saved 
# in a file named 'filtered_output.txt'.
#
#    perl gel_mobility_2-0.pl ./proteinGroups.txt 0.8 1.0 18 filtered_output.txt
#
# 
# Script version : 2.0
# Date: 27-May-2011
# Author: Manfredo Quadroni
# Contact: wwwpaf@unil.ch
#
#############################################################################################

#!/usr/bin/perl


my($tabfile, $less_value, $more_value, $slice_nb, $output_filename) = @ARGV;

#open input file
open( INPUTF, '<', $tabfile ) || die("Cannot open $tabfile in r mode: $! \n");
#open output list file (selected protein groups)
open( OUTPUTF, '>', $output_filename ) || die("Cannot open $output_filename in w mode: $! \n");
#open file for miscellaneous debug information
open( OUTPUTFC, '> debug.txt' ) || die("Cannot open debug.txt in w mode: $! \n");

$index = 0;

#read all input file
@roba = <INPUTF>;
print OUTPUTFC "Input file : \t","$tabfile\n";
print OUTPUTFC "Less value : \t","$less_value\n";
print OUTPUTFC "More value : \t","$more_value\n";
print OUTPUTFC "Number of slices : \t","$slice_nb\n";
print OUTPUTFC "Output file : \t","$output_filename\n";


foreach $line (@roba) {
	chomp($line);
	if ( $index eq 0 ) {
		
		print OUTPUTF "$line\n";
		
		@headercut = split( /\t/, $line );
		for ( $i = 0 ; $i < 900 ; $i += 1 )    #scan header
		{
			if ( $headercut[$i] =~ m/Normalized band/ ) {
				push( @colsratio, $i );
			}
			;    # @colsratio is list of column numbers containing ratios
			if ( $headercut[$i] =~ m/Count band/ ) { push( @colscount, $i ); }
			;    # @colscount is list of column numbers containing counts
		}
		print OUTPUTFC "\ncolumns w. ratios\n";
		print OUTPUTFC "@colsratio\n", "\ncolumns with counts\n";
		print OUTPUTFC "@colscount\n", "\nmmmmmmmmm\n";
		print OUTPUTFC "\ncolumns w. ratios\t", "@colsratio\n",
		  "\ncolumns with counts\n", "@colscount\n", "\nmmmmmmmmm\n";

	}

	else {
		$loratio = 0;
		$hiratio = 0;
		@cut     = split( /\t/, $line );
		for ( $k = 0 ; $k < $slice_nb ; $k += 1 ) {
			$coluindex = $colsratio[$k];
			$rr[$k] = $cut[$coluindex];
		};    # $rr[$k] is ratio for band $k
		for ( $k = 0 ; $k < $slice_nb ; $k += 1 ) {
			$countindex = $colscount[$k];
			$cc[$k] = $cut[$countindex];
		};    # $cc[$k] is count for band $k

		print OUTPUTFC "$cut[0]\t", "ratios\t", "@rr\n";
		print OUTPUTFC "$cut[0]\t", "counts\t", "@cc\n";

		for ( $k = 0 ; $k < $slice_nb ; $k += 1 ) {  # $k contains band number index -1
			if ( ( $rr[$k] <= $less_value ) and ( $cc[$k] >= 1 ) and ( $rr[$k] ne "" ) )
			{
				$loratio += 1;
				print OUTPUTFC "$cut[0]\t", "$colsratio[$k]\t", "$rr[$k]\t", "lora\t",
				  "$loratio\n";
			}
			;    #"$rr[$k]\t","$cc[$k]\t","oooo\t";};
			if ( ( $rr[$k] > $more_value ) and ( $cc[$k] >= 1 ) ) {
				$hiratio += 1;
				print OUTPUTFC "$cut[0]\t", "$colsratio[$k]\t", "$rr[$k]\t", "hira\t",
				  "$hiratio\n";
			}
			;    #"$rr[$k]\t","$cc[$k]\n";};
		}

# BEWARE : do not set if (($rr[$k] <= 0.8) and ($cc[$k] > 0)) : it will take all void cases "" as good for loratio
	}

	if ( ( $loratio > 0 ) and ( $hiratio > 0 ) ) {
		print OUTPUTF "$line\n";
		print OUTPUTFC "ID  ", "$cut[0]\t", "taken\n";
		#print "ID  ", "$cut[0]\t", "taken\n";
	}
	$loratio = 0;
	$hiratio = 0;

	$index++;
}

#close all files
close(OUTPUTF);
close(INPUTF);
close(OUTPUTFC);

