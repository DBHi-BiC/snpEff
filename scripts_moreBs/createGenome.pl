#!/usr/bin/perl

#------------------------------------------------------------------------------
#
# Create a fake genome 'methylated' genome
# 'C' bases at position multiple of 100 are methylated (in either strand)
#
#------------------------------------------------------------------------------

my($CHROMOS) = 1;
my($GENOME_LEN) = 10 * 1000;
my($REPEAT_LEN) = $GENOME_LEN / 10;

$chromo = "";

for( $c=1 ; $c <= $CHROMOS ; $c++ ) {
	print ">X$c\n";
	for( $i=1 ; $i <= $GENOME_LEN ; $i++ ) {
		my($r) = int(4 * rand());

		$base = "";

		#---
		# Create bases: {a, c|C, g|G, t} 
		# where 
		# 		'C' stand for "methylated C" 
		# 		'G' for "methylaterd C on minus strand"
		#---
		if( $r == 0 ) { $base = 'a'; }	
		elsif( $r == 1 ) { 
			if( $i % 100 == 0 ) { $base = 'C'; }
			else { $base = 'c'; }
		} elsif( $r == 2 ) { 
			if( $i % 100 == 0 ) { $base = 'G'; }
			else { $base = 'g'; }
		}	elsif( $r == 3 ) { $base = 't'; }	

		$chromo .= $base;
		print $base;
		if( ($i % 80) == 0 )	{ print "\n"; }
	}
}

# Create a 'Repeat' part
$repeat = substr($chromo, 0, $REPEAT_LEN);
print ">R\n";
print "$repeat\n";
