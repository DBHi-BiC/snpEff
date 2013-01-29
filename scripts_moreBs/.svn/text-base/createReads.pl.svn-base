#!/usr/bin/perl

use strict;

my($READ_LEN) = 76;

#------------------------------------------------------------------------------
# Reverse WC
#------------------------------------------------------------------------------
sub rwc($) {
	my($s) = @_;
	return reverse( wc($s) );
}

#------------------------------------------------------------------------------
# WC-complement
#------------------------------------------------------------------------------
sub wc($) {
	my($s) = @_;
	$s = uc $s;
	$s =~ tr/ACGT/TGCA/;
	return $s;
}

#------------------------------------------------------------------------------
# Main
#------------------------------------------------------------------------------

# Read reference sequence as 1 line
my($seq, $q);
$seq = <STDIN>;
chomp $seq;

# Create fake quality
$q = "Z" x $READ_LEN;

# Create fake reads
my($len, $i, $start, $strandOri, $strand, $read, $id);
$len = length($seq) - $READ_LEN -1;
for( $i=1 ; $i <= $len ; $i++ ) {
	$start = $i;

	$read = substr( $seq , $start, $READ_LEN );

	# Read strand
	$strand = rand();
	if( $strand >= 0.5 ) { 
		$id = "O_";
	} else {
		$read = rwc( $read );
		$id = "M_";
	}

	# Create an ID
	$id .= ($start+1);

	# Show read in fastq format
	print "\@$id\n$read\n+\n$q\n";
}
