#!/usr/bin/perl

use strict;

my($READ_LEN) = 76;
my($INSERT_SIZE) = 100;

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

# Parse command line
my($fq1, $fq2) = ($ARGV[0], $ARGV[1]);
if( $fq2 eq '' )	{ die "Usage: createReadsPe.pl outFastq1 outFastq2\n"; }
open FQ1, "> $fq1";
open FQ2, "> $fq2";

# Read reference sequence as 1 line
my($seq, $q);
$seq = <STDIN>;
chomp $seq;

# Create fake quality ('A' is Sanger 32)
$q = "A" x $READ_LEN;

# Create fake reads
my($len, $i, $start, $strandOri, $strand, $read, $read2, $id, $id2);
$len = length($seq) - $READ_LEN - $INSERT_SIZE - 1;
for( $i=1 ; $i <= $len ; $i++ ) {
	$start = $i;

	$read = substr( $seq , $start, $READ_LEN );
	$read2 = substr( $seq , $start + $INSERT_SIZE, $READ_LEN );
	$read2 = rwc($read2);

	# Read strand
	$strand = rand();
	if( $strand >= 0.5 ) { 
		$strand = 1;
		$id = "O_";
	} else {
		$strand = -1;
		$id = "M_";
		($read, $read2) = ($read2, $read);	# Reads are always 5' to 3'
	}

	# Create an ID
	$id2 = $id . ($start+1) . "_PE2";
	$id = $id . ($start+1) . "_PE1";

	# Show read in fastq format
	print FQ1 "\@$id\n$read\n+\n$q\n";
	print FQ2 "\@$id2\n$read2\n+\n$q\n";
}
