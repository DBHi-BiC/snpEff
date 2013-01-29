#!/usr/bin/perl

use strict;

my($READ_LEN) = 76;
my($INSERT_SIZE) = 100;
my($IDENTICAL_IDS) = 1;

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
my($seq, $seqCT, $seqGA, $q);
$seq = <STDIN>;
chomp $seq;

# Covert refrecen sequence
$seqCT = $seq;
$seqCT =~ tr/cC/TC/;
$seqGA = $seq;
$seqGA =~ tr/gG/AG/;

# Create fake quality
$q = "Z" x $READ_LEN;

# Create fake reads
my($len, $i, $start, $strandOri, $strand, $read, $read2, $id, $id2);
$len = length($seq) - 2 * $READ_LEN - $INSERT_SIZE - 1;
for( $i=0 ; $i <= $len ; $i++ ) {
	$start = $i;

	# Original strand
	if( ($i % 4) < 2 )	{ 
		$strandOri = 1;
		$read = substr( $seqCT , $start, $READ_LEN );
		$read2 = substr( $seqCT , $start + $READ_LEN + $INSERT_SIZE, $READ_LEN );
	} else {
		$strandOri = -1;
		$read = substr( $seqGA , $start, $READ_LEN );
		$read2 = substr( $seqGA , $start + $READ_LEN + $INSERT_SIZE, $READ_LEN );
	}
	$read2 = rwc( $read2 );

	# Read strand
	if( ($i % 2) == 0  ) { 
		$strand = 1;
	} else {
		$strand = -1;
		($read, $read2) = ($read2, $read);	# Reads are always 5' to 3'
	}

	# Create IDs
	if( ($strandOri == 1) && ($strand == 1) )	{ $id = "OT_"; }
	elsif( ($strandOri == 1) && ($strand == -1) )	{ $id = "MT_"; }
	elsif( ($strandOri == -1) && ($strand == 1) )	{ $id = "OB_"; }
	elsif( ($strandOri == -1) && ($strand == -1) )	{ $id = "MB_"; }

	$id = $id . ($start+1);
	if( $IDENTICAL_IDS ) { $id2 = $id ; }
	else {
		$id2 = $id . "_PE2";
		$id = $id . "_PE1";
	}

	# Show read in fastq format
	print FQ1 "\@$id\n$read\n+\n$q\n";
	print FQ2 "\@$id2\n$read2\n+\n$q\n";
}
