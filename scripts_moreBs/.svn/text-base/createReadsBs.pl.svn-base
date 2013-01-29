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
my($seq, $seqCT, $seqGA, $q);
$seq = <STDIN>;
chomp $seq;

# Covert reference sequence
$seqCT = $seq;
$seqCT =~ tr/cC/TC/;	# Methylated Cs are not converted to Ts
$seqGA = $seq;
$seqGA =~ tr/gG/AG/;	# Methylated Cs are not converted to Ts (opposite strand G not converted into A)

# Create fake quality
$q = "Z" x $READ_LEN;

# Create fake reads
my($len, $i, $start, $strandOri, $strand, $read, $id);
$len = length($seq) - $READ_LEN -1;
for( $i=0 ; $i <= $len ; $i++ ) {
	$start = $i;

	# Original strand
	if( ($i % 4) < 2 )	{ 
		$strandOri = 1;
		$read = substr( $seqCT , $start, $READ_LEN );
	} else {
		$strandOri = -1;
		$read = substr( $seqGA , $start, $READ_LEN );
	}

	# Read strand
	if( ($i % 2) == 0  ) { 
		$strand = 1;
	} else {
		$strand = -1;
		$read = rwc( $read );
	}

	# Create an ID
	if( ($strandOri == 1) && ($strand == 1) )	{ $id = "OT_"; }
	elsif( ($strandOri == 1) && ($strand == -1) )	{ $id = "MT_"; }
	elsif( ($strandOri == -1) && ($strand == 1) )	{ $id = "OB_"; }
	elsif( ($strandOri == -1) && ($strand == -1) )	{ $id = "MB_"; }
	$id .= ($start+1);

	# Show read in fastq format
	print "\@$id\n$read\n+\n$q\n";
}
