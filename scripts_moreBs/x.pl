#!/usr/bin/perl

use strict;

my( $l, $id);
while( $l = <STDIN> ) {
	chomp $l;
	if( $l =~ /^#/ ) {
		# Ignore line
	} elsif( $l =~ /^>.*\|(.+)$/ ) {
		$id = $1;
	} else {
		my($start, $meth, $unmeth, $p) = split(/\t/, $l);
		$start++;
		my($tot) = $meth + $unmeth;

		if( $p < 0.1)	{ print "$id\t$start\t+\t$tot\t$meth\tCpG\n"; }
		else { print STDERR "Ignored: $l\n"; }
	}
}
