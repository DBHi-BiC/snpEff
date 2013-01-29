#!/bin/sh

#-------------------------------------------------------------------------------
# From Alexander :
#	Currently I'm using MACS v.1.4 with params like:
#	--tsize 36
#	--mfold 4,30
#	--format=BED
#	-g 2.7e9 // as default
#	--diag --fe-min=0 --fe-max=30 --fe-step=2  // for diagnostics
#
# I didn't try new beta MACS 1.4rc2.
#-------------------------------------------------------------------------------

export PYTHONPATH=$HOME/tools/MACS-1.4.0beta/lib/python2.7/site-packages/
$HOME/tools/MACS-1.4.0beta/bin/macs14 $*

Alexander Mazur
 to me
	
show details 2:37 PM (6 hours ago)
	
Hi Pablo,

