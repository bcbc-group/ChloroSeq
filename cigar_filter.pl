#!/usr/bin/perl

use strict;
use warnings;

#Get file names from command line
if (!$ARGV[0]){
    print "Give an input sam file. \n";
    exit 1;
}

my $inlist = $ARGV[0];

#Open input file
open (INLIST, "<$inlist") || die "Cannot open list file.\n";

while (<INLIST>) {
    chomp($_);
    if ($_ =~ /^@/){
	print $_ . "\n";
	next;
    }
    my @feature = split (/\t/, $_);

    if ($feature[5] =~ /(.*M)([0-9]+)(N.*)/){

	if ($2 < 2500){ 
	    print join ("\t", @feature) . "\n";
	}
	else {next;}
    }
    else {print join ("\t", @feature) . "\n";}
}

  
