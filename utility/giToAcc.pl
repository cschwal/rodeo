#!/usr/bin/perl
use strict;
use warnings;
use LWP::Simple;

my $file = $ARGV[0];
my $dest = $ARGV[1];
open(INFO, $file);
open(FILE, '>>', $dest);
foreach my $line (<INFO>)
{
	chomp($line);
	my $gi = $line;

	my $base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
	my $url = $base .= "efetch.fcgi?db=protein&id=$gi&rettype=acc";

	my $output = get($url);
	chomp($output);
	print FILE "$gi,$output";

}


