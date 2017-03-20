#!/usr/bin/perl
use strict;
use warnings;
use Cwd;

my $toPeptides = 0;
my $rodeo = 0;
my $translate = 0;
my $hmmer = 0;
my $eutil = 0;
my $hmm = 0;
my $h3f = 0;
my $h3i = 0;
my $h3m = 0;
my $h3p = 0;
my $pfamA = 0;

my $dir = getcwd;
opendir(my $dh, $dir) || die;
while( my $curr_file = readdir($dh) ) {

	if ($curr_file eq "rodeo.pl") {
		$rodeo = 1;
		print "FOUND rodeo.pl\n";
		system("chmod +x rodeo.pl");
	}

	if ($curr_file eq "toPeptides.pl") {
		$toPeptides = 1;
		print "FOUND toPeptides.pl\n";
		system("chmod +x toPeptides.pl");
	}

	if ($curr_file eq "translate.pl") {
		$translate = 1;
		print "FOUND translate.pl\n";
		system("chmod +x translate.pl");
	}

	if ($curr_file eq "Pfam-A.hmm") {
		$hmm = 1;
		print "FOUND Pfam-A.hmm\n";
	}

	if ($curr_file eq "Pfam-A.hmm.h3f") {
		$h3f = 1;
		print "FOUND Pfam-A.hmm.h3f\n";
	}

	if ($curr_file eq "Pfam-A.hmm.h3i") {
		$h3i = 1;
		print "FOUND Pfam-A.hmm.h3i\n";
	}

	if ($curr_file eq "Pfam-A.hmm.h3m") {
		$h3m = 1;
		print "FOUND Pfam-A.hmm.h3m\n";
	}

	if ($curr_file eq "Pfam-A.hmm.h3p") {
		$h3p = 1;
		print "FOUND Pfam-A.hmm.h3p\n";
	}
}

if ($hmm == 1 || $h3f == 1 || $h3i == 1 || $h3m == 1 || $h3p == 1) {
	$pfamA = 1;
}

my $hmmer_test = qx(hmmscan);
my $eutil_test = qx(esearch -db protein -query 1);
my $perl_test = qx(perl -v);


if ( index($hmmer_test, "Incorrect number of command line arguments.") != -1 ){
	$hmmer = 1;
	print "HMMER INSTALLED\n";
}

if ( index($eutil_test, "<ENTREZ_DIRECT>") != -1){
	$eutil = 1;
	print "ENTREZ_DIRECT INSTALLED\n";
}

if ($toPeptides == 0 || $rodeo == 0 || $eutil == 0 || $hmmer == 0 || $translate == 0 || $pfamA == 0) {
	print "\nINSTALLATION INCOMPLETE. MISSING THE FOLLOWING REQUIREMENTS:\n\n";

	if ( $toPeptides == 0 ) {
		print "toPeptides.pl\n";
	}

	if ( $rodeo == 0 ) {
		print "rodeo.pl\n";
	}

	if ( $eutil == 0 ) {
		print "ENTREZ_DIRECT\n";
	}

	if ( $hmmer == 0 ) {
		print "HMMER\n";
	}

	if ( $translate == 0 ) {
		print "translate.pl\n";
	}

	if ( $pfamA == 0 ) {
		print "Pfam-A NOT FULLY INSTALLED OR IN THE WRONG LOCATION, MISSING THE FOLLOWING:\n";
		if ($hmm == 0) {
			print "Pfam-A.hmm\n";
		}
		if ($h3f == 0) {
			print "Pfam-A.hmm.h3f\n";
		}
		if ($h3i == 0) {
			print "Pfam-A.hmm.h3i\n";
		}
		if ($h3m == 0) {
			print "Pfam-A.hmm.h3m\n";
		}
		if ($h3p == 0) {
			print "Pfam-A.hmm.h3p\n";
		}
		print "THE ABOVE FILES NEED TO BE IN THE SAME DIRECTORY AS THIS INSTALLATION SCRIPT\n";
	}

}

else {
	print "\nINSTALLATION SUCCESSFUL\n\n";
}








