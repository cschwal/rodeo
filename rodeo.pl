#!/usr/bin/perl
#
# Copyright (C) 2015 Jonathan I. Tietz
# University of Illinois
# Department of Chemistry
#
# Copyright (C) 2015 Christopher J. Schwalen
# University of Illinois
# Department of Chemistry
#
# Copyright (C) 2015 Parth S. Patel
# University of Illinois
# Department of Chemistry
#
# Copyright (C) 2015 Douglas A. Mitchell
# University of Illinois
# Department of Chemistry
#
# License: GNU Affero General Public License v3 or later
# Complete license availabel in the accompanying LICENSE.txt.
# or <http://www.gnu.org/licenses/>.
#
# This file is part of RODEO.
#
# RODEO is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# RODEO is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#

$ENV{PERL_LWP_SSL_VERIFY_HOSTNAME}=0;

use strict;
#use warnings;
use POSIX;
use Time::Piece;
use List::Util qw[min max];
use LWP::Simple;
use LWP::UserAgent;
use File::Temp;
#for debugging
#use Data::Dumper;
 
my $version_string = "20170216b (CJS212a)";  # include date/version number as well as random key for verification
 
my $startTime = time; 
my $endTime;
my $numNeighbors = 6;
my $duration = 0;
my @gis = ();
my @GIarray;
my %fxnMatch;
my %PfamMatch;
my %secondPfamMatch;
my %thirdPfamMatch;
my %PfamName;
my %secondPfamName;
my %thirdPfamName;
my %thirdPfamEval;
my %PfamEval;
my %secondPfamEval;
my %geneINFO;
my @outputPrecs;
my @outputPrecsInfo;
my $giCount = 0;
my @fxns;
my @fxnsNoHyps;
my @scores;
my @precs;
my $peptideHTML = "";
my $species;
my $genus;
my @neighborGIs;
my $precMin = 30;
my $precMax = 100;
my $toc = "";
my $evaluateAllIDs = 0;
my %customPfams;
my %secondCustPfam;
my %customPfamEvals;
my %userHMMInfo;
my @allFoundPfams;
my %otherLassoPfamsEval;
my $evalAll = 0;
my %custHmmMatch;
my %translations;
my @accessions;
my $main_acc;
 
my $cPFam = "PF00733";
my $ePFam = "PF05402";
my $bPFam = "PF13471";
my $g1PFam = "PF00326";  
my $g2PFam = "PF00593";
 
# flags that are up when we should BLAST 
my $cFromHMMER = 0;
my $bFromHMMER = 0;
my $userPfam = 0;
 
my $firstGIinCluster;
my $intergenicHigh;
my $intergenicLow;
my $intergeneicPeptideLen;
my $intergenicGap;
my $overlap = 600;
 
my $cGI = "NONE";
my $bGI = "NONE";
my $eGI = "NONE";
my $dGI = "NONE";
my $fGI = "NONE";
my $gGI = "NONE";
my $mostLikelyPrec = "";
my $leaderCore = "";
my $coreStart = 0;
my $coreEnd = 1;
my $tie = 1;
 
my $mainGI;
my $architectureHTML;
my $isList = 0;
my %userPFamInfo;
my %giColors;
my $scaffHTML = "";
my @scaffs;
my @neighborINFO;
my $idGI = 0;
my $bypassScore = 0;
my @custHmms;
my $hmmdir = "hmm";
 
# get the complete command line for debugging or reproducibility purposes
my $command_line = join " ", $0, @ARGV;
 
my( $indexX )= grep { $ARGV[$_] eq "-x" } 0..$#ARGV;
if (defined($indexX))
{
    $bypassScore = 1;
}
 
my( $indexL )= grep { $ARGV[$_] eq "-l" } 0..$#ARGV;
my( $indexP )= grep { $ARGV[$_] eq "-p" } 0..$#ARGV;
my( $indexC )= grep { $ARGV[$_] eq "-c" } 0..$#ARGV;
my( $indexO )= grep { $ARGV[$_] eq "-o" } 0..$#ARGV;
my( $indexA )= grep { $ARGV[$_] eq "-a" } 0..$#ARGV;
my( $indexE )= grep { $ARGV[$_] eq "-e" } 0..$#ARGV;
my( $indexCSV )= grep { $ARGV[$_] eq "-csv" } 0..$#ARGV;
my( $indexCSVa )= grep { $ARGV[$_] eq "-csva" } 0..$#ARGV;
my( $indexGI )= grep { $ARGV[$_] =~ /[0-9]/ } 0..$#ARGV;
 
if (defined($indexP)) {
    if (defined($ARGV[$indexP+1])) {
            $hmmdir = $ARGV[$indexP+1];
    }
    @custHmms = getCustomHmms();
}
 
if (defined($indexCSV))
{
    if ($bypassScore == 0) 
    {
        open MEEJ, '>>', $ARGV[$indexCSV+1];
        print MEEJ "query acc,genus/species,nucleotide acc,leader,core,min,max,distance,coreMass,score,within 500 nuc?,within 150 nuc?,greater than 1000?,2/4 cysteines?,leader longer than core?,plausible lasso ring?,GxxxxxT?,does core start with G?,same direction?,ratio length leader/core < 2 and > 0.5,starts with C and even number of C?,has G in core?,core has at least 1 aromatic?,at least 2 aromatic?,core has odd number of C?,leader has W?,leader has K?,leader has C?,cluster has C?,E?,B?\n";
        print MEEJ ",,,,,,,,WEIGHT,1,1,-1,1,1,1,3,2,2,1,1,-2,1,-1,-2,1,-1,2\n";
        close MEEJ;
    }
    else
    {
        open MEEJ, '>>', $ARGV[$indexCSV+1];
        print MEEJ "query acc,genus/species,nucleotide acc,min,max,dir,peptide\n";
        close MEEJ;
    }
}
 
if (defined($indexCSVa))
{
    open MEEJ, '>>', $ARGV[$indexCSVa+1];
    print MEEJ "query accession,genus/species,nucleotide acc,protein acc,fxn,dir,min,max,PfamID1,description,E-value1,PfamID2,description,E-value2,PfamID3,description,E-value3,custom hmm,custom hmm E-value,custom hmm2,custom hmm2 E-value\n";
    close MEEJ;
}
 
my $printPrecs = 1;
 
if(!defined($indexL) && !defined($indexC))
{

    my $file = $ARGV[$indexO+1];

    printHeader($file);

    $mainGI = $ARGV[$indexGI];
    if ($mainGI =~ /[A-Z]+/) {

        my $base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
        my $url = $base .= "efetch.fcgi?db=protein&id=$mainGI&rettype=gi";
        my $ua = LWP::UserAgent->new;
        my $req = HTTP::Request->new(GET => $url);
        my $res = $ua->request($req);
        $mainGI = $res->as_string;
	my @resArray = split("\n",$mainGI);
        $mainGI = @resArray[-1];


       chomp($mainGI);

    }


    main($mainGI);
    if (defined($indexCSVa))
    {
        printCSVa();
    }
    if (defined($indexA)) {
        evalAllIDgis();
    }
     
    $endTime = time;
    $duration = ($endTime - $startTime);
    print FILE <<ENDHTML;
    <p>RUN TIME: $duration SECONDS</p></html>
ENDHTML
    close(FILE);
    if( defined($indexE) )
    {
        #emailFile();
        my $date = localtime->strftime('%y%m%d');
        my $dee = "RODEOoutput" . $date;
        if (defined($indexCSV) && defined($indexCSVa))
        {
            system("(uuencode $file $file;uuencode $ARGV[$indexCSV+1] $ARGV[$indexCSV+1];uuencode $ARGV[$indexCSVa+1] $ARGV[$indexCSVa+1];) | mail -s '$dee' $ARGV[$indexE+1]");
        }
        elsif (defined($indexCSV))
        {
            system("(uuencode $file $file;uuencode $ARGV[$indexCSV+1] $ARGV[$indexCSV+1];) | mail -s '$dee' $ARGV[$indexE+1]");
        }
        elsif (defined($indexCSVa))
        {
            system("(uuencode $file $file;uuencode $ARGV[$indexCSVa+1] $ARGV[$indexCSVa+1];) | mail -s '$dee' $ARGV[$indexE+1]");
        }
        else
        {
            system("uuencode $file $file | mail -s '$dee' $ARGV[$indexE+1]");
        }
    }
    exit();
}
 
if ( defined($indexL) && defined($indexC) )
{
    my @slurp;
    open (FH, $ARGV[$indexC+1]);
    @slurp = <FH>;
    close FH;
 
    $slurp[0] =~ /([0-9]+)/ ;
    $numNeighbors = $1;
 
    $slurp[1] =~ /([0-9]+)/ ;
    $overlap = $1;
 
    $slurp[2] =~ /([0-9]+)/ ;
    $precMin = $1;
     
    $slurp[3] =~ /([0-9]+)/ ;
    $precMax = $1;
 
    for(my $i = 5; $i < @slurp; $i++)
    {
        my @meej = split / /, $slurp[$i];
        if ( index($meej[1], "HMM") == -1 ) 
        {
            $userPFamInfo{$meej[3]}[0] = $meej[0]; # what the user named the pfam
            if( defined($meej[4]))
            {
                $userPfam = 1;
                chomp($meej[4]);
                $userPFamInfo{$meej[3]}[1] = $meej[4]; # what the user wants it colored 
            }
        }
 
        if ( index($meej[1], "HMM") != -1 ) 
        {
            my $name_len = scalar @meej;
            my $long_name = $meej[2];
            if (scalar @meej > 3) {
                for(my $j = 3; $j < scalar @meej - 1; $j++) {
                    $long_name .= " $meej[$j]";
                }
            }
            $userHMMInfo{$long_name}[0] = $meej[0]; # what the user named the hmm
            chomp($meej[3]);
            $userHMMInfo{$long_name}[1] = $meej[$name_len - 1]; # what the user wants it colored 
            $userPfam = 1;
             
        }
         
    }
 
    $printPrecs = index($slurp[4], "YES");
 
    $isList = 1;
    my $file = $ARGV[$indexO+1];
    printHeader($file);
    open (INFO, $ARGV[$indexL+1]) or die;
    my $gi_counter = 0;
    foreach my $line (<INFO>) 
    {             
        chomp($line);
        $mainGI = $line;
        if ($line =~ /[A-Z]+/) {
#            $mainGI = get($url);
        my $base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
        my $url = $base .= "efetch.fcgi?db=protein&id=$mainGI&rettype=gi";
        my $ua = LWP::UserAgent->new;
        my $req = HTTP::Request->new(GET => $url);
        my $res = $ua->request($req);
        $mainGI = $res->as_string;
        my @resArray = split("\n",$mainGI);
        $mainGI = @resArray[-1];

            chomp($mainGI);
            $gi_counter += 1;
            print "Processing $line ($gi_counter) ... ";
            main($mainGI);
            if (defined($indexCSVa))
        {
            printCSVa();
        }
            if (defined($indexA)) {
                evalAllIDgis();
                $evalAll = 0;
            }
            print "Done\n"; 
        }
        elsif ($line =~ /([0-9]+)/) {
            $mainGI = $1;
            $gi_counter += 1;
            print "Processing $mainGI ($gi_counter) ... ";
            main($mainGI);
            if (defined($indexCSVa))
        {
            printCSVa();
        }
            if (defined($indexA)) {
                evalAllIDgis();
                $evalAll = 0;
            }
             
            print "Done\n";         
        }
 
        else {
            next;
        }
         
    CONTINUE:   
        clear();
    }
    close(INFO);
    $endTime = time;
    $duration = ($endTime - $startTime)/60;
    $duration = substr($duration, 0, 4);
    print FILE <<ENDHTML;
    <p>RUN TIME: $duration MIN</p></html>
ENDHTML
    close(FILE);
    if( defined($indexE) )
    {
        my $date = localtime->strftime('%y%m%d');
        my $dee = "RODEOoutput" . $date;
        if (defined($indexCSV) && defined($indexCSVa))
        {
            system("(uuencode $file $file;uuencode $ARGV[$indexCSV+1] $ARGV[$indexCSV+1];uuencode $ARGV[$indexCSVa+1] $ARGV[$indexCSVa+1];) | mail -s '$dee' $ARGV[$indexE+1]");
        }
        elsif (defined($indexCSV))
        {
            system("(uuencode $file $file;uuencode $ARGV[$indexCSV+1] $ARGV[$indexCSV+1];) | mail -s '$dee' $ARGV[$indexE+1]");
        }
        elsif (defined($indexCSVa))
        {
            system("(uuencode $file $file;uuencode $ARGV[$indexCSVa+1] $ARGV[$indexCSVa+1];) | mail -s '$dee' $ARGV[$indexE+1]");
        }
        else
        {
            system("uuencode $file $file | mail -s '$dee' $ARGV[$indexE+1]");
        }
    }
    exit;
}
 
if ( defined($indexL) )
{
 
    $isList = 1;
    my $file = $ARGV[$indexO+1];
    printHeader($file);
    open (INFO, $ARGV[$indexL+1]) or die;
    my $gi_counter = 0;
    foreach my $line (<INFO>) 
    {               
        chomp($line);
        $mainGI = $line;
         
        if ($line =~ /[A-Z]+/) {
#       $mainGI = get($url);
        my $base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
        my $url = $base .= "efetch.fcgi?db=protein&id=$mainGI&rettype=gi";
	my $ua = LWP::UserAgent->new;
        my $req = HTTP::Request->new(GET => $url);
        my $res = $ua->request($req);
        $mainGI = $res->as_string;
        my @resArray = split("\n",$mainGI);
        $mainGI = @resArray[-1];


            chomp($mainGI);
            $gi_counter += 1;
            print "Processing $line ($gi_counter) ... ";
            main($mainGI);
            if ( defined($indexA) ) {
                evalAllIDgis();
                $evalAll = 0;
            }
            print "Done\n"; 
        }
        elsif ($line =~ /([0-9]+)/) {
            $mainGI = $1;
            $gi_counter += 1;
            print "Processing $mainGI ($gi_counter) ... ";
            main($mainGI);
            if ( defined($indexA) ) {
                evalAllIDgis();
                $evalAll = 0;
            }
            print "Done\n";         
        }
 
        else {
            next;
        }
 
        if (defined($indexCSVa))
        {
            printCSVa();
        }
    CONTINUE:   
        clear();
    }
    close(INFO);
    $endTime = time;
    $duration = ($endTime - $startTime)/60;
    $duration = substr($duration, 0, 4);
    print FILE <<ENDHTML;
    <p>RUN TIME: $duration MIN</p></html>
ENDHTML
    close(FILE);
    if( defined($indexE) )
    {
        my $date = localtime->strftime('%y%m%d');
        my $dee = "RODEOoutput" . $date;
        if (defined($indexCSV) && defined($indexCSVa))
        {
            system("(uuencode $file $file;uuencode $ARGV[$indexCSV+1] $ARGV[$indexCSV+1];uuencode $ARGV[$indexCSVa+1] $ARGV[$indexCSVa+1];) | mail -s '$dee' $ARGV[$indexE+1]");
        }
        elsif (defined($indexCSV))
        {
            system("(uuencode $file $file;uuencode $ARGV[$indexCSV+1] $ARGV[$indexCSV+1];) | mail -s '$dee' $ARGV[$indexE+1]");
        }
        elsif (defined($indexCSVa))
        {
            system("(uuencode $file $file;uuencode $ARGV[$indexCSVa+1] $ARGV[$indexCSVa+1];) | mail -s '$dee' $ARGV[$indexE+1]");
        }
        else
        {
            system("uuencode $file $file | mail -s '$dee' $ARGV[$indexE+1]");
        }
    }
}
 
if ( defined($indexC) )
{
    my @slurp;
    open (FH, $ARGV[$indexC+1]);
    @slurp = <FH>;
    close FH;
 
    $slurp[0] =~ /([0-9]+)/ ;
    $numNeighbors = $1;
 
    $slurp[1] =~ /([0-9]+)/ ;
    $overlap = $1;
 
    $slurp[2] =~ /([0-9]+)/ ;
    $precMin = $1;
     
    $slurp[3] =~ /([0-9]+)/ ;
    $precMax = $1;
 
    for(my $i = 5; $i < @slurp; $i++)
    {
        my @meej = split / /, $slurp[$i];
        if ( index($meej[1], "HMM") == -1 ) {
            $userPFamInfo{$meej[3]}[0] = $meej[0];
            if( defined($meej[4]))
            {
                $userPfam = 1;
                chomp($meej[4]);
                $userPFamInfo{$meej[3]}[1] = $meej[4]; # maps PFam ID -> color
            }
            else
            {
                $userPFamInfo{$meej[3]}[1] = "white";
            }
        }
     
 
        if ( index($meej[1], "HMM") != -1 ) 
        {
            my $name_len = scalar @meej;
            my $long_name = $meej[2];
            if (scalar @meej > 3) {
                for(my $j = 3; $j < scalar @meej - 1; $j++) {
                    $long_name .= " $meej[$j]";
                }
            }
            $userHMMInfo{$long_name}[0] = $meej[0]; # what the user named the hmm
            chomp($meej[3]);
            $userHMMInfo{$long_name}[1] = $meej[$name_len - 1]; # what the user wants it colored 
            $userPfam = 1;
            #print "$long_name and $meej[0] color it $meej[$name_len - 1]\n";
        }
         
    }
 
    $printPrecs = index($slurp[4], "YES");
 
    my $file = $ARGV[$indexO+1];
    printHeader($file);
    $mainGI = $ARGV[$indexGI];
    if ($mainGI =~ /[A-Z]+/) {
#        $mainGI = get($url);
        my $base = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
        my $url = $base .= "efetch.fcgi?db=protein&id=$mainGI&rettype=gi";
        my $ua = LWP::UserAgent->new;
        my $req = HTTP::Request->new(GET => $url);
        my $res = $ua->request($req);
        $mainGI = $res->as_string;
        my @resArray = split("\n",$mainGI);
        $mainGI = @resArray[-1];


        chomp($mainGI);
    }
 
    main($mainGI);  
    if (defined($indexCSVa))
    {
        printCSVa();
    }
    if (defined($indexA)) {
        evalAllIDgis();
                $evalAll = 0;
   }
     
    $endTime = time;
    $duration = ($endTime - $startTime);
    print FILE <<ENDHTML;
    <p>RUN TIME: $duration SECONDS</p></html>
ENDHTML
    if( defined($indexE) )
    {
        my $date = localtime->strftime('%y%m%d');
        my $dee = "RODEOoutput" . $date;
        if (defined($indexCSV) && defined($indexCSVa))
        {
            system("(uuencode $file $file;uuencode $ARGV[$indexCSV+1] $ARGV[$indexCSV+1];uuencode $ARGV[$indexCSVa+1] $ARGV[$indexCSVa+1];) | mail -s '$dee' $ARGV[$indexE+1]");
        }
        elsif (defined($indexCSV))
        {
            system("(uuencode $file $file;uuencode $ARGV[$indexCSV+1] $ARGV[$indexCSV+1];) | mail -s '$dee' $ARGV[$indexE+1]");
        }
        elsif (defined($indexCSVa))
        {
            system("(uuencode $file $file;uuencode $ARGV[$indexCSVa+1] $ARGV[$indexCSVa+1];) | mail -s '$dee' $ARGV[$indexE+1]");
        }
        else
        {
            system("uuencode $file $file | mail -s '$dee' $ARGV[$indexE+1]");
        }
    }
    close(FILE);
}
 
sub evalAllIDgis
{
    #print "running evalAllIDgis\n";
    my @tempScaffs = @scaffs;
    my $num_of_scaffs = scalar @scaffs;
    #print "running on $num_of_scaffs scaffs";
    for (my $r = 1; $r < @tempScaffs; $r++) {
        $idGI = $tempScaffs[$r];
        #print "    Evaluating nucleotide record with GI $idGI ... ";
        clearIDGI();
        main($mainGI);
        if (defined($indexCSVa))
    {
        printCSVa();
    }
        print "DONE\n";
    }
}
 
sub main 
{
    @gis = fetchNeighboringGIS($mainGI);
 
    # first use hmmer to see if cluster is interesting....
    for(my $i = 0; $i < @gis; $i++)
    {
        chomp($gis[$i]);
        hmmer($gis[$i], $cPFam, $ePFam, $bPFam, $g1PFam, $g2PFam);
    }
 
        printFxns();
        getIntergenicRegions();
 
    for(my $o = 0; $o < @gis-1; $o++)
    {
        if (defined($geneINFO{$gis[$o]}[0])) {
            if( $geneINFO{$gis[$o]}[0] eq "C")
            {
                $cGI = $gis[$o];
            }
            if($geneINFO{$gis[$o]}[0] eq "B")
            {
                $bGI = $gis[$o];
            }
            if($geneINFO{$gis[$o]}[0] eq "E")
            {
                $eGI = $gis[$o];
            }
            if($geneINFO{$gis[$o]}[0] eq "D")
            {
                if($dGI ne "NONE")
                {
                    $dGI .= "+$gis[$o]";
                }
                else
                {
                    $dGI = $gis[$o];
                }
            }
            if($geneINFO{$gis[$o]}[0] eq "F")
            {
                $fGI = $gis[$o];
            }
            if($geneINFO{$gis[$o]}[0] eq "G")
            {
                $gGI = $gis[$o];
            }
        }
    }
    getArcHTML();
    toHTML();
}
 
sub printCSVa
{
    my $lass = "";
    open MEEJ, '>>', $ARGV[$indexCSVa+1];
    for (my $k = 0; $k < @gis; $k++) 
    {
	chomp($main_acc);
	chomp($accessions[$k]);
	my $base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
	my $url = $base .= "efetch.fcgi?db=protein&id=$idGI&rettype=acc";
#	my $id_acc = get($url);

        my $ua = LWP::UserAgent->new;
        my $req = HTTP::Request->new(GET => $url);
        my $res = $ua->request($req);
	my $id_acc = $res->as_string;
        my @resArray = split("\n",$id_acc);
        $id_acc = @resArray[-1];


	chomp($id_acc);
	chomp($id_acc);
	
        print MEEJ "$main_acc,$genus,$id_acc,$accessions[$k],";
        for (my $l = 0; $l < 4; $l++) 
        {
            print MEEJ"$geneINFO{$gis[$k]}[$l],";
        }
        if ( index($PfamName{$gis[$k]} , ",") != -1)  {
            $PfamName{$gis[$k]} =~ tr/,/ /;
        }
        if ($PfamMatch{$gis[$k]} ne "NO PFAM MATCH")
        {
            print MEEJ "$PfamMatch{$gis[$k]},$PfamName{$gis[$k]},$PfamEval{$gis[$k]},";
            if ( exists($secondPfamMatch{$gis[$k]}))
            {
                if ( index($secondPfamName{$gis[$k]} , ",") != -1)  {
                    $secondPfamName{$gis[$k]} =~ tr/,/ /;
                }
                if (exists($thirdPfamName{$gis[$k]})) {
                    if ( index($thirdPfamName{$gis[$k]} , ",") != -1)  {
                        $thirdPfamName{$gis[$k]} =~ tr/,/ /;
                    }
                    print MEEJ "$secondPfamMatch{$gis[$k]},$secondPfamName{$gis[$k]},$secondPfamEval{$gis[$k]},$thirdPfamMatch{$gis[$k]},$thirdPfamName{$gis[$k]},$thirdPfamEval{$gis[$k]},";
                }
                else {
                    print MEEJ "$secondPfamMatch{$gis[$k]},$secondPfamName{$gis[$k]},$secondPfamEval{$gis[$k]},,,,";
                }
            }
            else {
                    print MEEJ ",,,,,,";
            }
        }
        else {
            print MEEJ "NO PFAM MATCH,,,,,,,,,";
        }
        

        if (defined($custHmmMatch{$gis[$k]})) 
        {
            my @hmmEvals;
            if( exists($custHmmMatch{$gis[$k]})) 
            {
                my @hmmMatchAndEval = split/,/,$custHmmMatch{$gis[$k]};
                #@hmmMatchAndEval = sort { $a <=> $b } @hmmMatchAndEval;
                for(my $i = 1; $i < @hmmMatchAndEval; $i+=2) {
                    push(@hmmEvals,$hmmMatchAndEval[$i]);
                }
                @hmmEvals = sort { $a <=> $b } @hmmEvals;
                for (my $i = 0; $i < @hmmEvals; $i++) {
                    for(my $j = 1; $j < @hmmMatchAndEval; $j+=2) {
                        if ($hmmMatchAndEval[$j] == $hmmEvals[$i]) {
                            print MEEJ "$hmmMatchAndEval[$j-1],$hmmMatchAndEval[$j],";
                        }
                    }
                }
 
            }
        }
 
        #else
        #{
            print MEEJ "\n";
        #}
 
        $lass = "";
    }
    close MEEJ;
}
 
sub printHeader
{
    if ( defined($indexL) )
    {
        open (INFO, $ARGV[$indexL+1]) or die;
        foreach my $line (<INFO>) 
        {               
            chomp($line);
            if ($line =~ /([A-Za-z0-9_.]+)/) {
                push(@GIarray, $line);              
            }
        }
        close(INFO);
 
        $toc = "<h3>Input Queries (click to navigate)</h3>";
        $toc .= '<ul style="list-style-type:none">';
        for (my $i = 0; $i < @GIarray; $i++) {
            $toc .= '<li><a href="#';
            $toc .= "$GIarray[$i]";
            $toc .= '"';
            $toc .= ">$GIarray[$i]";
            $toc .= "</a></li>";
        }
        $toc .= "</ul>";
    }
 
    open(FILE, '>', $_[0]);
 
my $time = localtime;
print FILE <<ENDHTML;
<html>
<head>
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap.min.css">
    <link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.6/css/bootstrap-theme.min.css">
</head>
<style media="screen" type="text/css">
 
.square {
  width: 54px;
  height: 14px;
  background-color: white;
  outline: #ffffff solid 1px;
  text-align: center;
  line-height: 14px;
  font-size: 12px;
}
 
table {
   font-size: 11px;
}
 
</style>
 
<script src='https://img.jgi.doe.gov//js/overlib.js'></script>
<div class="container">
<h1 align="center" id="header">RODEO</h1>
<div class="row">
     <div class="col-md-5">
        <h3>Parameters</h3>
             <table class="table table-condensed" style="width:100%;">
                    <tr><th scope="row">Version</th><td>$version_string</td></tr>
                    <tr><th scope="row">Command</th><td><code style="color:black;">$command_line</code></td></tr>
                    <tr><th scope="row">Run Time</th><td>$time</td></tr>
                    <tr><th scope="row">Gene Window</th><td>+/- $numNeighbors genes</td></tr>
                    <tr><th scope="row">Peptide Range</th><td>$precMin&ndash;$precMax aa</td></tr>
                    <tr><th scope="row">Fetch Distance</th><td>$overlap bp</td></tr>
              </table>
     </div>
     <div class="col-md-6">
                <div class="panel panel-default">
                 
ENDHTML
 
if ( defined($indexC) )
{
    print FILE <<ENDHTML;
<div class="panel-heading">Annotation Legend</div>
                <table class="table table-striped">
                        <tr>
                            <th>Appearance</th>
                            <th>Pfam or pHMM name</th>
                            <th>Annotation</th>
                             
                        </tr>
ENDHTML
    foreach my $group (keys %userPFamInfo) 
    {
        for(my $y = 0; $y < @{$userPFamInfo{$group}}; $y+=2) 
        {
        print FILE "<td><div class='square' style='outline-color:black;background-color:$userPFamInfo{$group}[1]; color:grey;'>$userPFamInfo{$group}[0]</div></td>";
        print FILE "<td>$group</td>";
        print FILE "<td>$userPFamInfo{$group}[0]</td>"; 
        print FILE '</tr>';
        }
    }
 
    foreach my $group (keys %userHMMInfo) 
    {
        for(my $y = 0; $y < @{$userHMMInfo{$group}}; $y+=2) 
        {
        print FILE "<td><div class='square' style='outline-color:black;background-color:$userHMMInfo{$group}[1]; color:grey;'>$userHMMInfo{$group}[0]</div></td>";
        print FILE "<td>$group</td>";
        print FILE "<td>$userHMMInfo{$group}[0]</td>"; 
        print FILE '</tr>';
        }
    }
}        
print FILE <<ENDHTML;           
 
</table>
</div></div>
</div>
$toc
ENDHTML
 
}
 
 
sub getArcHTML
{
    my $width;
    my $arrowWid;
    my $arrowWid2;
    my $arrowWid3;
    my $arcStart = 0;
    my $subBy = $geneINFO{$gis[0]}[2] - 500;
    my $scaleFactor = (660/( $geneINFO{$gis[@gis-1]}[3] - $geneINFO{$gis[0]}[2]) );
    my $letterX;
    $architectureHTML .= '<svg width="1060" height="53">';
 
    my $fillColor = "white";
     
    for(my $guy = 0; $guy <  @gis; $guy++)
    {
        if (defined($geneINFO{$gis[$guy]}[0])) {
            $fillColor = setFillColor($geneINFO{$gis[$guy]}[0], $gis[$guy]);
        } 
        if( $geneINFO{$gis[$guy]}[1] eq "fwd" )
        {
            $arrowWid = floor( ($geneINFO{$gis[$guy]}[2] - $subBy) * $scaleFactor );
            $arrowWid3 = floor( ($geneINFO{$gis[$guy]}[3] - $subBy) * $scaleFactor );
            if( ($arrowWid3-$arrowWid) < 40 )
            {
                $arrowWid2 = ($arrowWid + $arrowWid3)/2;
            }
            else
            {
                $arrowWid2 = $arrowWid3 - 20;
            }
            $letterX = floor( ($arrowWid+$arrowWid3)/2);
            $architectureHTML .= '<polygon points="';
            $architectureHTML .= "$arrowWid,10 $arrowWid,40 $arrowWid2,40 $arrowWid2,50 $arrowWid3,25 $arrowWid2,0 $arrowWid2,10 $arrowWid,10";
            if ($bypassScore == 0 || ($bypassScore == 1 && $userPfam == 1) )
            {
                $architectureHTML .= '" style="fill:';
                $architectureHTML .= "$fillColor";
                $architectureHTML .= ';stroke:black;stroke-width:.5" onMouseOver="return overlib(';
                #$architectureHTML .= "'$gis[$guy] - $PfamMatch{$gis[$guy]}  : $PfamName{$gis[$guy]}'";
                chomp($accessions[$guy]);
                $architectureHTML .= "'$accessions[$guy] - $PfamMatch{$gis[$guy]}  : $PfamName{$gis[$guy]}'";
                $architectureHTML .= ')" onMouseOut="return nd()"/>
                 <text x="';
                 $architectureHTML .= "$letterX";
                 $architectureHTML .= '"y="32"
                font-family="sans-serif"
               font-size="12px"
              text-anchor="middle"
              fill="grey">';
              $architectureHTML .= "$geneINFO{$gis[$guy]}[0]";
              $architectureHTML .= '</text>';
            }
            else
            {
                $architectureHTML .= '" style="fill:';
                $architectureHTML .= "$fillColor";
                $architectureHTML .= ';stroke:black;stroke-width:.5" onMouseOver="return overlib(';
                chomp($accessions[$guy]);
                $architectureHTML .= "'$accessions[$guy] - $PfamMatch{$gis[$guy]}  : $PfamName{$gis[$guy]}'";
                $architectureHTML .= ')" onMouseOut="return nd()"/>';
            }
        }
 
        else
        {
            $arrowWid3 = floor( ($geneINFO{$gis[$guy]}[2] - $subBy) * $scaleFactor );
            $arrowWid = floor( ($geneINFO{$gis[$guy]}[3] - $subBy) * $scaleFactor );
            if( ($arrowWid-$arrowWid3) < 40 )
            {
                $arrowWid2 = ($arrowWid + $arrowWid3)/2;
            }
            else
            {
                $arrowWid2 = $arrowWid3 + 20;
            }
            $letterX = floor( ($arrowWid+$arrowWid3)/2);
            $architectureHTML .= '<polygon points="';
            $architectureHTML .= "$arrowWid,10 $arrowWid,40 $arrowWid2,40 $arrowWid2,50 $arrowWid3,25 $arrowWid2,0 $arrowWid2,10 $arrowWid,10";
            if ($bypassScore == 0 || ($bypassScore == 1 && $userPfam == 1) )
            {
                $architectureHTML .= '" style="fill:';
                $architectureHTML .= "$fillColor";
                $architectureHTML .= ';stroke:black;stroke-width:.5" onMouseOver="return overlib(';
		chomp($accessions[$guy]);
                $architectureHTML .= "'$accessions[$guy] - $PfamMatch{$gis[$guy]}  : $PfamName{$gis[$guy]}'";
                $architectureHTML .= ')" onMouseOut="return nd()"/>
                 <text x="';
                 $architectureHTML .= "$letterX";
                 $architectureHTML .= '"y="32"
                font-family="sans-serif"
               font-size="12px"
              text-anchor="middle"
              fill="grey">';
              $architectureHTML .= "$geneINFO{$gis[$guy]}[0]";
              $architectureHTML .= '</text>';
            }
            else
            {
                $architectureHTML .= '" style="fill:';
                $architectureHTML .= "white";
                $architectureHTML .= ';stroke:black;stroke-width:.5" onMouseOver="return overlib(';
		chomp($accessions[$guy]);
                $architectureHTML .= "'accessions[$guy] - $PfamMatch{$gis[$guy]}  : $PfamName{$gis[$guy]}'";
                $architectureHTML .= ')" onMouseOut="return nd()"/>';
            }
        }
    }
    $architectureHTML .= '</svg>';
    addPrecursorHTML($scaleFactor, $subBy);
 
    my $barLength = $scaleFactor * 1000;
    my $barLegX = $barLength + 5;
    $architectureHTML .= '<svg width="500" height="23">';
    $architectureHTML .= '<polygon points="';
    $architectureHTML .= "0,10 $barLength,10";
    $architectureHTML .= '" style="fill:white;stroke:black;stroke-width:.5" />';
    $architectureHTML .= '<text x="';
    $architectureHTML .= "$barLegX";
    $architectureHTML .= '"y="12"
            font-family="sans-serif"
           font-size="10px"
          text-anchor="right"
          fill="black">1000 nucleotides</text>';
    $architectureHTML .= '<polygon points="';
    $architectureHTML .= "0,5 0,15";
    $architectureHTML .= '" style="fill:white;stroke:black;stroke-width:.5" />';
    $architectureHTML .= '<polygon points="';
    $architectureHTML .= "$barLength,5 $barLength,15";
    $architectureHTML .= '" style="fill:white;stroke:black;stroke-width:.5" />';
    $architectureHTML .= '</svg>';
}
 
sub clear
{
    @accessions = ();
    $idGI = 0;
    @allFoundPfams = ();
    $architectureHTML = "";
    @gis = ();
    @GIarray = ();
    %fxnMatch = ();
    %PfamMatch = ();
    %PfamName = ();
    %geneINFO = ();
    $giCount = 0;
    @fxns = ();
    @fxnsNoHyps = ();
    @scores = ();
    @precs = ();
    $peptideHTML = "";
    $cFromHMMER = 0;
    $bFromHMMER = 0;
    $firstGIinCluster = 0;
    $intergenicHigh = 0;
    $intergenicLow = 0;
    $intergeneicPeptideLen = 0;
    $intergenicGap = 0;
    $cGI = "NONE";
    $bGI = "NONE";
    $eGI = "NONE";
    $dGI = "NONE";
    $fGI = "NONE";
    $gGI = "NONE";
    $mostLikelyPrec = "";
    $leaderCore = "";
    $coreStart = 0;
    $coreEnd = 1;
    $tie = 1;
    @neighborGIs = ();
    @neighborINFO = ();
    @scaffs = ();
    $scaffHTML = "";
    @outputPrecsInfo = ();
    @outputPrecs = ();
    %PfamEval = ();
    %secondPfamEval = ();
}
 
sub clearIDGI
{
    @allFoundPfams = ();
    $architectureHTML = "";
    @gis = ();
    @GIarray = ();
    %fxnMatch = ();
    %PfamMatch = ();
    %PfamName = ();
    %geneINFO = ();
    $giCount = 0;
    @fxns = ();
    @fxnsNoHyps = ();
    @scores = ();
    @precs = ();
    $peptideHTML = "";
    $cFromHMMER = 0;
    $bFromHMMER = 0;
    $firstGIinCluster = 0;
    $intergenicHigh = 0;
    $intergenicLow = 0;
    $intergeneicPeptideLen = 0;
    $intergenicGap = 0;
    $cGI = "NONE";
    $bGI = "NONE";
    $eGI = "NONE";
    $dGI = "NONE";
    $fGI = "NONE";
    $gGI = "NONE";
    $mostLikelyPrec = "";
    $leaderCore = "";
    $coreStart = 0;
    $coreEnd = 1;
    $tie = 1;
    @neighborGIs = ();
    @neighborINFO = ();
    @scaffs = ();
    $scaffHTML = "";
    @outputPrecsInfo = ();
    @outputPrecs = ();
    %PfamEval = ();
    %secondPfamEval = ();
    %thirdPfamEval = ();
}
 
sub addPrecursorHTML
{
    my $width;
    my $letterX;
    my $arrowWid;
    my $arrowWid2;
    my $arrowWid3;
    my $scaleFactor = $_[0];
    my $subBy = $_[1];
    $architectureHTML .= '<svg width="1060" height="53">';
    for(my $q = 0; $q < @outputPrecsInfo; $q++)
    {
        if( $outputPrecsInfo[$q][2] eq "f")
        {
            $arrowWid = floor( ($outputPrecsInfo[$q][0] - $subBy) * $scaleFactor );
            $arrowWid3 = floor( ($outputPrecsInfo[$q][1] - $subBy) * $scaleFactor );
            $letterX = floor( ($arrowWid+$arrowWid3)/2);
            if( ($arrowWid3-$arrowWid) < 40 )
            {
                $arrowWid2 = ($arrowWid + $arrowWid3)/2;
            }
            else
            {
                $arrowWid2 = $arrowWid3 - 20;
            }
            $architectureHTML .= '<polygon points="';
            $architectureHTML .= "$arrowWid,10 $arrowWid,40 $arrowWid2,40 $arrowWid2,50 $arrowWid3,25 $arrowWid2,0 $arrowWid2,10 $arrowWid,10";
            $architectureHTML .= '" style="fill:';
            $architectureHTML .= "$outputPrecsInfo[$q][3]";
            $architectureHTML .= ';stroke:black;stroke-width:.5" />
             <text x="';
             $architectureHTML .= "$letterX";
             $architectureHTML .= '"y="32"
            font-family="sans-serif"
           font-size="12px"
          text-anchor="middle"
          fill="grey">';
          $architectureHTML .= "$outputPrecsInfo[$q][4]";
          $architectureHTML .= '</text>';
        }
 
        else
        {
            $arrowWid3 = floor( ($outputPrecsInfo[$q][0] - $subBy) * $scaleFactor );
            $arrowWid = floor( ($outputPrecsInfo[$q][1] - $subBy) * $scaleFactor );
            $letterX = floor( ($arrowWid+$arrowWid3)/2);
            if( ($arrowWid-$arrowWid3) < 40 )
            {
                $arrowWid2 = ($arrowWid + $arrowWid3)/2;
            }
            else
            {
                $arrowWid2 = $arrowWid3 + 20;
            }
            $architectureHTML .= '<polygon points="';
            $architectureHTML .= "$arrowWid,10 $arrowWid,40 $arrowWid2,40 $arrowWid2,50 $arrowWid3,25 $arrowWid2,0 $arrowWid2,10 $arrowWid,10";
            $architectureHTML .= '" style="fill:';
            $architectureHTML .= "$outputPrecsInfo[$q][3]";
            $architectureHTML .= ';stroke:black;stroke-width:.5" />
             <text x="';
             $architectureHTML .= "$letterX";
             $architectureHTML .= '"y="32"
            font-family="sans-serif"
           font-size="12px"
          text-anchor="middle"
          fill="grey">';
          $architectureHTML .= "$outputPrecsInfo[$q][4]";
          $architectureHTML .= '</text>';     
        }
    }
    $architectureHTML .= '</svg>';
}
 
sub setFillColor
{
    my $fxn = $_[0];
    my $scob = $_[1];
    if( exists($giColors{$scob}) )
    {
        return $giColors{$scob};
    }
    if( $fxn eq "C")
    {
        return "#377eb8";
    }
    if( $fxn eq "E")
    {
        return "#ffff33";
    }
    if( $fxn eq "B")
    {
        return "#ff7f00";
    }
    if( $fxn eq "D")
    {
        return "#984ea3";
    }
    if( $fxn eq "G1")
    {
        return "gray";
    }
    if( $fxn eq "G2")
    {
        return "gray";
    }
    return "white";
}
 
sub getIntergenicRegions
{
    my $start = 0;
    my $end = 0;
 
    my $clusterStart = 0;
    my $temp;
 
    for(my $j = 0; $j < @fxns-1; $j++)
    {
        if( index($fxns[$j], "?") == -1 )
        {
            # find when the core of the cluster starts
            if(index($fxns[$j], "C") != -1 || index($fxns[$j], "E") != -1 || index($fxns[$j], "B") != -1)
            {
                if($coreStart == 0)
                {
                    $coreStart = $j;
                }
            }
 
            if ($clusterStart == 0)
            {
                $clusterStart = $j;
            }
            $j++;
 
            while( index($fxns[$j], "?") != -1 && $j < @fxns-1)
            {
                $j++;
            }
 
            $temp = $j;
            while( (index($fxns[$temp], "C") != -1 || index($fxns[$temp], "E") != -1 || index($fxns[$temp], "B") != -1) && $temp < @fxns-1)
            {
                $temp++;
            }
            if( $coreEnd == 1)
            {
                $coreEnd = $temp-1;
            }
 
            if(index($fxns[$j], "?") == -1 )
            {
                $end = $j;
                if( index($fxns[$j], "C") != -1 || index($fxns[$j], "E") != -1 || index($fxns[$j], "B") != -1 )
                {
                    $coreEnd = $end;
                }
                $j--;
            }
        }
    }
 
     
    my $gene1start = $geneINFO{$gis[$clusterStart]}[2];
    my $gene1end = $geneINFO{$gis[$clusterStart]}[3];
 
    my $coreStartCoord = $geneINFO{$gis[$coreStart]}[2];
    my $coreEndCoord = $geneINFO{$gis[$coreEnd]}[3];
 
    my $fwdCount = 0;
    my $revCount = 0;
    my $clusterDir;
 
    for(my $c = $coreStart; $c <= $coreEnd; $c++)
    {
        if ($geneINFO{$gis[$c]}[1] eq "fwd")
        {
            $fwdCount++;
        }
 
        else
        {
            $revCount++;
        }
    }
 
    if($fwdCount > $revCount)
    {
        $clusterDir = "f";
    }
 
    else
    {
        $clusterDir = "r";
    }
 
    my $gene2start = $geneINFO{$gis[$end]}[2];
    my $gene2end = $geneINFO{$gis[$end]}[3];
 
    my $max1 = max($gene1start, $gene1end);
    my $max2 = max($gene2start, $gene2end);
 
    my $min1 = min($gene1start, $gene1end);
    my $min2 = min($gene2start, $gene2end);
 
    my $clusterMin = min($min1, $min2);
    my $clusterMax = max($max1, $max2);
 
    $clusterMax += $overlap;
    $clusterMin -= $overlap;
    if ($clusterMin < 0)
    {
        $clusterMin = 0;
    }
    my $distanceOverall = 1000000;
    my $distanceStartToStart = 0;
    my $distanceEndToStart = 0;
    my $distanceStartToEnd = 0;
    my $distanceEndToEnd = 0;
    my $peptideCoordMin = 0;
    my $peptideCoordMax = 0;
 
    my $curMin = 10000;
    my $score = 0;
    my $highScore = 0;
 
    system("esearch -db nuccore -query $idGI | efetch -format fasta -seq_start $clusterMin -seq_stop $clusterMax > geneCluster.fasta");

    system("./translate.pl geneCluster.fasta");
    system("perl toPeptides.pl $precMin $precMax $bypassScore >> peptides.txt");
 
    my $ct = 0;
 
    open (INFO, "<peptides.txt") or die "Could not open file\n";
    foreach my $line (<INFO>)
    {
        $score = 0;
        chomp($line);
        my $pepDir = $line =~ /([a-z]$)/;
        $pepDir = $1;
 
        my $pep = $line =~ m/^([A-Z]+)/;
        $pep = $1;
 
        my $peptideCoord = $line =~ /(\d+)/;
        $peptideCoord = $1;
 
        my @scoreArray;
 
        if ($pepDir eq "f")
        {
            $peptideCoordMin = $peptideCoord + $clusterMin - 1;
            $peptideCoordMax = $peptideCoordMin + ( (length($pep) + 1) * 3 ) - 1;
        }
 
        else
        {
            $peptideCoordMax = $peptideCoord + $clusterMin - 4;
            $peptideCoordMin = $peptideCoordMax - ( (length($pep) + 1) * 3 ) + 1;
        }
 
        for(my $p = 0; $p < @gis; $p++)
        {
            my $genLen = $geneINFO{$gis[$p]}[3] - $geneINFO{$gis[$p]}[2];
            my $genLenAA = $genLen/3;
            if ( $peptideCoordMin >= $geneINFO{$gis[$p]}[2] && $peptideCoordMax <= $geneINFO{$gis[$p]}[3] 
                && $genLenAA > $precMax)
            { 
                $distanceOverall = "x";
            }
 
            my $q = $p+1;
            if ($q < @gis && $peptideCoordMin < $geneINFO{$gis[$p]}[3] 
                && $peptideCoordMax > $geneINFO{$gis[$q]}[2] && $genLenAA > $precMax) {
                $distanceOverall = "x";
            }
        }
 
        if($bypassScore == 0)
        {
            for(my $p = $coreStart; $p <= $coreEnd; $p++)
            {
         
                if ( (($peptideCoordMin >= $geneINFO{$gis[$p]}[2] && $peptideCoordMin < $geneINFO{$gis[$p]}[3])
                 && $peptideCoordMax > $geneINFO{$gis[$p]}[3]) || ($peptideCoordMin < $geneINFO{$gis[$p]}[2] && $peptideCoordMax > $geneINFO{$gis[$p]}[2]) )
                { 
                    $distanceOverall = 0;
                }
 
                if($peptideCoordMin > $geneINFO{$gis[$p]}[2] && $peptideCoordMax < $geneINFO{$gis[$p]}[3] && $distanceOverall ne "x")
                {
                    $distanceOverall = 0;
                }
             
            }
 
            if ($distanceOverall ne "x" && $distanceOverall != 0)
            {
                my $minDistance = 10000;
 
                for (my $i = 0; $i < @gis; $i++) {
                    if ($geneINFO{$gis[$i]}[0] eq "C" || $geneINFO{$gis[$i]}[0] eq "E" || $geneINFO{$gis[$i]}[0] eq "B") {
 
                        my $geneStart = $geneINFO{$gis[$i]}[2];
                        my $geneEnd = $geneINFO{$gis[$i]}[3];
 
                        my $distanceTemp = min(abs($peptideCoordMin - $geneStart), abs($peptideCoordMin - $geneEnd), abs($peptideCoordMax - $geneStart), 
                            abs($peptideCoordMax - $geneEnd));
 
                        if ($distanceTemp == abs($peptideCoordMin - $geneStart)) {
                            $distanceTemp = $peptideCoordMin - $geneStart;
                        }
                        elsif ($distanceTemp == abs($peptideCoordMin - $geneEnd)) {
                            $distanceTemp = $peptideCoordMin - $geneEnd;
                        }
                        elsif ($distanceTemp == abs($peptideCoordMax - $geneStart)) {
                            $distanceTemp = $peptideCoordMax - $geneStart;
                        }
                        elsif ($distanceTemp == abs($peptideCoordMax - $geneEnd)) {
                            $distanceTemp = $peptideCoordMax - $geneEnd;
                        }
 
                        if ( abs($distanceTemp) < abs($minDistance) ) {
                            $minDistance = $distanceTemp;
                        }
                    }
                }
 
                if ($minDistance != 10000)  # meaning there is a C, E, or B in the cluster
                {
                    $distanceOverall = $minDistance;    
                }
 
                else {
                    my $inputStart = $geneINFO{$mainGI}[2];
                    my $inputEnd = $geneINFO{$mainGI}[3];
 
                    $distanceOverall = min( abs($peptideCoordMin - $inputStart), abs($peptideCoordMin - $inputEnd), abs($peptideCoordMax - $inputStart), 
                        abs($peptideCoordMax - $inputEnd));
 
                    if ($distanceOverall == abs($peptideCoordMin - $inputStart)) {
                        $distanceOverall = $peptideCoordMin - $inputStart;
                    }
                    elsif ($distanceOverall == abs($peptideCoordMin - $inputEnd)) {
                        $distanceOverall = $peptideCoordMin - $inputEnd;
                    }
                    elsif ($distanceOverall == abs($peptideCoordMax - $inputStart)) {
                        $distanceOverall = $peptideCoordMax - $inputStart;
                    }
                    elsif ($distanceOverall == abs($peptideCoordMax - $inputEnd)) {
                        $distanceOverall = $peptideCoordMax - $inputEnd;
                    }
                }
 
=pod                
                $distanceStartToStart = $peptideCoordMin - $coreStartCoord; # negative means outside of cluster
                $distanceEndToEnd = $coreEndCoord - $peptideCoordMax;
                $distanceStartToEnd = $peptideCoordMin - $coreEndCoord;
                $distanceEndToStart = $peptideCoordMax - $coreStartCoord;
                $distanceOverall = min( abs($distanceStartToStart), abs($distanceEndToEnd), abs($distanceStartToEnd), abs($distanceEndToStart) );
                if ($distanceOverall == abs($distanceStartToStart))
                {
                    $distanceOverall = $distanceStartToStart;
                }
                if ($distanceOverall == abs($distanceEndToEnd))
                {
                    $distanceOverall = $distanceEndToEnd;
                }
                if ($distanceOverall == abs($distanceStartToEnd))
                {
                    $distanceOverall = $distanceStartToEnd;
                }
                if ($distanceOverall == abs($distanceEndToStart))
                {
                    $distanceOverall = $distanceEndToStart;
                }
=cut                
            }
 
            if ($distanceOverall ne "x")
            {
                if( ($distanceOverall != 0) && abs($distanceOverall) < abs($curMin))
                {
                     
                        $curMin = $distanceOverall;
                        #$mostLikelyPrec = $pep;
                }
 
                if( abs($distanceOverall) <= 500 && $distanceOverall != 0)
                {
                    $score += 1;
                    push(@scoreArray, 1);
                    if( abs($distanceOverall) <= 150)
                    {
                        $score += 1;
                        push(@scoreArray, 1);
                    }
                    else
                    {
                        push(@scoreArray, 0);
                    }
                }
                else
                {
                    push(@scoreArray, 0);
                    push(@scoreArray, 0);
                }
 
 
                if( abs($distanceOverall) >= 1000 )
                {
                    $score -= 1;
                    push(@scoreArray, 1);
                }
                else
                {
                    push(@scoreArray, 0);
                }
 
                if( abs($distanceOverall) >= 12000 )
                {
                    next;
                }
 
            }
 
            else
            {
                $score -= 5;
                push(@scoreArray, 1);
                push(@scoreArray, 0);
                push(@scoreArray, 0);
                push(@scoreArray, 0);
            }
             
             
            # check if ring formed by core peptide is reasonable
            $leaderCore = extractLeaderCore($pep);
 
            my $coreIndex = index($leaderCore, ",");
            my $leader = substr($pep, 0, $coreIndex);
            my $core = substr($leaderCore, $coreIndex+2);
            my $plausibleRing = "NO";
 
            my $cCount = () = $core =~ /C/g;
 
            if( $cCount == 2 || $cCount == 4)
            {
                $score += 1;
                push(@scoreArray, 1);
            }
            else
            {
                push(@scoreArray, 0);
            }
 
            if( length($leader) > length($core) - 3 )
            {
                $score += 2;
                push(@scoreArray, 1);
            }
            else
            {
                push(@scoreArray, 0);
            }
 
            if( substr($core, 7, 1) eq "E" || substr($core, 8, 1) eq "E" || substr($core, 8, 1) eq "D" || substr($core, 9, 1) eq "D")
            {
                $plausibleRing = "YES";
                $score += 1;
                push(@scoreArray, 1);
            }
            else
            {
                push(@scoreArray, 0);
            }
 
            my $indexPenT = rindex($leader, "T");
            my $g = "";
            if ($indexPenT != -1)
            {
                $g = substr($leader, $indexPenT-6, 1 );
            }
            # checking for GxxxxxT motif
            if( $g eq "G")
            {
                $score += 3;
                push(@scoreArray, 1);
            }
            else
            {
                push(@scoreArray, 0);
            }
 
 
            if ( substr($core, 0, 1) eq "G" )
            {
                $score += 2;
                push(@scoreArray, 1);
            }
            else
            {
                push(@scoreArray, 0);
            }
 
            if( $clusterDir eq $pepDir)
            {
                $score += 1;
                push(@scoreArray, 1);
            }
            else
            {
                push(@scoreArray, 0);
            }
 
            my $ratio = length($leader)/length($core);
            if($ratio < 2 && $ratio > 0.5)
            {
                $score += 1;
                push(@scoreArray, 1);
            }
            else
            {
                push(@scoreArray, 0);
            }
 
            # core starts with C and has even number of C
            if ( substr($core, 0, 1) eq "C" && $cCount%2 == 0)
            {
                $score += 0;
                push(@scoreArray, 1);
            }
            else
            {
                push(@scoreArray, 0);
            }
 
            if ( index($core, "G") == -1)
            {
                $score -= 4;
                push(@scoreArray, 1);
            }
            else
            {
                push(@scoreArray, 0);
            }
 
            my $wCount = () = $core =~ /W/g;
            my $fCount = () = $core =~ /F/g;
            my $yCount = () = $core =~ /Y/g;
 
            # core has at least 1 aromatic residue
            if ( $wCount >= 1 || $fCount >= 1 || $yCount >= 1)
            {
                $score += 1;
                push(@scoreArray, 1);
            }
            else
            {
                push(@scoreArray, 0);
            }
 
            my $aromaticResInCore = $wCount + $fCount + $yCount;
            # core has at least 2 aromatic residues 
            if ($aromaticResInCore > 1)
            {
                $score += 2;
                push(@scoreArray, 1);
            }
            else
            {
                push(@scoreArray, 0);
            }
 
            # core has odd number of C
            if ( $cCount%2 == 1)
            {
                $score -= 2;
                push(@scoreArray, 1);
            }
            else
            {
                push(@scoreArray, 0);
            }
 
            $wCount = () = $leader =~ /W/g;
            if ($wCount > 0)
            {
                $score -= 1;
                push(@scoreArray, 1);
            }
            else
            {
                push(@scoreArray, 0);
            }
 
            my $kCount = () = $leader =~ /K/g;
            if ($kCount > 0)
            {
                $score += 1;
                push(@scoreArray, 1);
            }
            else
            {
                push(@scoreArray, 0);
            }
 
            $cCount = () = $leader =~ /C/g;
            if ($cCount > 0)
            {
                $score -= 2;
                push(@scoreArray, 1);
            }
            else
            {
                push(@scoreArray, 0);
            }
 
            my $cfam = 0;
            my $bfam = 0;
            my $efam = 0;
            for (my $i = 0; $i < @gis; $i++) {
                if ($PfamMatch{$gis[$i]} eq "PF00733") {
                    $cfam = 1;
                }
 
                if ($PfamMatch{$gis[$i]} eq "PF13471") {
                    $bfam = 1;
                }
 
                if ($PfamMatch{$gis[$i]} eq "PF05402") {
                    $efam = 1;
                }
            }
 
            push(@scoreArray, $cfam);
            push(@scoreArray, $efam);
            push(@scoreArray, $bfam);
 
            if ($bfam eq 0)
            {
                $score -= 2;
            }
 
            if ( $score >= $highScore )
            {
                push(@scores, $score);
                push(@precs, $pep);
            }
 
            if( $score > $highScore)
            {   
                $highScore = $score;
                $mostLikelyPrec = $pep;
            }
 
            $leaderCore =~ s/X/L/g;
            $leaderCore =~ s/Z/V/g;
            my $coreMass = calcMass($core);
 
            my @peebs = split(/,/, $leaderCore);
         
 
            if (defined($indexCSV) && $distanceOverall ne "x")
            {
		my $base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
		my $url = $base .= "efetch.fcgi?db=protein&id=$idGI&rettype=acc";
#		my $id_acc = get($url);
        my $ua = LWP::UserAgent->new;
        my $req = HTTP::Request->new(GET => $url);
        my $res = $ua->request($req);
        my $id_acc = $res->as_string;
        my @resArray = split("\n",$id_acc);
        $id_acc = @resArray[-1];

		chomp($id_acc);
		chomp($id_acc);
		
		$base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
		$url = $base .= "efetch.fcgi?db=protein&id=$mainGI&rettype=acc";
#		my $temp_acc = get($url);
        my $ua = LWP::UserAgent->new;
        my $req = HTTP::Request->new(GET => $url);
        my $res = $ua->request($req);
        my $temp_acc = $res->as_string;
        my @resArray = split("\n",$temp_acc);
        $temp_acc = @resArray[-1];


		chomp($temp_acc);
		chomp($temp_acc);

                open MEEJ, '>>', $ARGV[$indexCSV+1];
                print MEEJ "$temp_acc,$genus,$id_acc,$peebs[0],$peebs[1],$peptideCoordMin,$peptideCoordMax,$distanceOverall,$coreMass,$score,";
                for (my $p = 0; $p < @scoreArray-1; $p++) 
                {
                    print MEEJ "$scoreArray[$p],";
                }
                print MEEJ "$scoreArray[@scoreArray-1]\n";
                @scoreArray = ();
                close MEEJ;
            }
 
 
 
             
            if ($distanceOverall ne "x" && $bypassScore == 0)
            {
                $outputPrecsInfo[$ct][0] = $peptideCoordMin;
                $outputPrecsInfo[$ct][1] = $peptideCoordMax;
                $outputPrecsInfo[$ct][2] = $pepDir;
                $outputPrecsInfo[$ct][3] = generateRandomColor($score);
                $outputPrecsInfo[$ct][4] = $ct;
                $peptideHTML .= "<tr>\n"; 
                $peptideHTML .= "\t<td>$peebs[0]</td>
                <td>$peebs[1]</td>
                <td>$ct</td>
                <td>$peptideCoordMin</td>
                <td>$peptideCoordMax</td>
                <td>$distanceOverall</td> 
                <td>$pepDir</td>
                <td>$coreMass</td>
                <td>$score</td>";
                $peptideHTML .= "</tr>\n";
                $ct++;
            }
        }
 
     
 
        elsif ($bypassScore == 1 && $distanceOverall ne "x")
        {
            $pep =~ s/X/L/g;
            $pep =~ s/Z/V/g;
 
            if (defined($indexCSV))
            {

		my $base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
		my $url = $base .= "efetch.fcgi?db=protein&id=$idGI&rettype=acc";
#		my $id_acc = get($url);
        my $ua = LWP::UserAgent->new;
        my $req = HTTP::Request->new(GET => $url);
        my $res = $ua->request($req);
        my $id_acc = $res->as_string;
        my @resArray = split("\n",$id_acc);
        $id_acc = @resArray[-1];


		chomp($id_acc);
		chomp($id_acc);
		
		$base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
		$url = $base .= "efetch.fcgi?db=protein&id=$mainGI&rettype=acc";
#		my $temp_acc = get($url);
        my $ua = LWP::UserAgent->new;
        my $req = HTTP::Request->new(GET => $url);
        my $res = $ua->request($req);
        my $temp_acc = $res->as_string;
        my @resArray = split("\n",$temp_acc);
        $temp_acc = @resArray[-1];


		chomp($temp_acc);
		chomp($temp_acc);

                open MEEJ, '>>', $ARGV[$indexCSV+1];
                print MEEJ "$temp_acc,$species,$id_acc,$peptideCoordMin,$peptideCoordMax,$pepDir,$pep\n";
                close MEEJ;
            }
            $peptideHTML .= "<tr>\n"; 
            $peptideHTML .= "\t<td>$pep</td>
            <td>$peptideCoordMin</td>
            <td>$peptideCoordMax</td>
            <td>$pepDir</td>";
            $peptideHTML .= "</tr>\n";
        }
 
        $distanceOverall = 1000000;
         
    }
 
    close(INFO);
 
    $leaderCore = extractLeaderCore($mostLikelyPrec);
    $leaderCore =~ s/X/L/g;
    $leaderCore =~ s/Z/V/g;
 
    system("rm peptides.txt");
 
}
 
sub generateRandomColor
{
    my $score = $_[0];
    if ($score <= 4)
    {
        return "#FFFFFF";
    }
    if ($score == 5 || $score == 6)
    {
        return "#FFAAAA";
    }
    if ($score == 7 || $score == 8)
    {
        return "#FF5555";
    }
    if ($score >= 9)
    {
        return "#FF0000";
    }
}
 
sub calcMass
{
    my $a = $_[0];
    my @a = ();
    my $x = length $a;
    @a = split q{}, $a;
    my $b = 0;
    my %data = (
        A=>71.04,  R=>156.10,  D=>115.02,  N=>114.04,
        C=>103.01,  E=>129.04,  Q=>128.06,  G=>57.02,
        H=>137.06,  I=>113.08,  L=>113.08,  K=>128.09,
        M=>131.04,  F=>147.07,  P=>97.05,  S=>87.03,
        T=>101.05,  W=>186.08,  Y=>163.06,  V=>99.07,
        X=>113.08, Z=>99.07
 
    );
    for my $i( @a ) {
        $b += $data{$i};
    }
   
    return $b;
}
 
 
sub extractLeaderCore
{
    my $peptide = $_[0];
    my $core = $peptide =~ m/.*(T[A-Z]{7,10}(D|E)[A-Z]{5,20}$)/g;
    if ($1) 
    {
        $core = substr($1, 2);
    }
    my $coreIndex = index($peptide, $core);
    my $leader = substr($peptide, 0, $coreIndex);
 
    # if there's a methionine, start the leader there
    if ( index($leader, "M") != -1)
    {
        $leader = substr($leader, index($leader, "M") );
    }
 
    my $ret = $leader;
    $ret .= ", ";
    $ret .= $core;
    return $ret;
     
}
 
sub printFxns
{
    for(my $i = 0; $i < @gis; $i++)
    {
        chomp($gis[$i]);
        getFxn($gis[$i], $i);
    }
 
    if( $geneINFO{$gis[0]}[2] > $geneINFO{$gis[@gis-1]}[3] )
    {
        @gis = reverse(@gis);
        @fxns = ();
        for(my $i = 0; $i < @gis; $i++)
        {
            chomp($gis[$i]);
            getFxn($gis[$i], $i);
        }
    }
 
    @fxnsNoHyps = @fxns;
}
 
sub getFxn
{   
    my $gi = $_[0]; 
    my $cntr = $_[1];
    my $matchedFunction;
 
        $matchedFunction = $geneINFO{$gi}[0];
     
 
    toArray($gi);
 
     
    if(exists( $geneINFO{$gi} ) )
    {
        if($geneINFO{$gi}[1] eq "fwd")
        {
            $matchedFunction .= "> ";
        }
 
        elsif ($geneINFO{$gi}[1] eq "rev")
        {
            $matchedFunction .= "< ";
        }
    }
 
 
    push(@fxns, $matchedFunction);
 
}
 
 
sub fetchNeighboringGIS 
{


my @count_check = qx( (esearch -db protein -query $mainGI | elink -target nuccore) &> nul.out);



    open FOPEN, '<', "nul.out";
    foreach my $line (<FOPEN>)  
    {
        if (index($line, "not found") != -1 || index($line, "Empty result - nothing to do") != -1 ) 
        {
            print "QUERY NOT FOUND IN NCBI GENBANK\n\n";
            print FILE "<h2>$main_acc NOT FOUND IN NCBI GENBANK";

            if( $isList == 1)
            {
                print FILE '<a href="#header"><small><small>back to top</small></small></a></h2>' ;
                print "Done\n";
                goto CONTINUE;
            }
            else
            {
                print FILE '</h2>';
                exit;
            }
            system("rm nul.out");
            if ($isList == 1) 
            {
                goto CONTINUE;
            }
            else 
            {
                exit;
            }

        }
    }
    close FOPEN;

    my $prot_count = 0;

    foreach my $line (@count_check) {
        if ( index($line, "Count") != -1) {
            $line =~ /([0-9]+)/;
            $prot_count = $1;
            last;
        }
    }
 
    if ($prot_count > 100) {
        print FILE "<h2> RETURNED GENBANK FILE FOR $mainGI TOO LARGE TO EVALUATE";
        if( $isList == 1)
        {
            print FILE '<a href="#header"><small><small>back to top</small></small></a></h2>' ;
            print "Done\n";
            goto CONTINUE;
        }
        else
        {
            print FILE '</h2>';
            exit;
        }
    }
 
    #my $val = qx(esearch -db protein -query $mainGI | elink -target nuccore | efetch | xtract);
    my $val = qx(esearch -db protein -query $mainGI | elink -target nuccore | efetch |  tr -s [:space:] ' ');
    if ( index($val, "Seq-entry") == -1 )
    {
        print FILE "<h2>$mainGI NOT FOUND IN NCBI GENBANK";
        if( $isList == 1)
        {
            print FILE '<a href="#header"><small><small>back to top</small></small></a></h2>' ;
            print "Done\n";
            goto CONTINUE;
        }
        else
        {
            print FILE '</h2>';
            exit;
        }
    }
    $species = substr($val, index($val, 'taxname')+9);
    $species = substr($species, 0, index($species, '"'));
    $genus = substr($species, 0, index($species, ' '));
    #my @orfs = $val =~ /(product whole gi [0-9]* , location int { from [0-9]* , to [0-9]* .*?id gi [0-9]*)/g;



    #my @orfs = $val =~ /(product whole gi [0-9]* , location\s?[a-z]*?\s?[{]?\s? int { from [0-9]* , to [0-9]* .*?id gi [0-9]*)/g;
    my @orgs = $val =~ /(mRNA , ext name "[a-z]+ [a-z]+" } , partial [A-z]+ , )?(product whole gi [0-9]* , location\s?[a-z]*?\s?[{]?\s? int { from [0-9]* , to [0-9]* .*?id gi [0-9]*)/g;

    my @orfs = ();
    for (my $q = 0; $q < @orgs; $q++) 
    {
        if( index($orgs[$q], "mRNA") == -1 && index($orgs[$q+1], "mRNA") == -1)
        {
            push(@orfs, $orgs[$q+1])
        }

    }

    @neighborGIs = sort {
      my ($aa) = $a =~ /(?:from ([0-9]*))/;
      my ($bb) = $b =~ /(?:from ([0-9]*))/;
      $aa <=> $bb;
    } @orfs;
 
    @neighborGIs = sort {
      my ($aa) = $a =~ /(?:id gi ([0-9]*))/;
      my ($bb) = $b =~ /(?:id gi ([0-9]*))/;
      $aa <=> $bb;
    } @neighborGIs;
 
    my $foundFirstID = 0;
    for(my $j = 0; $j < @neighborGIs; $j++)
    {
        my @guy = split / /, $neighborGIs[$j];
        if (index($neighborGIs[$j],$mainGI) != -1)
        {
            if ($foundFirstID == 0 && $idGI == 0)
            {
                $idGI = $guy[@guy-1];
                $foundFirstID = 1;
            }
            push(@scaffs, $guy[@guy-1]);
        }
    }

=pod
    if ( $idGI == 0)
    {
        print FILE "<h2> NO RESULTS FOR $mainGI";
        if( $isList == 1)
        {
            print FILE '<a href="#header"><small><small>back to top</small></small></a></h2>' ;
            goto CONTINUE;
        }
        else
        {
            print FILE '</h2>';
            exit;
        }
    }
=cut
 
    my $scaffsSize = @scaffs;
    if($scaffsSize > 1)
    {

        my $base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
        my $url = $base .= "efetch.fcgi?db=protein&id=$idGI&rettype=acc";

#        my $temp_acc = get($url);
        my $ua = LWP::UserAgent->new;
        my $req = HTTP::Request->new(GET => $url);
        my $res = $ua->request($req);
        my $temp_acc = $res->as_string;
        my @resArray = split("\n",$temp_acc);
        $temp_acc = @resArray[-1];



        if ($evalAll == 0) {
            $scaffHTML .= "<h4>MULTIPLE NUCLEOTIDE RECORDS FOUND, EVALUATING <a href='https://www.ncbi.nlm.nih.gov/protein/$temp_acc'>$temp_acc</a></h4><p><small>ALL FOUND ENTRIES:<br>";
            for(my $j = 0; $j < @scaffs; $j++)
            {
                $base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
                $url = $base .= "efetch.fcgi?db=protein&id=$scaffs[$j]&rettype=acc";

#                $temp_acc = get($url);
        my $ua = LWP::UserAgent->new;
        my $req = HTTP::Request->new(GET => $url);
        my $res = $ua->request($req);
        $temp_acc = $res->as_string;
        my @resArray = split("\n",$temp_acc);
        $temp_acc = @resArray[-1];


		chomp($temp_acc);
                $scaffHTML .= "<a href='https://www.ncbi.nlm.nih.gov/protein/$temp_acc'>$temp_acc</a><br>";
            }
            $scaffHTML .= "</small></p>";
 
            if (defined($indexA)) {
                $evalAll = 1;
            }
        }
        else {
            $scaffHTML .= "<h4>EVALUATING <a href='https://www.ncbi.nlm.nih.gov/protein/$temp_acc'>$temp_acc</a></h4>";
        }
    }
     
    my @neighboringEvals;
    my $count = 0;
    # find GIs that have a matching idGI and extract all of them
    foreach my $line (@neighborGIs) {
        if( index($line, $idGI) != -1 )
        {
            $line =~ /(\d+)/;
            push(@neighboringEvals, $1);
            push(@neighborINFO, $line);
            $neighboringEvals[$count] .= "\n";
            $count++;
        }
    }
 
    my $indexOfInput = 0;
    # now find input gi and get it's neighbors
    for (my $i=0; $i < @neighboringEvals; $i++)
    {
        if ($neighboringEvals[$i] == $mainGI)
        {
            $indexOfInput = $i;
            last; # end loop once we've found index
        }
    }
    # now narrow down neighbors to +/- 6
    my @conciseEvals;
    my $start = $indexOfInput - $numNeighbors;
    my $end = $indexOfInput + $numNeighbors;
    my $length = @neighboringEvals;
 
    if ($start < 0)
    {
        $start = 0;
    }
 
    if ($end > $length-1)
    {
        $end = $length-1;
    }
 
    for(my $j = 0; $start <= $end; $start++)
    {
        $conciseEvals[$j] = $neighboringEvals[$start];
        $j++;
    }

    foreach my $line (@conciseEvals) 
    {
        chomp($line);
        my $gi = $line;

        my $base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
        my $url = $base .= "efetch.fcgi?db=protein&id=$gi&rettype=acc";

#        my $output = get($url);
        my $ua = LWP::UserAgent->new;
        my $req = HTTP::Request->new(GET => $url);
        my $res = $ua->request($req);
        my $output = $res->as_string;
        my @resArray = split("\n",$output);
        $output = @resArray[-1];


        chomp($output);
        push(@accessions, $output);
        #print "$output\n";
    }

    return @conciseEvals;
}
 
 
sub toArray
{
 
    my $gi = $_[0];
    my $beg;
    my $end;
    my $dir = "fwd"; 
 
    # find our gi, then extract relevant info from it
    foreach my $line (@neighborINFO)             # iterate through formatted list of GIs
    {
        my $gid = $line =~ m/($gi)/;            # if GI number of line matches desired
        if ($gid == 1 && !exists($geneINFO{$gi}[2]) ) 
        {
            my @result = $line =~ /(\d+)/g;        # create an array of three numbers
             
            $beg = $result[1] + 1;        # set beginning and end of GENE
            $end = $result[2] + 1;
 
            my $strand = $line =~ m/(minus)/; 
 
            if($strand == 1)
            {
                $dir = "rev";
            }
            if( !exists( $geneINFO{$gi}[1] ) )
            {
                $geneINFO{$gi}[1] = $dir;
                $geneINFO{$gi}[2] = $beg;
                $geneINFO{$gi}[3] = $end;
            }
        }
    }
 
    close(INFO);
} 
 
sub hmmer
{
    my $cPFamID = $_[1];
    my $ePFamID = $_[2];
    my $bPFamID = $_[3];
    my $g1PFamID = $_[4];  
    my $g2PFamID = $_[5];
    my $IDtoCompare = "";
    my $pfamName = "";
    my $inputGI = $_[0];
 
    # first get fasta for given gi 
 
    system( "esearch -db protein -query $inputGI | efetch -format fasta > temp.fasta");

    #my $tmp = File::Temp->new();
    #$tmp = qx(esearch -db protein -query $inputGI | efetch -format $tmp);


    # now use hmmscan 
    system( "hmmscan -o trash.txt --noali --domtblout pfamStuff.tab Pfam-A.hmm temp.fasta");
    #qx(hmmscan -o trash.txt --noali --domtblout pfamStuff.tab Pfam-A.hmm $tmp);
 
    # grep to get the best matched PFam ID
    system( 'grep -o "PF[0-9]*" pfamStuff.tab > PFamID.txt');
    system( 'grep -w ">>" trash.txt > PFamName.txt');
 
    # extract the Evals
    my $flag = 0;
    open FOPEN, '<', "trash.txt";
 
    foreach my $line (<FOPEN>)
    {
        if( index($line, "inclusion threshold") != -1 )
        {
            next;
        }
 
        if( index($line, "E-value") != -1 )
        {
            $flag = 1;
            next;
        }
        if ($flag == 1) 
        {
            $flag = 2;
            next;
        }
        if ($flag == 2)
        {
            $flag = 3;
            #$line =~ /([0-9][.]?[0-9e-]+)/;
            if( $line =~ /([0-9][.]?[0-9e-]+)/ )
            {
                $PfamEval{$inputGI} = $1;
            }
            next;
        }
        if ($flag == 3)
        {
            $flag = 4;
            #$line =~ /([0-9][.]?[0-9e-]+)/;        
            if( $line =~ /([0-9][.]?[0-9e-]+)/ )
            {
                $secondPfamEval{$inputGI} = $1;
            }
            next;
        }
        if ($flag == 4) 
        {
            #$line =~ /([0-9][.]?[0-9e-]+)/;        
            if( $line =~ /([0-9][.]?[0-9e-]+)/ )
            {
                $thirdPfamEval{$inputGI} = $1;
            }
            last;
        }
    } 
 
    close FOPEN;
 
    my $gotFirstPFam = 0;
    open (INFO, "<PFamID.txt");
 
    foreach my $line (<INFO>)
    {
        chomp($line);
 
        my $yoyo = $line;
        $yoyo .= $inputGI;
         
        push(@allFoundPfams, $yoyo);
 
        if ($gotFirstPFam == 2 && $line ne $IDtoCompare && $line ne "PF"  && $line ne $secondPfamMatch{$inputGI})
        {
            $thirdPfamMatch{$inputGI} = $line;
            $gotFirstPFam = 3;
        }

        if ($gotFirstPFam == 1 && $line ne $IDtoCompare && $line ne "PF" )
        {
            $secondPfamMatch{$inputGI} = $line;
            $gotFirstPFam = 2;
        }
 
        if ($gotFirstPFam == 0)
        {
            $IDtoCompare = $line;
            $gotFirstPFam = 1;
        }
 
    }
 
    close(INFO);
 
    my $gotFirstName = 0;
 
    open (INFO, "<PFamName.txt") ;
 
    foreach my $line (<INFO>)
    {
        chomp($line);
        $line = substr($line, 3);
        if($gotFirstName == 1)
        {
            my $pfamName2 = $line;
            #$pfamName2 = substr($pfamName2, 3);
            $secondPfamName{$inputGI} = $pfamName2;
            $gotFirstName = 2;
        }
        if($gotFirstName == 2 && $line ne $secondPfamName{$inputGI})
        {
            my $pfamName3 = $line;
            #$pfamName3 = substr($pfamName3, 3);
            $thirdPfamName{$inputGI} = $pfamName3;
            $gotFirstName = 3;
        }


        if($gotFirstName == 0)
        {
            $pfamName = $line;
            #$pfamName = substr($pfamName, 3);
            $gotFirstName = 1;
        }
    }
 
    my $currEval = "1e0";
    my $currEval2 = "1e0";
 
    my $hmmString = "";
    my $bestHMM = "_none_";
    my %hmm_list_for_sorting;
    if ( defined($indexP) ) 
    {
        foreach my $custHmm(@custHmms) 
        {
            chomp($custHmm);
            if (evalCustHmm($custHmm, $inputGI, 0) ne "") {
                my $custHmmHit = evalCustHmm($custHmm, $inputGI, 0);
                my $custHmmEval = getCustomEval();
                $hmm_list_for_sorting{$custHmmHit} = $custHmmEval;
                if ($custHmmEval < $currEval) {
                    $customPfams{$inputGI} = $custHmmHit; 
                    $customPfamEvals{$inputGI} = $custHmmEval;  
                    $currEval = $customPfamEvals{$inputGI};
                }
            }
        }
         
        foreach my $name (sort { $hmm_list_for_sorting{$b} <=> $hmm_list_for_sorting{$a} } keys %hmm_list_for_sorting) {
            $hmmString .= "$name,$hmm_list_for_sorting{$name},";
            $bestHMM = $name;
        }
        $custHmmMatch{$inputGI} = $hmmString;
    }
 
    close(INFO);
 
    $PfamName{$inputGI} = $pfamName;
     
    if( defined($userPFamInfo{$IDtoCompare}) && !defined($geneINFO{$inputGI}[0]) )
    {
        $geneINFO{$inputGI}[0] = $userPFamInfo{$IDtoCompare}[0];
        $giColors{$inputGI} = $userPFamInfo{$IDtoCompare}[1];
        $PfamMatch{$inputGI} = $IDtoCompare;
        return;
    }
 
     
    if($IDtoCompare eq $cPFamID)
    {
        $cFromHMMER = 1;
        if ( !exists($geneINFO{$inputGI}[0]) && $bypassScore == 0 ) {
            $geneINFO{$inputGI}[0] = "C";
        }
        $PfamMatch{$inputGI} = $IDtoCompare;
    }
 
    elsif($IDtoCompare eq $ePFamID)
    {
        if ( !exists($geneINFO{$inputGI}[0]) && $bypassScore == 0) {
            $geneINFO{$inputGI}[0] = "E";
        }
        $PfamMatch{$inputGI} = $IDtoCompare;
    }
 
    elsif($IDtoCompare eq $g1PFamID)
    {
        if ( !exists($geneINFO{$inputGI}[0]) && $bypassScore ==0 ) {
            $geneINFO{$inputGI}[0] = "G1";
        }
        $PfamMatch{$inputGI} = $IDtoCompare;
    }
 
    elsif($IDtoCompare eq $g2PFamID)
    {
        if ( !exists($geneINFO{$inputGI}[0]) && $bypassScore == 0) {
            $geneINFO{$inputGI}[0] = "G2";
        }
        $PfamMatch{$inputGI} = $IDtoCompare;
    }
 
    elsif($IDtoCompare eq $bPFamID)
    {
        if ( !exists($geneINFO{$inputGI}[0]) && $bypassScore == 0) {
            $geneINFO{$inputGI}[0] = "B";
        }
        $bFromHMMER = 1;
        $PfamMatch{$inputGI} = $IDtoCompare;
    }
 
    else
    {
        $geneINFO{$inputGI}[0] = "";
        $PfamMatch{$inputGI} = $IDtoCompare;
        $giColors{$inputGI} = "white";
         
     if ( $bestHMM ne "_none_"     ) {
            $bestHMM =~ s/ +/ /g;
             
            if ( defined($userHMMInfo{$bestHMM})) {
                $geneINFO{$inputGI}[0] = $userHMMInfo{$bestHMM}[0];
             
             
                $giColors{$inputGI} = $userHMMInfo{$bestHMM}[1];
            } 
            }
         
    }
     
    if( index($PfamMatch{$inputGI}, "P") == -1 )
    {
        $PfamMatch{$inputGI} = "NO PFAM MATCH";
        $PfamName{$inputGI} = "NO PFAM MATCH";
    }
     
    system("rm trash.txt temp.fasta pfamStuff.tab");
}
 
sub getCustomEval {
    my $flag = 0;
    open FOPEN, '<', "custPfam.txt";
 
    foreach my $line (<FOPEN>)
    {
        if( index($line, "E-value") != -1 )
        {
            $flag = 1;
            next;
        }
        if ($flag == 1) 
        {
            $flag = 2;
            next;
        }
        if ($flag == 2)
        {
            $flag = 3;
            #$line =~ /([0-9][.]?[0-9e-]*)/; # use regex below for superior performance
            $line =~ /\s+(\S+)\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+\S+/g;
            #if( $1 )
            #{
                system("rm custPfam.txt custPFamName.txt");
                return $1;
            #}
            next;
        }
} 
    close FOPEN;
}
 
 
sub getCustomHmms {
    my @files;
    opendir(DIR,$hmmdir) or die "SPECIFIED HMM DIRECTORY NOT FOUND";
    while( (my $filename = readdir(DIR))){
        if ( index($filename,"Pfam-A.hmm") == -1 &&index(lc($filename),".hmm") != -1 && index($filename, ".h3") == -1) {
            push(@files, $filename);
        }
    }
    closedir(DIR);
 
    return @files;
}
 
 
sub evalCustHmm {
    my $cycName;
    my $hmm = $_[0];
    my $inputGI = $_[1];
    my $flag = 0;
 
    system( "hmmscan -o custPfam.txt --noali --domtblout tab.tab $hmmdir/$hmm temp.fasta");
    system( "sed 's/,/ /g' custPfam.txt > custPfamC.txt && mv custPfamC.txt custPfam.txt");
    system( 'grep -w ">>" custPfam.txt > custPFamName.txt');
 
    open (INFO, "<custPFamName.txt") ;
 
    foreach my $line (<INFO>)
    {
        chomp($line);
        if( index($line, ">>") != -1)
        {
            if ( defined($userHMMInfo{$hmm}) ) {
                $geneINFO{$inputGI}[0] = $userHMMInfo{$hmm}[0];
                #$giColors{$inputGI} = $userHMMInfo{$hmm}[1];
            }
            $cycName = $line;
            close(INFO);
            $cycName = substr($cycName, 3);
            $cycName =~ s/^\s+|\s+$//g ;
            if ($_[2] == 0) {
                return $cycName;
            }
            elsif ($flag == 0){
                $flag = 1;
                next;
            }
            if($flag == 1) {
                return $cycName;
            }
 
        }
        else
        {
            close(INFO);
            return "";
        }
    }
}
 
sub toHTML
{


    my $base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/';
    my $url = $base .= "efetch.fcgi?db=protein&id=$mainGI&rettype=acc";

#    $main_acc = get($url);
        my $ua = LWP::UserAgent->new;
        my $req = HTTP::Request->new(GET => $url);
        my $res = $ua->request($req);
        $main_acc = $res->as_string;
        my @resArray = split("\n",$main_acc);
        $main_acc = @resArray[-1];


    chomp($main_acc);
    chomp($main_acc);
print FILE <<ENDHTML;
<h2 id="$main_acc">results for $main_acc [$species] 
ENDHTML
if(defined($indexL))
{
    print FILE '<a href="#header"><small><small>back to top</small></small></a>';
}
print FILE <<ENDHTML;
</h2>
<p></p>
$scaffHTML
<h3>architecture</h3>
$architectureHTML
<br><br>
ENDHTML
# check to see if whole genome is available
my $url = "https://www.ncbi.nlm.nih.gov/nuccore/$idGI";
        my $ua = LWP::UserAgent->new;
        my $req = HTTP::Request->new(GET => $url);
        my $res = $ua->request($req);
        my $content = $res->as_string;

#my $content = get($url);
my $whole_gen = 0;

if ( index($content, "whole genome") != -1) 
{
    $whole_gen = 1;
    print FILE <<ENDHTML;
    <a href='https://www.ncbi.nlm.nih.gov/nuccore/$idGI'>Link to nucleotide sequence</a>
ENDHTML
}

else 
{
        print FILE <<ENDHTML;
    <a href='https://www.ncbi.nlm.nih.gov/nuccore/$idGI'>Link to nucleotide sequence</a>
ENDHTML
}

print FILE <<ENDHTML;
<br><br><table class="table table-condensed">
  <tbody>
    <tr>
      <th scope="col">Accession</th>
      <th scope="col">min</th>
      <th scope="col">max</th>
      <th scope="col">direction</th>
      <th scope="col">length (aa)</th>
ENDHTML
 
print FILE '<th scope="col">Pfam/HMM</th>';
print FILE '<th scope="col">E-value</th>';     
 
print FILE <<ENDHTML;      
      <th scope="col">description</th>
    </tr>
ENDHTML
for(my $giCount = 0; $giCount < @gis; $giCount++)
{
    my $aaLen = floor(($geneINFO{$gis[$giCount]}[3]-$geneINFO{$gis[$giCount]}[2])/3);
    print FILE ("<tr>\n");    
    if ($PfamMatch{$gis[$giCount]} eq "NO PFAM MATCH")
    {

            print FILE ("\t<td><a href='https://www.ncbi.nlm.nih.gov/protein/$accessions[$giCount]'>$accessions[$giCount]</a></td>
            <td>$geneINFO{$gis[$giCount]}[2]</td> 
            <td>$geneINFO{$gis[$giCount]}[3]</td>
            <td>$geneINFO{$gis[$giCount]}[1]</td>
            <td>$aaLen</td>
            <td>$PfamMatch{$gis[$giCount]}");
    }
    else 
    {
 
            print FILE ("\t<td><a href='https://www.ncbi.nlm.nih.gov/protein/$accessions[$giCount]'>$accessions[$giCount]</a></td>
            <td>$geneINFO{$gis[$giCount]}[2]</td> 
            <td>$geneINFO{$gis[$giCount]}[3]</td>
            <td>$geneINFO{$gis[$giCount]}[1]</td>
            <td>$aaLen</td>
           <td><a href='http://pfam.xfam.org/family/$PfamMatch{$gis[$giCount]}'>$PfamMatch{$gis[$giCount]}</a>");
        
    }
 
    if ( exists($secondPfamMatch{$gis[$giCount]}) )
    {
        print FILE "<br><a href='http://pfam.xfam.org/family/$secondPfamMatch{$gis[$giCount]}'>$secondPfamMatch{$gis[$giCount]}</a>";
        #print "$gis[$giCount] : 2 - $secondPfamMatch{$gis[$giCount]}\n";
        if ( exists($thirdPfamMatch{$gis[$giCount]}) ) {
            #print "$gis[$giCount] : 3 - $thirdPfamMatch{$gis[$giCount]}\n";
            print FILE "<br><a href='http://pfam.xfam.org/family/$thirdPfamMatch{$gis[$giCount]}'>$thirdPfamMatch{$gis[$giCount]}</a>";
        }
    }
 
    my @hmmEvals;
    my @hmmMatches;
    if( exists($custHmmMatch{$gis[$giCount]})) {
        my @hmmMatchAndEval = split/,/,$custHmmMatch{$gis[$giCount]};
        #@hmmMatchAndEval = sort { $a <=> $b } @hmmMatchAndEval;
        for(my $i = 1; $i < @hmmMatchAndEval; $i+=2) {
            push(@hmmEvals,$hmmMatchAndEval[$i]);
        }
        @hmmEvals = sort { $a <=> $b } @hmmEvals;
        for (my $i = 0; $i < @hmmEvals; $i++) {
            for(my $j = 1; $j < @hmmMatchAndEval; $j+=2) {
                if ($hmmMatchAndEval[$j] == $hmmEvals[$i]) {
                    print FILE "<br>$hmmMatchAndEval[$j-1]";
                }
            }
        }
 
    }
 
    if( index($PfamMatch{$gis[$giCount]}, "NO") == -1 && exists($PfamEval{$gis[$giCount]}))
    {
        print FILE ("</td><td>$PfamEval{$gis[$giCount]}");  
    } 
 
    if ( exists($secondPfamMatch{$gis[$giCount]}) && defined($secondPfamMatch{$gis[$giCount]}) )
    {
        print FILE "<br>$secondPfamEval{$gis[$giCount]}";
        if ( exists($thirdPfamMatch{$gis[$giCount]}) && defined($thirdPfamMatch{$gis[$giCount]}) ) {
             print FILE "<br>$thirdPfamEval{$gis[$giCount]}";
        }
    }
 
    if ( exists($customPfams{$gis[$giCount]}) )
    {
        if (index($PfamMatch{$gis[$giCount]}, "NO") == -1) {
            for (my $i = 0; $i < @hmmEvals; $i++) {
                print FILE "<br>$hmmEvals[$i]";
            }
            #print FILE "<br>$customPfamEvals{$gis[$giCount]}";
        }
        if (index($PfamMatch{$gis[$giCount]}, "NO") != -1) {
            print FILE "</td><td>";
            for (my $i = 0; $i < @hmmEvals; $i++) {
                print FILE "<br>$hmmEvals[$i]";
            }
        }
 
    }
 
 
    if ($PfamName{$gis[$giCount]} eq "NO PFAM MATCH") {
        if ( exists($customPfams{$gis[$giCount]}) ){
            print FILE"</td><td>$PfamName{$gis[$giCount]}</td>";
        }
        else {
            print FILE"</td><td></td><td>$PfamName{$gis[$giCount]}</td>";         
        }
 
    }
    else {
        print FILE ("</td><td>$PfamName{$gis[$giCount]}");
    }
    if ( exists($secondPfamName{$gis[$giCount]}) )
    {
        print FILE "<br>$secondPfamName{$gis[$giCount]}";
        if ( exists($thirdPfamName{$gis[$giCount]}) ) {
            print FILE "<br>$thirdPfamName{$gis[$giCount]}";
        }
 
    }
    print FILE ("</td>");
    print FILE ("</tr>\n");
}
print FILE <<ENDHTML;  
  </tbody>
</table>
<p></p>
ENDHTML
 
if( $printPrecs != -1 && $bypassScore == 0)
{   
print FILE <<ENDHTML;
<table class="table table-bordered">
  <tbody>
    <tr>
      <th scope="col">leader</th>
      <th scope="col">core</th>
      <th scope="col">id</th>
      <th scope="col">min</th>
      <th scope="col">max</th>
      <th scope="col">dist.</th>
      <th scope="col">dir</th>
      <th scope="col">core mass</th>
      <th scope="col">score</th>
    </tr>
    $peptideHTML
    </tbody>
</table>
ENDHTML
}
 
elsif( $bypassScore == 1 && $printPrecs != -1)
{
    print FILE <<ENDHTML;
<table class="table table-bordered">
  <tbody>
    <tr>
      <th scope="col">peptide</th>
      <th scope="col">min coord</th>
      <th scope="col">max coord</th>
      <th scope="col">dir</th>
    </tr>
    $peptideHTML
    </tbody>
</table>
ENDHTML
}
 
else
{
    print FILE <<ENDHTML; 
ENDHTML
}
 
}
