#!/usr/bin/perl

# find all possible peptides in out.txt

use strict;
use warnings;

my $file = 'out.txt';
open INFO, $file or die;
my $frameNum = 1; 
my $count = 0;
my @pepIndices;
my @coords;
my $fwdOrRev = "f";
my $min = $ARGV[0];
my $max = $ARGV[1];
my $bypassScore = $ARGV[2];
my $stuff;
my @peptides;

if ($bypassScore == 0)
{
	foreach my $line (<INFO>)  {   
	   my $string= $line;
	   @peptides = $line =~/((M|Z|X)[^!\n]+)/sg;  #Look for an M, followed by 1 or more of any character except an exclamation point.
	   if( $1 && $peptides[0]) 
	   {
		   	while ($peptides[$count])
		   	{
	   		
		   		# these next few lines truncate the match to the smallest possible peptide that starts with M 
		   		my $stuff = $peptides[$count];
		   		while ( index($stuff, "M") != -1 && length(substr($stuff, index($stuff, "M"))) >= 30 )
				{
					$peptides[$count] = substr($stuff, index($stuff, "M") );
					$stuff = substr($stuff, index($stuff, "M")+1 );
				}
							
		   		my $pepLength = length($peptides[$count]);

		   		# only want to process peptide if the length is within
		   		# the boundaries 
		   		if ($pepLength >= $min && $pepLength <= $max )
		   		{
		   			# calculate coordinates and assign forward or reverse
			   		$pepIndices[$count] = index($line, $peptides[$count]);
			   		if ($frameNum > 3)
			   		{
			   			$frameNum = $frameNum - 3;
			   			$fwdOrRev = "r";
			   		}
				   	if($fwdOrRev eq "f")
				   	{
				   		$coords[$count] = ($pepIndices[$count] * 3) + ($frameNum - 1);
				   	}
				   	else 
				   	{
				   		$coords[$count] = (length($line)*3) - ($pepIndices[$count] * 3) - ($frameNum - 1);
				   	}

				   	# looking for a T follwed by 7-10 of anything ending with a D or E
				   	my $match = $peptides[$count] =~ m/T[A-Z]{7,10}(D|E)[A-Z]{5,20}$/g;
				 
				   	# if we found the match, add it to the file in the specified format 
				   	if($match)
				   	{	
					   			my $deebs = length($line);
					   			print "$peptides[$count],$coords[$count],$fwdOrRev\n";
				   		
				   	}
				}

			   	$count++;
		   	}
		   	$count = 0;
	   }
	   if($line =~ /^$/) 
	  	{
	  		$frameNum++;
	  	}
	}
}
else
{
foreach my $line (<INFO>)  {   
	   my $string= $line;
	   @peptides = $line=~/((M|Z|X)[^!\n]+)/sg;  #Look for an M, followed by 1 or more of any character except an exclamation point.
	   if( $1 && $peptides[0]) 
	   {
		   	while ($peptides[$count])
		   	{
	   		
		   		# these next few lines truncate the match to the smallest possible peptide that starts with M 
		   		my $stuff = $peptides[$count];
		   		while ( index($stuff, "M") != -1 && length(substr($stuff, index($stuff, "M"))) >= 30 )
				{
					$peptides[$count] = substr($stuff, index($stuff, "M") );
					$stuff = substr($stuff, index($stuff, "M")+1 );
				}
							
		   		my $pepLength = length($peptides[$count]);

		   		# only want to process peptide if the length is within
		   		# the boundaries 
		   		if ($pepLength >= $min && $pepLength <= $max )
		   		{
		   			# calculate coordinates and assign forward or reverse
			   		$pepIndices[$count] = index($line, $peptides[$count]);
			   		if ($frameNum > 3)
			   		{
			   			$frameNum = $frameNum - 3;
			   			$fwdOrRev = "r";
			   		}
				   	if($fwdOrRev eq "f")
				   	{
				   		$coords[$count] = ($pepIndices[$count] * 3) + ($frameNum - 1);
				   	}
				   	else 
				   	{
				   		$coords[$count] = (length($line)*3) - ($pepIndices[$count] * 3) - ($frameNum - 1);
				   	}

		   			print "$peptides[$count],$coords[$count],$fwdOrRev\n";

				}

			   	$count++;
		   	}
		   	$count = 0;
	   }
	   if($line =~ /^$/) 
	  	{
	  		$frameNum++;
	  	}
	}	
}


close INFO;

system("rm out.txt");
