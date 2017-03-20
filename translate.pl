#!/usr/bin/perl

$filename = $ARGV[0];

# Remove the newline from the filename
chomp $filename;

# SUBROUTINE 1
# Get the contents of the file
(@fasta_file_data) = get_file_data($filename);

# SUBROUTINE 2
# Substract the clean DNA sequence
($sequence) = extract_sequence_from_fasta_data(@fasta_file_data);


# Translate the DNA to protein in six reading frames
#   and print the protein in lines 70 characters long
$protein = translate($sequence);
print_sequence($protein, 70000);

$sequence_frame = frame($sequence, 2);
$protein = translate($sequence_frame);
print_sequence($protein, 70000);

$sequence_frame = frame($sequence, 3);
$protein = translate($sequence_frame);
print_sequence($protein, 70000);

# Calculate reverse complement
$revcom = complement($sequence);

$protein = translate($revcom);
print_sequence($protein, 70000);

$sequence_frame = frame($revcom, 2);
$protein = translate($sequence_frame);
print_sequence($protein, 70000);

$sequence_frame = frame($revcom, 3);
$protein = translate($sequence_frame);
print_sequence($protein, 70000);


exit;


################################################################################
# Subroutines for this example
################################################################################

# SUBROUTINE 1
# get_file_data
#
# A subroutine to get data from a file given its filename

sub get_file_data {

    my($filename) = @_;

    # Initialize variables
    my @fasta_file_data = (  );

    unless( open(FILEDATA, $filename) ) {
        exit;
    }

    @fasta_file_data = <FILEDATA>;

    close FILEDATA;
    my @stuff;

    for(my $i = 0; $i < @fasta_file_data; $i++)
    {
        if($fasta_file_data[$i] ne "\n")
        {
            push(@stuff, $fasta_file_data[$i]);
        }
        
        if($fasta_file_data[$i] eq "\n")
        {
            last;
        }

    }

    return @stuff;

}


# SUBROUTINE 2
# extract_sequence_from_fasta_data
#
# A subroutine to extract FASTA sequence data from an array

sub extract_sequence_from_fasta_data {

    my(@fasta_file_data) = @_;

    # Declare and initialize variables
    my $sequence = '';

    foreach my $line (@fasta_file_data) {

        # discard fasta header line
        if($line =~ /^>/) {
            next;

        # keep line, add to sequence string
        } else {
            $sequence .= $line;
        }
    }

    # remove non-sequence data (in this case, whitespace) from $sequence string
    $sequence =~ s/\s//g;

    return $sequence;
}


# SUBROUTINE 3
sub translate {
    my($dna) = @_;

    $dna = uc $dna;
	my $protein = '';

    my(%genetic_code) = (
    
    'TCA' => 'S',    # Serine
    'TCC' => 'S',    # Serine
    'TCG' => 'S',    # Serine
    'TCT' => 'S',    # Serine
    'TTC' => 'F',    # Phenylalanine
    'TTT' => 'F',    # Phenylalanine
    'TTA' => 'L',    # Leucine
    'TTG' => 'X',    # Leucine
    'TAC' => 'Y',    # Tyrosine
    'TAT' => 'Y',    # Tyrosine
    'TAA' => '!',    # Stop
    'TAG' => '!',    # Stop
    'TGC' => 'C',    # Cysteine
    'TGT' => 'C',    # Cysteine
    'TGA' => '!',    # Stop
    'TGG' => 'W',    # Tryptophan
    'CTA' => 'L',    # Leucine
    'CTC' => 'L',    # Leucine
    'CTG' => 'L',    # Leucine
    'CTT' => 'L',    # Leucine
    'CCA' => 'P',    # Proline
    'CCC' => 'P',    # Proline
    'CCG' => 'P',    # Proline
    'CCT' => 'P',    # Proline
    'CAC' => 'H',    # Histidine
    'CAT' => 'H',    # Histidine
    'CAA' => 'Q',    # Glutamine
    'CAG' => 'Q',    # Glutamine
    'CGA' => 'R',    # Arginine
    'CGC' => 'R',    # Arginine
    'CGG' => 'R',    # Arginine
    'CGT' => 'R',    # Arginine
    'ATA' => 'I',    # Isoleucine
    'ATC' => 'I',    # Isoleucine
    'ATT' => 'I',    # Isoleucine
    'ATG' => 'M',    # Methionine
    'ACA' => 'T',    # Threonine
    'ACC' => 'T',    # Threonine
    'ACG' => 'T',    # Threonine
    'ACT' => 'T',    # Threonine
    'AAC' => 'N',    # Asparagine
    'AAT' => 'N',    # Asparagine
    'AAA' => 'K',    # Lysine
    'AAG' => 'K',    # Lysine
    'AGC' => 'S',    # Serine
    'AGT' => 'S',    # Serine
    'AGA' => 'R',    # Arginine
    'AGG' => 'R',    # Arginine
    'GTA' => 'V',    # Valine
    'GTC' => 'V',    # Valine
    'GTG' => 'Z',    # Valine
    'GTT' => 'V',    # Valine
    'GCA' => 'A',    # Alanine
    'GCC' => 'A',    # Alanine
    'GCG' => 'A',    # Alanine
    'GCT' => 'A',    # Alanine
    'GAC' => 'D',    # Aspartic Acid
    'GAT' => 'D',    # Aspartic Acid
    'GAA' => 'E',    # Glutamic Acid
    'GAG' => 'E',    # Glutamic Acid
    'GGA' => 'G',    # Glycine
    'GGC' => 'G',    # Glycine
    'GGG' => 'G',    # Glycine
    'GGT' => 'G',    # Glycine
    );

   
	# Translate each three-base codon into an amino acid, and append to a protein 
	for(my $i=0; $i < (length($dna) - 2) ; $i += 3) 
		{
		$codon = substr($dna,$i,3);

		if(exists $genetic_code{$codon}) 
			{
			$aminoacid = $genetic_code{$codon};
			$protein .= $aminoacid;
			}
		else
			{
            $protein .= "X";
			#print "Bad codon \"$codon\"!!\n";
			#exit;
			}
		}
	
	return($protein);
}


# SUBROUTINE 4
# print_sequence
#
# A subroutine to format and print sequence data 

sub print_sequence {

    my($sequence, $X) = @_;

	my $line;

    # Print sequence in lines of $length
    for ( my $pos = 0 ; $pos < length($sequence) ; $pos += $X ) {
		$line = substr($sequence, $pos, $X);
        #print $line . "\n";
        system("echo '$line\n' >> out.txt")
    }
}


sub complement {

	my($DNA) = @_;

	# Make a new copy of the DNA (see why we saved the original?)
	$revcom = reverse $DNA;

	# See the text for a discussion of tr///
	$revcom =~ tr/ACGTacgt/TGCAtgca/;

	return($revcom);

}


# frame
#
# A subroutine to translate a frame of DNA

sub frame {

    my($sequence, $start) = @_;

	# Subtract $sequence from ($start - 1) to the end of the sequence
	# $start - 1 because in Perl the first posicion is 0, not 1
    my $sequence_frame = substr($sequence, $start - 1);

	return $sequence_frame;

}
