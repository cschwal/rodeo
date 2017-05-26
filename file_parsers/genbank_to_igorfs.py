#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed May 24 21:53:38 2017

@author: bryce
"""

import argparse
from Bio import Entrez
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import SeqIO

class My_record(object):
    """Contains member functions and variables for handline *.gbk files"""
    def __init__(self, record, min_aa_length, max_aa_length):
        """record is the output of Bio.SeqIO.parse()"""
        self.id = record.id
        self.record = record
        self.sequence = record.seq
        self.CDSs = []
        self.intergenic_seqs = []
        self.intergenic_orfs = []
        self.start_codons = ['ATG','GTG', 'TTG']
        self.stop_codons = ['TAA','TAG','TGA']
        self.set_CDSs()
        self.set_intergenic_seqs()
        self.min_aa_length = min_aa_length
        self.set_intergenic_orfs(min_aa_length, max_aa_length)
            
        
    class Sub_seq(object):
        """Useful for storing subsequences and their coordinates"""
        def __init__(self, seq, start, end):
            self.start = start
            self.end = end
            self.sequence = seq
            
                
    def set_CDSs(self):
        """Get CDSs from record"""
        for feature in self.record.features:
            if feature.type == 'CDS':
                start = int(feature.location.start)
                end = int(feature.location.end)
                cds = self.Sub_seq(seq = "", start=start, end=end) #Don't need CDS seq...
                self.CDSs.append(cds)
       
    def set_intergenic_seqs(self):
        start = 0
        for cds in self.CDSs:
            end = cds.start
            nt_seq = self.sequence[start:end]
            intergenic_sequence = self.Sub_seq(seq=nt_seq, 
                                               start = start,
                                               end=end)
            self.intergenic_seqs.append(intergenic_sequence)
            start = cds.end
    
    def set_intergenic_orfs(self, min_aa_seq_length, max_aa_seq_length=1000000):
        for intergenic_seq in self.intergenic_seqs:
            #OFFSET will adjust for the main sequence
            OFFSET = intergenic_seq.start 
            for strand, sequence in [(1, intergenic_seq.sequence),
                                     (-1, intergenic_seq.sequence.reverse_complement())]:
                start = 0
                for start_codon in self.start_codons:
                    start = sequence.find(start_codon, start)
                    if start == -1:
                        continue
                    #Start searching right before our threshold
                    end = start + (self.min_aa_length-1)*3
                    found_stop = False
                    while end < len(sequence):
                        codon = sequence[end:end+3]
                        if str(codon) in self.stop_codons:
                            found_stop = True
                            break
                        else:
                            end = end + 3
                    
                    if not found_stop:
                        continue
                    
                    if (end-start)/3 > max_aa_seq_length or \
                        (end-start)/3 < min_aa_seq_length:
                        continue
                    
                    nt_subsequence = sequence[start:end+3]
                    aa_sequence = nt_subsequence.translate(11)
                    
                    #get nucleotide coords for original strand
                    if strand == -1:
                        old_end = end
                        end = len(sequence) - start
                        start = len(sequence) - old_end - 3
                        potential_orf = self.Sub_seq(aa_sequence, end+OFFSET, start+OFFSET)
                    else:
                        end = end + 3
                        potential_orf = self.Sub_seq(aa_sequence, start+OFFSET, end+OFFSET)
                    
                    
                    self.intergenic_orfs.append((strand, potential_orf))




#==============================================================================
# For testing!
#==============================================================================
#input_file = "chris_file.gbk"
#upper_lim = 10000
#lower_lim = 50

def __main__():
    parser = argparse.ArgumentParser("Get intergenic potential ORFs from *.gbk file")
    parser.add_argument('file', type=str)  #input file
    parser.add_argument('-u', type=int, default=100, help='Maximum size of potential ORF') 
    parser.add_argument('-l', type=int, default=20, help='Minimum size of potential ORF') 
    ARGS, unparsed = parser.parse_known_args()
    
    input_file = ARGS.file
    upper_lim = ARGS.u
    lower_lim = ARGS.l

    records = SeqIO.parse(input_file,"genbank") 
    for record in records: #may have multiple records in file
        main_record = My_record(record, lower_lim, upper_lim)
        print("Record " + main_record.id +'\n')
        for strand, orf in main_record.intergenic_orfs:
            print("Potential ORF of length " + str(len(orf.sequence)) + 
                  " found at " + str(orf.start) + ":" + str(orf.end) + 
                  " on strand " + str(strand))
            print(orf.sequence + '\n')
        print("*"*60)


__main__()