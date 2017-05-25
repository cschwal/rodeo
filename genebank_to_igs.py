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
    
    def __init__(self, record):
        self.record = record
        self.sequence = record.seq
        self.CDSs = []
        self.intergenic_seqs = []
        self.intergenic_orfs = []
        self.start_codons = ['ATG','GTG', 'TTG']
        self.stop_codons = ['TAA','TAG','TGA']
            
        
    class Sub_seq(object):
        def __init__(self, seq, start, end):
            self.start = start
            self.end = end
            self.sequence = seq
            
                
    def get_CDSs(self):
        for feature in self.record.features:
            if feature.type == 'CDS':
                start = int(feature.location.start)
                end = int(feature.location.end)
                orf = self.Sub_seq(seq = "", start=start, end=end) #Don't need CDS seq...
                self.CDSs.append(orf)
       
    def get_intergenic_seqs(self):
        start = 0
        for orf in self.CDSs:
            end = orf.start
            nn_seq = self.sequence[start:end]
            intergenic_sequence = self.Sub_seq(seq=nn_seq, 
                                               start = start,
                                               end=end)
            self.intergenic_seqs.append(intergenic_sequence)
            start = orf.end
    
    def get_intergenic_orfs(self, min_aa_seq_length, max_aa_seq_length=1000000):
        for intergenic_seq in self.intergenic_seqs:
            OFFSET = intergenic_seq.start #OFFSET will adjust for the main sequence
            for strand, sequence in [(1, intergenic_seq.sequence),
                                     (-1, intergenic_seq.sequence.reverse_complement())]:
                start = 0
                for start_codon in self.start_codons:
                    start = sequence.find(start_codon, start)
                    if start == -1:
                        continue
                    end = start + 3
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
    parser.add_argument('-u', type=int, default=1000000, help='Maximum size of potential ORF') 
    parser.add_argument('-l', type=int, default=50, help='Minimum size of potential ORF') 
    ARGS, unparsed = parser.parse_known_args()
    
    input_file = ARGS.file
    upper_lim = ARGS.u
    lower_lim = ARGS.l

    records = SeqIO.parse(input_file,"genbank")
    for record in records:
        main_record = My_record(record)
        main_record.get_CDSs()
        main_record.get_intergenic_seqs()
        main_record.get_intergenic_orfs(lower_lim, upper_lim)
    
    for strand, orf in main_record.intergenic_orfs:
        print("Potential ORF of length " + str(len(orf.sequence)) + 
              " found at " + str(orf.start) + ":" + str(orf.end) + 
              " on strand " + str(strand))
        print(orf.sequence + '\n')


__main__()