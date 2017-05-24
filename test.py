# -*- coding: utf-8 -*-

#==============================================================================
# #main TODOS:
#   Get multiple query responses working (WPs)
#   Make sure to account for different dictionary keys (brute force search?)
#==============================================================================

import RODEO_utils as R_utils
acc_seq = []
acc_seq1 = "CCM09445.1"
acc_seq2 = "BAL72546.1"
acc_seq3 = "WP_050838132.1"
#records_save = R_utils.get_efetch_results(query=acc_seq2)
results = R_utils.get_efetch_results(query=acc_seq2)
results.print_info()
print("*"*20)

results.trim_to_n_nucleotides(N=2000)
#results.trim_to_n_neighbors(3)

results.print_info()

from Bio import SeqIO
from Bio import Seq
from Bio.Alphabet import generic_dna
table = 11
min_pro_len = 30

def find_orfs_with_trans(seq, trans_table, min_protein_length):
    answer = []
    seq_len = len(seq)
    for strand, nuc in [(+1, seq), (-1, seq.reverse_complement())]:
        for frame in range(3):
            trans = str(nuc[frame:].translate(trans_table))
            trans_len = len(trans)
            aa_start = 0
            aa_end = 0
            while aa_start < trans_len:
                aa_end = trans.find("*", aa_start)
                if aa_end == -1:
                    aa_end = trans_len
                if aa_end-aa_start >= min_protein_length:
                    if strand == 1:
                        start = frame+aa_start*3
                        end = min(seq_len,frame+aa_end*3+3)
                    else:
                        start = seq_len-frame-aa_end*3-3
                        end = seq_len-frame-aa_start*3
                    answer.append((start, end, strand,
                                   trans[aa_start:aa_end]))
                aa_start = aa_end+1
    answer.sort()
    return answer


seq  = Seq.Seq(results.cluster_sequence, generic_dna)
orf_list = find_orfs_with_trans(seq, table, min_pro_len)
for start, end, strand, pro in orf_list:
    print("%s...%s - length %i, strand %i, %i:%i" \
          % (pro[:30], pro[-3:], len(pro), strand, start, end))
