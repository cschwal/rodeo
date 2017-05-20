# -*- coding: utf-8 -*-

import RODEO_utils as R_utils
acc_seq = []
acc_seq1 = "CCM09445.1"
acc_seq2 = "BAL72546.1"
acc_seq3 = "WP_050838132.1"
#records_save = R_utils.get_efetch_results(query=acc_seq2)
results = R_utils.get_efetch_results(query=acc_seq1)

for orf in results.orfs:
    print("Accession:  " + orf.accession)
    print("Coords:  " + str(orf.start) + " to " +  str(orf.end))