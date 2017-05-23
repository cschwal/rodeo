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
