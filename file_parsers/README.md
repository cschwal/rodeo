## Genbank_to_igorfs.py
Genbank file parser for extracting potential intergenic ORFs from  intergenic sequences. 

### Usage
```python genebank_to_igorfs.py file [-l LOWER_LIMIT ] [-u UPPER_LIMIT] ```
The limits represent the length constraints on the potential ORFs. The default range is 40 to 1000000.  

When ran with ```test_file.gbk```, the output is as appears below:
~~~~
user@machine:/~path$ python genebank_to_igorfs.py test_file.gbk

Record c00021_N912Ps

Potential ORF of length 85 found at 325:70 on strand -1
MASTSALAVGPVTSWGSGTVRSVSSQAGTPSTARTPSDTSRIPRTVEPGRNTGISVPVAGPTSSARPSSNTRMSAWPSGSTRCG*

Potential ORF of length 139 found at 428:11 on strand -1
VAPQHPGRLGHPSVRLGAGQPARDGGLRAGVPGADGVHQRPGRRPGHQLGQRHRAVGQLPGRHPEHRAHSVGHQPHPAHGGTGAQHRYLRPGRRPDQLGAAVVEHEDVGMAVRQHALRVARTAVGSPLLPHPVDRVPP*

Potential ORF of length 91 found at 284:11 on strand -1
LGQRHRAVGQLPGRHPEHRAHSVGHQPHPAHGGTGAQHRYLRPGRRPDQLGAAVVEHEDVGMAVRQHALRVARTAVGSPLLPHPVDRVPP*

Potential ORF of length 75 found at 1538:1763 on strand 1
MVPEVRAGRRAYTRAGGRAGRAAASAWTVATPAAPPRKVRTPQGRVVANGNPGRLAGQCHRKQTAPARERAEQG*

Potential ORF of length 59 found at 2144:1967 on strand -1
MAFRRGQAWSTSSCSRRRAAVDTALRPNRVPKLAAEDPNQDRTSVRAGHRPKSGRAGL*

Potential ORF of length 104 found at 3596:3284 on strand -1
MPPGGCWTPPATAPSSASAQPIPKVMPPAAHSRWRRTLGVPLQATFRLIYSAKSPRPGSGGRTESLSSIGTGIRLRPETIPRSSRRDATRGFYRPNVLDCRHT*

************************************************************
~~~~

### Requirements
 -[Biopython](http://biopython.org/wiki/Biopython)
