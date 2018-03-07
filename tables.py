"""
The set of genetic tables.
"""
# Common genetic table for bacteria.
table = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
         'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
         'TAT': 'Y', 'TAC': 'Y', 'TGT': 'C', 'TGC': 'C',
         'TGG': 'W', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L',
         'CTG': 'L', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P',
         'CCG': 'P', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q',
         'CAG': 'Q', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R',
         'CGG': 'R', 'ATT': 'I', 'ATC': 'I', 'ATA': 'I',
         'ATG': 'M', 'ACT': 'T', 'ACC': 'T', 'ACA': 'T',
         'ACG': 'T', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K',
         'AAG': 'K', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R',
         'AGG': 'R', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V',
         'GTG': 'V', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A',
         'GCG': 'A', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E',
         'GAG': 'E', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G',
         'GGG': 'G', 'TAA': 0, 'TAG': 0, 'TGA': 0}
# Codons to their bit representation in decimal values.
codon2bit = {'TTT': 63, 'TTC': 62, 'TTA': 60, 'TTG': 61, 'TCT': 59, 'TCC': 58,
             'TCA': 56, 'TCG': 57, 'TAT': 51, 'TAC': 50, 'TGT': 55, 'TGC': 54,
             'TGG': 53, 'CTT': 47, 'CTC': 46, 'CTA': 44, 'CTG': 45, 'CCT': 43,
             'CCC': 42, 'CCA': 40, 'CCG': 41, 'CAT': 35, 'CAC': 34, 'CAA': 32,
             'CAG': 33, 'CGT': 39, 'CGC': 38, 'CGA': 36, 'CGG': 37, 'ATT': 15,
             'ATC': 14, 'ATA': 12, 'ATG': 13, 'ACT': 11, 'ACC': 10, 'ACA': 8,
             'ACG': 9, 'AAT': 3, 'AAC': 2, 'AAA': 0, 'AAG': 1, 'AGT': 7,
             'AGC': 6, 'AGA': 4, 'AGG': 5, 'GTT': 31, 'GTC': 30, 'GTA': 28,
             'GTG': 29, 'GCT': 27, 'GCC': 26, 'GCA': 24, 'GCG': 25, 'GAT': 19,
             'GAC': 18, 'GAA': 16, 'GAG': 17, 'GGT': 23, 'GGC': 22, 'GGA': 20,
             'GGG': 21, 'TAA': 48, 'TAG': 49, 'TGA': 52}
# Bit to letter representation of codons.
bit2codon = {63: 'TTT', 62: 'TTC', 60: 'TTA', 61: 'TTG', 59: 'TCT', 58: 'TCC',
             56: 'TCA', 57: 'TCG', 51: 'TAT', 50: 'TAC', 55: 'TGT', 54: 'TGC',
             53: 'TGG', 47: 'CTT', 46: 'CTC', 44: 'CTA', 45: 'CTG', 43: 'CCT',
             42: 'CCC', 40: 'CCA', 41: 'CCG', 35: 'CAT', 34: 'CAC', 32: 'CAA',
             33: 'CAG', 39: 'CGT', 38: 'CGC', 36: 'CGA', 37: 'CGG', 15: 'ATT',
             14: 'ATC', 12: 'ATA', 13: 'ATG', 11: 'ACT', 10: 'ACC', 8: 'ACA',
             9: 'ACG', 3: 'AAT', 2: 'AAC', 0: 'AAA', 1: 'AAG', 7: 'AGT',
             6: 'AGC', 4: 'AGA', 5: 'AGG', 31: 'GTT', 30: 'GTC', 28: 'GTA',
             29: 'GTG', 27: 'GCT', 26: 'GCC', 24: 'GCA', 25: 'GCG', 19: 'GAT',
             18: 'GAC', 16: 'GAA', 17: 'GAG', 23: 'GGT', 22: 'GGC', 20: 'GGA',
             21: 'GGG', 48: 'TAA', 49: 'TAG', 52: 'TGA'}
# Bit codons with their aminoacids.
bit2aa = {63: 'F', 62: 'F', 60: 'L', 61: 'L', 59: 'S', 58: 'S', 56: 'S',
          57: 'S', 51: 'Y', 50: 'Y', 55: 'C', 54: 'C', 53: 'W', 47: 'L',
          46: 'L', 44: 'L', 45: 'L', 43: 'P', 42: 'P', 40: 'P', 41: 'P',
          35: 'H', 34: 'H', 32: 'Q', 33: 'Q', 39: 'R', 38: 'R', 36: 'R',
          37: 'R', 15: 'I', 14: 'I', 12: 'I', 13: 'M', 11: 'T', 10: 'T',
          8: 'T', 9: 'T', 3: 'N', 2: 'N', 0: 'K', 1: 'K', 7: 'S', 6: 'S',
          4: 'R', 5: 'R', 31: 'V', 30: 'V', 28: 'V', 29: 'V', 27: 'A',
          26: 'A', 24: 'A', 25: 'A', 19: 'D', 18: 'D', 16: 'E', 17: 'E',
          23: 'G', 22: 'G', 20: 'G', 21: 'G', 48: 0, 49: 0, 52: 0}
# Reverse-complement bit codons with their aminoacids.
revcomplbit2aa = {0: 'F', 16: 'F', 48: 'L', 32: 'L', 4: 'S', 20: 'S', 52: 'S',
                  36: 'S', 12: 'Y', 28: 'Y', 8: 'C', 24: 'C', 40: 'W', 1: 'L',
                  17: 'L', 49: 'L', 33: 'L', 5: 'P', 21: 'P', 53: 'P', 37: 'P',
                  13: 'H', 29: 'H', 61: 'Q', 45: 'Q', 9: 'R', 25: 'R', 57: 'R',
                  41: 'R', 3: 'I', 19: 'I', 51: 'I', 35: 'M', 7: 'T', 23: 'T',
                  55: 'T', 39: 'T', 15: 'N', 31: 'N', 63: 'K', 47: 'K',
                  11: 'S', 27: 'S', 59: 'R', 43: 'R', 2: 'V', 18: 'V', 50: 'V',
                  34: 'V', 6: 'A', 22: 'A', 54: 'A', 38: 'A', 14: 'D', 30: 'D',
                  62: 'E', 46: 'E', 10: 'G', 26: 'G', 58: 'G', 42: 'G', 60: 0,
                  44: 0, 56: 0}
