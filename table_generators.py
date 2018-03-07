"""
Several dictionaries for handling genetic tables.
"""


def char2bit(nucl):
    """
    Returns a bit value of any base according to
    given dictionary. Process unknown nucleotide
    and empty string.
    """
    dictmap = {'A': 0, 'a': 0, 'T': 3, 't': 3, 'G': 1,
               'g': 1, 'C': 2, 'c': 2, '': None, "N": 4, 'n': 4}
    return dictmap[nucl]


def string2bit(string):
    """
    Transforms given nucleotide sequence
    to an integer according to char2bit
    dictionary.
    """
    value = 0
    for i in string:
        value <<= 2
        value += char2bit(i)
    return value


def revcompl(value1, length):
    """
    Returns a revcompl nucl seq of given length in value1
    in bit format.
    """
    value2 = 0
    value1 = ~value1 & 63
    for _ in range(length):
        buff = value1 & 3
        value2 <<= 2
        value2 += buff
        value1 >>= 2
    return value2


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
stop_codons = ['TAA', 'TAG', 'TGA']
start_codon = 'ATG'
buff = dict.fromkeys(table.keys())
for i in buff.keys():
    buff[i] = string2bit(i)
print("codon2bit = ", buff, sep='')
translate = dict((v, k) for k, v in buff.items())
print("bit2codon = ", translate, sep='')
bittable = dict()
for i in translate.keys():
    bittable[i] = table[translate[i]]
print("bit2aa = ", bittable, sep='')
revbittable = dict()
for i in bittable.keys():
    revbittable[revcompl(i, 3)] = bittable[i]
print("revcomplbit2aa = ", revbittable, sep='')
