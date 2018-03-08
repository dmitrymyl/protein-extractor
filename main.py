from sys import argv
from collections import Counter, deque
from itertools import cycle
from functools import wraps
import datetime


log_file = open("log.txt", "w")


def log_wrapper(name, file):
    """
    Writes time of calling the specified function in given file.
    """
    def decorator(func):
        @wraps(func)
        def wrapper(*args, **kwargs):
            current = datetime.datetime.now()
            file.write("{}\t{}\tfunc started\n".format(name, current))
            result = func(*args, **kwargs)
            current = datetime.datetime.now()
            file.write("{}\t{}\tfunc ended\n".format(name, current))
            return result
        return wrapper
    return decorator


@log_wrapper("argv_parser", log_file)
def argv_parser(argv_list, qualifiers):
    """
    Parses CLI argument line. If no qualifiers specified, prints
    description and table of qualifiers. If not empty, checks whether
    all requiered and no unknown qualifiers are specified with correct values.
    """
    if len(argv_list) < 2:
        print("The program extracts all possible protein sequences "
              "from given nucleotide fasta file.")
        print("Obligatory qualifiers:")
        print("-input\t\tstring\tName of nucleotide file in .fasta format.\n"
              "-output\t\tstring\tname of file output.\n"
              "-minlen\t\tint\tMinimal length of protein to output.\n"
              "-binsize\tint\tLength of bin.")
        exit(0)
    else:
        qualifiers.update(dict(zip(argv_list[1::2], argv_list[2::2])))
        standard = ["-input", "-output", "-binsize", "-minlen"]
        for i in standard:
            if i not in qualifiers.keys():
                print("Died:")
                print("{} not found.".format(i))
                exit(1)
        for i in qualifiers.keys():
            if i not in standard:
                print("Died:")
                print("Unknown qualifier: {}".format(i))
                exit(1)
        try:
            qualifiers["-binsize"] = int(qualifiers["-binsize"])
        except ValueError:
            print("Died:")
            print("-binsize must be an integer value.")
            exit(1)
        try:
            qualifiers["-minlen"] = int(qualifiers["-minlen"])
        except ValueError:
            print("Died:")
            print("-minlen must be an integer value.")
            exit(1)


# @log_wrapper("char2bit", log_file)
# No logging as it is frequently used function.
def char2bit(nucl):
    """
    Returns a bit value of any base according to
    given dictionary. Process unknown nucleotide
    and empty string.
    """
    dictmap = {'A': 0, 'a': 0, 'T': 3, 't': 3, 'G': 1, 'g': 1,
               'C': 2, 'c': 2, '': None, "N": 4, 'n': 4, '\n': 5}
    return dictmap[nucl]


# @log_wrapper("bit2char", log_file)
# No logging as it is frequently used function.
def bit2char(value):
    """
    Returns a base for given int value 0..3 according
    to given dictmap.
    """
    dictmap = {0: 'A', 1: 'G', 2: 'C', 3: 'T'}
    return dictmap[value]


# @log_wrapper("revcompl", log_file)
# No logging as it is frequently used function.
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


@log_wrapper("binnate", log_file)
def binnate(input_iterable, binsize):
    """
    v1.2. Returns generator with bins as lists.
    """
    input_list = sorted(input_iterable)
    bound = binsize
    for value in input_list:
        while value >= bound:
            yield buff
            buff.clear()
            bound += binsize
        buff.append(value)
    yield buff


@log_wrapper("fasta_writer", log_file)
def fasta_writer(seqiter, file, start='', end='', number='', reverse=False):
    """
    Writes any sequence in seqiter into file in fasta format.
    Optional arguments are coordinates in GenBank notation and
    number of sequence.
    """
    rev = ''
    if reverse:
        rev = "(reverse)"
    file.write(">{}_{}..{}{}\n".format(number, start, end, rev))
    for i in range(0, len(seqiter), 60):
        file.write("{}\n".format(''.join(seqiter[i:i + 59])))


@log_wrapper("counter_writer", log_file)
def counter_writer(counter, file, number=''):
    """
    Writes any counter in file in fasta format.
    Optional argument is number of counter.
    """
    file.write(">{}\n".format(number))
    for k, v in counter.items():
        file.write("{}\t{}\n".format(k, v))


table_forward = {63: 'F', 62: 'F', 60: 'L', 61: 'L', 59: 'S', 58: 'S', 56: 'S',
                 57: 'S', 51: 'Y', 50: 'Y', 55: 'C', 54: 'C', 53: 'W', 47: 'L',
                 46: 'L', 44: 'L', 45: 'L', 43: 'P', 42: 'P', 40: 'P', 41: 'P',
                 35: 'H', 34: 'H', 32: 'Q', 33: 'Q', 39: 'R', 38: 'R', 36: 'R',
                 37: 'R', 15: 'I', 14: 'I', 12: 'I', 13: 'M', 11: 'T', 10: 'T',
                 8: 'T', 9: 'T', 3: 'N', 2: 'N', 0: 'K', 1: 'K', 7: 'S',
                 6: 'S', 4: 'R', 5: 'R', 31: 'V', 30: 'V', 28: 'V', 29: 'V',
                 27: 'A', 26: 'A', 24: 'A', 25: 'A', 19: 'D', 18: 'D', 16: 'E',
                 17: 'E', 23: 'G', 22: 'G', 20: 'G', 21: 'G', 48: 0, 49: 0,
                 52: 0}
table_reverse = {0: 'F', 16: 'F', 48: 'L', 32: 'L', 4: 'S', 20: 'S', 52: 'S',
                 36: 'S', 12: 'Y', 28: 'Y', 8: 'C', 24: 'C', 40: 'W', 1: 'L',
                 17: 'L', 49: 'L', 33: 'L', 5: 'P', 21: 'P', 53: 'P', 37: 'P',
                 13: 'H', 29: 'H', 61: 'Q', 45: 'Q', 9: 'R', 25: 'R', 57: 'R',
                 41: 'R', 3: 'I', 19: 'I', 51: 'I', 35: 'M', 7: 'T', 23: 'T',
                 55: 'T', 39: 'T', 15: 'N', 31: 'N', 63: 'K', 47: 'K', 11: 'S',
                 27: 'S', 59: 'R', 43: 'R', 2: 'V', 18: 'V', 50: 'V', 34: 'V',
                 6: 'A', 22: 'A', 54: 'A', 38: 'A', 14: 'D', 30: 'D', 62: 'E',
                 46: 'E', 10: 'G', 26: 'G', 58: 'G', 42: 'G', 60: 0, 44: 0,
                 56: 0}
codon_counter = Counter({63: 0, 62: 0, 60: 0, 61: 0, 59: 0, 58: 0,
                         56: 0, 57: 0, 51: 0, 50: 0, 55: 0, 54: 0,
                         53: 0, 47: 0, 46: 0, 44: 0, 45: 0, 43: 0,
                         42: 0, 40: 0, 41: 0, 35: 0, 34: 0, 32: 0,
                         33: 0, 39: 0, 38: 0, 36: 0, 37: 0, 15: 0,
                         14: 0, 12: 0, 13: 0, 11: 0, 10: 0, 8: 0,
                         9: 0, 3: 0, 2: 0, 0: 0, 1: 0, 7: 0, 6: 0,
                         4: 0, 5: 0, 31: 0, 30: 0, 28: 0, 29: 0, 27: 0,
                         26: 0, 24: 0, 25: 0, 19: 0, 18: 0, 16: 0,
                         17: 0, 23: 0, 22: 0, 20: 0, 21: 0})
bitcodon = {63: 'TTT', 62: 'TTC', 60: 'TTA', 61: 'TTG', 59: 'TCT', 58: 'TCC',
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
qual = dict()
argv_parser(argv, qual)
filename = qual["-input"]
minlen = qual["-minlen"]
binsize = qual["-binsize"]
mnemonic = qual["-output"]
try:  # Can it open the input file?
    fasta = open(filename, "r")
except IOError as e:
    print(e)
    exit(1)
protein_file = open("{}.fasta".format(mnemonic), "w")
aa_file = open("{}.aacontent".format(mnemonic), "w")
title = fasta.readline()
print("Parsing {}".format(title[1:].strip()))
orf = cycle([0, 1, 2])
base = 5  # Like a newline.
cur_pos = 0
cur_orf = 0
triplet_size = 0
nucl_counter = {0: 0, 1: 0, 2: 0, 3: 0}
handler = [[] for _ in range(6)]  # A structure for all growing proteins in six RFs.
triplet = 0
length_list = []
protein_amount = 0
while base is not None:  # Not EOF
    if base < 4:  # Base is a regular one
        triplet &= 15
        triplet <<= 2
        triplet += base
        nucl_counter[base] += 1
        cur_pos += 1
        triplet_size += 1
        cur_orf = next(orf)
    elif base == 4:  # base is N, hence any peptide sequence disrupts.
        seed_size = 0
        triplet = 0
        cur_pos += 1
        for i in handler:
            i.clear()
        cur_orf = next(orf)
    else:  # Got newline.
        base = char2bit(fasta.read(1))
        continue
    if triplet_size >= 3:  # If codon is ready.
        # Check forward string.
        if table_forward[triplet] == 'M':  # Start codon: start a new seq.
            handler[cur_orf].append([cur_pos - 2, deque()])
            for protein in handler[cur_orf]:
                protein[1].append(triplet)
        elif table_forward[triplet] == 0:  # Stop codon: terminate.
            for protein in handler[cur_orf]:
                protein.append(cur_pos)
            while handler[cur_orf]:
                buff = handler[cur_orf].pop()
                if len(buff[1]) >= minlen:  # Good protein, write it.
                    protein_amount += 1
                    codon_counter += Counter(buff[1])
                    protseq = list()
                    length_list.append(len(buff[1]))
                    while buff[1]:
                        protseq.append(table_forward[buff[1].popleft()])  # Like queue.
                    fasta_writer(protseq, protein_file, start=buff[0],
                                 end=buff[2], number=protein_amount)
                    counter_writer(Counter(protseq), aa_file,
                                   number=protein_amount)
        else:  # Any other codon.
            for protein in handler[cur_orf]:
                protein[1].append(triplet)
        # Check reverse string.
        rev_triplet = revcompl(triplet, 3)
        if table_reverse[rev_triplet] == 'M':  # Start codon: terminate seq.
            for protein in handler[cur_orf + 3]:
                protein[1].append(rev_triplet)
                protein.append(cur_pos)
            while handler[cur_orf + 3]:
                buff = handler[cur_orf + 3].pop()
                if len(buff[1]) >= minlen:  # Good protein, write it.
                    protein_amount += 1
                    codon_counter += Counter(buff[1])
                    protseq = list()
                    length_list.append(len(buff[1]))
                    while buff[1]:
                        protseq.append(table_reverse[buff[1].pop()])  # Like stack.
                    fasta_writer(protseq, protein_file, start=buff[0],
                                 end=buff[2], number=protein_amount,
                                 reverse=True)
                    counter_writer(Counter(protseq), aa_file,
                                   number=protein_amount)
        elif table_reverse[rev_triplet] == 0:  # Stop codon: start a new seq.
            handler[cur_orf + 3].clear()
            handler[cur_orf + 3].append([cur_pos - 2, deque()])
        else:  # Any other codon: append.
            for protein in handler[cur_orf + 3]:
                protein[1].append(rev_triplet)
    base = char2bit(fasta.read(1))
# Do output.
print("Found {} proteins".format(protein_amount))
print("Length distribution across proteins:")
start, end = 0, binsize
for i in binnate(length_list, binsize):
    print("{}-{}\t{}".format(start, end, len(i)))
    start, end = end, end + binsize
binfile = open("bin.txt", "w")
binfile.write(str(length_list))
binfile.close()
print("Total codon content:")
for k, v in codon_counter.items():
    print(bitcodon[k], v, sep="\t")
print("Nucleotide content:")
for k, v in nucl_counter.items():
    print(bit2char(k), v, sep="\t")
# Close all files.
fasta.close()
protein_file.close()
log_file.close()
aa_file.close()
