# protein-extractor
A python 3.6 program that finds all possible open reading frames in whole genome.

## Contents of repo
File | Description
-----|------------
`main.py` | The main file with working program.
`tables.py` | Sample genetic tables in format of python dictionaries.
`table_generators.py` | The program that generates tables in tables.py
`random_n_gen.py` | The program that modifies given fasta sequence. Insert N in random places with an average defined by user (some integer out of 10000 input bases).
`ecoli.fasta` | Sample bacterium genome.
`garbage.fasta` | Sample bacterium genome modified with `random_n_gen.py`.
`proteins.png` | An image for README.md
`deque.png` | An image for README.md

## Main program
The program takes only one file in `.fasta` format with only one record in it. It takes into account the reverse strand and unknown nucleotides as N. 
### Qualifiers
Qualifier | type | Description
----------|------|------------
-input | str | Filename in `.fasta` format.
-minlen | int | Minimal length of proteins to be reported.
-binsize | int | Size of bins for protein length report.
-output | str | `mnemonics` for output files.
### Input
Only one file in `.fasta` containing one genome of any creature.
### Output
File | Description
-----|------------
`mnemonics.fasta` | `.fasta` file with all proteins which length is greater than minimal.
`mnemonics.aacontent` | The file with aminoacid content of each reported protein.
`stdout` | Number of reported proteins;
`stdout` | distribution of protein lengths across the bins of given length;
`stdout` | codon content;
`stdout` | given DNA nucleotide content.
`log.txt` | The file with timestamps of calling all the functions in `main.py` program except for exceptionally frequently called ones.

### Benchmarks
The program works around 65 sec with `-minlen 0` and 35 sec with `-minlen 100` on a sample genome. The difference is possibly due to the fact that the bulk of found proteins are shorter the 100 aa, so that it takes a lot of time to write them in `.fasta` file.

## Algorithm
### Bit representation
To enhance speed, all characters are translated into integers in the following manner:

Char | A | G | C | T | N | n
-----|---|---|---|---|---|---
decimal | 0 | 1 | 2 | 3 | 4 | 5
binary | 00 | 01 | 10 | 11 | 100 | 101

The main point is notating ATGC in two-bit alphabet such that the complementary sequence is returned with NOT operation and two sequences can be checked on complementarity with XOR operation. N and \n are notated for convenience. Integer representation allows for presenting any nucleotide sequences of defined length as a unique integer value and easily taking slices with bit shifts and AND operations. Hence, reverse complement string of specified length can be derieved with simple bitwise operations.

### Flow genome parsing

The program `main.py` works with whole genome in a flow manner reading the sequence per byte after parsing the header. Each new character is checked wether it is a base (thus completing the codon), N (thus discrupting the codon and all growing protein sequences in all ORFs) or newline (thus continuing reading). The completed codon is subsequently checked If it is start codon on forward string or stop on reverse, a new protein is started to grow in a particular ORF. If it is a regular codon, it is added to growing proteins in particular ORF regarding to strandness. Finally, if it is a stop codon on forward string or start codon on reverse, all protein sequences in particular ORF are checked individually whether they are long enough to be reported. If so, the translation tooks place, aminoacid and codon content are count and protein sequence is written in file in 'fasta' format as well as its length is added to a list of reported proteins' lengths.

Thus, two strands are processed with only one file read. This algorithm implies the definition of protein on forward strand as any sequence that starts with start codon and ends with the closest stop codon and on reverse strand as the shortest sequence from stop codon to the closest start codon.

![alt text](https://github.com/dmitrymyl/protein-extractor/blob/master/proteins.png "Proteins")

### Implementation

To get ORF each time a new character is read, the cyclic iterator upon (0, 1, 2) is used. All ORFs are stored in list of 6 lists, which indexes are number of current ORF for forward string and increased by 3 for reverse string.

Protein sequences are stored in double-ended queues (deques) as codons are being read from the file. This allows to differentially access forward and reverse sequences with the same datatype. Being read forward, codons are placed in deque in a stack manner (LIFO). Protein sequences in forward ORFs are coded from start to end of deque, whereas those in reverse ORFs are coded from end to start. So, in case of forward ORFs, aminoacids are accessed in queue manner (FIFO) and in stack manner in case of reverse ORFs. The deque type provide left and right `pop()` method with a constant time of execution.

![alt text](https://github.com/dmitrymyl/protein-extractor/blob/master/deque.png "Deque")
