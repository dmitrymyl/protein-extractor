# protein-extractor
A python 3.6 program that finds all possible open reading frames in whole genome.

## Contents of repo
1. `main.py` -- the main file with working program.
2. `tables.py` -- used genetic tables in format of python dictionaries.
3. `table_generators.py` -- the program that generates tables in tables.py
4. `ecoli.fasta` -- sample bacterium genome.

## main.py
The program takes only one file in `.fasta` format with only one record in it. It takes into account the reverse strand and unknown nucleotides as N. 
### Qualifiers
* -input    str Filename in `.fasta` format.
* -minlen   int Minimal length of proteins to be reported.
* -binsize  int Size of bins for protein length report.
* -output   str `mnemonics` for output files.
### Input
Only one file in `.fasta` containing one genome of any creature.
### Output
* `mnemonics.fasta`     `.fasta` file with all proteins which length is greater than minimal.
* `mnemonics.aacontent` The file with aminoacid content of each reported protein.
* `stdout`              Number of reported proteins, distribution of protein lengths across the bins of given length, codon content, given DNA nucleotide content.
* `log.txt`             The file with timestamps of calling all the functions in `main.py` program except for exceptionally frequently called ones.
