"""
Puts "N" bases in random place of given sequence.
The first qualifier is an input sequence filename,
the second qualifier is an ouput sequnce filename,
the third qualifier (int) is an average number of
"N" out of 10000 input bases.
"""
from sys import argv
import random

foutput = open(argv[2], "w")
sixty_counter = 0
with open(argv[1], "r") as finput:
    header = finput.readline()
    foutput.write(header)
    base = finput.read(1)
    while base != "":
        if base != "\n":
            happen = random.randint(0, 10000)
            if happen > 10000 - int(argv[3]):
                foutput.write("N")
                sixty_counter += 1
                if sixty_counter % 60 == 0:
                    foutput.write("\n")
            foutput.write(base)
            sixty_counter += 1
            if sixty_counter % 60 == 0:
                foutput.write("\n")
        base = finput.read(1)
foutput.close()
