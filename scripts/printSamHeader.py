"""
This script takes as input a reference genome fasta file, 
and prints (to stdout) SAM-format header lines,
containing information about IDs and lengths of sequences in the fasta file .
Usage:
printSamHeader.py reference.fa 
"""
import sys

fasta_in = sys.argv[1]

with open(fasta_in) as f:
    while True:
        line = f.readline()
        if line[0] == ">":
            break

    while True:
        header = line[1:].rstrip()
        lines = []
        line = f.readline()
        while True:
            if not line or line[0] == ">":
                break
            lines.append(line.rstrip())
            line = f.readline()
        print('@SQ\tSN:{}\tLN:{}'.format(header, len("".join(lines))))
        if not line:
            break
print("@PG\tID:LAST-SPLIT-PE\tPN:LAST-SPLIT-PE\tVN:0.00")
