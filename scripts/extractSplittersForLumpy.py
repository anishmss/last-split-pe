#!/usr/bin/env python

import argparse
from itertools import groupby

def samReadAlnsGen(fileHandle): #returns all alignments of one read
    for k,samLines in groupby(fileHandle, lambda x : (x.split()[0],int(x.split()[1])&64)):
        yield samLines #creates empty list whenever there is a switch


if __name__ == "__main__":
	
    parser = argparse.ArgumentParser()
    parser.add_argument("SAM", action = "store", help="alignments reported in the SAM format")
    args = parser.parse_args() 

    f = open(args.SAM)
    
    for samAlns in samReadAlnsGen(f):
        samAlns = list(samAlns)
        if len(samAlns) == 2 : # Lumpy can only handle splits with two flanks
            for line in samAlns:
                contents = line.split()
                if int(contents[1])&64:
                    contents[0] += "_1"
                else:
                    contents[0] += "_2"
                    
                print "\t".join(contents)
