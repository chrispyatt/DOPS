#!/usr/bin/env python3

import sys
import argparse
from Bio import SeqIO, Seq




parser = argparse.ArgumentParser(description='Do reverse complementing')
parser.add_argument('inputSeq', type=str, help='the DNA sequence to be reverse complemented, as a string or in a FASTA formatted file (specify with --format)')
parser.add_argument('--format', dest='format', choices=['string', 'fasta'], default='string')
parser.add_argument('--cheat', dest='cheat', action='store_true')
parser.add_argument('--caps', dest='caps', action='store_true')
parser.add_argument('--out', dest='out', default=None)
parser.set_defaults(cheat=False, inputSeq="", caps=False)

args = parser.parse_args()


def checkDNA(sequence):
    """
    Checks that the input sequence contains only contains A, T, G, or C 
    (or N, or - for gapped sequences). Character check is case-insensitive.
    """
    allowed = ["a", "t", "g", "c", "n", "-"]
    match = [characters.lower() in allowed for characters in sequence.lower()]
    
    return all(match)

def revComp(sequence):
    outputSeq = ""
    for character in sequence.lower():
        if character == "a":
            outputSeq = outputSeq + "t"
        elif character == "g":
            outputSeq = outputSeq + "c"
        elif character == "t":
            outputSeq = outputSeq + "a"
        elif character == "c":
            outputSeq = outputSeq + "g"
        elif character == "n":
            outputSeq = outputSeq + "n"
        elif character == "-":
            outputSeq = outputSeq + "-"
        else:
            sys.exit("Something's gone wrong! I found a \"{}\"!".format(character))
    if caps:
        return outputSeq.upper()
    else:
        return outputSeq


def main():
    if cheat:
        # use biopython to reverse complement the sequence
        if format=='string':
            seq = Seq.Seq(inputSeq)
            output = seq.reverse_complement()
        elif format=='fasta':
            with open(inputSeq, "r") as infile:
            outputFasta = 
                for seq_record in SeqIO.parse(infile, "fasta"):
                   new_seq_record = seq_record.reverse_complement()
    else:
        # check sequence is DNA - if it is, reverse complement it - if not, exit program
        if checkDNA(inputSeq):
            revComp(inputSeq)
        else:
            sys.exit("Non compliant sequence provided. Is it DNA? I can handle \"N\" and \"-\" but nowt else")
        
    

if __name__ == "__main__":
    main()
