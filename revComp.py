#!/usr/bin/env python3

import sys
import argparse
from Bio import SeqIO, Seq


parser = argparse.ArgumentParser(description='Do reverse complementing')
parser.add_argument('inputSeq', type=str, help='the DNA sequence to be reverse complemented, as a string or in a FASTA formatted file (specify with --format)')
parser.add_argument('--format', dest='formt', choices=['string', 'fasta'], default='string')
parser.add_argument('--cheat', dest='cheat', action='store_true')
parser.add_argument('--caps', dest='caps', action='store_true')
parser.add_argument('--out', dest='out', default=None, help='Specify the name for your output file, without file extension.')
parser.set_defaults(cheat=False, inputSeq="", caps=False)

args = parser.parse_args()


def checkformat():
    """
    Check whether format matches flag.
    """
    if formt == 'fasta':
        assert os.path.isfile(inputSeq), "Oh no! I expected a FASTA file but all I got was a string!"
    elif formt == 'string':
        assert checkDNA(inputSeq), "If --format is set to \'string\', please enter just the sequence, without any FASTA header."


def checkDNA(sequence):
    """
    Checks that the input sequence contains only contains A, T, G, or C 
    (or N, or - for gapped sequences). Character check is case-insensitive.
    """
    allowed = ["a", "t", "g", "c", "n", "-"]
    match = [characters.lower() in allowed for characters in sequence.lower()]
    return all(match)


def revComp(sequence):
    """
    Reverse complements the provided sequence by substituting the appropriate character and 
    prepending this to a new string (the output).
    """
    outputSeq = ""
    for character in sequence.lower():
        if character == "a":
            outputSeq = "t" + outputSeq
        elif character == "g":
            outputSeq = "c" + outputSeq
        elif character == "t":
            outputSeq = "a" + outputSeq
        elif character == "c":
            outputSeq = "g" + outputSeq
        elif character == "n":
            outputSeq = "n" + outputSeq
        elif character == "-":
            outputSeq = "-" + outputSeq
        else:
            sys.exit("Something's gone wrong! I found a \"{}\"!".format(character))
    if caps:
        return outputSeq.upper()
    else:
        return outputSeq


def revCompCheat(sequence):
    """
    Reverse complements using biopython. If input is a simple string, the Seq library
    is used, else the SeqIO library handles full FASTA files.
    """
    if formt == 'string':
        seq = Seq.Seq(inputSeq)
        return = str(seq.reverse_complement())
    elif formt == 'fasta':
        seq_objects = []
        with open(inputSeq, "r") as infile:
            for seq_record in SeqIO.parse(infile, "fasta"):
                ID = seq_record.id
                new_seq_record = seq_record.reverse_complement()
                new_seq_record.id = ID
                seq_objects.append(new_seq_record)
        #convert list of sequence objects back to string
        output = ''
        for obj in seq_objects:
            output = output + obj.id + '\n' + obj.seq + '\n'
        #return the output string
        return output


def main():
    checkformat()
    if cheat:
        # use biopython to reverse complement the sequence
        cheated = revCompCheat(inputSeq)
        if out:
            out = out + '.fasta'
            with open(out, 'w') as outfile:
                outfile.write(cheated)
            return 'Output written to {}'.format(out)
        else:
            return output
    else:
        # check sequence is DNA - if it is, reverse complement it - if not, exit program
        if checkDNA(inputSeq):
            revComp(inputSeq)
        else:
            sys.exit("Non compliant sequence provided. Is it DNA? I can handle \"N\" and \"-\" but nowt else")
        
    

if __name__ == "__main__":
    main()
