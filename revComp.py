#!/usr/bin/env python3

import sys
import os
import argparse
from Bio import SeqIO, Seq


parser = argparse.ArgumentParser(description='Do reverse complementing')
parser.add_argument('inputSeq', type=str, help='the DNA sequence to be reverse complemented, as a string or in a FASTA formatted file (specify with --fasta)')
parser.add_argument('--fasta', dest='fasta', action='store_true')
parser.add_argument('--cheat', dest='cheat', action='store_true')
parser.add_argument('--caps', dest='caps', action='store_true')
parser.add_argument('--out', dest='out', default=None, help='Specify the name for your output file, without file extension. If no name entered, output will be printed to the console.')
parser.set_defaults(cheat=False, fasta=False, inputSeq="", caps=False)

args = parser.parse_args()


def checkformat(fasta,inputSeq):
    """
    Check whether format matches flag.
    """
    if fasta:
        assert os.path.isfile(inputSeq), "Oh no! I expected a FASTA file but all I got was a string!"
        with open(inputSeq, 'r') as handle:
            fas = SeqIO.parse(handle, 'fasta')
            return any(fas)
    else:
        assert checkDNA(inputSeq), "Please enter just the sequence, without any FASTA header (or invoke --fasta)."
        return True


def checkDNA(sequence):
    """
    Checks that the input sequence contains only contains A, T, G, or C 
    (or N, or - for gapped sequences). Character check is case-insensitive.
    """
    allowed = ["a", "t", "g", "c", "n", "-"]
    match = [characters.lower() in allowed for characters in sequence.lower()]
    return all(match)


def revComp(sequence,caps=False):
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


def revCompCheat(sequence,fasta):
    """
    Reverse complements using biopython. If input is a simple string, the Seq library
    is used, else the SeqIO library handles full FASTA files.
    """
    if not fasta:
        seq = Seq.Seq(sequence)
        return str(seq.reverse_complement())
    else:
        seq_objects = []
        with open(sequence, "r") as infile:
            for seq_record in SeqIO.parse(infile, "fasta"):
                ID = seq_record.id
                new_seq_record = seq_record.reverse_complement()
                new_seq_record.id = ID
                seq_objects.append(new_seq_record)
        #convert list of sequence objects back to string
        output = ''
        for obj in seq_objects:
            output = output + ">" + str(obj.id) + '\n' + str(obj.seq) + '\n'
        #return the output string
        return output


def main():
    checkformat(args.fasta,args.inputSeq)
    output = ''
    if args.cheat or args.fasta:
        # use biopython to reverse complement the sequence
        output = revCompCheat(args.inputSeq,args.fasta)
    else:
        # check sequence is DNA - if it is, reverse complement it - if not, exit program
        if checkDNA(args.inputSeq):
            output = revComp(args.inputSeq)
        else:
            sys.exit("Non compliant sequence provided. Is it DNA? I can handle \"N\" and \"-\" but nowt else")
    if args.out:
        outfilename = str(args.out) + '.fasta'
        with open(outfilename, 'w') as outfile:
            outfile.write(output)
        print('Output written to {}.fasta'.format(args.out))
    else:
        print(output)
        
    

if __name__ == "__main__":
    main()
