from revComp import *
import unittest

class TestRevComp(unittest.TestCase):
    def test_revcomp(self):
        self.assert(revComp("accgttaattgccgt"),"acggcaattaacggt)
        self.assert(revComp("ACCGTTAATTGCCGT",True),"ACGGGCAATTAACGGT")
        self.assert(revComp("agtcgahgattc"),1)
        self.assert(revComp("agtcgtagcnnn---taagct"),"agctta---nnngctacgact")

    def test_revcompcheat(self):
        self.assert(revCompCheat("fasta.fasta","fasta"),">label1\ntttttcaagaagaccc")
        self.assert(revCompCheat("multifasta.fasta","fasta"),">label2\nctttgtataatccc\n>label3\ngtcctctaaat\nctttgtataatccc>label4\ntagcgatcgggaattcat")
        self.assert(revCompCheat("agggcttca"),"tgaagccct")

    def test_checkformat(self):
        self.assert(checkformat("fasta","fasta.fasta"),True)
        self.assert(checkformat("string","aagggttac"),True)
        self.assert(checkformat("fasta","attaggsc"),False)
        self.assert(checkformat("string","fasta.fasta"),False)

    def test_checkdna(self):
        self.assert(checkDNA("agggcttca"),True)
        self.assert(checkDNA("agggcnA-ttca"),True)
        self.assert(checkDNA("agdggcfttca"),False)
        self.assert(checkDNA("agggRYcttca"),False)



def test():
    sequence = "accgttaattgccgt"
    sequence2 = "ACCGTTAATTGCCGT"
    sequence3 = "agtcgahgattc"
    sequence4 = "agtcgtagcnnn---taagct"
    fasta = "fasta.fasta"
    multifasta = "multifasta.fasta"
    

