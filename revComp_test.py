from revComp import *
import unittest


class TestRevComp(unittest.TestCase):
    def test_revcomp(self):
        self.assertEqual(revComp("accgttaattgccgt"),"acggcaattaacggt")
        self.assertEqual(revComp("ACCGTTAATTGCCGT",True),"ACGGGCAATTAACGGT")
        self.assertEqual(revComp("agtcgahgattc"),1)
        self.assertEqual(revComp("agtcgtagcnnn---taagct"),"agctta---nnngctacgact")


    def test_revcompcheat(self):
        self.assertEqual(revCompCheat("fasta.fasta","fasta"),">label1\ntttttcaagaagaccc")
        self.assertEqual(revCompCheat("multifasta.fasta","fasta"),">label2\nctttgtataatccc\n>label3\ngtcctctaaat\nctttgtataatccc>label4\ntagcgatcgggaattcat")
        self.assertEqual(revCompCheat("agggcttca","string"),"tgaagccct")


    def test_checkformat(self):
        self.assertEqual(checkformat("fasta","fasta.fasta"),True)
        self.assertEqual(checkformat("string","aagggttac"),True)
        self.assertEqual(checkformat("fasta","attaggsc"),False)
        self.assertEqual(checkformat("string","fasta.fasta"),False)


    def test_checkdna(self):
        self.assertEqual(checkDNA("agggcttca"),True)
        self.assertEqual(checkDNA("agggcnA-ttca"),True)
        self.assertEqual(checkDNA("agdggcfttca"),False)
        self.assertEqual(checkDNA("agggRYcttca"),False)


if __name__ == '__main__':
    unittest.main()

