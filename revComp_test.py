from revComp import *
import unittest


class TestRevComp(unittest.TestCase):
    def test_revcomp(self):
        self.assertEqual(revComp("accgttaattgccgt"),"acggcaattaacggt")
        self.assertEqual(revComp("ACCGTTAATTGCCGT",True),"ACGGCAATTAACGGT")
        self.assertRaises(SystemExit, revComp, "agtcgahgattc")
        self.assertEqual(revComp("agtcgtagcnnn---taagct"),"agctta---nnngctacgact")


    def test_revcompcheat(self):
        self.assertEqual(revCompCheat("testFiles/fasta.fasta",True),">label1\ntttttcaagaagaccc\n")
        self.assertEqual(revCompCheat("testFiles/multifasta.fasta",True),">label2\nctttgtataatccc\n>label3\ngtcctctaaat\n>label4\ntagcgatcgggaattcat\n")
        self.assertEqual(revCompCheat("agggcttca",False),"tgaagccct")


    def test_checkformat(self):
        self.assertTrue(checkformat(True,"testFiles/fasta.fasta"))
        self.assertTrue(checkformat(False,"aagggttac"))
        self.assertRaises(AssertionError, checkformat, True,"attaggsc")
        self.assertRaises(AssertionError, checkformat, False,"testFiles/fasta.fasta")
        self.assertFalse(checkformat(True,"testFiles/notfasta.csv"))


    def test_checkdna(self):
        self.assertTrue(checkDNA("agggcttca"))
        self.assertTrue(checkDNA("agggcnA-ttca"))
        self.assertFalse(checkDNA("agdggcfttca"))
        self.assertFalse(checkDNA("agggRYcttca"))


if __name__ == '__main__':
    unittest.main()

