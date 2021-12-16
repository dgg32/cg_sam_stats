import unittest
#import pandas as pd
#from pandas.testing import assert_series_equal
import numpy as np
import helper_function


class Test_Function(unittest.TestCase):


    def test_get_version(self):
        ist= helper_function.get_version()
        soll = "0.1-2-ge241d54"

        self.assertEqual(ist, soll)

    def test_read_sam(self):
        test_sam_file = "./test/SEnoBarcode/Intensities/FQMAP/B1234567_L02_read.sam"

        sam_object = helper_function.read_sam(test_sam_file)

        soll = sam_object["B1234567L2C001R003_366917"]["RNAME"] 
        ist = "gi|49175990|ref|NC_000913.2|"

        self.assertEqual(ist, soll)

    def test_read_sam_annotation(self):
        test_sam_file = "./test/SEnoBarcode/Intensities/FQMAP/B1234567_L02_read.sam"

        sam_object = helper_function.read_sam(test_sam_file)

        soll = sam_object["B1234567L2C001R003_366917"]["ANNOTATION"]["MD:Z"] 
        ist = "0A39"

        self.assertEqual(ist, soll)
    
    def test_get_ref_from_sam(self):
        test_sam_file = "./test/SEnoBarcode/Intensities/FQMAP/B1234567_L02_read.sam"

        sam_object = helper_function.read_sam(test_sam_file)

        soll = helper_function.get_ref(sam_object["B1234567L2C001R003_366917"])
        ist = "AAAAAAGATTTAAGCAAATATAAAAAAAGACAATGGTTTC"

        self.assertEqual(ist, soll)


    def test_get_ref_from_case_1(self):
        
        soll = helper_function.get_ref({"SEQ": "CGATACGGGGACATCCGGCCTGCTCCTTCTCACATG", "CIGAR": "36M", "ANNOTATION": {"MD:Z": "1A0C0C0C1T0C0T27"}})
        ist = "CACCCCTCTGACATCCGGCCTGCTCCTTCTCACATG"

        self.assertEqual(ist, soll)

    def test_get_ref_from_insertion(self):
        
        soll = helper_function.get_ref({"SEQ": "GAGACGGGGTGACATCCGGCCTGCTCCTTCTCACAT", "CIGAR": "6M1I29M", "ANNOTATION": {"MD:Z": "0C1C0C1C0T0C27"}})
        ist = "CACCCC-TCTGACATCCGGCCTGCTCCTTCTCACAT"

        self.assertEqual(ist, soll)

    def test_get_ref_from_deletion(self):
        
        soll = helper_function.get_ref({"SEQ": "AGTGATGGGGGGGTTCCAGGTGGAGACGAGGACTCC", "CIGAR": "9M9D27M", "ANNOTATION": {"MD:Z": "2G0A5^ATGATGTCA27"}})
        ist = "AGGAATGGGATGATGTCAGGGGTTCCAGGTGGAGACGAGGACTCC"

        self.assertEqual(ist, soll)

    def test_get_ref_from_mix(self):
        
        soll = helper_function.get_ref({"SEQ": "AGTGATGGGAGGATGTCTCGTCTGTGAGTTACAGCA", "CIGAR": "2M1I7M6D26M", "ANNOTATION": {"MD:Z": "3C3T1^GCTCAG26"}})
        ist = "AG-GCTGGTAGCTCAGGGATGTCTCGTCTGTGAGTTACAGCA"

        self.assertEqual(ist, soll)

    def test_get_ref_from_simple(self):
    
        soll = helper_function.get_ref({"SEQ": "CAAAAAGATTTAAGCAAATATAAAAAAAGACAATGGTTTC", "CIGAR": "40M", "ANNOTATION": {"MD:Z": "0A39"}})
        ist = "AAAAAAGATTTAAGCAAATATAAAAAAAGACAATGGTTTC"

        self.assertEqual(ist, soll)

    def test_get_triplets(self):
        soll = helper_function.get_triplets("ATCGT")
        ist = ["ATC", "TCG", "CGT"]

        self.assertEqual(ist, soll)
    
    def test_get_triplets_mutations(self):
        soll = helper_function.get_triplets_mutations({"SEQ": "CAAA", "CIGAR": "4M", "ANNOTATION": {"MD:Z": "0A3"}})
        #ist = "AAAA"
        ist = {"AAA": {"CAA": 1, "AAA": 1}} 

        self.assertEqual(ist, soll)

    def test_get_triplets_mutations_smart_1(self):
        soll = helper_function.get_triplets_mutations_smart({"SEQ": "CAAA", "CIGAR": "4M", "ANNOTATION": {"MD:Z": "4"}})
        #ist = "AAAA"
        ist = {"CAA": {"CAA": 1}, "AAA": {"AAA": 1}} 

        self.assertEqual(ist, soll)

    def test_get_triplets_mutations_smart_2(self):
        soll = helper_function.get_triplets_mutations_smart({"SEQ": "CAAA", "CIGAR": "4M", "ANNOTATION": {"MD:Z": "0A3"}})
        #ist = "AAAA"
        ist = {"AAA": {"CAA": 1, "AAA": 1}} 

        self.assertEqual(ist, soll)

    def test_get_triplets_mutations_smart_3(self):
        sam_object = {"SEQ": "CAAAAAGATTTAAGCAAATATAAAAAAAGACAATGGTTTC", "CIGAR": "40M", "ANNOTATION": {"MD:Z": "0A39"}}
        soll = helper_function.get_triplets_mutations_smart(sam_object)
        #ist = "AAAA"
        ist = helper_function.get_triplets_mutations(sam_object) 

        self.assertEqual(ist, soll)

    def test_get_triplets_mutations_smart_4(self):
        sam_object = {"SEQ": "GAGACGGGGTGACATCCGGCCTGCTCCTTCTCACAT", "CIGAR": "6M1I29M", "ANNOTATION": {"MD:Z": "0C1C0C1C0T0C27"}}
        soll = helper_function.get_triplets_mutations_smart(sam_object)
        #ist = "AAAA"
        ist = {} 

        self.assertEqual(ist, soll)
    
    def test_get_triplets_mutations_smart_format_1(self):
        soll_triplets = helper_function.get_triplets_mutations_smart({"SEQ": "CAAA", "CIGAR": "4M", "ANNOTATION": {"MD:Z": "0A3"}})
        soll = helper_function.format_triplet_mutations(soll_triplets)
        #ist = "AAAA"
        
        ist = "X--	-X-	--X	Subset	A--	C--	G--	T--	N--	-A-	-C-	-G-	-T-	-N-	--A	--C	--G	--T	--N\nA\tA\tA\tRead 1\t1\t1\t0\t0\t0\t2\t0\t0\t0\t0\t2\t0\t0\t0\t0\n"

        #print ("ist\n",ist)

        self.assertEqual(ist, soll)


if __name__ == '__main__':
    unittest.main()