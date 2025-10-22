"""Simple test suite for code to find STRs in a sequence.

Can see test cases coded for, and test new ones to help
adjust code when needed.

"""

import re
import numpy as np
import sys

# get parent directory
# and append to path

current_path = sys.path[0]
# get parent directory
parent_directory = current_path[:current_path.rfind("/")]
sys.path.append(parent_directory)

from utils.FetchVariantUMIs import determine_STR_status


def run_tests():
    """Run tests for the determine_STR_status function."""
    
    test_cases = [
        # Format: (name, start, allele, sequence, expected_result)
        
        # Single base repeat cases
        ("Simple A repeat", 10, "A", "CTGATAAAAAAGCT", True),
        ("No repeat", 5, "C", "CCATGCTGACTA", False),
        ("Borderline repeat", 5, "G", "ATCGGGACT", True),
        ("Start edge case", 0, "T", "TTTACGTA", True),
        ("End edge case", 7, "A", "TCGATAAAA", True),
        
        # Multi-base repeat cases
        ("Direct dinucleotide repeat", 2, "AT", "GCATATATCG", True),
        ("Separated repeat", 2, "GC", "TAGCATGCTT", False),
        ("Dinucleotide pattern", 3, "CG", "ATACGCGCGTA", True),
        
        # Complex cases
        ("Complex repeat", 5, "ATTG", "CGAAAATTGATTGCC", True),
        ("Partial repeat", 3, "GCTA", "ATGCTAGCGCTAGC", False),
        ("Long sequence repeat", 10, "ACGTACGT", "TGACACGTACGTACGTACGTTC", True),
        ("Half repeat", 3, "ACAC", "TGACACACACGT", True),
        ("No half repeat", 3, "ACGT", "TGAACGTACGTT", True),
        ("Empty allele", 5, "", "ACGTACGT", False),
        ("Special chars in allele", 3, "A+T", "CGAAA+TA+TACGT", np.nan),  # Now expects False due to skipping
        ("Trinucleotide repeat", 5, "CAG", "ATGGCAGCAGCAGCAGTA", True),
        ("Complex STR", 10, "ATATAT", "GCGCATATATATATATGCGC", True),
        ("Long repeat",5, "AAAAAA", "GCGCAAAAAAAAGCGC", True), # homopolymer
        ("short reference allele homoplymer", 5, "AA", "GCGCAAAAAAAAGCGC", True), # homopolymer
        ("repeated allele, but not STR", 2, "AA", "GCAATTGAAATGCA", False), # homopolymer

    ]
    
    passed = 0
    failed = 0
    
    print("\nRunning STR Detection Tests...\n")
    
    for name, start, allele, seq, expected in test_cases:
        result = determine_STR_status(start, allele, seq)
        if result == expected:
            passed += 1
        elif np.isnan(result) and np.isnan(expected):
            passed += 1
        else:
            failed += 1
            print(f"FAILED: {name} - Got {result}, Expected {expected}")
            print(f"  Allele: '{allele}', Start: {start}, Sequence: '{seq}'\n")
    # color text green if all tests passed, red if any failed
    if failed == 0:
        print("\033[92mAll tests passed!\033[0m")
    else:
        print(f"\033[91m{failed} tests failed!\033[0m")
    
    return passed == len(test_cases)

def main():
    run_tests()

if __name__ == "__main__":
    main()