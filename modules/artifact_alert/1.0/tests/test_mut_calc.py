import pandas as pd
import os
import sys
sys.path.append("/home/kyakimovich/repos/kurt/modules/artifact_alert/1.0/scripts")

import calculate_mutation_rates as mut_calc

TEST_OUT = "/home/kyakimovich/repos/kurt/modules/artifact_alert/1.0/tests/test_out"
TEST_DATA = "/home/kyakimovich/repos/kurt/modules/artifact_alert/1.0/tests/test_data"

def debug_position(chromosome, position, ref_base, depth, bases):
    """Debug function to see what's happening with base parsing"""
    parser = mut_calc.PileupParser()
    
    print(f"\nDEBUG - {chromosome}:{position} (ref={ref_base}, depth={depth})")
    print(f"Raw bases: {bases[:100]}{'...' if len(bases) > 100 else ''}")  # Truncate for readability
    
    # Show indel counting
    indel_counts = parser.count_indels(bases)
    print(f"Indels found: {indel_counts}")
    
    # Show cleaned string
    cleaned = parser.clean_bases_string(bases)
    print(f"Cleaned bases: {cleaned[:100]}{'...' if len(cleaned) > 100 else ''}")
    print(f"Raw length: {len(bases)}, Cleaned length: {len(cleaned)}")
    
    # Show final counts
    base_counts = parser.parse_position(bases, ref_base)
    print(f"Final counts: A={base_counts.A}, T={base_counts.T}, G={base_counts.G}, C={base_counts.C}")
    print(f"Ref={base_counts.ref}, Del={base_counts.deletions}, Ins={base_counts.insertions}")
    print(f"Total calculated: {base_counts.get_total()}")
    print(f"Expected depth: {depth}, Got total: {base_counts.get_total()}, Diff: {depth - base_counts.get_total()}")

def run_debug_test():
    """Run debug analysis on real pileup positions"""
    
    input_pileup = os.path.join(TEST_DATA, "real_test.pileup")
    
    if not os.path.exists(input_pileup):
        print(f"ERROR: Input file not found: {input_pileup}")
        return
    
    print("DEBUGGING REAL PILEUP POSITIONS:")
    print("="*80)
    
    # Read the pileup file, skipping comment lines
    columns = ['chromosome', 'position', 'ref_base', 'depth', 'bases', 'qualities']
    df = pd.read_csv(input_pileup, sep='\t', names=columns, header=None, comment='#')
    
    # Debug all positions to see the patterns
    for _, row in df.iterrows():
        debug_position(
            row['chromosome'], 
            row['position'], 
            row['ref_base'], 
            row['depth'], 
            row['bases']
        )

def run_test():
    """Run the mutation rate calculation test"""
    
    # Input and output paths
    input_pileup = os.path.join(TEST_DATA, "real_test.pileup")
    output_file = os.path.join(TEST_OUT, "real_mutation_rates.tsv")
    
    # Ensure output directory exists
    os.makedirs(TEST_OUT, exist_ok=True)
    
    # Check input file exists
    if not os.path.exists(input_pileup):
        print(f"ERROR: Input file not found: {input_pileup}")
        return
    
    # Set test parameters
    min_depth = 30
    min_alt_count = 1
    
    print(f"Reading test data from: {input_pileup}")
    print(f"Writing results to: {output_file}")
    print(f"Parameters: min_depth={min_depth}, min_alt_count={min_alt_count}")
    
    try:
        # Run the analysis - this SHOULD create the output file
        mut_calc.process_pileup_file(input_pileup, output_file, min_depth, min_alt_count)
        
        # Verify the file was created
        if os.path.exists(output_file):
            print(f"\nSUCCESS: Output file created at {output_file}")
            
            # Read and display results
            results = pd.read_csv(output_file, sep='\t')
            print("\n" + "="*80)
            print("TEST RESULTS:")
            print("="*80)
            print(results.to_string(index=False))
            
            # Summary statistics
            print(f"\nSUMMARY:")
            print(f"Total positions analyzed: {len(results)}")
            if len(results) > 0:
                print(f"Positions with SNVs: {(results['snv_alt_count'] > 0).sum()}")
                print(f"Positions with indels: {(results['indel_alt_count'] > 0).sum()}")
                print(f"Mean total mutation rate: {results['total_mutation_rate'].mean():.4f}")
        else:
            print(f"ERROR: Output file was not created at {output_file}")
            
    except Exception as e:
        print(f"ERROR during processing: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    # Run both the regular test and debug
    run_test()
    
    print("\n" + "="*80)
    
    # Run debug to see what's happening with real positions
    run_debug_test()
