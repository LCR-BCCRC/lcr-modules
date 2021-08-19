import argparse
from SigProfilerExtractor import sigpro as sig
import os

os.environ["OPENBLAS_NUM_THREADS"] = "1"

def main():
    args = parse_arguments()
    sig.sigProfilerExtractor(input_type = 'matrix', \
                             output = args.output, \
                             input_data = args.input_data, \
                             reference_genome = args.ref, \
                             opportunity_genome = args.ref, \
                             context_type = args.context_type, \
                             exome = args.exome, \
                             minimum_signatures = args.minimum_signatures, \
                             maximum_signatures = args.maximum_signatures, \
                             nmf_replicates = args.nmf_replicates, \
                             matrix_normalization = args.matrix_normalization, \
                             nmf_init = args.nmf_init, \
                             cpu = args.cpu, \
                             make_decomposition_plots = False)
    return()

def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_data')
    parser.add_argument('output')
    parser.add_argument('ref')
    parser.add_argument('context_type', default = '96,DINUC,ID')
    parser.add_argument('exome', default = False, type = bool)
    parser.add_argument('minimum_signatures', default = 1, type = int)
    parser.add_argument('maximum_signatures', default = 25, type = int)
    parser.add_argument('nmf_replicates', default = 100 , type = int)
    parser.add_argument('matrix_normalization', default = 'gmm')
    parser.add_argument('nmf_init', default = 'nndsvd_min')
    parser.add_argument('cpu', default = -1, type = int)
    
    return(parser.parse_args())

if __name__ == '__main__':
    main()
