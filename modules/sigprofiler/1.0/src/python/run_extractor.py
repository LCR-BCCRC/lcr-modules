import argparse
from SigProfilerExtractor import sigpro as sig
import os

os.environ["OPENBLAS_NUM_THREADS"] = "1" # ulimit fix

def main():
    """ Performs NMF decomposition with provided matrix """
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
    parser.add_argument('input_data', help = 'Path to matrix that will undergo NMF.')
    parser.add_argument('output', help = 'Path to output directory.')
    parser.add_argument('ref', help = 'Reference build of matrix.')
    parser.add_argument('context_type', default = '96,DINUC,ID', help = 'A string of mutaion context name/names separated by comma (","). The items in the list defines the mutational contexts to be considered to extract the signatures. The default value is "96,DINUC,ID", where "96" is the SBS96 context, "DINUC" is the DINUCLEOTIDE context and ID is INDEL context.')
    parser.add_argument('exome', default = False, type = bool, = 'Defines if the exomes will be extracted. The default value is "False".')
    parser.add_argument('minimum_signatures', default = 1, type = int, help = "The minimum number of signatures to be extracted. The default value is 1.")
    parser.add_argument('maximum_signatures', default = 25, type = int, help = "The maximum number of signatures to be extracted. The default value is 25.")
    parser.add_argument('nmf_replicates', default = 100 , type = int, help = "The number of iteration to be performed to extract each number signature. The default value is 100.")
    parser.add_argument('matrix_normalization', default = 'gmm', help = 'Method of normalizing the genome matrix before it is analyzed by NMF. Default is value is "gmm". Other options are, "log2", "custom" or "none".')
    parser.add_argument('nmf_init', default = 'nndsvd_min', help = "The initialization algorithm for W and H matrix of NMF. Options are 'random', 'nndsvd', 'nndsvda', 'nndsvdar' and 'nndsvd_min'. Default is 'nndsvd_min'.")
    parser.add_argument('cpu', default = -1, type = int, help = 'Number of threads to use. Default use all threads.')
    
    return(parser.parse_args())

if __name__ == '__main__':
    main()
