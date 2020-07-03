"""
This script will convert segmentation files into BED format for conversion of genomic coordinates by liftOver tool.
In addition, this script can generate segmentation file from BED file after liftOver conversion.
"""


#!/usr/bin/python

import pandas as pd
import argparse
import simplejson
import sys



def main():
    # initiate the parser and handle arguments from command line
    args = parse_args()

    # read required first argument from shell command as path to .seg file to be converted
    input_file = args.input
    # read required second argument from shell command as path to resulting .bed file
    output_file = args.output
    # determine the format of input file to understand which format to convert
    input_format = input_file[-3:]
    # if input is .seg file, convert it to the .bed file and have separate file created with a header
    if input_format == "seg":
        check_arguments(args, input_format)
        chromosome = args.chromColnum-1
        start = args.startColnum-1
        end = args.endColnum-1
        seg_to_bed (input_file, output_file, chromosome, start, end)
    # if input is .bed file, load appropriate header and convert it to the .seg file
    elif input_format == "bed":
        check_arguments(args, input_format)
        column_names_path=args.column_header
        with open(column_names_path) as f:
            column_names=simplejson.load(f)
        bed_to_seg (input_file, column_names, output_file)
    else:
        sys.exit('Please provide input file in either .seg or .bed format')


def seg_to_bed(input_file, output_file, chromosome, start, end):
    # import .seg file
    seg = pd.read_table(input_file)
    seg.rename(columns={seg.columns[chromosome]: "chrom", seg.columns[start]: "start", seg.columns[end]: "end"}, inplace=True)
    seg.fillna('NA', inplace=True)

    # rearrange columns order to have first 3 cols according the BED format
    bed_required = seg[['chrom', 'start', 'end']]
    bed_other = seg.drop(['chrom', 'start', 'end'], axis=1)
    bed = bed_required.join(bed_other)
    # shift start position by 1 to the left
    bed['start'] = bed['start'].apply(lambda x: x-1)
    # check that chromosomes are prefixed and prefix if they are not
    chrom = list(bed['chrom'])
    for i in range(len(chrom)):
        if 'chr' not in chrom[i]:
            chrom[i]='chr'+str(chrom[i])
    bed['chrom']=chrom

    # concatenate all information not required by BED into one column
    # first concatenate column names
    names = list(bed_other.columns)
    col_names=''
    for name in names:
        col_names += str(name) + str('|')
    # then concatenate values of each feature
    bed[col_names] = ''
    for name in names:
        bed[col_names] = bed[col_names].astype(str)+bed_other[name].astype(str)+'|'

    # remove all columns with extra information that was just concatenated into a single column
    bed = bed[['chrom', 'start', 'end', col_names]]

    # write resulting data frame to the output file
    bed.to_csv(output_file, header=False, index=False, sep="\t")

    # save column names in a separate file with the same name and .header
    # create a list of column names
    col_names = list(bed.columns.values)
    # write to a file
    output_col_names=output_file+'.header'
    outF=open(output_col_names, "w")
    simplejson.dump(col_names, outF)
    outF.close()


def bed_to_seg(input_file, column_names, output_file):
    # import .seg file
    seg = pd.read_table(input_file, index_col=None, header=None, names=column_names)
    seg.fillna('NA', inplace=True)
    # shift start position by 1 to the right
    seg['start'] = seg['start'].apply(lambda x: x+1)

    # split information in the last column into separate columns to match original segmentation file
    # first get the names of the columns
    other_columns = seg.iloc[:,-1]
    names = str(other_columns.name)
    col_names = list(names.split('|')[:-1])

    # then split values for each feature in a column
    seg = seg.join(seg[names].str.split('|', expand=True))
    seg = seg.drop(names, axis=1)
    seg = seg.iloc[:, :-1]
    # rename the columns to match original segmentation file
    seg.columns = seg.columns[:3].tolist()+col_names

    # rearrange columns order to match it original .seg file
    ID_column = seg[['ID']]
    other_columns = seg.drop(['ID'], axis=1)
    seg = ID_column.join(other_columns)

    # write resulting data frame to the output file
    seg.to_csv(output_file, header=True, index=False, sep="\t")


def check_arguments(args, input_format):
    if input_format == 'seg' and not all([args.chromColnum, args.startColnum, args.endColnum]):
        raise ValueError ('Must specify number of columns in segmentation file containing name of chromosome, starting and ending position of the feature')

    if input_format == 'bed' and not args.column_header:
        raise ValueError ('Must specify file containing header of the BED file used for liftOver conversion')


def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--input",
                        help="Initial file to convert to/from .bed format", required=True)
    parser.add_argument("--output",
                        help="Resulting file after convertion", required=True)
    parser.add_argument("--column-header",
                        help="When converting from .bed to .seg format, provide file containing column headers")
    parser.add_argument("--chromColnum", type=int,
                        help="number of column in a segmentation file that contains information about name of chromosome")
    parser.add_argument("--startColnum", type=int,
                        help="number of column in a segmentation file that contains information about start position of the feature")
    parser.add_argument("--endColnum", type=int,
                        help="number of column in a segmentation file that contains information about end position of the feature")

    args, unknown = parser.parse_known_args()

    return args



if __name__ == '__main__':
    main()
