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
    input_format = args.inType
    # if input is .seg file, convert it to the .bed file and have separate file created with a header
    if input_format in ["seg", "bedpeA", "bedpeB"]:
        check_arguments(args, input_format)
        chromosome = args.chromColnum-1
        start = args.startColnum-1
        end = args.endColnum-1
        tsv_to_bed(input_file, output_file, chromosome, start, end, input_format)
    # if input is .bed file, load appropriate header and convert it to the .seg file
    elif input_format == "bed":
        check_arguments(args, input_format)
        column_names_path=args.column_header
        with open(column_names_path) as f:
            column_names=simplejson.load(f)
        bed_to_seg(input_file, column_names, output_file)
    elif input_format == "bedpeA_bedpeB": 
        check_arguments(args, input_format)
        bedpeA = input_file
        bedpeB = input_file.replace("bedpeA", "bedpeB")
        column_names_path=args.column_header
        with open(column_names_path) as f:
            column_names=simplejson.load(f)
        bed_to_bedpe(bedpeA, bedpeB, column_names, output_file)
    else:
        sys.exit('Please provide input file in either .seg, .bedpe, or .bed format')

# Counts VCF header lines for bedpe file parsing
def count_header(fname):
    i=0
    with open(fname, 'r') as f:
        for line in f:
            if "##" in line:
                i += 1
    return i


def tsv_to_bed(input_file, output_file, chromosome, start, end, inType):
    
    # Some bedpe files start with a VCF header
    # Since the header is distinguished by "##" and Pandas can only use a single-character comment, 
    # we must parse the file and count the comment lines and indicate to read_table how many 
    # lines to skip. 

    if inType in ["bedpeA", "bedpeB"]: 
        skip_lines = count_header(input_file)
        print("Skipping " + str(skip_lines) + " lines from the bedpe VCF header. ")
    else:
        skip_lines = 0
    
    # import .seg file
    seg = pd.read_table(input_file, skiprows=skip_lines)

    seg.rename(columns={seg.columns[chromosome]: "chrom", seg.columns[start]: "start", seg.columns[end]: "end"}, inplace=True)
    seg.fillna('NA', inplace=True) 

    # Drop partner coordinates for bedpe so column 4 is identical in both files
    if(inType == "bedpeA"):
        seg.drop(seg.columns[[3, 4, 5]], axis = 1, inplace=True)
    if(inType == "bedpeB"): 
        seg.drop(seg.columns[[0, 1, 2]], axis = 1, inplace=True)

    # rearrange columns order to have first 3 cols according the BED format
    bed = seg.loc[:, ['chrom', 'start', 'end']]
    bed_other = seg.drop(['chrom', 'start', 'end'], axis=1)

    # Create collapsed column name from all non-coordinate colnames
    other_colnames = "|".join(list(bed_other.columns))

    # Create a new df storing all non-coordinate column values collapsed
    # other_collapsed = pd.DataFrame()
    # other_collapsed = other_collapsed.append({other_colnames: bed_other.apply(lambda x: '|'.join(x.astype(str).values), axis=1)}, ignore_index = True)

    bed.loc[:, other_colnames] = bed_other.apply(lambda x: '|'.join(x.astype(str).values), axis=1)
    # shift start position by 1 to the left
    bed.loc[:, 'start'] = bed['start'].apply(lambda x: int(x-1))
    # check that chromosomes are prefixed and prefix if they are not
    chrom = list(bed['chrom'])
    for i in range(len(chrom)):
        if 'chr' not in chrom[i]:
            chrom[i]='chr'+str(chrom[i])
    bed.loc[:, 'chrom']=chrom

    # remove all columns with extra information that was just concatenated into a single column
    bed = bed[['chrom', 'start', 'end', other_colnames]]

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


def restore_columns(df): 
    
    # Extract the coordinate columns (everything but the concatenated columns)
    coord_columns = df.iloc[:,:-1]

    # Get the remaining non-coordinate column
    other_columns = df.iloc[:,-1:].reset_index()
    # Create a list of column names from the split colname
    names = other_columns.columns.values[-1]
    col_names = list(str(names).split('|'))

    # then split values for each feature in a column
    other_columns = other_columns.join(other_columns[names].str.split('|', expand=True))
    other_columns = other_columns.iloc[:, 2:]
    other_columns.columns = col_names
    # Join the coordinate columns with the other columns
    filled = coord_columns.join(other_columns)

    return(filled)

    
def bed_to_seg(input_file, column_names, output_file):
    # import .seg file
    seg = pd.read_table(input_file, index_col=None, header=None, names=column_names)
    seg.fillna('NA', inplace=True)
    # shift start position by 1 to the right
    seg['start'] = seg['start'].apply(lambda x: x+1)

    seg = restore_columns(seg)

    # rearrange columns order to match it original .seg file
    ID_column = seg[['ID']]
    other_columns = seg.drop(['ID'], axis=1)
    seg = ID_column.join(other_columns)

    # write resulting data frame to the output file
    seg.to_csv(output_file, header=True, index=False, sep="\t")




def bed_to_bedpe(bedpeA, bedpeB, column_names, output_file):
    # import .bed files
    bedpeA = pd.read_table(bedpeA, index_col=None, header=None, names=column_names)
    bedpeA.fillna('NA', inplace=True)
    # Rename columns to 
    bedpeA.rename(columns={"chrom": "CHROM_A", "start": "START_A", "end": "END_A"}, inplace=True)

    bedpeB = pd.read_table(bedpeB, index_col=None, header=None, names=column_names)
    bedpeB.fillna('NA', inplace=True)
    # shift start position by 1 to the right
    bedpeB.rename(columns={"chrom": "CHROM_B", "start": "START_B", "end": "END_B"}, inplace=True)

    # Join bed files to make bedpe

    bedpe = pd.merge(bedpeA, bedpeB, on = column_names[-1])
    bedpe = bedpe.iloc[:, [0, 1, 2, 4, 5, 6, 3]]

    bedpe = restore_columns(bedpe)

    bedpe.to_csv(output_file, header=True, index=False, sep="\t")


def check_arguments(args, input_format):
    if input_format in ['seg', "bedpeA", "bedpeB"] and not all([args.chromColnum, args.startColnum, args.endColnum]):
        raise ValueError ('Must specify number of columns in segmentation file containing name of chromosome, starting and ending position of the feature')

    if input_format in ['bed', 'bedpeA_bedpeB'] and not args.column_header:
        raise ValueError ('Must specify file containing header of the BED file used for liftOver conversion')



def parse_args():
    parser = argparse.ArgumentParser()

    parser.add_argument("--input",
                        help="Initial file to convert to/from .bed format", required=True)
    parser.add_argument("--output",
                        help="Resulting file after conversion", required=True)
    parser.add_argument("--column-header",
                        help="When converting from .bed to .seg format, provide file containing column headers")
    parser.add_argument("--chromColnum", type=int,
                        help="number of column in a segmentation file that contains information about name of chromosome")
    parser.add_argument("--startColnum", type=int,
                        help="number of column in a segmentation file that contains information about start position of the feature")
    parser.add_argument("--endColnum", type=int,
                        help="number of column in a segmentation file that contains information about end position of the feature")
    parser.add_argument("--inType", 
                        help="Type of input file (bedpeA, bedpeB, seg, bed, or bedpeA_bedpeB)", type=str, required=True)

    args, unknown = parser.parse_known_args()

    return args



if __name__ == '__main__':
    main()
