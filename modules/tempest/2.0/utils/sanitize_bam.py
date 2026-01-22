import pysam
import argparse


def get_args():
    parser = argparse.ArgumentParser()
    # config
    parser.add_argument('--tagstoremove',required=True,type=str,nargs='+',help='')
    parser.add_argument('--max_base_qual',required=True,type=int,help='')
#    parser.add_argument('--min_base_qual',required=True,type=int,help='')
    parser.add_argument('--threads',required=True,type=int,help='')
    # IO
    parser.add_argument('--in_bam',required=True,type=str,help='')
    parser.add_argument('--out_bam',required=True,type=str,help='')

    return parser.parse_args()

def san_bams(input_bam: str, threads: int, tagstoremove: list, max_base_qual: int, output_bam: str) -> None:
    inFile = pysam.AlignmentFile(input_bam, check_sq=False, mode = "rb", threads=threads)  # We have to provide check_sq=False in case this is an unaligned BAM
    tagstoremove = set(tagstoremove)
    # Remove the entries for these tags from the BAM header
    inHeader = inFile.header.to_dict()  # Returns a multi-level dictionary
    outHeader = {}
    for level, attributes in inHeader.items():
        # If we are removing the RG tag, remove this level
        if "RG" in tagstoremove and level == "RG":
            continue
        outHeader[level] = attributes

    outFile = pysam.AlignmentFile(output_bam, header=outHeader, mode = "wb")

    # Process reads
    for read in inFile.fetch(until_eof=True):
        # Cap the upper limit of base qualities
        outqual = list(qual if qual <= max_base_qual else max_base_qual for qual in read.query_qualities)
        # Mask bases with quality scores below the specified threshold
        masked_seq = "".join(read.query_sequence[i] if read.query_qualities[i] >= 20 else "N" for i in range(0, read.query_length))
        read.query_sequence = masked_seq
        read.query_qualities = outqual

        # Remove the unwanted tags
        # pysam returns read tags as a list of tuples
        # ex. [(NM, 2), (RG, "GJP00TM04")]
        outtags = list(tag for tag in read.get_tags() if tag[0] not in tagstoremove)
        read.set_tags(outtags)

        outFile.write(read)

    # For safety, manually close files
    inFile.close()
    outFile.close()


def main():
    args = get_args()
    san_bams(args.in_bam, args.threads, args.tagstoremove, args.max_base_qual, args.out_bam)

if __name__ == '__main__':
    main()

