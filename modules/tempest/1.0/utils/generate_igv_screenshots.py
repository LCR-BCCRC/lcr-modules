#!/usr/bin/env python

import os
import argparse
import subprocess
import tempfile

class mutation():

    def __init__(self, sample, chrom, pos, ref, alt):
        self.sample = sample
        self.chrom = chrom
        self.pos = pos
        self.ref = ref  # Ref allele
        self.alt = alt  # Alt allele


def is_valid_file(filepath, parser):
    if not os.path.exists(filepath):
        raise parser.error("Unable to locate \'%s\': No such file or directory" % filepath)
    else:
        return filepath


def get_args():
    parser = argparse.ArgumentParser(description="Generates a set of IGV screenshots of variants in a specified sample")
    parser.add_argument("-m", "--maf", metavar="MAF", required=True, type=lambda x: is_valid_file(x, parser), help="Input MAF file specifying the variants of interest")
    parser.add_argument("-b", "--bam", metavar="BAM/CRAM", required=True, type=lambda x: is_valid_file(x, parser), nargs="+", help="One or more BAM files which correspond to the sample of interest")
    parser.add_argument("--genome_build", metavar="hg19,hg38", default="hg38", choices=["hg38", "hg19"], help="Reference genome build [Default: hg38]")
    parser.add_argument("-o", "--outdir", metavar="/path/to/output/directory/", required=True, type=str, help="Output directory to place IGV screenshots")
    return parser.parse_args()


def load_mutations(mutfile, samplecol = "Tumor_Sample_Barcode", chromcol = "Chromosome", poscol = "Start_Position", refcol = "Reference_Allele", altcol = "Tumor_Seq_Allele2"):
    """
    Load annotated mutations from a TSV file

   @params:
        :input mutfile: A string containing a filepath to the mutations of interest
    """

    muts = []  # List containing mutation() objects representing the mutations
    header = {samplecol: None, chromcol: None, poscol: None, refcol: None, altcol: None}
    header_proc = None
    with open(mutfile) as f:
        for line in f:
            line = line.rstrip("\n").rstrip("\r")

            if line.startswith("#"):
                continue  # Skip comment lines

            cols = line.split("\t")
            if not header_proc:  #i.e we haven't seen the header yet
                # Find out which column index corresponds with which attribute
                i = 0
                for col in cols:
                    if col in header:
                        header[col] = i
                    i += 1

                # Check that all the required columns have been found
                for column, index in header.items():
                    if index is None:
                        raise AttributeError(f"Unable to locate column \'{column}\' in the file header")
                header_proc = True  # Header has been successfully processed
                continue

            # Process variant
            sample_id = cols[header[samplecol]]
            chrom = cols[header[chromcol]]
            pos = int(cols[header[poscol]])
            ref = cols[header[refcol]]
            alt = cols[header[altcol]]

            # Store variant
            mut = mutation(sample_id, chrom, pos, ref, alt)
            muts.append(mut)

    return muts

def generate_igv_screenshots(mutations: iter, bams: iter, outdir: str, genome_build: str):
    """
    Using the provided mutation location and BAM files, generate an IGV screenshot
    for each mutation for each sample
    """

    igv_batch = tempfile.mkstemp(text=True, suffix=".igv")[1]

    with open(igv_batch, "w") as o:
        o.write("new" + os.linesep)
        o.write("genome " + genome_build + os.linesep)

        # Specify BAM files to load
        for bam in bams:
            o.write("load " + bam + os.linesep)

        o.write("snapshotDirectory " + outdir + os.linesep)

        # Now this is where the real magic happens. Go to the genomic coordinates of interest, and generate IGV screenshots
        for mut in mutations:
            # Add a buffer on either side of the position (ex. 100bp)
            start = mut.pos - 100
            end = mut.pos + 100
            coord = mut.chrom + ":" + str(start) + "-" + str(end)
            o.write("goto " + coord + os.linesep)
            o.write("sort BASE " + mut.chrom + ":" + str(mut.pos) + os.linesep)
            # Set output image name
            outname = mut.sample + "_" + mut.chrom + "-" + str(mut.pos) + "_" + mut.ref + "-" + mut.alt + ".png"
            o.write("snapshot " + outname + os.linesep)

        o.write("exit" + os.linesep)

    # Now that we have created the IGV batch file, actually run IGV
    xvfb_com = "xvfb-run -a -n 0 -e /dev/stderr --server-args '-screen 0 1920x1080x24 -terminate' igv --batch " + igv_batch
    subprocess.check_call(xvfb_com, shell=True)
    os.remove(igv_batch)

    # Screenshots generated

def main(args=None):
    if args is None:
        args = get_args()

    # Step 1. Load mutation from MAF file
    muts = load_mutations(args.maf)

    # Generate an IGV batch file and run IGV
    generate_igv_screenshots(muts, args.bam, args.outdir, args.genome_build)


if __name__ == "__main__":
    main()
