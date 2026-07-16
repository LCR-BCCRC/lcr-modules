#!/usr/bin/env python

import argparse
import os
import sys

# Converts and merges IgBLAST results with MiXCR results

# Set arguments:
# Input file: Igblast database results file
# Input file: Original mixcr results
# Output file: Updated mixcr results

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d','--database',required=True)
    parser.add_argument('-m','--mixcr',required=True)
    parser.add_argument('-o','--output',required=True)
    parser.add_argument('-s','--sequence',required=True)
    args = parser.parse_args()

    return args

def main(args):
    db_file = args.database
    mixcr_results = args.mixcr
    output_file = args.output
    seq_info_file = args.sequence

    # lists to make sure that mixcr results and igblastn results have same amount of clones
    mixcr_clones = []
    igblastn_clones = []

    # get header line to extract column positions
    mixcr_count = 0
    mixcr_positions = {}

    # open original mixcr file in order to append results to each clone
    # open output file -- not sure if  rewriting file or making new file is better? less redundant mayb?
    mixcr = open(mixcr_results, 'r')
    out = open(output_file, 'w')

    for mixcr_line in mixcr:
        mixcr_count += 1
        mixcr_line = mixcr_line.strip('\n')
        if mixcr_count == 1:
            # add new columns to mixcr header
            new_header = mixcr_line + "\tProductive\tigblastnTopVAllele\tigblastnTopSequenceIdentity\tmutatedStatus\tigblastnTopBitScore\tigblastnTopEVal\tigblastnVAllele\tsequenceIdentity\tallMutatedStatus\tallBitScores\tallEVals\tigblastMutations\tigblastnJAllele\tigblastnTopVAlleles\tnumMissingRegions\tmissingRegions\n"  
            out.write(new_header)

            # get column positions  
            columns = mixcr_line.split('\t')
            pos = 0
            for feature in columns:
                mixcr_positions[feature] = pos
                pos += 1
        # extract results for each clone
        if mixcr_count > 1:
            data = mixcr_line.split('\t')
            m_cloneId = data[mixcr_positions["cloneId"]]
            m_cloneFraction = data[mixcr_positions["cloneFraction"]]
            m_cloneCount = data[mixcr_positions["cloneCount"]]

            # add clones to list to make sure ending up with same # of results
            mixcr_clones.append(m_cloneId)

            # Parse seq_info file to get missing regions for each clone
            with open(seq_info_file,'r') as seq_handle:
                seq_info_counter = 0
                for seq_info_line in seq_handle:
                    seq_info_counter += 1
                    seq_info_line = seq_info_line.strip("\n")
                    if seq_info_counter == 1:
                        seq_info_fields_positions = 0
                        seq_info_columns = seq_info_line.split("\t")
                        seq_info_fields = {}
                        for field in seq_info_columns:
                            seq_info_fields[field] = seq_info_fields_positions
                            seq_info_fields_positions += 1
                    if seq_info_counter > 1:
                        seq_info_line_values = seq_info_line.split("\t")
                        s_cloneId = seq_info_line_values[seq_info_fields["cloneId"]]
                        s_cloneFraction = seq_info_line_values[seq_info_fields["cloneFraction"]]
                        s_cloneCount = seq_info_line_values[seq_info_fields["cloneCount"]]
                        s_numMissing = seq_info_line_values[seq_info_fields["numMissing"]]
                        s_regionsMissing = seq_info_line_values[seq_info_fields["regionsMissing"]]
                        if s_cloneId == m_cloneId and s_cloneFraction == m_cloneFraction and s_cloneCount == m_cloneCount:
                            break

            # Extract alignment data for each clone
            with open(db_file, 'r') as handle:
                blast_lines = 0
                for line in handle:
                    blast_lines += 1
                    #print(f"lines : {blast_lines}")
                    if blast_lines == 1:
                        vdj_check = 'no'
                        continue

                    # The query line will have cloneId, cloneFraction, and cloneCount to make sure it matches up with MiXCR results
                    if line.startswith("# Query:"):
                        line = line.strip('\n')
                        line = line[9:]
                        clone_info = line.split('_')
                        b_cloneId = clone_info[1]
                        b_cloneFraction = clone_info[3]
                        b_cloneCount = clone_info[5]

                        # Set up variables used downstream
                        cloneMatch = 'F'
                        v_identity_line = None
                        vdj_check = 'no'
                        #v_allele_pos = None
                        #productive_pos = None

                    # Check if igblastn query results are for corresponding MiXCR clonotype
                    if b_cloneId == m_cloneId and b_cloneFraction == m_cloneFraction and b_cloneCount == m_cloneCount:
                        cloneMatch = 'T'
                    else:
                        cloneMatch = 'F'

                    # VDJ rearrangement summary will have top V match, top D match, and top J match plus whether in frame or stop codon
                    # Save the line number of the next line as a variable because it will contain the V gene information
                    if 'rearrangement summary' in line and cloneMatch == 'T':
                        vdj_summary = blast_lines + 1
                        vdj_check = 'yes'
                        info_line = line.strip('\n')
                        info_line = info_line.strip(".  Multiple equivalent top matches, if present, are separated by a comma.").split(' (')
                        info_split = info_line[1].split(', ')
                        info_position = 0
                        for info in info_split:
                            if info == 'Top V gene match' or info == 'Top V gene match)':
                                v_allele_pos = info_position
                            if info == 'Productive' or info == 'Productive)':
                                productive_pos = info_position
                            if info == 'Top J gene match' or info == 'Top J gene match)':
                                j_allele_pos = info_position
                            info_position += 1

                    if vdj_check == 'yes':
                        if blast_lines == vdj_summary:
                            vdj_line = line.split('\t')
                            # Top V genes returned by igblastn
                            b_v_allele = vdj_line[v_allele_pos]
                            b_productive = vdj_line[productive_pos]
                            b_j_allele = vdj_line[j_allele_pos]
                            # turn into list to double check when extracting sequence identity %
                            if ',' in b_v_allele:
                                b_v_alleles_list = b_v_allele.split(',')
                            else: 
                                b_v_alleles_list = b_v_allele

                        elif line.startswith("# Fields:"):
                            field_columns = line.strip("# Fields: ").strip("\n").split(", ")
                            field_position = 0
                            hit_positions = {}
                            for info in field_columns:
                                field_position += 1
                                hit_positions[info] = field_position

                        # Get sequence alignment, productive status, and mutations for each allele
                        elif line.startswith('V\t') and m_cloneFraction in line:
                            info_line = line.strip('\n').split('\t')
                            
                            # check that v allele matched top v hits
                            #if info_line[hit_positions["subject id"]] in b_v_alleles_list:

                            # get relevant alignment info
                            v_subject_id =  info_line[hit_positions["subject id"]]
                            v_identity = info_line[hit_positions["% identity"]]
                            v_mutations = info_line[hit_positions["BTOP"]]
                            v_bit_score = info_line[hit_positions["bit score"]]
                            v_eval = info_line[hit_positions["evalue"]]

                            # set mutation status
                            if float(v_identity) < 98:
                                mutation_status = "MUTATED"
                            else:
                                mutation_status = "UNMUTATED"

                            if v_identity_line is not None:
                                v_subject_id_line = v_subject_id_line + f",{v_subject_id}"
                                v_identity_line = v_identity_line + f",{v_identity}"
                                v_mutations_line = v_mutations_line + f",{v_mutations}"
                                mutation_status_line = mutation_status_line + f",{mutation_status}"
                                v_bit_score_line = v_bit_score_line + f",{v_bit_score}"
                                v_eval_line = v_eval_line + f",{v_eval_line}"
                            else:
                                v_subject_id_line = f"{v_subject_id}"
                                v_identity_line = f"{v_identity}"
                                v_mutations_line = f"{v_mutations}"
                                mutation_status_line = f"{mutation_status}"
                                v_bit_score_line = f"{v_bit_score}"
                                v_eval_line = f"{v_eval}"

                        if not blast_lines == 1 and (line.startswith("# IGBLASTN") or line.startswith("Total queries")):

                            v_subject_id_line_list = v_subject_id_line.split(",")
                            bit_scores = v_bit_score_line.split(",")
                            v_identity_line_list = v_identity_line.split(",")
                            evalues = v_eval_line.split(",")

                            if "," in b_v_allele:
                                # Extract top V allele and alignment using bit score and eval if there are multiple "top" hits

                                bit_scores_numeric = []
                                evalues_numeric = []
                                for score in bit_scores:
                                    bit_scores_numeric.append(float(score))
                                for eval_score in evalues:
                                    evalues_numeric.append(float(eval_score))
                                
                                # Sort positions based on descending bitscore and then based on ascending e-value
                                top_hits = []
                                # order of bit scores should be in same order as other relevant info
                                for position in range(len(bit_scores)):
                                    top_hits.append([position, bit_scores_numeric[position], evalues_numeric[position]])
                                # sort by descending bit score and ascending e-value
                                sorted_scores = sorted(top_hits, key = lambda x: (-x[1], x[2]))
                                final_bit_position = sorted_scores[0][0]

                                igblast_top_bit_score = bit_scores[final_bit_position]
                                igblast_top_e_val = evalues[final_bit_position]
                                igblast_top_v_gene = v_subject_id_line_list[final_bit_position]
                                igblast_top_v_alignment = v_identity_line_list[final_bit_position]

                            else:
                                top_v_allele_position = v_subject_id_line_list.index(b_v_allele)

                                igblast_top_v_gene = v_subject_id_line_list[top_v_allele_position]
                                igblast_top_bit_score = bit_scores[top_v_allele_position]
                                igblast_top_e_val = evalues[top_v_allele_position]
                                igblast_top_v_alignment = v_identity_line_list[top_v_allele_position]

                            if float(igblast_top_v_alignment) < 98:
                                igblast_top_mutated_status = "MUTATED"
                            else:
                                igblast_top_mutated_status = "UNMUTATED"

                            new_line = mixcr_line + f"\t{b_productive}\t{igblast_top_v_gene}\t{igblast_top_v_alignment}\t{igblast_top_mutated_status}\t{igblast_top_bit_score}\t{igblast_top_e_val}\t{v_subject_id_line}\t{v_identity_line}\t{mutation_status_line}\t{v_bit_score_line}\t{v_eval_line}\t{v_mutations_line}\t{b_j_allele}\t{b_v_allele}\t{s_numMissing}\t{s_regionsMissing}\n"
                            out.write(new_line)

                            igblastn_clones.append(b_cloneId)
                    
                    elif m_cloneId == b_cloneId and vdj_check =='no' and line.startswith("# 0 hits found"):
                        igblastn_clones.append(b_cloneId)
                        new_line = mixcr_line + f"\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\t{s_numMissing}\t{s_regionsMissing}\n"
                        out.write(new_line)

    # Sanity checks to make sure same # of clones in results
    if len(igblastn_clones) != len(mixcr_clones):
        print(f"Different numbers of MiXCR clones and igblast clones for " + str(mixcr_results) + f"... double check results\nigblastn clones: {igblastn_clones}\nmixcr clones: {mixcr_clones}")
    else:
        # Makes sure that lists are not empty in the case of empty MiXCR results/empty igblastn results
        if mixcr_clones != [] and igblastn_clones != []:
            for item in igblastn_clones:
                if item in mixcr_clones:
                    check = 't'
                else:
                    check = 'n'
            if check == 'n':
                print("Different clone numbers being used by MiXCR and Igblast... double check results")
    
    out.close()
    mixcr.close()

if __name__ == "__main__":
    args = parse_args()
    main(args)