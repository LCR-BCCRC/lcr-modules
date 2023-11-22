#!/usr/bin/env python

import argparse
import csv
from Bio.PDB import PDBParser
import operator
import itertools as it
import gzip
import os
import sys
import numpy as np
from math import sqrt
import copy

def parse_arguments():
    parser = argparse.ArgumentParser()
    
    parser.add_argument('-s', '--hotspots',
                        type=str, required=True,
                        help='Hotspot results file')
    parser.add_argument('-m', '--mupit-annotation',
                        type=str, required=True,
                        help='Mupit annotation file')
    parser.add_argument('-t', '--output-merged',
                        type=str, required=True,
                        help='Output merged file created after running hotspot.py')
    parser.add_argument('-c', '--mtc-file',
                        type=str, required=True,
                        help='Multiple test correction file created by multiple_testing_correction.py')
    parser.add_argument('-q', '--q-value',
                        type=float, required=True,
                        help='Q-value used to create hotspot output')
    parser.add_argument('-p','--pdb-info',
                        type=str, required=True,
                        help='File containing PDB info and paths to PDB structures')
    parser.add_argument('-a','--angstroms',
                        type=float, required=True,
                        help='Radius used to find neighbor residues.')
    parser.add_argument('-x', '--metadata-out',
                        type=str, required=True,
                        help='Directory to output results file.')
    parser.add_argument('-g', '--gene',
                        type=str, default=None,
                        help='Verbose output for selected gene')
    parser.add_argument('-y','--maf-mode',
                        action='store_true',
                        help='Write output format in MAF mode')
    args = parser.parse_args()
    args = vars(args)
    return args


# Read in hotspot results
def read_hotspot_file(hotspots,gene_filter=None):
    hotspot_results = []
    with open(hotspots) as handle:
        results = csv.reader(handle, delimiter='\t')
        for line in results:
            hotspots_count = 0
            ttype = line[1]
            if ttype == "REF":
                continue
            gene = line[0]
            if gene_filter != None:
                if gene != gene_filter:
                    continue
                else:
                    print(f"Found entry match for gene {gene_filter}")
            residues_long = line[2:]
            for rs in residues_long:
                hotspots_count += 1
                hotspot = [gene, str(hotspots_count)]
                residues_medium = rs.split(';')
                for transcript_residue in residues_medium:
                    #residue_pos = residue.split(':')[1]
                    #transcript = residue.split(':')[0]
                    #if not transcript in hotspot:
                    #    if len(hotspot) == 2:
                    #        hotspot.append(transcript)
                    #    else if len(hotspot) > 2:
                    #        old_transcript = hotspot[2]
                    #        new_transcript = old_transcript + ":" + transcript
                    #        hotspot[2] = new_transcript
                    hotspot.append(transcript_residue)
                hotspot_results.append(hotspot)

    return hotspot_results


def read_mupit_file(annotations):
    mupit_reverse_dict = {}
    mupit_dict = {}

    with open(annotations) as handle:
        myreader = csv.reader(handle, delimiter='\t')
        header = next(myreader)
        pdb_ix = header.index('pdb_id')
        chain_ix = header.index('chain')
        residue_ix = header.index('residue')
        transcript_ix = header.index('Reference Transcript')
        chromosome_ix = header.index('Chromosome')
        sample_ix = header.index('Sample ID')
        gene_ix = header.index('HUGO symbol')
        ref_res_ix = header.index('Reference Codon Position')
        ref_gen_ix = header.index('Reference Genomic Position')
        reference_ix = header.index('Reference base(s)')

        for line in myreader:
            if not (line[pdb_ix], line[gene_ix], line[chain_ix], line[residue_ix]) in mupit_reverse_dict.keys():
                mupit_reverse_dict[(line[pdb_ix], line[gene_ix], line[chain_ix], line[residue_ix])] = [line[ref_res_ix], line]
            pdb = line[pdb_ix]
            gene = line[gene_ix]
            chain = line[chain_ix]
            ref_codon = str(line[ref_res_ix])
            position = line[ref_gen_ix]
            key = (pdb, gene, chain, ref_codon, position)
            if not key in mupit_dict.keys():
                mupit_dict[key] = [line]
            else:
                mupit_dict[key].append(line)
    # Check that residues are the same
    residue_dicts = {k: [mupit_dict[k][x][ref_res_ix] for x in range(len(mupit_dict[k]))] for k in mupit_dict.keys()}
    for k, v in residue_dicts.items():
        if len(set(v)) > 1:
            key_pdb = k[0]
            key_gene = k[1]
            key_chain = k[2]
            key_codon = k[3]
            key_position = k[4]
            print(f"Issue reading MuPIT file: There is conflicting structure residue positions at lines matching the following column values:\nPDB: {key_pdb}\tGENE: {key_gene}\tCHAIN: {key_chain}\tCODON: {key_codon}\tGENOMIC_POS: {key_position}")
            print(f"Sets of residue_ix: {set(v)}")
    return (mupit_dict, mupit_reverse_dict, pdb_ix, chain_ix, residue_ix, sample_ix, gene_ix, ref_res_ix, reference_ix, ref_gen_ix, chromosome_ix)


def read_merged_output(outmerged, sig_level):
    #long_output = [['Structure', 'Model', 'Chain', 'Residue Position', 'p-value']]
    # Make dictionary with format: {(residue, p-value): [all columns]}
    long_output = {}
    input_path = os.path.abspath(outmerged)
    output_dir = os.path.dirname(input_path)
    output_file = output_dir + "/long_output_merged.txt"
    outfile = open(output_file, 'w')
    header = ['Structure', 'Model', 'Chain', 'Structure Residue', 'Residue Mutation Count', 'Mutation Density', 'Hotspot P-value']
    outfile.write('\t'.join(header) + "\n")

    with open(outmerged) as handle:
        myreader = csv.reader(handle, delimiter='\t')
        header = next(myreader)
        pdb_ix = header.index('Structure')
        model_ix = header.index('Model')
        chain_ix = header.index('Chain')
        res_ix = header.index('Mutation Residues')
        count_ix = header.index('Residue Mutation Count')
        density_ix = header.index('Mutation Density')
        p_ix = header.index('Hotspot P-value')
        for line in myreader:
            if not line[p_ix]:
                continue
            structure = line[pdb_ix]
            chains = line[chain_ix].split(',')
            models = line[model_ix].split(',')
            res_pos = line[res_ix].split(',')
            res_pval = map(float, line[p_ix].split(','))
            count = line[count_ix].split(',')
            mut_density = line[density_ix].split(',')
            #if structure=="1w72":
                #print(f"Structure: {structure}")
                #print(f"Number of chains: {len(chains)}")
                #print(f"Number of residue positions: {len(res_pos)}")
                #print(f"Residue positions: {res_pos}")
            for i, p in enumerate(res_pval):
                if p <= 1.1:
                    out_key = (structure, chains[i], res_pos[i], p)
                    out_list = [line[pdb_ix], models[i], chains[i], res_pos[i], count[i], mut_density[i], p]
                    long_output[out_key] = out_list
                    #if structure=="1w72":
                        #print(out_key)
                        #print(long_output[out_key])
                    out_file_list = [line[pdb_ix], models[i], chains[i], str(res_pos[i]), str(count[i]), str(mut_density[i]), str(p)]
                    outfile.write('\t'.join(out_file_list) + "\n")
    outfile.close()
    return long_output


def read_mtc_output(mtc):
    output = {}
    with open(mtc) as handle:
        myreader = csv.reader(handle, delimiter='\t')
        header = next(myreader)
        gene_ix = header.index('HUGO Symbol')
        transcript_ix = header.index('Sequence Ontology Transcript')
        res_ix = header.index('CRAVAT Res')
        min_p_ix = header.index('Min p-value')
        q_ix = header.index('q-value')
        g_pos_ix = header.index('genomic position')
        for line in myreader:
            out_key = (line[gene_ix], line[transcript_ix], line[res_ix])
            if out_key not in output.keys():
                output[out_key] = {}
            output[out_key]["min_p"] = str(line[min_p_ix])
            output[out_key]["q_val"] = str(line[q_ix])
            output[out_key]["genomic_position"] = line[g_pos_ix]
            #outline = [line[gene_ix], line[transcript_ix], line[res_ix], line[min_p_ix], line[q_ix]]
            #output.append(outline)
    return output


def distance(Coordinates):
    x, y, z = zip(*Coordinates)
    return sqrt((x[0] - x[1])**2 + (y[0] - y[1])**2 + (z[0] - z[1])**2)


def read_pdb_info(pdb_file):
    pdb_dict = {}
    pdb_gene_dict = {}
    pdb_to_gene_dict = {}

    with open(pdb_file) as handle:
        handle.readline()
        myreader = csv.reader(handle, delimiter='\t')
        for pdbid, lines in it.groupby(myreader, lambda x: x[0]):
            lines = list(lines)
            gene2chain = {}

            # Get list of pdbids associated with different genes to help match values from mtc_output_min and merged_output
            
            #gene = lines[0][2]
            #if not gene in pdb_gene_dict.keys():
            #    pdb_gene_dict[gene] = []
            #    pdb_gene_dict[gene].append(pdbid)
            #else:
            #    pdb_gene_dict[gene].append(pdbid)

            for chain_description, lines_subset in it.groupby(lines, lambda x: x[5]):
                gene2chain[chain_description] = [l[1] for l in lines_subset]
            gene2chain['path'] = lines[0][4]
            pdb_dict[pdbid] = gene2chain
    with open(pdb_file) as handle:
        handle.readline()
        myreader = csv.reader(handle, delimiter='\t')
        for line in myreader:
            pdb = line[0]
            chain = line[1]
            gene = line[2]
            pdb_chain_key = (pdb, chain)

            if not gene in pdb_gene_dict.keys():
                pdb_gene_dict[gene] = []
                pdb_gene_dict[gene].append(pdb)
            else:
                pdb_gene_dict[gene].append(pdb)

            if not pdb_chain_key in pdb_to_gene_dict.keys():
                pdb_to_gene_dict[pdb_chain_key] = gene
            else:
                print(f"More than two genes correspond to pdb {pdb} and chain {chain}!")

    return pdb_dict, pdb_gene_dict, pdb_to_gene_dict


def get_cog(pdb_dict, pdbid):
    pdb_parser = PDBParser(QUIET=True)
    struct_info = pdb_dict[pdbid]

    struct_temp = copy.deepcopy(struct_info)
    pdb_path = struct_temp.pop("path")

    struct_chains = []
    for k in struct_temp.keys():
        struct_chains.extend(struct_temp[k])
    struct_temp = ''

    cog_dict = {}
    if not pdb_path:
        print(f"No pdb_path available for {pdbid}. . . Skipping . . .")
        return None
    try: 
        if pdb_path.endswith(".gz"):
            with gzip.open(pdb_path, 'rt') as handle:
                structure = pdb_parser.get_structure(pdbid, handle)
        else:
            structure = pdb_parser.get_structure(pdbid, pdb_path)
        for model in structure:
            for chain in model:
                if chain.id == " ":
                    chain.id = "A"
                    if not pdbid.startswith("ENSP") and not pdbid.startswith("NP"):
                        del model.child_dict[' ']
                        model.child_dict['A'] = chain
        # Generate center of geometries for atoms
        models = []
        for model in structure:
            models.append(model.get_id())
            for chain in model:
                if chain.get_id() == ' ':
                    chain.id = 'A'
                chain_id = chain.get_id()
                for residue in chain:
                    # This is to ignore hetero residues
                    if not (residue.get_full_id()[3][0] == ' '):
                        continue
                    full_id = residue.get_full_id()
                    center_of_geometry = np.sum(atom.coord for atom in residue) / len(residue)
                    if (full_id[2] in struct_chains):
                        cog_dict[full_id] = center_of_geometry
        return cog_dict
    except KeyboardInterrupt:
        raise
    except:
        print(f"Failed to read structure {pdbid} and generate center of geometries or neighbors.")
        return None


def get_neighbors(cog_dict, residue_id, angstroms):
    neighbors_dict = {}
    neighbors = []
    for residue_2 in cog_dict.keys():
        if (residue_2 == residue_id):
            continue
        tmp_dist = distance([cog_dict[residue_id], cog_dict[residue_2]])
        if tmp_dist <= angstroms:
            neighbors.append(residue_2)
    neighbors_dict[residue_id] = neighbors
    #print("Neighbors are:")
    #print(neighbors_dict)
    return neighbors_dict


def create_neighborhood_dict(pdb, hotspot_gene, residue, residue_neighbors, pdb_to_gene_dict, mupit_dict, mupit_reverse_dict):
    neighborhood = {}
    
    residue_id = list(residue_neighbors)[0]

    chain = residue_id[2]
    chain_pos = str(residue_id[3][1])
    pdb_pos = f"{chain}_{chain_pos}"

    #neighborhood["pdb"] = [pdb]
    #neighborhood["gene"] = [hotspot_gene]
    #neighborhood["residue"] = [residue]
    #neighborhood["chain_pos"] = [pdb_pos]

    # Create multi level dictionary

    neighborhood["pdb"] = pdb
    neighborhood["mutated_residues"] = {}
    neighborhood["mutated_residues"][residue] = {}
    neighborhood["mutated_residues"][residue]["pdb"] = pdb
    neighborhood["mutated_residues"][residue]["gene"] = hotspot_gene
    neighborhood["mutated_residues"][residue]["chain_pos"] = pdb_pos

    samples = []
    for k, v in mupit_dict.items():
        occurences = len(v)
        for i in range(occurences):
            s = v[i]
            #print(s)
            if s[0]==pdb and s[6]==hotspot_gene and s[1]==chain and s[2]==chain_pos and s[8]==residue:
                sample = s[5]
                samples.append(sample)
    if samples != []:
        neighborhood["mutated_residues"][residue]["samples"] = samples

    neighbors = residue_neighbors[residue_id]
    for neighbor in neighbors:
        neighbor_samples = []

        pdb = neighbor[0]
        chain = neighbor[2]
        chain_pos = str(neighbor[3][1])
        pdb_pos = f"{chain}_{chain_pos}"
        gene = pdb_to_gene_dict[(pdb, chain)]

        mupit_reverse_key = (pdb, gene, chain, chain_pos)
        if mupit_reverse_key in mupit_reverse_dict.keys():
            neighbor_res = mupit_reverse_dict[mupit_reverse_key][0]

            #neighborhood["residue"].append(neighbor_res)

            neighborhood["mutated_residues"][neighbor_res] = {}
            neighborhood["mutated_residues"][neighbor_res]["pdb"] = pdb
            neighborhood["mutated_residues"][neighbor_res]["gene"] = gene
            neighborhood["mutated_residues"][neighbor_res]["chain_pos"] = pdb_pos
        else:
            #print(f"No mutations for neighbor residue: {mupit_reverse_key}") 
            continue
        
        # get neighbor samples
        for k, v in mupit_dict.items():
            occurrences = len(v)
            for i in range(occurrences):
                s = v[i]
                if s[0]==pdb and s[6]==gene and s[1]==chain and s[2]==chain_pos and s[8]==neighbor_res:
                    sample = s[5]
                    neighbor_samples.append(sample)
        if neighbor_samples != []:
            neighborhood["mutated_residues"][neighbor_res]["samples"] = neighbor_samples
    
   #print(neighborhood)
    return neighborhood

def write_metadata(neighborhood_dict, metadata_output, verbose=False):

    #mupit_abs_file = os.path.abspath(mupit_mutations)
    #mupit_file = mupit_abs_file + "_detailed"
    #mupit_out = open(mupit_file, 'w')

    with open(metadata_output, 'a') as writer:
        for gene, hotspots in neighborhood_dict.items():
            for hotspot, residues in hotspots.items():
                out_line = []
                out_line.append(gene)
                out_line.append(hotspot)

                # Get list of all residues HotMAPS said were in hotspot
                hotmaps_residues = []
                if verbose:
                    hotmaps_residues_verbose = []
                    gene_position_verbose_col = []
                    point_samples = []

                # Contains all residues involved in hotspot in GENE:POSITION format
                gene_position_col = []
                #gene_position_verbose_col = []

                # All samples with mutations in hotspot
                samples_counter = []
                #point_samples = []

                # List of unique genes that contribute to hostpot
                unique_genes = []

                for residue, description in residues.items():
                    hotspot_pdb = description["pdb"]
                    hotmaps_residues.append(str(residue))

                    if verbose:
                        hotmaps_residues_verbose.append(str(residue) + ":" + hotspot_pdb)

                    for hotspot_residue, residue_info in description["mutated_residues"].items():
                        residue_gene = residue_info["gene"]
                        residue_position = str(hotspot_residue)
                        residue_pdb = residue_info["pdb"]
                        residue_chain_pos = residue_info["chain_pos"]

                        gene_position = residue_position + ":" + residue_gene

                        if verbose:
                            gene_position_verbose = residue_position + ":" + residue_gene + ":" + residue_chain_pos + ":" + hotspot_pdb
                            gene_position_verbose_col.append(gene_position_verbose)

                        gene_position_col.append(gene_position)

                        unique_genes.append(residue_gene)

                        # Get samples that correspond to specific residue associated with hotspot
                        residue_samples = residue_info["samples"]
                        samples_counter.extend(residue_samples)
                        if verbose:
                            point_samples_line = residue + ":" + residue_position + ":" + ','.join(residue_samples)
                            point_samples.append(point_samples_line)

                hotmaps_residues = ','.join(hotmaps_residues)
                unique_genes_col = ','.join(list(set(unique_genes)))
                hotmaps_residue_col = ','.join(list(set(hotmaps_residues)))
                gene_position_col = ','.join(list(set(gene_position_col)))
                if verbose:
                    hotmaps_residue_v_col = ','.join(list(set(hotmaps_residues_verbose)))
                    gene_position_v_col = ','.join(list(set(gene_position_verbose_col)))
                    point_samples = ';'.join(list(set(point_samples)))
                unique_samples = ','.join(list(set(samples_counter)))
                #point_samples = ';'.join(list(set(point_samples)))
                #n_samples = str(len(samples_counter))
                n_samples = str(len(list(set(samples_counter))))

                out_line.append(n_samples)
                out_line.append(unique_genes_col)
                out_line.append(hotmaps_residues)
                out_line.append(gene_position_col)

                if verbose:
                    out_line.append(hotmaps_residue_v_col)
                    out_line.append(gene_position_v_col)
                    out_line.append(point_samples)
                    #print('\t'.join(out_line))
                writer.write('\t'.join(out_line) + "\n")

def write_maf(neighborhood_dict, metadata_output, mupit_reverse_dict, maf_format, mupit_gene_pos_ix, mupit_chromosome_ix):

    #mupit_abs_file = os.path.abspath(mupit_mutations)
    #mupit_file = mupit_abs_file + "_detailed"
    #mupit_out = open(mupit_file, 'w')
    if maf_format:
        with open(metadata_output, 'a') as writer:
            for gene, hotspots in neighborhood_dict.items():
                for hotspot, residues in hotspots.items():
                    #out_line = []
                    #out_line.append(gene)
                    #out_line.append(hotspot)
                    hotspot_gene_combo = gene + "_" + hotspot

                    hotmaps_residues = []
                    gene_position_col = []
                    
                    for residue, description in residues.items():
                        hotspot_pdb = description["pdb"]
                        hotmaps_residue = residue
                        for hotspot_residue, residue_info in description["mutated_residues"].items():
                            residue_gene = residue_info["gene"]
                            residue_position = str(hotspot_residue)
                            residue_pdb = residue_info["pdb"]
                            residue_chain_pos = residue_info["chain_pos"].split("_")[1]
                            residue_chain = residue_info["chain_pos"].split("_")[0]
                            mupit_key = (residue_pdb, residue_gene, residue_chain, residue_chain_pos)
                            mupit_line = mupit_reverse_dict[mupit_key][1]
                            mupit_chromosome = mupit_line[mupit_chromosome_ix]
                            mupit_genomic_pos = mupit_line[mupit_gene_pos_ix]

                            mupit_genomic_position = mupit_genomic_pos.replace("b'", "").replace("'", "").split(",")
                            for genomic_pos in mupit_genomic_position:
                                out_line = [residue_gene, mupit_chromosome, genomic_pos, hotspot_gene_combo, residue_position, residue_pdb, residue_chain, residue_chain_pos, hotmaps_residue]
                                writer.write('\t'.join(out_line) + "\n")


def main(args):
    # Returns hotspots in list within lists format ['gene', 'hotspot_num', 'transcript', 'residues']
    if args['gene'] != None:
        gene_filter = args['gene']
        print(f"Verbose option selected for gene {gene_filter}")
    elif args['gene'] == None:
        gene_filter = None

    print("Loading HotMAPS results...")
    hotspot_results = read_hotspot_file(args['hotspots'], gene_filter)

    # Returns a dictionary of {(pdb, gene, chain, ref codon) : line} each line of mupit file as a SORTED list within a list, and indexes of relevant columns
    print("Loading MuPIT annotations...")
    mupit_stuff = read_mupit_file(args['mupit_annotation'])
    (mupit_dict, mupit_reverse_dict, mupit_pdb_ix, mupit_chain_ix, mupit_residue_ix, mupit_sample_ix, mupit_gene_ix, mupit_ref_res_ix, mupit_ref_aa_ix, mupit_gen_pos_ix, mupit_chromosome_ix) = mupit_stuff

    # Returns long version of output merged as a dictionary, format :
    # {(residue, p-value) : [pdb, model, chain, res_pos, count, mut_density, p]}
    print("Loading output_merged file...")
    long_output_dict = read_merged_output(args['output_merged'], args['q_value'])

    # Returns multiple testing correction file as dictionary with format :
    # [(gene, transcript, residue): {"min_p": x, "q_val":y}]
    print("Loading multiple testing correction file...")
    mtc_dict = read_mtc_output(args['mtc_file'])

    # Read in PDB data as a dictionary with format:
    # '{'pdb_id' : {'chain_description': x, 'path': y}, ...}
    print("Loading PDB structure info . . .")
    pdb_dict, pdb_gene_dict, pdb_to_gene_dict = read_pdb_info(args['pdb_info'])

    if args['gene'] == None:
        verbose = False

    if args['gene'] != None:
        verbose = True

    maf_metadata_out = os.path.abspath(args['metadata_out'])
    exp_metadata_out = maf_metadata_out.replace("maf","expanded")
    
    if verbose == True:
        verbose_metadata_out = maf_metadata_out.replace("maf",f"verbose_{args['gene']}")

    assert not os.path.exists(maf_metadata_out), f"{maf_metadata_out} file already exists. Exiting."
    assert not os.path.exists(exp_metadata_out), f"{exp_metadata_out} file already exists. Exiting." 

    if not verbose:
        with open(exp_metadata_out, 'w') as handle:
                output_header = ['GENE', 'HOTSPOT_NUM', 'N_SAMPLES', 'ALL_GENES', 'HOTMAPS_RES', 'MUTATED_RES']
                handle.write('\t'.join(output_header) + "\n")
        if args['maf_mode']:
            with open(maf_metadata_out, 'w') as handle:
                    output_header = ['Hugo_Symbol', 'Chromosome', 'Start_Position', 'Hotspot_ID', 'Protein_Residue', 'PDB_Structure', 'PDB_Chain', 'PDB_Chain_Position', 'Central_Residue']
                    handle.write('\t'.join(output_header) + "\n")
    # verbose version
    if verbose:
        print("Proceeding with verbose output...")
        with open(verbose_metadata_out, 'w') as handle:
                output_header = ['GENE', 'HOTSPOT_NUM', 'N_SAMPLES', 'ALL_GENES', 'HOTMAPS_RES', 'MUTATED_RES', 'HOTMAPS_RES_V', 'MUTATED_RES_V', 'SAMPLES_V']
                handle.write('\t'.join(output_header) + "\n")

    # Create a hotspot dict to keep track of origin and metadata of hotspot residues
    # hotspot_dict = {'gene': 
    #                       {'hotspot_1': 
    #                           {83: 
    #                               {"pdb": 'pdb', 
    #                                "mutated_residues":
    #                                           {residue:
    #                                             {"pdb":pdb_res,
    #                                              "gene":gene,
    #                                              "chain_pos":
    #                                              "samples":[samples]
    #                                               },
    #                                             {residue:
        #                                          "pdb":pdb_res,
    #                                              "gene":gene,
    #                                              "chain_pos":
    #                                              "samples":[samples]
    #                                               }
    #                                           }
    #                                }
    #                              }
    #                         }, 
    #                       {'hotspot_2': ...}
    #                 }

   #hotspot_dict = {}

    print("Initializing...")
    # Get q-values for hotspots
    for hotspot in hotspot_results:
        hotspot_dict = {}
        hotspot_gene = hotspot[0]
        hotspot_num = hotspot[1]
        #transcript_id = hotspot[2]
        hotspot_residues = hotspot[2:]

        #if hotspot_gene == "IKZF3":
        print(f"---------*----*---*--*-*\nWORKING ON HOTSPOT:\nGENE:     {hotspot_gene}\nHOTSPOT:  {str(hotspot_num)}\nRESIDUES:  {','.join(hotspot_residues)}")

        # fill out dict

        if hotspot_gene not in hotspot_dict.keys():
            hotspot_dict[hotspot_gene] = {}
            hotspot_dict[hotspot_gene][hotspot_num] = {}
        else:
            hotspot_dict[hotspot_gene][hotspot_num] = {}

        # for each residue in hotspot, get min p-value to find what structure the hotspot is based off in the mtc_output results file
        for transcript_residue in hotspot_residues:
            transcript_id = transcript_residue.split(':')[0]
            residue = transcript_residue.split(':')[1]

            #### SKIP RESIDUES THAT HAVE VALUE -1
            if residue=="-1":
                continue

            hotspot_neighbors = None

            # set up residue level in hotspot dict
            hotspot_dict[hotspot_gene][hotspot_num][residue] = {}

            hotspot_key = (hotspot_gene, transcript_id, residue)

            print(f"...\nSearching for hotspot key {hotspot_key} in multiple test correction results... ")

            #if ":" in transcript_id:
            #    transcript_1 = transcript_id.split(":")[0]
            #    transcript_2 = transcript_id.split(":")[1]
            #    potential_hotspot_key_1 = (hotspot_gene, transcript_1, residue)
            #    potential_hotspot_key_2 = (hotspot_gene, transcript_2, residue)
            #    if potential_hotspot_key_1 in mtc_dict.keys():
            #        hotspot_key = potential_hotspot_key_1
            #    elif potential_hotspot_key_2 in mtc_dict.keys():
            #        hotspot_key = potential_hotspot_key_2
            #    else:
            #        print(f"Could not find hotspot key for hotspot {gene}, {residue}")
            if hotspot_key in mtc_dict.keys():

                print(f"Found hotspot key {hotspot_key} in multile test correction results! Value is:\n{mtc_dict[hotspot_key]}")
                res_p = None

                min_p = mtc_dict[hotspot_key]["min_p"]
                q_val = mtc_dict[hotspot_key]["q_val"]
                genomic_pos = mtc_dict[hotspot_key]["genomic_position"]
                #for k in long_output_dict.keys():
                #    if k[0] == '1w72':
                #        print(k)
                #        print(long_output_dict[k])

                # get chain residue
                #(pdb, gene, chain, ref_codon, position)

                print(f"....\nSearching for PDB structures at {hotspot_gene}:{residue} with p-value matching {min_p}...")
                for pdbid in pdb_gene_dict[hotspot_gene]:
                    structure_dict = pdb_dict[pdbid]
                    structure_temp = copy.deepcopy(structure_dict)
                    pdb_path = structure_temp.pop('path')

                    pos_chains = {"chains": []}
                    for k, v in structure_temp.items():
                        chains = structure_temp[k]
                        pos_chains["chains"] = pos_chains["chains"] + chains

                    for chain in pos_chains["chains"]:
                        chain = chain[0]
                        mupit_chain_key = (pdbid, hotspot_gene, chain, residue, genomic_pos)
                        if mupit_chain_key in mupit_dict.keys():
                            #print(f"Mupit Chain Key: {mupit_chain_key}")
                            #print(f"Mupit results: {mupit_dict[mupit_chain_key]}")

                            #chain_residue = (mupit_dict[mupit_chain_key][0][mupit_residue_ix])

                            #### TO FIX THE ISSUE OF MULTIPLE ENTRIES HAVING SAME PROTEIN RESIDUE BUT DIFFERENT CHAIN RESIDUE VALUES
                            possible_chain_residues = []
                            for mupit_line in mupit_dict[mupit_chain_key]:
                                possible_chain_residue = mupit_line[mupit_residue_ix]
                                possible_chain_residues.append(possible_chain_residue)
                            possible_chain_residues = set(possible_chain_residues)

                            #possible_res_p = (pdbid, chain, chain_residue, float(min_p))

                            for pcr in possible_chain_residues:
                                possible_res_p = (pdbid, chain, pcr, float(min_p))
                                ###### THIS IS A TEMPORARY FIX TO THE DNP MUTATION ISSUE WHEN THE AMINO ACID CHANGE OCCURS DUE TO THE SECOND NUCLEOTIDE CHANGE (END POSITION COLUMN)
                                #if hotspot_gene=="GRHPR" and residue=="41" and pdbid =="NP_036335.1_1":
                                    #print(f"Manually overriding chain residue for GRHPR and pdbid {pdbid}")
                                    #possible_res_p = (pdbid, chain, "40", float(min_p))

                                #if hotspot_gene=="BTG2" and residue=="44" and pdbid=="3dju":
                                    #print(f"Manually overriding chain residue for BTG2 and pdbid {pdbid}")
                                    #possible_res_p = (pdbid, chain, "43", float(min_p))

                                #if hotspot_gene=="BTG1" and residue=="43" and pdbid=="NP_001722.1_1":
                                    #print(f"Manually overriding chain residue for BTG1 and pbdid {pdbid}")
                                    #possible_res_p = (pdbid, chain, "44", float(min_p))

                                #print(possible_res_p)
                    #possibl    e_res_p = (pdbid, str(residue), float(min_p))
                                if possible_res_p in long_output_dict.keys():
                                    chain_residue = pcr
                                    print(f"Found PDB match for {hotspot_gene}:{residue}:\tpdb: {pdbid} chain: {chain} chain_residue: {chain_residue} p_value: {min_p}")
                                    res_p = possible_res_p
                                    chain_winner = chain
                                    residue_winner = chain_residue
                                    # added a break to reduce time since there are so many for loops now
                                    break
                        #else:
                            #print(f"Mupit chain key is NOT in mupit_dict: {mupit_chain_key}")
                            #print(f"{list(mupit_dict)[0]}")
                            #print("Testing to see if key exists at all..")
                            #test_key = ('1w72', 'B2M', 'B', '73', "b'45007770,45007771,45007772'")
                            #print(test_key)
                            #print(f"{mupit_dict[test_key]}")
                if res_p is None:
                    print(f"ERROR: Could not find PDB structure match in output_merged for GENE: {hotspot_gene}\tP_VAL: {min_p}\tRESIDUE: {residue} in the structures:\n{pdb_gene_dict[hotspot_gene]}")
                    sys.exit()
                # Find structure corresponding to pdbid-residue-p value pair
                if not res_p is None and res_p in long_output_dict.keys():
                    #print(f"Using long_output_dict values for key: {res_p}\nLong_output_value is: {long_output_dict[res_p]}")
                    output_line = long_output_dict[res_p]
                    pdb = output_line[0]
                    model = output_line[1]
                    chain = output_line[2]
                    mutation_count = output_line[4]
                    density = output_line[5]
                    mupit_key = (pdb, hotspot_gene, chain, residue, genomic_pos)

                    # set up pdb structure that hotspot is based off of
                    hotspot_dict[hotspot_gene][hotspot_num][residue]["pdb"] = pdb

                    print(f".....\nSearching for mupit_key {mupit_key} in mupit dict in order to retrieve chain and chain residue...")

                    if mupit_key in mupit_dict.keys():
                        # Begin filling out neighborhood level of hotspot_dict
                        hotspot_dict[hotspot_gene][hotspot_num][residue]["neighborhood"] = {}

                        # Get the corresponding structure residue from mupit annotation file
                        chain_residue = mupit_dict[mupit_key][0][mupit_residue_ix]
                        mupit_lines = mupit_dict[mupit_key]

                        # Get structure
                        structure_cog = get_cog(pdb_dict, pdb)

                        # Get density and neighbor information
                        residue_full_id = (pdb, int(model), chain, (' ', int(chain_residue), ' '))
                        assert residue_full_id in structure_cog.keys(), f"Full id {residue_full_id} does not match center of geometry dictionary keys."

                        hotspot_neighbors = get_neighbors(structure_cog, residue_full_id, args['angstroms'])
                        print("Identifying residue neighbors...")
                    else:
                        print(f"ERROR: Could not find key {mupit_key} in mupit annotation file...")

                if hotspot_neighbors is not None:
                    assert (len(hotspot_neighbors.keys()) == 1), f"More than one key in hotspot neighbors hotspot at step: {hotspot_gene} ; {hotspot_residue}. Keys: {hotpot_neighbors.keys()}"

                    # Create the neighborhood dictionary
                    print("Finding mutations in neighborhood...")
                    neighborhood_dict = create_neighborhood_dict(pdb, hotspot_gene, residue, hotspot_neighbors, pdb_to_gene_dict, mupit_dict, mupit_reverse_dict)
                else:
                    print(f"ERROR: Could not find key {res_p} in output_merged file...")
                    continue
            else:
                print(f"ERROR: Could not find hotspot {hotspot_key} in multiple testing correction results...")
                continue

            hotspot_dict[hotspot_gene][hotspot_num][residue] = neighborhood_dict

        if args['maf_mode']:
            write_maf(hotspot_dict, maf_metadata_out, mupit_reverse_dict, args['maf_mode'], mupit_gen_pos_ix, mupit_chromosome_ix)
        #if not args['maf_mode']:
        write_metadata(hotspot_dict, exp_metadata_out, verbose)


if __name__ == '__main__':
    args = parse_arguments()
    main(args)

