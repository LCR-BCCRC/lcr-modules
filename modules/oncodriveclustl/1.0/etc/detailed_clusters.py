#!/usr/bin/env python

import argparse
import gzip
import csv
import numpy as np
import pandas as pd
import copy

def parse_arguments():
    parser = argparse.ArgumentParser()

    parser.add_argument('-e', '--elements',
                        type=str, required=True,
                        help='OncodriveCLUSTL elements results file')
    parser.add_argument('-c', '--clusters',
                        type=str, required=True,
                        help='OncodriveCLUSTL clusters results file')
    parser.add_argument('-q', '--q_value',
                        type=float, default=0.01,
                        help='Analytical q-value threshold for OncodriveCLUSTL elements (not clusters)')
    parser.add_argument('-n', '--samples',
                        type=int, default=5,
                        help='Minimum sample threshold for OncodriveCLUSTL clusters')
    parser.add_argument('-s', '--score',
                        type=int, required=False,
                        default=None,
                        help='Minimum score threshold for OncodriveCLUSTL clusters')
    parser.add_argument('-o', '--output',
                        type=str, required=True,
                        help='Output file name for OncodriveCLUSTL genomic coordinates file')

def filter_results(elements, clusters, q_val, samples, score=None):
    # Filter the files
    filtered_elements = elements[~elements['Q_ANALYTICAL'].isna()]
    filtered_elements = filtered_elements[filtered_elements['Q_ANALYTICAL'] <= q_val]

    # Only include clusters from the filtered elements
    filtered_clusters = clusters[clusters['SYMBOL'].isin(filtered_elements['SYMBOL'])]
    filtered_clusters = filtered_clusters[filtered_clusters['N_SAMPLES'] >= samples]

    if score is not None:
        filtered_clusters = filtered_clusters[filtered_clusters['SCORE'] >= score]

    return filtered_clusters

def reformat_clusters(clusters):
    clusters = clusters[["SYMBOL","CHROMOSOME","COORDINATES"]]

    # Number the clusters beforehand
    cluster_genes = {}
    hotspot_ids = []
    for index, row in clusters.iterrows():
        gene = row['SYMBOL']
        if gene not in list(cluster_genes):
            cluster_genes[gene] = 1
        else:
            cluster_genes[gene] = cluster_genes[gene] + 1
        hotspot_id = f'{gene}_{str(cluster_genes[gene])}'
        hotspot_ids.append(hotspot_id)

    clusters['HOTSPOT_ID'] = hotspot_ids

    # Expand clusters that span over introns
    expanded_clusters = clusters.copy()
    expanded_clusters = expanded_clusters['COORDINATES'].str.split(';',expand=True).stack().reset_index(level=1, drop=True).reset_index()
    expanded_clusters.columns = ['index','COORDINATES']

    # Merge back into original clusters df
    merged_clusters = pd.merge(clusters.drop(['COORDINATES'], axis=1), expanded_clusters, left_index=True, right_on='index')
    merged_clusters = merged_clusters[['SYMBOL','CHROMOSOME','COORDINATES']]
    merged_clusters.reset_index(drop=True, inplace=True)

    # Expand COORDINATES column to include all coordinates in the hotspot
    expanded_coordinates['COORDINATES'] = merged_clusters['COORDINATES'].apply(lambda x: [int(i) for i in x.split(',')])
    expanded_coordinates['COORDINATES'] = expanded_coordinates['COORDINATES'].apply(lambda x: list([range(x[0], x[1]+1)) if len(x) > 1 else [x[0]])

    clusters_expanded = expanded_coordinates.explode('COORDINATES')
    clusters_expanded.reset_index(drop=True, inplace=True)

    return clusters_expanded

def main(args):
    # Read in the files
    if args['elements'].endswith(".gz"):
        elements_df = pd.read_csv(args['elements'], compression='gzip', sep='\t', comment='#')
    else:
        elements_df = pd.read_csv(args['elements'], sep='\t', comment='#')
    
    if args['clusters'].endswith(".gz"):
        clusters_df = pd.read_csv(args['clusters'], compression='gzip', sep='\t', comment='#')
    else:
        clusters_df = pd.read_csv(args['clusters'], sep='\t', comment='#')

    # Filter elements and clusters by q-value and samples
    clusters_filtered = filter_results(elements=elements_df, clusters=clusters_df, q_val=args['q_value'], samples=args['samples'], score=args['score'])

    # Convert clusters dataframe to genomic coordinates style
    clusters_expanded = reformat_clusters(clusters_filtered)

    # Write output file
    clusters_expanded.to_csv(args['output'], na_rep="NA", index=False)

if __init__ == '__main__':
    args = parse_arguments()
    main(args)
