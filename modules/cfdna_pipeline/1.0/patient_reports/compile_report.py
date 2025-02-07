import os
import argparse
import papermill as pm

def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_notebook',required=True,type=str,help='name of template notebook'),
    parser.add_argument('--out_notebook',required=True,type=str,help='name of compiled notebook'),
    parser.add_argument('--maf_files',required=True,type=str,nargs='+',help='list of paths to maf files of sage calls for patient'),
    parser.add_argument('--samplesheet_path',required=True,type=str,help='path to latest samplesheet'),
    parser.add_argument('--patient_id',required=True,type=str,help='patient ID'),
    parser.add_argument('--hs_metrics',required=True,type=str,nargs="+",help='path to file with target HS metrics'),
    parser.add_argument('--targ_cov',required=True,type=str,nargs="+",help='file path to target coverage file'),
    parser.add_argument('--repo_path',required=True,type=str,help='path to local pipeline repo'),
    parser.add_argument('--completion_files',required=False,type=str,nargs='+',help='paths to files containing completion times for each sample'),
    parser.add_argument('--lymphgen_output',required=False,type=str,nargs='+',help='output from lymphgen'),

    return parser.parse_args()

def call_papermill(in_notebook:str, out_notebook:str, args: argparse.Namespace) -> None:
    """Call papermill to run the notebook"""
    pm.execute_notebook(in_notebook, out_notebook,
    parameters=vars(args)) # turn to dict

def main():
    args = get_args()
    call_papermill(args.in_notebook, args.out_notebook, args)


if __name__ == '__main__':
    main()
