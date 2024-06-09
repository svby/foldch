#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path
from typing import List

import pandas as pd
import matplotlib.pyplot as plt

class Arguments:
    input_files: List[str]
    input_paths: List[Path]

    # reference_sample: str
    # reference_target: str


parser = argparse.ArgumentParser(
    prog='barchart',
    description='''
        Generates charts for foldch output
    '''
)
parser.add_argument(
    'input_files',
    type=str,
    nargs='+',
    help='input data (foldch output)'
)
# parser.add_argument('--force', '-f', action='store_true', help='whether to overwrite existing output file')
# parser.add_argument('--reference_sample', '--rs', type=str, default='H1975 par', help='biological control sample (e.g. parental cell line)')
# parser.add_argument('--reference_target', '--rt', type=str, default='28S', help='control target (e.g. housekeeping gene)')


def read_input(path: Path) -> pd.DataFrame:
    if path.suffix == '.xlsx' or path.suffix == '.xls':
        return pd.read_excel(path)
    elif path.suffix == '.csv':
        return pd.read_csv(path)
    else:
        print(f'Input file ({path.name}) is of unknown type (neither .xlsx/.xls nor .csv)', file=sys.stderr)
        sys.exit(1)
        

def concat_inputs(paths: List[Path]) -> pd.DataFrame:
    return pd.concat([read_input(path) for path in paths], axis=0)
    

def main(args: Arguments) -> None:
    input_df = concat_inputs(args.input_paths)
    
    for (sample,), group in input_df.groupby(['Sample']):
        group = group.reset_index()
        
        errors = group.apply(lambda entry: [entry['Fold'] - entry['Fold CI 68 Lower'], entry['Fold CI 68 Upper'] - entry['Fold']], axis=1, result_type='expand').transpose()
        
        plt.figure()
        plt.suptitle('Relative gene expression')
        plt.title(f'Sample: {sample}')
        plt.bar(group['Target'], group['Fold'], edgecolor='black', color='none')
        plt.errorbar(group['Target'], group['Fold'], yerr=errors, fmt='o', color='red', capsize=3, markersize=5)
        plt.axhline(y=1.0, color='gray', linestyle=(0, (5, 5))).set_linewidth(0.5)
    
    plt.show()


if __name__ == "__main__":
    args = parser.parse_args(None, Arguments)
    
    args.input_paths = [Path(file).resolve() for file in args.input_files]

    nonexistent = [path for path in args.input_paths if not path.exists()]
    if any(nonexistent):
        print(f'Some input files do not exist ({', '.join(str(path) for path in nonexistent)})', file=sys.stderr)
        sys.exit(1)
    
    main(args)
