#!/usr/bin/env python3

# TODO: only one reference sample, so PCR plates with
# different cell lines aren't particularly ideal currently

import argparse
import sys
from pathlib import Path
from typing import List

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

class Arguments:
    input_files: List[str]
    output_file: str | None
    force: bool
    input_paths: List[Path]
    output_path: Path
    reference_sample: str
    reference_target: str


parser = argparse.ArgumentParser(
    prog='foldch',
    description='''
        Fold gene expression analysis, for Bio-Rad CFX qPCR output.
    '''
)
parser.add_argument(
    'input_files',
    type=str,
    nargs='+',
    help='''
        input data, either CSV or Excel spreadsheets (.xls/.xlsx);
        if input is an Excel spreadsheet, the first sheet is used.
    '''
)
parser.add_argument('--force', '-f', action='store_true', help='whether to overwrite existing output file')
parser.add_argument('--output_file', '--out', '-o', type=str, required=False, help='path to output file (defaults to "Analysis - [input_file].xlsx")')
parser.add_argument('--reference_sample', '--rs', type=str, default='H1975 par', help='biological control sample (e.g. parental cell line)')
parser.add_argument('--reference_target', '--rt', type=str, default='28S', help='control target (e.g. housekeeping gene)')


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


    def get_group(df: pd.DataFrame, sample: str, target: str) -> pd.DataFrame:
        return df.groupby(['Sample', 'Target']).get_group((sample, target))

    
    output_df = input_df[['Sample', 'Target']].drop_duplicates()
    
    # Calculated from input
    output_df['Ct Avg'] = output_df.apply(lambda row: get_group(input_df, row['Sample'], row['Target'])['Cq'].mean(), axis=1)
    output_df['Ct Std'] = output_df.apply(lambda row: get_group(input_df, row['Sample'], row['Target'])['Cq'].std(), axis=1)
    output_df['Ct Avg Ctrl'] = output_df.apply(lambda row: get_group(input_df, row['Sample'], args.reference_target)['Cq'].mean(), axis=1)
    output_df['Ct Std Ctrl'] = output_df.apply(lambda row: get_group(input_df, row['Sample'], args.reference_target)['Cq'].std(), axis=1)

    # Calculated from intermediate values
    output_df['Delta Ct'] = output_df['Ct Avg'] - output_df['Ct Avg Ctrl']
    output_df['Delta Ct Std'] = np.sqrt((output_df['Ct Std'] ** 2) + (output_df['Ct Std Ctrl'] ** 2))
    
    output_df['Delta Ct Avg Ctrl Sample'] = output_df.apply(lambda row: get_group(output_df, args.reference_sample, row['Target'])['Delta Ct'].mean(), axis=1)
    output_df['Delta Delta Ct'] = output_df['Delta Ct'] - output_df['Delta Ct Avg Ctrl Sample']
    
    output_df['Fold'] = 2 ** (-output_df['Delta Delta Ct'])
    output_df['Fold CI 68 Lower'] = 2 ** (-output_df['Delta Delta Ct'] - output_df['Delta Ct Std'])
    output_df['Fold CI 68 Upper'] = 2 ** (-output_df['Delta Delta Ct'] + output_df['Delta Ct Std'])
    output_df['Fold CI 68'] = output_df.apply(lambda row: f'[{row['Fold CI 68 Lower']:.2g}, {row['Fold CI 68 Upper']:.2g}]', axis=1)

    output_df['Log Fold'] = np.log2(output_df['Fold'])
    output_df['Log Fold CI 68 Lower'] = np.log2(output_df['Fold CI 68 Lower'])
    output_df['Log Fold CI 68 Upper'] = np.log2(output_df['Fold CI 68 Upper'])
    output_df['Log Fold CI 68'] = output_df.apply(lambda row: f'[{row['Log Fold CI 68 Lower']:.2g}, {row['Log Fold CI 68 Upper']:.2g}]', axis=1)
    
    simple_output = output_df[['Sample', 'Target', 'Fold', 'Fold CI 68', 'Log Fold', 'Log Fold CI 68']]
    simple_output = simple_output.rename(columns={'Fold': 'Rel Expr', 'Log Fold': 'Log Rel Expr'})
    simple_output = simple_output.loc[simple_output['Sample'] != args.reference_sample]
    
    output_df.to_excel(args.output_path)
    
    print(simple_output)
    print('Reference sample:', args.reference_sample, '; reference target:', args.reference_target)
    
    for (sample,), group in output_df.groupby(['Sample']):
        if sample == args.reference_sample: continue
        
        group = group.reset_index()
        group = group.loc[group['Target'] != args.reference_target]
        
        errors = group.apply(lambda entry: [entry['Fold'] - entry['Fold CI 68 Lower'], entry['Fold CI 68 Upper'] - entry['Fold']], axis=1, result_type='expand').transpose()
        
        plt.figure()
        plt.suptitle('Relative expression')
        plt.title(sample)
        plt.bar(group['Target'], group['Fold'], edgecolor='black', color='none')
        plt.errorbar(group['Target'], group['Fold'], yerr=errors, fmt='o', color='red', capsize=3, markersize=5)
        plt.axhline(y=1.0, color='gray', linestyle=(0, (5, 5))).set_linewidth(0.5)
    
    plt.show()


if __name__ == "__main__":
    args = parser.parse_args(None, Arguments)
    
    args.input_paths = [Path(file).resolve() for file in args.input_files]
    if args.output_file is not None:
        args.output_path = Path(args.output_file)
    else:
        if len(args.input_paths) == 1:
            first_input = args.input_paths[0]
            args.output_path = Path(f'Analysis - {first_input.stem}.xlsx')
        else:
            args.output_path = Path('Analysis.xlsx')

        print(f'No output file specified, defaulting to {args.output_path}')
    
    nonexistent = list(filter(lambda path: not path.exists(), args.input_paths))
    if any(nonexistent):
        print(f'Some input files do not exist ({', '.join(str(path) for path in nonexistent)})', file=sys.stderr)
        sys.exit(1)
    
    if args.output_path.exists():
        if args.force:
            args.output_path.unlink()
        else:
            print(f'Output file already exists ({args.output_path})', file=sys.stderr)
            sys.exit(1)
    
    main(args)
