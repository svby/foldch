#!/usr/bin/env python3

# TODO: only one reference sample, so PCR plates with
# different cell lines aren't particularly ideal currently

import argparse
import sys
from pathlib import Path
from typing import List

import pandas as pd
import numpy as np

class Arguments:
    input_paths: List[Path]
    output_path: Path

    input_files: List[str]
    force: bool
    output_file: str | None
    reference_sample: str
    reference_target: str
    concat: str


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
parser.add_argument(
    '--concat', '-c',
    type=str,
    default='input',
    choices=['input', 'output'],
    help='''
        specify how to process multiple input files.
        if argument is 'input', inputs will be concatenated and treated as a single experiment;
        if argument is 'output', input files will be treated as separate experiments and the results aggregated for chart generation.
    '''
)


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


def foldch(input_df: pd.DataFrame) -> pd.DataFrame:
    output_df = input_df[['Sample', 'Target']].drop_duplicates()


    def get_group(df: pd.DataFrame, sample: str, target: str) -> pd.DataFrame:
        return df.groupby(['Sample', 'Target']).get_group((sample, target))

    
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
    
    output_df = output_df.loc[output_df['Sample'] != args.reference_sample]
    output_df = output_df.loc[output_df['Target'] != args.reference_target]
    
    return output_df
    

def main(args: Arguments) -> None:
    output_df: pd.DataFrame
    
    if args.concat == 'input':
        output_df = foldch(concat_inputs(args.input_paths))
    elif args.concat == 'output':
        output_df = pd.concat([foldch(read_input(input)) for input in args.input_paths], axis=0)
    else:
        raise ValueError(f'invalid concat strategy {args.concat}')
    
    simple_output = output_df[['Sample', 'Target', 'Fold', 'Fold CI 68', 'Log Fold', 'Log Fold CI 68']]
    simple_output = simple_output.rename(columns={'Fold': 'Rel Expr', 'Log Fold': 'Log Rel Expr'})
    simple_output = simple_output.loc[simple_output['Sample'] != args.reference_sample]
    
    print(simple_output)
    print('Reference sample:', args.reference_sample, '; reference target:', args.reference_target)
    
    if args.output_path.suffix == '.xlsx' or args.output_path.suffix == '.xls':    
        output_df.to_excel(args.output_path)
    else:
        output_df.to_csv(args.output_path)
        
    print(f'Results saved to {args.output_path}')


if __name__ == "__main__":
    args = parser.parse_args(None, Arguments)
    
    args.input_paths = [Path(file).resolve() for file in args.input_files]
    if args.output_file is not None:
        args.output_path = Path(args.output_file)
    else:
        if len(args.input_paths) == 1:
            first_input = args.input_paths[0]
            args.output_path = Path(f'Analysis - {first_input.stem}.csv')
        else:
            args.output_path = Path('Analysis.csv')

        print(f'No output file specified, defaulting to {args.output_path}')
    
    nonexistent = [path for path in args.input_paths if not path.exists()]
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
