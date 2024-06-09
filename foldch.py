#!/usr/bin/env python3

import argparse
import sys
import re
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd
import numpy as np

class Arguments:
    input_paths: List[Path]
    output_path: Path

    input_files: List[str]
    force: bool
    output_file: str | None
    reference_samples: List[List[str]] = []
    reference_targets: List[List[str]] = []
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
parser.add_argument('--output-file', '--out', '-o', type=str, required=False, help='path to output file (defaults to "Analysis - [input_file].xlsx")')

parser.add_argument(
    '--reference-sample', '--rs',
    dest='reference_samples',
    type=str,
    action='append',
    metavar=('CONDITION', 'SAMPLE_NAME'),
    nargs=2,
    help='''
        specify an associated biological control sample for input samples (e.g. parental cell lines, or untreated sample) based on a specified condition.
        The condition can either be 'all', a /regular expression/, or a string.
        When the condition is 'all', it matches all samples that did not match earlier conditions.
        When the condition is a regular expression or a string, the name of each sample is matched with it.
        When no condition matches the name of a sample, the first sample in the input file is used as a fallback.
    '''
)
parser.add_argument(
    '--reference_target', '--rt',
    dest='reference_targets',
    type=str,
    action='append',
    metavar=('CONDITION', 'SAMPLE_NAME'),
    nargs=2,
    help='''
        specify an associated control target for input samples (i.e. housekeeping gene) based on a specified condition.
        See also help text for --reference-sample.
    '''
)

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


def get_reference_condition(reference_list: List[Tuple[str, str]], value: str) -> str | None:
    for (condition, ref) in reference_list:
        if condition == 'all':
            return ref
        elif condition.startswith('/') and condition.endswith('/'):
            if re.match(condition[1:-1], value) is not None:
                return ref
        elif condition == value:
            return ref
    return None

def foldch(args: Arguments, input_df: pd.DataFrame) -> pd.DataFrame:
    output_df = input_df[['Sample', 'Target']].drop_duplicates()
    
    reference_samples: Dict[str, str] = {}
    reference_targets: Dict[str, str] = {}
    
    first_row = output_df.iloc[0]
    default_sample = first_row['Sample']
    default_target = first_row['Target']
    
    for (_, row) in output_df.iterrows():
        sample = row['Sample']
        reference_sample = get_reference_condition(args.reference_samples, sample)
        reference_target = get_reference_condition(args.reference_targets, sample)
        
        if reference_sample is None:
            print(f'No reference sample specified for sample {sample}, defaulting to {default_sample}')
            reference_sample = default_sample
        if reference_target is None:
            print(f'No reference target specified for sample {sample}, defaulting to {default_target}')
            reference_target = default_target
        
        reference_samples[sample] = reference_sample
        reference_targets[sample] = reference_target


    def get_group(df: pd.DataFrame, sample: str, target: str) -> pd.DataFrame:
        return df.groupby(['Sample', 'Target']).get_group((sample, target))

    
    output_df['Ref Sample'] = output_df.apply(lambda row: reference_samples[row['Sample']], axis=1)
    output_df['Ref Target'] = output_df.apply(lambda row: reference_targets[row['Sample']], axis=1)
    
    # Calculated from input
    output_df['Ct Avg'] = output_df.apply(lambda row: get_group(input_df, row['Sample'], row['Target'])['Cq'].mean(), axis=1)
    output_df['Ct Std'] = output_df.apply(lambda row: get_group(input_df, row['Sample'], row['Target'])['Cq'].std(), axis=1)
    output_df['Ct Avg Ctrl'] = output_df.apply(lambda row: get_group(input_df, row['Sample'], row['Ref Target'])['Cq'].mean(), axis=1)
    output_df['Ct Std Ctrl'] = output_df.apply(lambda row: get_group(input_df, row['Sample'], row['Ref Target'])['Cq'].std(), axis=1)

    # Calculated from intermediate values
    output_df['Delta Ct'] = output_df['Ct Avg'] - output_df['Ct Avg Ctrl']
    output_df['Delta Ct Std'] = np.sqrt((output_df['Ct Std'] ** 2) + (output_df['Ct Std Ctrl'] ** 2))
    
    output_df['Delta Ct Avg Ctrl Sample'] = output_df.apply(lambda row: get_group(output_df, row['Ref Sample'], row['Target'])['Delta Ct'].mean(), axis=1)
    output_df['Delta Delta Ct'] = output_df['Delta Ct'] - output_df['Delta Ct Avg Ctrl Sample']
    
    output_df['Fold'] = 2 ** (-output_df['Delta Delta Ct'])
    output_df['Fold CI 68 Lower'] = 2 ** (-output_df['Delta Delta Ct'] - output_df['Delta Ct Std'])
    output_df['Fold CI 68 Upper'] = 2 ** (-output_df['Delta Delta Ct'] + output_df['Delta Ct Std'])
    output_df['Fold CI 68'] = output_df.apply(lambda row: f'[{row['Fold CI 68 Lower']:.2g}, {row['Fold CI 68 Upper']:.2g}]', axis=1)

    output_df['Log Fold'] = np.log2(output_df['Fold'])
    output_df['Log Fold CI 68 Lower'] = np.log2(output_df['Fold CI 68 Lower'])
    output_df['Log Fold CI 68 Upper'] = np.log2(output_df['Fold CI 68 Upper'])
    output_df['Log Fold CI 68'] = output_df.apply(lambda row: f'[{row['Log Fold CI 68 Lower']:.2g}, {row['Log Fold CI 68 Upper']:.2g}]', axis=1)
    
    for (_, reference) in reference_samples.items():
        output_df = output_df.loc[output_df['Sample'] != reference]
    for (_, reference) in reference_targets.items():
        output_df = output_df.loc[output_df['Target'] != reference]
    
    return output_df
    

def main(args: Arguments) -> None:
    output_df: pd.DataFrame
    
    if args.concat == 'input':
        output_df = foldch(args, concat_inputs(args.input_paths))
    elif args.concat == 'output':
        output_df = pd.concat([foldch(args, read_input(input)) for input in args.input_paths], axis=0)
    else:
        raise ValueError(f'invalid concat strategy {args.concat}')
    
    simple_output = output_df[['Sample', 'Target', 'Ref Sample', 'Ref Target', 'Fold', 'Fold CI 68', 'Log Fold', 'Log Fold CI 68']]
    simple_output = simple_output.rename(columns={'Fold': 'Rel Expr', 'Log Fold': 'Log Rel Expr'})
    
    print(simple_output)
    
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
