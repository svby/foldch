#!/usr/bin/env python3

# TODO: only one reference sample, so PCR plates with
# different cell lines aren't particularly ideal currently

import argparse
import sys
import pandas as pd
import numpy as np
from pathlib import Path


class Arguments:
    input_file: str
    output_file: str | None
    input_path: Path
    output_path: Path
    reference_sample: str
    reference_target: str


parser = argparse.ArgumentParser(
    prog='foldch'
)
parser.add_argument('input_file', type=str)
parser.add_argument('--output_file', '--out', '-o', type=str, required=False)
parser.add_argument('--reference_sample', '--rs', type=str, default='H1975 par')
parser.add_argument('--reference_target', '--rt', type=str, default='28S')


def main(args: Arguments) -> None:
    input_df: pd.DataFrame
    if args.input_path.suffix == '.xlsx' or args.input_path.suffix == '.xls': input_df = pd.read_excel(args.input_path)
    elif args.input_path.suffix == '.csv': input_df = pd.read_csv(args.input_path)
    else:
        print('Input file is of unknown type (neither .xlsx nor .csv)', file=sys.stderr)
        sys.exit(1)

    targets = input_df['Target'].unique()
    samples = input_df['Sample'].unique()

    combinations = np.meshgrid(samples, targets)
    combinations = np.vstack(np.swapaxes(np.stack(np.swapaxes(combinations, 1, 2), 1), 1, 2))
    
    output_df = pd.DataFrame(combinations, columns=['Sample', 'Target'])
    
    def get_group(df: pd.DataFrame, sample: str, target: str) -> pd.DataFrame:
        return df.groupby(['Sample', 'Target']).get_group((sample, target))
    
    # Calculated from input
    output_df['Ct Avg'] = output_df.apply(lambda x: get_group(input_df, x['Sample'], x['Target'])['Cq'].mean(), axis=1)
    output_df['Ct Std'] = output_df.apply(lambda x: get_group(input_df, x['Sample'], x['Target'])['Cq'].std(), axis=1)
    output_df['Ct Avg Ctrl'] = output_df.apply(lambda x: get_group(input_df, x['Sample'], args.reference_target)['Cq'].mean(), axis=1)
    output_df['Ct Std Ctrl'] = output_df.apply(lambda x: get_group(input_df, x['Sample'], args.reference_target)['Cq'].std(), axis=1)

    # Calculated from intermediate values
    output_df['Delta Ct'] = output_df['Ct Avg'] - output_df['Ct Avg Ctrl']
    output_df['Delta Ct Std'] = np.sqrt((output_df['Ct Std'] ** 2) + (output_df['Ct Std Ctrl'] ** 2))
    
    output_df['Delta Ct Avg Ctrl Sample'] = output_df.apply(lambda x: get_group(output_df, args.reference_sample, x['Target'])['Delta Ct'].mean(), axis=1)
    output_df['Delta Delta Ct'] = output_df['Delta Ct'] - output_df['Delta Ct Avg Ctrl Sample']
    
    output_df['Fold'] = 2 ** (-output_df['Delta Delta Ct'])
    output_df['Fold CI 68 Lower'] = 2 ** (-output_df['Delta Delta Ct'] - output_df['Delta Ct Std'])
    output_df['Fold CI 68 Upper'] = 2 ** (-output_df['Delta Delta Ct'] + output_df['Delta Ct Std'])
    output_df['Fold CI 68'] = output_df.apply(lambda x: f'[{x['Fold CI 68 Lower']:.2g}, {x['Fold CI 68 Upper']:.2g}]', axis=1)

    output_df['Log Fold'] = np.log2(output_df['Fold'])
    output_df['Log Fold CI 68 Lower'] = np.log2(output_df['Fold CI 68 Lower'])
    output_df['Log Fold CI 68 Upper'] = np.log2(output_df['Fold CI 68 Upper'])
    output_df['Log Fold CI 68'] = output_df.apply(lambda x: f'[{x['Log Fold CI 68 Lower']:.2g}, {x['Log Fold CI 68 Upper']:.2g}]', axis=1)
    
    simple_output = output_df[['Sample', 'Target', 'Fold', 'Fold CI 68', 'Log Fold', 'Log Fold CI 68']]
    
    output_df.to_excel(args.output_path)
    
    print(simple_output)
    print('Reference sample:', args.reference_sample, '; reference target:', args.reference_target)


if __name__ == "__main__":
    args = parser.parse_args(None, Arguments)
    args.input_path = Path(args.input_file).resolve()
    args.output_path = Path(args.output_file) if args.output_file is not None else args.input_path.parent / f'Analysis - {args.input_path.stem}.xlsx'
    
    if args.output_path.exists():
        print('Output file already exists', file=sys.stderr)
        sys.exit(1)
    
    main(args)
