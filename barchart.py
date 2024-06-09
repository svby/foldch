#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path
from typing import List

from matplotlib.axes import Axes
import pandas as pd
import matplotlib.pyplot as plt

class Arguments:
    input_files: List[str]
    input_paths: List[Path]
    groups: List[List[str]] = []
    
    linear: bool


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
parser.add_argument('--linear', type=bool, default=False, help='plot non-log fold values')
parser.add_argument('--group', dest='groups', type=str, nargs='+', action='append', help='group samples in one plot')


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
    print(args.groups)
    
    grouping = input_df.groupby(['Sample'])
    unplotted_samples = [group[0][0] for group in grouping]
    
    def plot(axes: Axes, group: pd.DataFrame) -> None:
        values: pd.DataFrame
        errors: pd.DataFrame
        if args.linear:
            values = group['Fold']
            errors = group.apply(lambda entry: [entry['Fold'] - entry['Fold CI 68 Lower'], entry['Fold CI 68 Upper'] - entry['Fold']], axis=1, result_type='expand').transpose()
            axes.axhline(y=1.0, color='gray', linestyle=(0, (5, 5))).set_linewidth(0.5)
        else:
            values = group['Log Fold']
            errors = group.apply(lambda entry: [entry['Log Fold'] - entry['Log Fold CI 68 Lower'], entry['Log Fold CI 68 Upper'] - entry['Log Fold']], axis=1, result_type='expand').transpose()
            axes.axhline(y=0.0, color='gray', linestyle='--').set_linewidth(0.5)
        
        axes.bar(group['Target'], values, edgecolor='black', color='none')
        axes.errorbar(group['Target'], values, yerr=errors, fmt='o', color='red', capsize=3, markersize=5)
        
        axes.title.set_text(f'Sample: {sample}')

    for plot_group in args.groups:
        (fig, axes_list) = plt.subplots(1, len(plot_group), sharex=True, sharey=True)
        fig.suptitle('Relative gene expression')
        fig.set_size_inches(5 * len(plot_group), 5)
        
        controls = []
        for i in range(len(plot_group)):
            sample = plot_group[i]
            axes: Axes = axes_list[i]
            
            group = grouping.get_group((sample,))
            (reference_sample,) = group['Ref Sample'].unique()
            (reference_target,) = group['Ref Target'].unique()
            controls.append(f'{reference_sample}/{reference_target}')
            
            unplotted_samples.remove(sample)
            plot(axes, group)

            if i == 0: 
                axes.set_ylabel('fold change' if args.linear else 'log2(fold change)')
        
        if len(set(controls)) == 1: controls = [controls[0]]
        fig.text(0.015, 0.015, f'Control: {', '.join(controls)}').set_fontsize('medium')
        
    for sample in unplotted_samples:
        (fig, axes) = plt.subplots(1)
        fig.suptitle('Relative gene expression')
        fig.set_size_inches(5, 5)
        
        group = grouping.get_group((sample,))
        (reference_sample,) = group['Ref Sample'].unique()
        (reference_target,) = group['Ref Target'].unique()
        plot(axes, group)
        
        fig.text(0.015, 0.015, f'Control: {reference_sample}/{reference_target}').set_fontsize('medium')
        
    plt.show()


if __name__ == "__main__":
    args = parser.parse_args(None, Arguments)
    
    args.input_paths = [Path(file).resolve() for file in args.input_files]

    nonexistent = [path for path in args.input_paths if not path.exists()]
    if any(nonexistent):
        print(f'Some input files do not exist ({', '.join(str(path) for path in nonexistent)})', file=sys.stderr)
        sys.exit(1)
    
    main(args)
