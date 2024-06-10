#!/usr/bin/env python3

import argparse
import sys
from pathlib import Path
from typing import List

from matplotlib.axes import Axes
from matplotlib.figure import Figure
import pandas as pd
import matplotlib.pyplot as plt

class Arguments:
    input_files: List[str]
    input_paths: List[Path]
    groups: List[List[str]] = []
    sample_renames: List[List[str]] = []
    output_dir: str
    file_suffix: str
    
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
parser.add_argument('--linear', action='store_true', default=False, help='plot non-log fold values')
parser.add_argument('--group', dest='groups', type=str, nargs='+', action='append', help='group samples in one plot')
parser.add_argument('--output', '-o', dest='output_dir', type=str, default='.', help='output directory for figures')
parser.add_argument('--figsuffix', dest='file_suffix', type=str, default='', help='suffix for file names')
parser.add_argument(
    '--rename-sample', '-S',
    dest='sample_renames',
    type=str,
    metavar=('SAMPLE', 'NEW_NAME'),
    nargs=2,
    action='append',
    help='samples to rename'
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


def get_sample_name(args: Arguments, sample: str) -> str:
    for (old, new) in args.sample_renames:
        if old == sample: return new
    return sample


def main(args: Arguments) -> None:
    output_path = Path(args.output_dir).resolve()
    
    input_df = concat_inputs(args.input_paths)
    
    grouping = input_df.groupby(['Sample'])
    unplotted_samples = [group[0][0] for group in grouping]
    
    def plot(axes: Axes, group: pd.DataFrame) -> None:
        values: pd.DataFrame
        errors: pd.DataFrame
        if args.linear:
            values = group['Fold']
            errors = group.apply(lambda entry: [entry['Fold'] - entry['Fold CI 68 Lower'], entry['Fold CI 68 Upper'] - entry['Fold']], axis=1, result_type='expand').transpose()
            axes.axhline(y=1.0, color='gray', linestyle='--').set_linewidth(0.5)
        else:
            values = group['Log Fold']
            errors = group.apply(lambda entry: [entry['Log Fold'] - entry['Log Fold CI 68 Lower'], entry['Log Fold CI 68 Upper'] - entry['Log Fold']], axis=1, result_type='expand').transpose()
            axes.axhline(y=0.0, color='gray', linestyle=(0, (5, 5))).set_linewidth(0.5)
        
        axes.bar(group['Target'], values, edgecolor='black', color='none')
        axes.errorbar(group['Target'], values, yerr=errors, fmt='o', color='red', capsize=3, markersize=5)
        axes.tick_params(axis='x', labelrotation=45)
        
        axes.title.set_text(get_sample_name(args, sample))

    figcount = 1
    def new_figname() -> str:
        nonlocal figcount
        name = f'figure{figcount}'
        if len(args.file_suffix) != 0: name = name + f'-{args.file_suffix}'
        figcount = figcount + 1
        return name
    
    
    def save_fig(fig: Figure) -> None:
        figure_path = output_path / f'{new_figname()}.pdf'
        fig.savefig(figure_path, format='pdf')
        print(f'Saved figure to {figure_path}')
    

    for plot_group in args.groups:
        (fig, axes_list) = plt.subplots(1, len(plot_group), sharex=True, sharey=True)
        fig.suptitle('Relative gene expression', weight='bold').set_fontsize('x-large')
        fig.set_size_inches(6 * len(plot_group), 5)
        
        for i in range(len(plot_group)):
            sample = plot_group[i]
            axes: Axes = axes_list[i]
            
            group = grouping.get_group((sample,))
            (reference_sample,) = group['Ref Sample'].unique()
            (reference_target,) = group['Ref Target'].unique()
            
            if sample in unplotted_samples:
                unplotted_samples.remove(sample)

            plot(axes, group)
            axes.text(0.01, 0.99, f'Control: {get_sample_name(args, reference_sample)}/{reference_target}', ha='left', va='top', transform=axes.transAxes).set_fontsize('small')

            if i == 0: 
                axes.set_ylabel('fold change' if args.linear else 'log2(fold change)')
        
        fig.tight_layout()
        save_fig(fig)
        
    for sample in unplotted_samples:
        (fig, axes) = plt.subplots(1)
        fig.suptitle('Relative gene expression', weight='bold').set_fontsize('x-large')
        fig.set_size_inches(6, 5)
        
        group = grouping.get_group((sample,))
        (reference_sample,) = group['Ref Sample'].unique()
        (reference_target,) = group['Ref Target'].unique()
        plot(axes, group)
        
        axes.text(0.01, 0.99, f'Control: {get_sample_name(args, reference_sample)}/{reference_target}', ha='left', va='top', transform=axes.transAxes).set_fontsize('small')
        
        fig.tight_layout()
        save_fig(fig)
        
    plt.show()


if __name__ == "__main__":
    args = parser.parse_args(None, Arguments)
    
    args.input_paths = [Path(file).resolve() for file in args.input_files]

    nonexistent = [path for path in args.input_paths if not path.exists()]
    if any(nonexistent):
        print(f'Some input files do not exist ({', '.join(str(path) for path in nonexistent)})', file=sys.stderr)
        sys.exit(1)
    
    main(args)
