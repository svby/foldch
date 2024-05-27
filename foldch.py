#!/usr/bin/env python3

# TODO: only one reference sample, so PCR plates with
# different cell lines aren't particularly ideal currently

import math
import argparse
import pandas as pd
import numpy as np

class Arguments:
    input_file: str
    reference_sample: str
    reference_target: str

parser = argparse.ArgumentParser(
    prog='foldch'
)
parser.add_argument('input_file', type=str)
parser.add_argument('--reference_sample', '--rs', type=str, default='H1975 par')
parser.add_argument('--reference_target', '--rt', type=str, default='28S')

def main(args: Arguments):
    df = pd.read_excel(args.input_file)

    targets = df['Target'].unique()
    samples = df['Sample'].unique()

    combinations = np.meshgrid(samples, targets)
    combinations = np.vstack(np.swapaxes(np.stack(np.swapaxes(combinations, 1, 2), 1), 1, 2))
    
    print('Reference sample:', args.reference_sample, '; reference target:', args.reference_target)

    # This is not idiomatic, probably.
    # Improve this later, I've spent too much time here moving this from Excel
    grouped_combinations = df.groupby(['Sample', 'Target'])

    intermediate = []
    for sample, target in combinations:
        group = grouped_combinations['Cq'].get_group((sample, target))
        ctrl_target_group = grouped_combinations['Cq'].get_group((sample, args.reference_target))

        ct_avg = group.mean()
        ct_avg_ctrl = ctrl_target_group.mean()
        
        ct_std = group.std()
        ct_std_ctrl = ctrl_target_group.std()

        delta_ct = ct_avg - ct_avg_ctrl
        delta_ct_std = math.sqrt((ct_std ** 2) + (ct_std_ctrl ** 2))

        intermediate.append(
            np.array([
                sample,
                target,
                delta_ct,
                delta_ct_std
            ])
        )

    output = pd.DataFrame(np.asarray(intermediate), columns=['Sample', 'Target', 'Delta Ct', 'Std Delta Ct'])
    output[['Delta Ct', 'Std Delta Ct']] = output[['Delta Ct', 'Std Delta Ct']].apply(pd.to_numeric)
    # output = output.convert_dtypes()

    results = []
    for sample, target in combinations:
        intermediate_entry = output.groupby(['Sample', 'Target']).get_group((sample, target))
        
        # Somewhat pointless with one control cell line, but just out of principle
        # To do: 3D df probably would remove the need for this

        avg_delta_ct_bio_ctrl = output.groupby(['Sample', 'Target']).get_group((args.reference_sample, target))['Delta Ct'].mean()
        delta_delta_ct = (intermediate_entry['Delta Ct'] - avg_delta_ct_bio_ctrl).mean()
        delta_ct_std = intermediate_entry['Std Delta Ct'].mean()
        
        fold = 2 ** (-delta_delta_ct)
        fold_ci68 = (2 ** (-delta_delta_ct - delta_ct_std), 2 ** (-delta_delta_ct + delta_ct_std))

        log_fold = math.log2(fold)
        log_fold_ci68 = (math.log2(fold_ci68[0]), math.log2(fold_ci68[1]))

        results.append([
            sample, target, fold, f'[{fold_ci68[0]:.2g}, {fold_ci68[1]:.2g}]', log_fold, f'[{log_fold_ci68[0]:.2g}, {log_fold_ci68[1]:.2g}]'
        ])

    print(pd.DataFrame(results, columns=['Sample', 'Target', 'Fold Expr', 'Fold 68% CI', 'Log Fold', 'Log Fold 68% CI']))


if __name__ == "__main__":
    main(parser.parse_args(None, Arguments))
