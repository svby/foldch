# foldch

Fold gene expression analysis, for Bio-Rad CFX qPCR output.

Takes a path to an Excel spreadsheet, a reference sample/cell line and a reference target (i.e. housekeeping gene) as input.

## Dependencies

- numpy
- pandas
- matplotlib

## Example usage

### Basic usage

```
# Will create a file titled
# "Analysis - HEK293T - Quantification Summary.csv"
# Treats "HEK293T parental" as the reference sample,
# and "18S" as the reference target, for all samples
./foldch.py \
    --rs all "HEK293T parental" \
    --rt all "18S" -- \
    "HEK293T - Quantification Summary.csv"
```

### PCR plates with separate experimental setups

```
# Will use "Untreated 4h" as the reference sample
# for all samples ending in 4h, and similarly for 12h
./foldch.py \
    --rt all GAPDH \
    --rs /.+ 4h/ "Untreated 4h" \
    --rs /.+ 12h/ "Untreated 12h" \
    -o analysis.csv -- \
    *Quantification\ Summary*.csv
```

### Bar chart creation

```
./foldch.py ...
./barchart.py "Analysis - HEK293T - Quantification Summary.csv"
```
