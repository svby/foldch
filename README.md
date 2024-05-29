# foldch

Fold gene expression analysis, for Bio-Rad CFX qPCR output.

Takes a path to an Excel spreadsheet, a reference sample/cell line and a reference target (i.e. housekeeping gene) as input.

## Dependencies

- numpy
- pandas

## Example usage

```
# Will create a file titled "Analysis - H1975 - Quantification Summary.xlsx"
./foldch.py "H1975 - Quantification Summary.xlsx" --rs "H1975" --rt "Mdh1"
```