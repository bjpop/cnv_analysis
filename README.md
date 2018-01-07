# Overview 

CNV analysis tools

# Licence

This program is released as open source software under the terms of [MIT License](https://raw.githubusercontent.com/cnv_analysis-team/cnv_analysis/master/LICENSE)

# Installing

cnv_analysis can be installed using `pip` in a variety of ways (`%` indicates the command line prompt):

1. Inside a virtual environment: 
```
% python3 -m venv cnv_analysis_dev 
% source cnv_analysis_dev/bin/activate
% pip install -U /path/to/cnv_analysis
```
2. Into the global package database for all users:
```
% pip install -U /path/to/cnv_analysis
```
3. Into the user package database (for the current user only):
```
% pip install -U --user /path/to/cnv_analysis
```

# Usage

Assuming all cnv data is in the file `data_all.tsv`

Merge all CNVs per family:
```
merge_cnvs data_all.tsv > merged_cnvs_by_family.tsv

``` 

Case-control analysis of all CNVs. Duplicates will be output to `data_all.dups.csv`.
```
case_control_cnvs --merged merged_cnvs_by_family.tsv --all data_all.tsv
```

Generate graphs of CNV relations:
```
graph_cnvs case_control_by_family.tsv 

```
