# SQANTI3 OUTPUT COMPARATOR

The script **sqanti3_output_comparison** is a plugin for the SQANTI3 tool ([publication](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5848618/) and [code repository](https://github.com/ConesaLab/SQANTI)). This tool aims to compare different SQANTI3 outputs in order to explore the similarities and differences from different samples, replicates, pipelines or technologies.

## Running

This tool is divided in 2 scripts:

- **sqanti3_output_report.R**: loads the input data and runs the R Markdown script.
- **SQANTI3_output_report.Rmd**: generates a report with the comparison of the given information

### sqanti3_output_report.R

The input data are the classification and junction (*\*.txt*) files returned by the SQANTI3 quality control software. It is possible to add as many samples as you desire but it is mandatory to introduce for each sample both the classification and junction file with the same prefix and containing the default sufix given by SQANTI3 ( \*\_classification.txt and \*\_junctions.txt). All this files must be in the same directory which relative path will be given as an argument to the script.

```
$ ls ~/home/dir_in

sample1_classification.txt        sample1_junctions.txt

sample2_classification.txt        sample2_junctions.txt

sample3_classification.txt        sample3_junctions.txt
```

#### Usage

The -d argument in mandatory.

- **-d**: Name of the directory containing the SQANTI3 output files (classification and junction files are required)

- **-o**: Output directory for the report and CSV file (working directory as default)

- **-n**: Output name for the HTML report (without extension)

```
#!bash
$ Rscript sqanti3_comparison_report.R -h

Usage: sqanti3_comparison_report.R [-i DIRIN] [-o DIROUT] [-n OUTNAME]


Options:
	-d DIRIN, --dir=DIRIN
		directory with input files (classification and junction files)

	-o DIROUT, --outdir=DIROUT
		Output directory for the report and CSV file [default= .]

	-n OUTNAME, --name=OUTNAME
		Output name for the HTML report and CSV file (without extension) [default= comparison_output]

	-h, --help
		Show this help message and exit
```

#### Example run in bash

Using only the -d argument:

`Rscript sqanti3_comparison_report.R -d exampledata`

Using all arguments:

`Rscript sqanti3_comparison_report.R -d exampledata -o exampleout -n nameout`
