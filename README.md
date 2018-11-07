# dbgapXML
Parse dbGaP XML data dictionaries and variable reports to tab-delimited text files.

## Parameters
`study` phs study accession number

`--table` pht table accession number. No `--table` extracts all variable reports and data dictionaries in study.

`--dir` directory to which files are written.

## Example

To parse the data dictionaries and variable reports from the Framingham Heart Study (study accession number phs000007) dataset "Clinic Exam, Original Cohort Exam 22", which has the table accession number pht000024, you would use the following line:

```
Rscript dbgap-variable-summary-extraction-v3.R phs000007 --table pht000024
```
