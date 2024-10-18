### Protein and PTM Quantitative Analysis Workflow.

### The workflow

The present workflow accepts PSM reports of multiple MS runs. Merges the PSM reports by the functions of meRgreeMS to create a single peptide report, as it can be seen in the meRgreeMS_schema.pdf file. The peptide report is then imported in VIQoR for protein quantification analysis, while PTMomics data are separated and normalized by the protein expressions. Protein and PTM feautures are statistically inferred by PolySTest and clustered by VSClust. The schematic of the complete workflow is illustrated in Analysis_schema.pdf file.


### Installation
Download repository to your computer. Keep the folders for the 3 steps of the analysis but remove the existing data files. Place the PSM reports and the corresponding .fasta protein sequence database in "1. Quality Control" folder and follow the analysis provided by comments in PTMs_Analysis.R script.

Install the following R packages:

```R
install.packages(c("readxl", "ggplot2", "reshape", "reshape2", "grid", 
                   "gridExtra", "scales", "protr", "stringi"))
```
The dependences are also included in PTMs_Analysis.R and meRgreeMS.R files.

Load the PTMs_Analysis.R file in [Rstudio](http://rstudio.com) and follow the step by step instructions provided in the script.


### Tools web-services
The tools needed for the analysis are hosted by the [Computational Proteomics group](http://computproteomics.bmb.sdu.dk:). The tools online services can be accessed by any web browser (preferably Chrome or Mozilla).

* VIQoR can be found [here](http://computproteomics.bmb.sdu.dk:8192/app_direct/VIQoR/).
* PolySTest can be found [here](http://computproteomics.bmb.sdu.dk:8192/app_direct/PolySTest/).
* VSClust can be found [here](http://computproteomics.bmb.sdu.dk:8192/app_direct/VSClust/).

### Contact
For software issues and general questions, please submit an issue.

### License
Apache-2.0