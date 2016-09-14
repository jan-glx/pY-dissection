# Dissecting phosphotyrosine dependent protein-protein interactions
A set of R scripts used to obtain the results presented in the same-named report on an internship of mine. Likely not interesting for anyone. More a reference for myself and maybe for someone in the group of Rob Russell working on the project.

## Requirements
* a `Data` folder parallel to the repository clone containing:
  * [IM-22632.txt](http://www.ebi.ac.uk/intact/export?format=mitab_27&query=IM-22632&negative=false&spoke=false&ontology=false&sort=intact-miscore&asc=false)
  * a folder `Out` containing the csv export of Mechismo (see step ?"
  * `naccess_human_uniprot.tsv.gz` and `pfam_instances_uniprot.tsv.gz` available upon request to me or maybe Matthew Betts
  * `Phosphosite_seq.txt`, `Kinase_Substrate_Dataset.tsv` and `Phosphorylation_site_dataset.tsv` available from [PhosphoSitePlus®](http://www.phosphosite.org/staticDownloads.action)
* R packages
 * `install.packages(c("magrittr", "tidyr", "dplyr", "stringr", "ggplot2")`
 * `install.packages("data.table", type = "source", repos = "http://Rdatatable.github.io/data.table")`
 * a wrapper R package for amodified version of ANCHOR available only with written permission of Zsuzsanna Dosztányi upon request

## Usage
Run the R scripts in the order specified by the initial digit of their file name. Intermediate results are saved in the `../Data` directory. Step 3 requires you to use [Mechismo](http://mechismo.russelllab.org/) with the file `mechismoinput.txt` generated in step 2 and download the results as `.csv` in the `../Data/Out` folder. The files starting with digit `9` are the final analysis scripts producing the reported results and can be run in any order.

Please file an issue if you encounter any difficulties.

