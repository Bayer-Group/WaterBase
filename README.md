# WaterBase



This is an R package for reproducing the Waterbase data analysis and the results for manuscript.

## Data

- Raw data handling is located in ./data-raw
- User would need to set up a folder containing downloaded data files from the water base website:  `Waterbase_v2019_1_S_WISE6_SpatialObject_DerivedData_repaired.csv` and `Waterbase_v2020_1_T_WISE6_AggregatedData.csv`, those can be downloaded from the **Water Framework Directive** website.
- The primary source for toxicity benchmarks (BMs) was the set of species sensitivity distributions (SSDs) from Posthuma et al. (2019). This dataset consisted of parameters (median and standard deviation of log-transformed toxicity data) of log-normal SSDs for 12386 chemicals based on chronic (NOEC or EC10) endpoints (Posthuma et al., 2019b). Details about the SSDs, underlying ecotoxicological data, and quality assurance can be found elsewhere (Posthuma et al., 2019b). If you would like to obtain the raw SSD data file, please either contact Posthuma or the authors of the article. 


## Reproduce the Results in the Manuscript

### Installing


```
install_github("Bayer-Group/WaterBase")
```

### Data Selections and Curation

Note that this step has been done inside the package and the pre-processed water base data is a data object `Data.AggBySiteID.2` included in the WaterBase package. You can use `data("Data.AggBySiteID.2")` to load the data object. Also, the chemical classes are defined in the data object `chemclass`.

The detailed pre-processing of the data including data selection and curation are included in the file 

```
vignettes/articles/Stepwise-Curation_WB_v2.Rmd`. 
```

The file below include processing the the chemicalclass information.

```
data-raw/DATASET_class.R
```


### Tables and Figures in the Manuscript

The figures and tables in the manuscript can be reproduced using the following R markdown file. The generated figures are saved in `inst/manuscript_2022`

```
vignettes/articles/Manuscript_Prep_2022_article.Rmd
```

Sensitivity analysis is documented in 

```
vignettes/articles/Sensitivity.Rmd
```

### Tables and Figures in the Supplemental Materials

The figures and tables in the supplemental information can be reproduced using the following R markdown files. The generated figures are saved in `inst/manuscript_2022` and `inst/manuscript_maps`, tables are in `inst/SM`

```
vignettes/Manuscript_Supplemental.Rmd
vignettes/articles/Manuscript_supp_summaries.Rmd
```

## Authors

* **Ismael M. Rodea Palomares** - [Ismael M. Rodea Palomares](mailto:ismaelm.rodeapalomares@bayer.com)
* **Zhenglei Gao** - [Zhenglei Gao](mailto:zhenglei.gao@bayer.com)


## Acknowledgments

* Thanks to Arnd Weyers & Markus Ebeling for their valuable inputs.


