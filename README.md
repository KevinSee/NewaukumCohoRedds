
<!-- README.md is generated from README.Rmd. Please edit that file -->

# NewaukumCohoRedds

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/KevinSee/NewaukumCohoRedds/master?urlpath=rstudio)

This repository contains the data and code for simulating GRTS surveys
for Coho redds on the Newaukum river. The initial goal is to test
various GRTS study designs and test them against several years of census
count redd surveys.

## Contents

The **analysis** directory contains:

-   [:file_folder: data](/analysis/data): Data used in the analysis.
    -   [:file_folder: raw_data](/analysis/data/raw_data) contains tag
        lists and data downloaded from PTAGIS or other sources
    -   [:file_folder: derived_data](/analysis/data/derived_data)
        contains analysis products that have been saved along the way,
        including model fits and cleaned up habitat and fish data.
-   [:file_folder: figures](/analysis/figures): Plots and other
    illustrations
-   [:file_folder: paper](/analysis/paper): R Markdown source documents.
-   [:file_folder: scripts](/analysis/scripts): R scripts  
-   [:file_folder:
    supplementary-materials](/analysis/supplementary-materials):
    Supplementary materials including notes and other documents prepared
    and collected during the analysis.

## How to run in your browser or download and run locally

This research compendium has been developed using the statistical
programming language R. To work with the compendium, you will need
installed on your computer the [R
software](https://cloud.r-project.org/) itself and optionally [RStudio
Desktop](https://rstudio.com/products/rstudio/download/).

You can download the compendium as a zip from from this URL:
[master.zip](/archive/master.zip). After unzipping: - open the `.Rproj`
file in RStudio - run `devtools::install()` to ensure you have the
packages this analysis depends on (also listed in the
[DESCRIPTION](/DESCRIPTION) file).

### Licenses

**Text and figures :**
[CC-BY-4.0](http://creativecommons.org/licenses/by/4.0/)

**Code :** See the [DESCRIPTION](DESCRIPTION) file

**Data :** [CC-0](http://creativecommons.org/publicdomain/zero/1.0/)
attribution requested in reuse

### Contributions

We welcome contributions from everyone. Before you get started, please
see our [contributor guidelines](CONTRIBUTING.md). Please note that this
project is released with a [Contributor Code of Conduct](CONDUCT.md). By
participating in this project you agree to abide by its terms.
