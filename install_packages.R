# Installed packages: shiny, tidyverse, devtools
# Taken from here: https://github.com/xmc811/ShinyProxy-template/blob/master/optional/install_packages.R
# some necessities
install.packages(c('devtools',
                   'dplR',
                   'corrplot',
                   'dbscan',
                   'pryr',
                   'magrittr',
                   'R.utils',
                   'stringr',
                   'precrec',
                   'optparse',
                   'e1071',
                   'readr',
                   'dplyr',
                   'magrittr',
                   'stringr',
                   'vcfR',
                   'reticulate'
                   ), repos = 'http://cran.rstudio.com/')

# Install packages for VARPP app
devtools::install_github("Hobbeist/varppRule@v1.0")

