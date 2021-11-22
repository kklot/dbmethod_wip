#!/bin/bash

Rscript -e 'rmarkdown::render("fit_dhs_4cc_tmb.Rmd", "html_document",
"Uganda_male", params = list(country = "Uganda", sex = 1))'

Rscript -e 'rmarkdown::render("fit_dhs_4cc_tmb.Rmd", "html_document",
"Guinea_male", params = list(country = "Guinea", sex = 1))'

Rscript -e 'rmarkdown::render("fit_dhs_4cc_tmb.Rmd", "html_document",
"Rwanda_female", params = list(country = "Rwanda", sex = 2))'

Rscript -e 'rmarkdown::render("fit_dhs_4cc_tmb.Rmd", "html_document",
"Mozambique_female", params = list(country = "Mozambique", sex = 2))'

