#!/bin/bash

Rscript -e 'rmarkdown::render("sim_fit.Rmd","html_document","o21000nonenone",params=list(nsv =2,sample_size=1000,bias="none",trend="none"))'
Rscript -e 'rmarkdown::render("sim_fit.Rmd","html_document","o21000noneincrease",params=list(nsv =2,sample_size=1000,bias="none",trend="increase"))'
Rscript -e 'rmarkdown::render("sim_fit.Rmd","html_document","o21000nonedecrease",params=list(nsv =2,sample_size=1000,bias="none",trend="decrease"))'

Rscript -e 'rmarkdown::render("sim_fit.Rmd","html_document","o31000nonenone",params=list(nsv =3,sample_size=1000,bias="none",trend="none"))'
Rscript -e 'rmarkdown::render("sim_fit.Rmd","html_document","o31000noneincrease",params=list(nsv =3,sample_size=1000,bias="none",trend="increase"))'
Rscript -e 'rmarkdown::render("sim_fit.Rmd","html_document","o31000nonedecrease",params=list(nsv =3,sample_size=1000,bias="none",trend="decrease"))'

Rscript -e 'rmarkdown::render("sim_fit.Rmd","html_document","o21000logisnone",params=list(nsv =2,sample_size=1000,bias="logis",trend="none"))'
Rscript -e 'rmarkdown::render("sim_fit.Rmd","html_document","o21000logisincrease",params=list(nsv =2,sample_size=1000,bias="logis",trend="increase"))'
Rscript -e 'rmarkdown::render("sim_fit.Rmd","html_document","o21000logisdecrease",params=list(nsv =2,sample_size=1000,bias="logis",trend="decrease"))'

Rscript -e 'rmarkdown::render("sim_fit.Rmd","html_document","o31000logisnone",params=list(nsv =3,sample_size=1000,bias="logis",trend="none"))'
Rscript -e 'rmarkdown::render("sim_fit.Rmd","html_document","o31000logisincrease",params=list(nsv =3,sample_size=1000,bias="logis",trend="increase"))'
Rscript -e 'rmarkdown::render("sim_fit.Rmd","html_document","o31000logisdecrease",params=list(nsv =3,sample_size=1000,bias="logis",trend="decrease"))'

# 4 surveys
Rscript -e 'rmarkdown::render("sim_fit.Rmd","html_document","o41000logisnone",params=list(nsv =4,sample_size=1000,bias="logis",trend="none"))'
Rscript -e 'rmarkdown::render("sim_fit.Rmd","html_document","o41000logisincrease",params=list(nsv =4,sample_size=1000,bias="logis",trend="increase"))'
Rscript -e 'rmarkdown::render("sim_fit.Rmd","html_document","o41000logisdecrease",params=list(nsv =4,sample_size=1000,bias="logis",trend="decrease"))'

