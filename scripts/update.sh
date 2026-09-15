#!/bin/sh

git pull -Xours

Rscript data-raw/update_cis.R

echo National
git status | grep -q cis.csv && Rscript scripts/estimate.R -d 1 && Rscript scripts/analyse.R

echo Regional
git status | grep -q cis.csv && Rscript scripts/estimate.R -r -n -d 1 && Rscript scripts/analyse.R -r

echo Age
git status | grep -q cis_age.csv && Rscript scripts/estimate.R -g -d 1 && Rscript scripts/analyse.R -g

echo Variants
git status | grep -q cis_variants.csv && Rscript scripts/estimate.R -v -d 1 && Rscript scripts/analyse.R -v

Rscript scripts/plot_estimates.r
Rscript -e 'rmarkdown::render("writeup/report.Rmd", output_dir = "pkgdown/assets")'
Rscript -e 'rmarkdown::render("writeup/paper.Rmd", output_format = "all", output_dir = "pkgdown/assets")'
