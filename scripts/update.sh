#!/bin/sh

git pull -Xours

Rscript data-raw/update_cis.R

echo Space
git status | grep -q cis.csv && Rscript scripts/estimate.R -n -d 1 && Rscript scripts/analyse.R

echo Age
git status | grep -q cis_age.csv && Rscript scripts/estimate.R -g -d 1 && Rscript scripts/analyse.R -g

echo Variants
git status | grep -q cis_variants.csv && Rscript scripts/estimate.R -v -d 1 && Rscript scripts/analyse.R -v

Rscript scripts/plot_estimates.r
cp figures/additional/national_infections_all.svg figures
cp figures/additional/national_r_3months.svg figures
cp figures/additional/national_Rt_3months.svg figures

git add data/*.csv
git add outputs/*.csv
git diff-index --quiet HEAD || git commit -m"Updated estimates $(date)"
git push -v
