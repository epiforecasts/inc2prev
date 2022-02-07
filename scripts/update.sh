#!/bin/sh

Rscript data-raw/update_cis.R

echo Space
Rscript scripts/estimate.R -n -d 1
Rscript scripts/analyse.R

echo Age
Rscript scripts/estimate.R -g -d 1
Rscript scripts/analyse.R -g

echo Variants
Rscript scripts/estimate.R -v -d 1
Rscript scripts/analyse.R -v
