#!/bin/sh

git pull -Xours

Rscript data-raw/update_cis.R
Rscript data-raw/update_ab.R

echo Local
Rscript scripts/estimate.R -l --chains=4 --warmup=1000 $* && Rscript scripts/analyse.R -l $*

echo Regional
Rscript scripts/estimate.R -r -i -a --chains=4 --warmup=1000 $* && Rscript scripts/analyse.R -r -a -i $*
Rscript scripts/estimate.R -r -a --chains=4 --warmup=1000 $* && Rscript scripts/analyse.R -r -a $*

Rscript scripts/generate_local_ab.r $*
