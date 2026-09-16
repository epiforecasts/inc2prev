#!/bin/sh

# Age-stratified estimates using antibody and vaccination data.
# Additional arguments (e.g. -m <date>) are passed to both scripts.
Rscript scripts/estimate.R -g -a -d 1 --chains=4 --warmup=1000 $* && Rscript scripts/analyse.R -g -a $*
