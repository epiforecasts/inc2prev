#!/bin/bash

while read line; do
  echo standard $line
  Rscript scripts/estimate.R -d 1 -g -a -m $line
done < scripts/dates

