#!/bin/sh

git pull -Xours
git add data-processed/*.csv
git add outputs/*.csv
git add pkgdown/assets
git diff-index --quiet HEAD || git commit -m"Updated estimates $(date)"
git push -v
