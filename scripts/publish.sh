#!/bin/sh

git pull -Xours
git add data/*.csv
git add outputs/*.csv
git diff-index --quiet HEAD || git commit -m"Updated estimates $(date)"
git push -v
