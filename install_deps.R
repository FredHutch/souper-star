#!/usr/bin/env Rscript

# Following instructions on https://www.archrproject.com/

install.packages("ragg")
library(ragg)
install.packages("pkgdown")
library(pkgdown)
install.packages("devtools")
library(devtools)
install.packages("BiocManager")
BiocManager::install("chromVAR")
BiocManager::install("motifmatchr")

devtools::install("/opt/r_deps/ArchR-1.0.2", repos = BiocManager::repositories())
library(ArchR)
# ArchR::installExtraPackages()