#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA

clusters = pd.read_csv("souporcell/clusters.tsv")

barcode_counts = pd.read_csv("barcodes.tsv.gz", sep="\\t")

sample_assignments = pd.read_csv("sample_manifest.csv", header=None, names=["index", "sample"])

print(clusters.head())
print(barcode_counts.head())
print(sample_assignments.head())
assert False
