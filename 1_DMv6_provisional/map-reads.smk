import os, re
import pandas as pd
from Bio import SeqIO
shell.executable("bash")

units = pd.read_csv(config["units"], index_col=["sample", "unit"], dtype=str, sep = "\t")
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])

short = units.loc[units["readlen"] == "short"]
short.index = short.index.remove_unused_levels()

long  = units.loc[units["readlen"] == "long"]
long.index = long.index.remove_unused_levels()

short_samples = {}
for i in short.index.levels[0]:
    short_samples[i] = ["data/clipOverlap/{}-{}.bam".format(i, j) for j in short.loc[i].index]

long_samples = {}
for i in long.index.levels[0]:
    long_samples[i] = ["{}-{}".format(i, j) for j in long.loc[i].index]

dedup_samples = {}
for i in short.index.levels[0]:
    dedup_samples[i] = ["data/dedup/{}-{}.bam".format(i, j) for j in short.loc[i].index]

def is_single_end(sample,unit):
    return pd.isnull(short.loc[(sample, unit), "fq2"])

def get_fastq(wildcards):
    return "data/reads/"+short.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()

def get_trimmed(wildcards):
    if not is_single_end(**wildcards):
        return expand("data/trimmed/{sample}-{unit}-{group}.fastq.gz",
            group=[1,2], **wildcards)
    return "data/trimmed/{sample}-{unit}.fastq.gz".format(**wildcards)

def get_long(wildcards):
    return "data/reads/"+long.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    
rule all:
    input:
        ["data/merged/{}.bam".format(x) for x in units.index.levels[0]],
        ["data/read_depth/{}.tsv".format(x) for x in units.index.levels[0]]
        
include: "rules/merge_long.rules"
include: "rules/merge_dedup.rules"
include: "rules/window_read_depth.rules"
include: "rules/samtools_index_pe.rules"
include: "rules/samtools_merge.rules"
include: "rules/clip_overlap.rules"
include: "rules/filter_good_pairs.rules"
include: "rules/mark_duplicates.rules"
include: "rules/filter_mapped_long.rules"
include: "rules/map_long.rules"
include: "rules/align.rules"
include: "rules/index_asm.rules"
include: "rules/cutadapt_pe.rules"
include: "rules/cutadapt.rules"
include: "rules/fastqc.rules"