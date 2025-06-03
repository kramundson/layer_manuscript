import os, re
import pandas as pd
from Bio import SeqIO
shell.executable("bash")

units = pd.read_csv(config["units"], index_col=["sample", "unit"], dtype=str, sep = "\t")
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])

short = units.loc[units["hap8_discovery"] == "TRUE"]
short.index = short.index.remove_unused_levels()

geno = units.loc[units["hap8_genotyping"] == "TRUE"]
geno.units = geno.index.remove_unused_levels()

# short_samples = {}
# for i in short.index.levels[0]:
#     short[i] = ["data/merged/{}-{}.bam".format(i, j) for j in short.loc[i].index]

targets = glob_wildcards("ref/ref_windows/windows-{id}").id

rule all:
    input:
        # "data/calls/firstpass-filter.vcf.gz",
        "data/calls/bb-mf93s1-targets.vcf"
        
# variant discovery, parallelized across regions of the genome
rule freebayes_interval:
    input:
        ref=config["genome"],
        short_bam=["data/merged/{}.bam".format(x) for x in short.index.levels[0]],
        short_bai=["data/merged/{}.bam.bai".format(x) for x in short.index.levels[0]],
        target="ref/ref_windows/windows-{chunk}"
    output:
        "data/calls/chunk/{chunk}-calls.vcf"
    params:
        options=config["params"]["freebayes"]
    log:
        "log/freebayes/{chunk}.log"
    shell: """
        freebayes \
            -t {input.target} \
            --fasta-reference {input.ref} \
            --bam {input.short_bam} \\
            {params.options} \
            > {output} \
            2> {log}
    """
    
# gather region-specific genotypes
rule merge_vcfs:
    input:
        # vcfs=["data/calls/interval/{}-calls.vcf".format(re.sub("\t", "_", x)) for x in intervals],
        # intervals=config["intervals"]
        vcfs=["data/calls/chunk/{}-calls.vcf".format(x) for x in targets]
    output:
        lst=temp("tmp_vcf_files.txt"),
        vcf="data/calls/all-calls.vcf.gz"
    params:
        # jarpath=config["params"]["mark_duplicates"]["jarpath"],
        # java_heap=config["params"]["mark_duplicates"]["java_heap"]
    log:
        "log/cat_vcfs.log"
    # conda: "../env/bcftools.yaml"
    shell: """
        module load bcftools/1.12
        echo {input.vcfs} | tr ' ' '\n' | sort > {output.lst}
        bcftools concat --file-list {output.lst} > {output.vcf} 2> {log}
    """
    
# apply first-pass filters, to remove low-quality variants
# mostly done to reduce file size for loading into R
rule first_pass_filter:
    input:
        "data/calls/all-calls.vcf.gz"
    output:
        "data/calls/firstpass-filter.vcf.gz"
    shell: """
        module load bcftools/1.12
        bcftools view -i "NUMALT==1" {input} | \
        bcftools view -i "CIGAR=='1X'" | \
        bcftools view -i "QUAL>=20" | \
        bcftools view -i "MQM>=50" | \
        bcftools view -i "MQMR>=50" | \
        bcftools view -i "abs(MQM-MQMR)<10" | \
        bcftools view -i "RPPR<=20" | \
        bcftools view -i "RPP<=20" | \
        bcftools view -i "EPPR<=20" | \
        bcftools view -i "EPP<=20" | \
        bcftools view -i "SAP<=20" | \
        bcftools view -i "SRP<=20" | \
        bgzip > {output}
    """
    
# apply sample-specific filters
rule hap8_bed_vcf:
    input:
        "data/calls/firstpass-filter.vcf.gz"
    output:
        "data/calls/hap8.bed"
    shell: """
        Rscript scripts/hap8-filter.R {input} {output}
    """
    
rule hap8_vcf:
    input:
        "data/calls/firstpass-filter.vcf.gz",
        "data/calls/hap8.bed"
    output:
        "data/calls/hap8.vcf"
    shell: """
        bedtools intersect -header -a {input[0]} -b {input[1]} > {output}
    """
    
rule genotype_targets:
    input:
        vcf="data/calls/hap8.vcf",
        bed="data/calls/hap8.bed",
        bam=["data/merged/{}.bam".format(x) for x in geno.index.levels[0]],
        bai=["data/merged/{}.bam.bai".format(x) for x in geno.index.levels[0]],
        ref=config["genome"]
    output:
        "data/calls/bb-mf93s1-targets.vcf"
    log:
    params: "--hwe-priors--off --min-mapping-quality 20 --min-base-qality 20 --report-monomorphic --ploidy 4 --use-best-n-alleles 6 --pooled-continuous --only-use-input-alleles"
    shell: """
        freebayes \
            --fasta-reference {input.ref} \
            --bam {input.bam} \
            --targets {input.bed} \
            --input-alleles {input.vcf} \
            {params} \
            --vcf {output}
    """
