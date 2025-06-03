import os, re
import pandas as pd
from Bio import SeqIO
shell.executable("bash")

configfile: "config.yaml"
units = pd.read_csv(config["units"], index_col =["sample", "unit"], dtype=str, sep="\t")
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])

short = units.loc[units["readlen"] == "short"]
short.index = short.index.remove_unused_levels()

long  = units.loc[units["readlen"] == "long"]
long.index = long.index.remove_unused_levels()

long_samples = {}
for i in long.index.levels[0]:
    long_samples[i] = ["{}-{}".format(i, j) for j in long.loc[i].index]


def get_long(wildcards):
    return "data/reads/"+long.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()

rule all:
    input:
        "ref/unpurged/RP_tetra_homcov140.p_ctg.fasta",
        "ref/RP_tetra_homcov140_purged.p_ctg.fasta" # this makes asm "RP_tetra_homcov140" throughout
    
rule hifiasm:
    input:
        ["data/reads/{}".format(x) for x in long.loc[long["readlen"] == "long", "fq1"]]
    output:
        "data/hifiasm/{asm}.p_ctg.gfa",
        "ref/unpurged/{asm}.p_ctg.fasta"
    log:
        "log/hifiasm/{asm}.log"
    threads: 
        config["params"]["hifiasm"]["threads"] # 32
    params:
        "--hom-cov 140 --n-hap 4 --primary"
    shell: """
        hifiasm -t {threads} {params} -o data/hifiasm/{wildcards.asm} {input}
        grep '^S' {output[0]} | cut -f 2,3 | awk '{{print ">"$1"\n"$2}}' > {output[1]}
    """
    
# hifiasm -o RP_tetra_homcov140 --primary -t 32 --n-hap 4 --hom-cov 140 m64069_211118_031642.hifi_reads.fastq.gz m64069_211125_193709.hifi_reads.fastq.gz m84066_230512_192902_s4.hifi_reads.default.fastq.gz
    
rule index_asm:
    input:
        "ref/unpurged/{asm}.p_ctg.fasta"
    output:
        "ref/unpurged/{asm}.p_ctg.fasta.bwt"
        "ref/unpurged/{asm}.p_ctg.fasta.fai"
        "ref/unpurged/{asm}.p_ctg.fasta.genome"
    log:
        "log/index_asm/{asm}.log"
    shell: """
        bwa index {input} 2> {log}
        samtools faidx {input}
        cut -f 1-2 {output[1]} > {output[2]}
    """
    
rule map_long_to_asm:
    input:
        reads=get_long,
        asm="ref/unpurged/RP_tetra_homcov140.p_ctg.fasta"
    output:
        "data/aligned_long/{sample}-{unit}_raw.bam",
        "data/aligned_long/{sample}-{unit}.bam"
    log:
        "log/minimap2/{sample}-{unit}.log"
    threads: 20
    params:
        rg="'@RG\\tID:{unit}\\tSM:{sample}'"
    shell: """
        module load minimap2 2> {log}
        minimap2 -R {params.rg} -t {threads} -ax map-hifi {input.asm} {input.reads} 2>> {log} | samtools sort -@ {threads} -o {output[0]} - 2>> {log}
        samtools view -q 1 -F 3840 -b {output[0]} > {output[1]}
    """
    
rule samtools_merge_long:
    input:
        lambda x: ["data/aligned_long/{}.bam".format(x) for x in long_samples[x.sample]]
    output:
        "data/merge_long/{sample}.bam"
    shell:
        "samtools merge {output} {input}"
    
rule purge_haplotigs_hist:
    input:
        "data/merge_long/{sample}.bam",
        "ref/unpurged/RP_tetra_homcov140.p_ctg.fasta"
    output:
        "{sample}.bam.gencov"
    log:
    threads: 12
    params:
    shell: """
        module load purge_haplotigs
        source activate purge_haplotigs
        purge_haplotigs hist -b {input[0]} -g {input[1]} -t {threads}
    """
    
rule purge_haplotigs_cov:
    input:
        "{sample}.bam.gencov"
    output:
        "{sample}_coverage_stats.csv"
    log:
    threads: 4
    params:
        "-l 15 -m 52 -h 190"
    shell: """
        module load purge_haplotigs
        source activate purge_haplotigs
        purge_haplotigs cov -i {input} {params} -o {output}
    """
    
# LEFT OFF HERE - see 2023-07-05 notes for code to finish this part off
# this does a whole minimap2 alignment
rule purge_haplotigs_purge:
    input:
        asm="ref/unpurged/{asm}.p_ctg.fasta",
        stat="PI310467_coverage_stats.csv"
    output:
        "ref/{asm}_purged.p_ctg.fasta"
    log:
    threads: 12
    params:
    shell: """
        purge_haplotigs purge -g {input.asm} -c {input.stat} -o {wildcards.asm}_purged -t {threads}
    """
