shell.executable("bash")
import re
import pandas as pd

units = pd.read_csv(config["units"], index_col=["sample", "unit"], dtype=str, sep = "\t")
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])

samples = {}
for i in units.index.levels[0]:
    samples[i] = ["data/dedup/{}-{}.bam".format(i, j) for j in units.loc[i].index]

def is_single_end(sample,unit):
    return pd.isnull(units.loc[(sample, unit), "fq2"])

def get_fastq(wildcards):
    return "data/rnaseq/"+units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()

def get_trimmed(wildcards):
    if not is_single_end(**wildcards):
        return expand("data/trimmed/{sample}-{unit}-{group}.fastq.gz",
            group=[1,2], **wildcards)
    return "data/trimmed/{sample}-{unit}.fastq.gz".format(**wildcards)

rule all:
    input:
        re.sub(".fa", ".1.ht2", config["genome"]),
        ["data/merged/{}.bam".format(x) for x in units.index.levels[0]],
        # "RP_braker3/braker.gtf"
        "braker3.done"

rule cutadapt_pe:
    input:
        get_fastq
    output:
        fastq1="data/trimmed/{sample}-{unit}-1.fastq.gz",
        fastq2="data/trimmed/{sample}-{unit}-2.fastq.gz",
        qc="data/trimmed/{sample}-{unit}.qc.txt"
    threads: config["params"]["cutadapt-pe"]["threads"]
    params:
        "-a {} -A {} -q {} -m {}".format(config["adapter"], config["adapter"], config["params"]["cutadapt-pe"]["qual"], config["params"]["cutadapt-pe"]["minlength"])
    log:
        "log/cutadapt/{sample}-{unit}.log"
    shell: """
        cutadapt \
            {params} \
            -j {threads} \
            -o {output.fastq1} \
            -p {output.fastq2} \
            {input} \
            > {output.qc}
            2> {log}
    """

rule hisat2_index:
    input:
        config["genome"]
    output:
        re.sub(".fa", ".1.ht2", config["genome"])
    log:
        "log/hisat2/index.log"
    params:
        outpref=re.sub(".fa", "", config["genome"])
    conda:
        "env/hisat2.yaml"
    shell: """
        module load hisat2 2> {log}
        hisat2-build {input} {params.outpref} 2>> {log}
    """
    
rule hisat2:
    input:
        idx=re.sub(".fa", ".1.ht2", config["genome"]),
        reads=get_trimmed
    output:
        "data/aligned_reads/{sample}-{unit}.bam"
    log:
        "log/hisat2/{sample}-{unit}.log"
    params:
        pref=re.sub(".fa", "", config["genome"]),
        etc=config["params"]["hisat2"]["etc"]
    threads:
        config["params"]["hisat2"]["threads"]
    conda:
        "env/hisat2.yaml"
    shell: """
        module load hisat2
        hisat2 -x {params.pref} {params.etc} -1 {input.reads[0]} -2 {input.reads[1]} --threads {threads}
 2>> {log} | samtools sort -@ {threads} -o {output} - 2>> {log}
    """

rule mark_duplicates:
    input:
        "data/aligned_reads/{sample}-{unit}.bam"
    output:
        bam="data/dedup/{sample}-{unit}.bam",
        metrics="data/dedup/{sample}-{unit}-metrics.txt"
    log:
        "logs/picard/{sample}-{unit}.log"
    params:
        jarpath=config["params"]["mark_duplicates"]["jarpath"],
        java_heap=config["params"]["mark_duplicates"]["java_heap"],
        opt=config["params"]["mark_duplicates"]["opt"]
    threads: 2
    shell:
        "java {params.java_heap} -jar {params.jarpath} MarkDuplicates INPUT={input} "
        "OUTPUT={output.bam} METRICS_FILE={output.metrics} {params.opt} 2> {log}"
        
rule samtools_merge:
    input:
        lambda x: samples[x.sample]
    output:
        "data/merged/{sample}.bam"
    conda:
        "env/samtools.yaml"
    shell: """
        samtools merge {output} {input}
    """
    
rule edta_ltr:
    input:
        asm=config["genome"],
        cds="ref/DM_1-3_516_R44_potato.v6.1.hc_gene_models.cds.fa"
    output:
        "../EDTA/RP_tetra_homcov140.p_ctg.fa.mod.EDTA.raw/RP_tetra_homcov140.p_ctg.fa.mod.LTR.intact.fa",
        "../EDTA/RP_tetra_homcov140.p_ctg.fa.mod.EDTA.raw/RP_tetra_homcov140.p_ctg.fa.mod.LTR.intact.gff3",
        "../EDTA/RP_tetra_homcov140.p_ctg.fa.mod.EDTA.raw/RP_tetra_homcov140.p_ctg.fa.mod.LTR.raw.fa"
    log:
        "RP_tetra_LTR.err"
    params: "ltr"
    threads: 32
    shell: """
        perl ../EDTA/EDTA_raw.pl \
        --genome {input.asm} \
        --species others \
        --type {params} \
        --threads {threads} \
        --cds {input.cds} \
        > {output}
        2> {log}
    """
    
rule edta_tir:
    input:
        asm=config["genome"],
        cds="ref/DM_1-3_516_R44_potato.v6.1.hc_gene_models.cds.fa"
    output:
        "../EDTA/RP_tetra_homcov140.p_ctg.fa.mod.EDTA.raw/RP_tetra_homcov140.p_ctg.fa.mod.TIR.intact.fa",
        "../EDTA/RP_tetra_homcov140.p_ctg.fa.mod.EDTA.raw/RP_tetra_homcov140.p_ctg.fa.mod.TIR.intact.gff3",
        "../EDTA/RP_tetra_homcov140.p_ctg.fa.mod.EDTA.raw/RP_tetra_homcov140.p_ctg.fa.mod.TIR.raw.fa"
    log:
        "RP_tetra_TIR.err"
    params: "tir"
    threads: 32
    shell: """
        perl ../EDTA/EDTA_raw.pl \
        --genome {input.asm} \
        --species others \
        --type {params} \
        --threads {threads} \
        --cds {input.cds} \
        > {output}
        2> {log}
    """
    
rule edta_helitron:
    input:
        asm=config["genome"],
        cds="ref/DM_1-3_516_R44_potato.v6.1.hc_gene_models.cds.fa"
    output:
        "../EDTA/RP_tetra_homcov140.p_ctg.fa.mod.EDTA.raw/RP_tetra_homcov140.p_ctg.fa.mod.Helitron.intact.fa",
        "../EDTA/RP_tetra_homcov140.p_ctg.fa.mod.EDTA.raw/RP_tetra_homcov140.p_ctg.fa.mod.Helitron.intact.gff3",
        "../EDTA/RP_tetra_homcov140.p_ctg.fa.mod.EDTA.raw/RP_tetra_homcov140.p_ctg.fa.mod.Helitron.raw.fa"
    log:
        "RP_tetra_helitron.err"
    params: "helitron"
    threads: 32
    shell: """
        perl ../EDTA/EDTA_raw.pl \
        --genome {input.asm} \
        --species others \
        --type {params} \
        --threads {threads} \
        --cds {input.cds} \
        > {output}
        2> {log}
    """
    
rule edta_combine:
    input:
        asm=config["genome"],
        cds="ref/DM_1-3_516_R44_potato.v6.1.hc_gene_models.cds.fa",
        ltr="../EDTA/RP_tetra_homcov140.p_ctg.fa.mod.EDTA.raw/RP_tetra_homcov140.p_ctg.fa.mod.LTR.raw.fa",
        tir="../EDTA/RP_tetra_homcov140.p_ctg.fa.mod.EDTA.raw/RP_tetra_homcov140.p_ctg.fa.mod.TIR.raw.fa",
        hel="../EDTA/RP_tetra_homcov140.p_ctg.fa.mod.EDTA.raw/RP_tetra_homcov140.p_ctg.fa.mod.Helitron.raw.fa"
    output:
        "../EDTA/RP_tetra_homcov140.p_ctg.fa.mod.EDTA.TElib.fa"
    log:
        "EDTA_combine.log"
    params:
    threads: 12
    shell: """
        perl ../EDTA/EDTA.pl \
        --overwrite 0 \
        --genome {input.asm} \
        --species others \
        --cds {input.cds} \
        2> {log}
    """
    
rule RepeatMasker:
    input:
        asm=config["genome"],
        lib="../EDTA/RP_tetra_homcov140.p_ctg.fa.mod.EDTA.TElib.fa"
    output:
        directory("RP_tetra_homcov140_rmout"),
        "repeatmasker.done"
        # "RP_tetra_homcov140_rmout/RP_tetra_homcov140_purged.p_ctg.fa.masked"
    log:
        "RP_rm.err"
    threads: 32
    shell: """
        module load repeatmasker/4.0.9.p2
        RepeatMasker \
        -lib {input.lib} \
        -pa {threads} \
        -gff \
        -no_is \
        -xsmall \
        -dir {output[0]} \
        {input.asm}
        2> {log}
        touch {output[1]}
    """
    
# todo modify paths in RP_braker3.sh
rule braker3:
    input:
        bam=["data/merged/{}.bam".format(x) for x in units.index.levels[0]],
        rep="repeatmasker.done",
        pro="ref/DM_M82_ODB10.fa"
    
    output:
        "braker3.done"
    params:
    threads: 36
    shell: """
        module load singularity
        bash scripts/RP_braker3.sh
        touch {output}
    """