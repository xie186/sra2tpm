import os
sample_list = config["sample_list"]
outdir = config["outdir"]

def read_meta_info(meta_file):
    dict_sam_meta = {}
    with open(meta_file) as meta:
        next(meta) ## skip header line
        for line in meta:
            line = line.rstrip("\n")
            ele = line.split("\t")
            dict_sam_meta[ele[0]] = {}
            #Alias   Sample  RUNs    LibType 
            srr = ele[2].split(",")
            #print(srr)
            dict_sam_meta[ele[0]]["RUNs"] = srr
            dict_sam_meta[ele[0]]["LibType"] = ele[3]
            dict_sam_meta[ele[0]]["SE"] = []
            dict_sam_meta[ele[0]]["PE1"] = []
            dict_sam_meta[ele[0]]["PE2"] = []
            for run in srr:
                if(ele[3] == "SE"):
                    dict_sam_meta[ele[0]]["SE"].append(os.path.join(config["inputdir"], run+".fastq.gz"))
                if(ele[3] == "PE"):
                    dict_sam_meta[ele[0]]["SE"].append("I_dont_exist.fq.gz")
                    dict_sam_meta[ele[0]]["PE1"].append(os.path.join(config["inputdir"], run+"_R1.fastq.gz"))
                    dict_sam_meta[ele[0]]["PE2"].append(os.path.join(config["inputdir"], run+"_R2.fastq.gz"))
    return dict_sam_meta
dict_sam_meta = read_meta_info(sample_list)
print(dict_sam_meta)


SAMPLE = list(dict_sam_meta.keys())
print(SAMPLE)
rule all:
    input:
        expand(os.path.join(outdir, "merged_fastq/{sample}.done"), sample=SAMPLE),
        expand(os.path.join(outdir, "trim_galore/{sample}.log"), sample=SAMPLE),
        expand(os.path.join(outdir, "kallisto_quant/{sample}.log"), sample=SAMPLE),
        expand(os.path.join(outdir, "hisat2_align/{sample}.log"), sample=SAMPLE),
        expand(os.path.join(outdir, "string_expr/{sample}.gtf"), sample=SAMPLE),
        expand(os.path.join(outdir, "string_expr/{sample}.tab"), sample=SAMPLE),

rule hisat2_index:
    input:
        fa = config["reference_fasta"],
    params:
        index_prefix = os.path.join(outdir, "hisat2_index/index/hisat2.index"),
    output:
        log = os.path.join(outdir, "hisat2_index/index/hisat2.index.log"),
    shell:
        """
hisat2-build {input.fa} {params.index_prefix} > {output.log}
"""


ruleorder: merge_readSE > merge_readPE
rule merge_readSE:
    input:
        fastq = lambda wildcards: dict_sam_meta[wildcards.sample]["SE"]
    params:
        base_name = "{sample}"
    output:
        fq = os.path.join(outdir, "merged_fastq/{sample}.fastq.gz"),
        done = os.path.join(outdir, "merged_fastq/{sample}.done"),
    shell:
        """
zcat {input} |gzip - > {output.fq}
touch {output.done}
"""

rule merge_readPE:
    input:
        fq1 = lambda wildcards: dict_sam_meta[wildcards.sample]["PE1"],
        fq2 = lambda wildcards: dict_sam_meta[wildcards.sample]["PE2"],
    params:
        base_name = "{sample}"
    output:
        fq1 = os.path.join(outdir, "merged_fastq/{sample}_R1.fastq.gz"),
        fq2 = os.path.join(outdir, "merged_fastq/{sample}_R2.fastq.gz"),
        done = os.path.join(outdir, "merged_fastq/{sample}.done"),
    shell:
        """
zcat {input.fq1} |gzip - > {output.fq1}
zcat {input.fq2} |gzip - > {output.fq2}
touch {output.done}
"""

ruleorder: trimSE > trimPE
rule trimSE:
    input:
        fastq = rules.merge_readSE.output.fq
    params:
        length = config["trim_length"],
        sample = "{sample}"
    output:
        tmp_outdir = directory(os.path.join(outdir, "trim_galore/{sample}/")),
        log = os.path.join(outdir, "trim_galore/{sample}.log"),
        fq = os.path.join(outdir, "trim_galore/{sample}/{sample}_trimmed.fq.gz"),
    shell:
        """
mkdir -p {output.tmp_outdir}
trim_galore --output_dir {output.tmp_outdir} --basename {params.sample} {input.fastq} > {output.log}
touch {output.fq}
"""

rule trimPE:
    input:
        fq1 = rules.merge_readPE.output.fq1,
        fq2 = rules.merge_readPE.output.fq2,
    params:
        length = config["trim_length"],
        sample = "{sample}"
    output:
        tmp_outdir = directory(os.path.join(outdir, "trim_galore/{sample}/")),
        log = os.path.join(outdir, "trim_galore/{sample}.log"),
        fq1 = os.path.join(outdir, "trim_galore/{sample}/{sample}_val_1.fq.gz"),
        fq2 = os.path.join(outdir, "trim_galore/{sample}/{sample}_val_2.fq.gz"),
    shell:
        """
mkdir -p {output.tmp_outdir}
trim_galore --paired --length {params.length} --output_dir {output.tmp_outdir} --basename {params.sample} {input.fq1} {input.fq2} > {output.log}
"""

rule kallisto_index:
    input:
        cdna = config["cDNA"]
    output:
        index = os.path.join(outdir, "trim_galore/index/kallisto.index"),
        log = os.path.join(outdir, "trim_galore/index/kallisto.index.log"),
    shell:
        """
kallisto index -i {output.index} {input.cdna} > {output.log}
"""

ruleorder: kallisto_quant_SE > kallisto_quant_PE
rule kallisto_quant_SE:
    input:
        fq=rules.trimSE.output.fq,
        idx=rules.kallisto_index.output.index
    output:
        tmp_outdir = directory(os.path.join(outdir, "kallisto_quant/{sample}/")),
        log = os.path.join(outdir, "kallisto_quant/{sample}.log")
    params:
        extra=config["kallisto_params"],
        extraSE=config["kallisto_params4SE"]
    shell:
        """
kallisto quant --single {params.extraSE} -i {input.idx} -o {output.tmp_outdir} {params.extra} {input.fq} 2> {output.log}
"""

rule kallisto_quant_PE:
    input:
        fq1=rules.trimPE.output.fq1,
        fq2=rules.trimPE.output.fq2,
        idx=rules.kallisto_index.output.index
    output:
        tmp_outdir = directory(os.path.join(outdir, "kallisto_quant/{sample}/")),
        log = os.path.join(outdir, "kallisto_quant/{sample}.log")
    params:
        extra=config["kallisto_params"]
    shell:
        """
kallisto quant -i {input.idx} -o {output.tmp_outdir} {params.extra} {input.fq1} {input.fq2} 2> {output.log}
"""

ruleorder: hisat2_align_SE > hisat2_align_PE
rule hisat2_align_SE:
    input:
        fq=rules.trimSE.output.fq,
        idx=rules.hisat2_index.output.log
    params:
        index_prefix = os.path.join(outdir, "hisat2_index/index/hisat2.index"),
    output:
        bam = os.path.join(outdir, "hisat2_align/{sample}.bam"),
        log = os.path.join(outdir, "hisat2_align/{sample}.log")
    shell:
        """
touch {input.idx}
hisat2 --new-summary --summary-file {output.log} -x {params.index_prefix} -U {input.fq} |samtools view -bS - | samtools sort - -o {output.bam}
samtools index {output.bam}
"""

rule hisat2_align_PE:
    input:
        fq1=rules.trimPE.output.fq1,
        fq2=rules.trimPE.output.fq2,
        idx=rules.hisat2_index.output.log
    output:
        bam = os.path.join(outdir, "hisat2_align/{sample}.bam"),
        log = os.path.join(outdir, "hisat2_align/{sample}.log")
    params:
        index_prefix = os.path.join(outdir, "hisat2_index/index/hisat2.index"),
    shell:
        """
touch {input.idx}
hisat2 --new-summary --summary-file {output.log} -x {params.index_prefix} -1 {input.fq1} -2 {input.fq2} |samtools view -bS - | samtools sort - -o {output.bam}
samtools index {output.bam}
"""


rule stringtie_expr:
    input:
        bam = os.path.join(outdir, "hisat2_align/{sample}.bam"),
        gtf = config["gtf"]
    threads: config["stringtie_thread_number"]
    output:
        gtf = os.path.join(outdir, "string_expr/{sample}.gtf"),
        tab = os.path.join(outdir, "string_expr/{sample}.tab"),
    shell:
        """
stringtie {input.bam} -p {threads} -G {input.gtf} -A {output.tab} -B -e -o {output.gtf}
"""
