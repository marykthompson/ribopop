def get_fastq(wildcards):
    return units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()


rule cutadapt_pe:
    input:
        get_fastq
    output:
        fastq1="trimmed/{sample}-{unit}.1.fastq.gz",
        fastq2="trimmed/{sample}-{unit}.2.fastq.gz",
        qc="trimmed/{sample}-{unit}.qc.txt"
    params:
        "-a {} {}".format(config["trimming"]["adapter"], config["params"]["cutadapt-pe"])
    log:
        "logs/cutadapt/{sample}-{unit}.log"
    wrapper:
        "0.17.4/bio/cutadapt/pe"


rule cutadapt:
    input:
        get_fastq
    output:
        fastq="trimmed/{sample}-{unit}.fastq.gz",
        qc="trimmed/{sample}-{unit}.qc.txt"
    params:
        "-a {} {}".format(config["trimming"]["adapter"], config["params"]["cutadapt-se"])
    log:
        "logs/cutadapt/{sample}-{unit}.log"
    wrapper:
        "0.17.4/bio/cutadapt/se"


rule bbtrim:
    input:
        get_fastq
    output:
        fastq="bb_trimmed/{sample}-{unit}.fastq.gz"
    params:
        extra="ref={} {}".format(','.join(config["trimming"]["contaminant_files"]), config["params"]["bbtrim"])
        #contaminant_files="{}".format(','.join(config["trimming"]["contaminant_files"])),
        #extra="ref={} {}".format(contaminant_files, config["params"]["bbtrim"])
    log:
        "logs/bbtrim/{sample}-{unit}.log"
    wrapper:
        "file:///Users/maryk.thompson/Desktop/Davislab/comp_labbook_backup/data/computational_projects/rateseq_pipelines/rna-quantseq-deseq/wrappers/wrapper_bbtrim"
