def get_fq(wildcards):
    if config["trimming"]["skip"]:
        # no trimming, use raw reads
        return units.loc[(wildcards.sample, wildcards.unit), ["fq1", "fq2"]].dropna()
    else:
        # yes trimming, use trimmed data
        trimmed_dir = config["trimming"]["trimmed_dir"]
        if not is_single_end(**wildcards):
            # paired-end sample
            return expand("{}/{sample}-{unit}.{group}.fastq.gz".format(trimmed_dir),
                          group=[1, 2], **wildcards)
        # single end sample
        return "{}/{sample}-{unit}.fastq.gz".format(trimmed_dir, **wildcards)


rule align:
    input:
        sample=get_fq
    output:
        # see STAR manual for additional output files
        "star/{sample}-{unit}/Aligned.out.bam",
        "star/{sample}-{unit}/ReadsPerGene.out.tab"
    log:
        "logs/star/{sample}-{unit}.log"
    params:
        # path to STAR reference genome index
        index=config["ref"]["index"],
        # optional parameters
        extra="--quantMode GeneCounts --sjdbGTFfile {} {}".format(
              config["ref"]["annotation"], config["params"]["star"])
    threads: 24
    wrapper:
        "file:///Users/maryk.thompson/Desktop/Davislab/comp_labbook_backup/data/computational_projects/rateseq_pipelines/rna-quantseq-deseq/wrappers/wrapper_star"
