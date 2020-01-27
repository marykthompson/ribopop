# Snakemake workflow: ribopop_differential_expression

This workflow performs a differential expression analysis between Ribo-pop
treated and untreated RNA-Seq libraries. Transcript abundance is estimated with
Kallisto and gene-level differential expression analysis is performed with DESeq2.
Reads are also aligned with STAR and subjected to QC analysis with RSeQC.
## Authors

* Mary Thompson (@marykthompson)
