This workflow performs differential expression analysis on single- or paired-end RNA-seq data.
After adapter removal with `Cutadapt <http://cutadapt.readthedocs.io>`_, or BBtrim,
transcript abundance was estimated with Kallisto. Differential expression analysis
was conducted with `DESeq2 <https://bioconductor.org/packages/release/bioc/html/DESeq2.html>`_
using the transcript import function.
