__author__ = "Mary Thompson"
__copyright__ = "Copyright 2020, Mary Thompson"
__email__ = "mary.thompson@bioch.ox.ac.uk"
__license__ = "MIT"
#Modified by MKT to use the gunzip -c command for readfiles

import os
from snakemake.shell import shell

extra = snakemake.params.get("extra", "")
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

sample = [snakemake.input] if isinstance(snakemake.input, str) else snakemake.input
n = len(sample)
assert n == 1 or n == 2, "input->sample must have 1 (single-end) or 2 (paired-end) elements."

outprefix = os.path.dirname(snakemake.output[0]) + "/"

shell(
    "bbduk.sh "
    "in={snakemake.input} "
    "out={snakemake.output.fastq} "
    "{snakemake.params.extra} "
    "{log}")
