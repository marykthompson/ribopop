__author__ = 'Mary Thompson'
__copyright__ = 'Copyright 2020, Mary Thompson'
__email__ = 'mary.thompson@bioch.ox.ac.uk'
__license__ = 'MIT'
#Modified by MKT from wrapper repo to use pre-built kallisto index
from snakemake.shell import shell

# Creating log
log = snakemake.log_fmt_shell(stdout=True, stderr=True)

# Placeholder for optional parameters
extra = snakemake.params.get('extra', '')

# Allowing for multiple FASTQ files
fastq = snakemake.input.get('fastq')
assert fastq is not None, 'input-> a FASTQ-file is required'
fastq = ' '.join(fastq) if isinstance(fastq, list) else fastq

shell(
    'kallisto quant '
    '{extra} '  # Optional parameters
    '--threads={snakemake.threads} '  # Number of threads
    '--index={snakemake.params.index} '  # Input file
    '--output-dir={snakemake.params.outdir} '  # Output directory
    '{fastq} '  # Input FASTQ files
    '{log}'  # Logging
)
