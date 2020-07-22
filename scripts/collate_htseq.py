'''
collate_htseq.py
- Collect the counts for all experiments and output in a .csv
'''

import pandas as pd
from collections import defaultdict

txt_2_gene_file = snakemake.input['txt_2_gene_file']
with open(txt_2_gene_file, 'r') as f:
    txt_2_gene = dict(line.strip().split('\t')[0:2] for line in f)

#gene may not have a symbol, for example if a synthetic spike-in RNA
with open(txt_2_gene_file, 'r') as f:
    gene_2_symbol = {}
    for line in f:
        fields = line.strip().split('\t')
        if len(fields) == 3:
            gene_2_symbol[fields[1]] = fields[2]
        else:
            gene_2_symbol[fields[1]] = fields[1]

units_file = snakemake.params['units_file']
units = pd.read_table(units_file, dtype=str).set_index(['sample', 'unit'], drop=False)
# enforce str in index
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])

exps = units.index.levels[0]
reps = units.index.levels[1]

abundance_files = snakemake.input['infiles']
sample_tuples = [(i.sample, i.unit) for i in units.itertuples()]

counts = [pd.read_csv(f, sep = '\t', names = ['gene', 'counts'], index_col = 'gene') for f in abundance_files]

df_dict = defaultdict(list)
for c, info in zip(counts, sample_tuples):
    experiment, rep = info
    df_dict[experiment].append(c)

exp_dfs = [pd.concat(df_dict[i], keys = reps, names = ['replicate']) for i in exps]
df = pd.concat(exp_dfs, keys = exps, names = ['experiment'])
df['symbol'] = df.index.get_level_values('gene').map(gene_2_symbol)
df.set_index('symbol', append = True, inplace = True)
df = df.reorder_levels(['gene', 'symbol', 'experiment', 'replicate'])
df.to_csv(snakemake.output['gene_table'])
