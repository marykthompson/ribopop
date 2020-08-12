'''
MKT 200715
- Write the reverse complement (i.e. target seqs) of the selected Dmel probes.
- Also write the reverse complement of the rRNA sequences.
- Use these sequences as input for Stellar local alignment.
'''

from Bio.Seq import Seq
import pandas as pd
import sys
from pyfaidx import Fasta

def write_probe_rcs(probe_file, outprefix):
    '''Write reverse complements --i.e. target seqs of the probes to fasta file'''

    probe_rc_outfile = f'{outprefix}_proberc.fa'
    probe_df = pd.read_csv(probe_file)
    probe_df['target_sequence'] = probe_df['sequence'].apply(lambda x: str(Seq(x).reverse_complement()))
    with open(probe_rc_outfile, 'w') as g:
        for i in probe_df.itertuples():
            g.write(f'>probe_{i.probe_num}\n{i.target_sequence}\n')

def write_rrna_rcs(txt_fasta, outprefix):
    '''Extract rRNA seqs and write reverse complements to a fasta file.'''

    outfile = f'{outprefix}_txtrc.fa'
    txt_seqs = Fasta(txt_fasta)
    ids = ['FBtr0346882', 'FBtr0346885']
    with open(outfile, 'w') as g:
        for i in ids:
            g.write(f'>{i}\n{str(Seq(str(txt_seqs[i])).reverse_complement())}\n')

def main():
    probe_file, txt_fasta, outprefix = sys.argv[1:]
    write_probe_rcs(probe_file, outprefix)
    write_rrna_rcs(txt_fasta, outprefix)

if __name__ == '__main__':
    main()
