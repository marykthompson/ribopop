#MKT 200715
#create a fasta file containing only a subset of transcripts
import sys
import gffutils
from pyfaidx import Fasta

def write_txts(db_file, gene_file, txt_file, outname):
    '''write fasta file with the transcripts in the  gene file'''

    db = gffutils.FeatureDB(db_file)
    txt_seqs = Fasta(txt_file)

    with open(gene_file, 'r') as f:
        genes = set([i.strip('\n') for i in f.readlines()])

    with open(f'{outname}.fa', 'w') as g:
        for record in txt_seqs:
            gene = next(db.parents(record.name, featuretype = 'gene')).id
            if gene in genes:
                g.write('>%s\n%s\n' % (record.name, str(record)))

def main():
    db_file, gene_file, txt_file, outname = sys.argv[1:]
    write_txts(db_file, gene_file, txt_file, outname)

if __name__ == '__main__':
    main()
