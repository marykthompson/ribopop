#MKT 200715
#create a fasta file containing only the longest transcript of each gene
import sys
import gffutils
from pyfaidx import Fasta

def find_longest_txts(db_file):
    '''return a set of the longest transcripts for each gene'''

    db = gffutils.FeatureDB(db_file)
    longest_txts = set()
    genes = db.features_of_type('gene')
    for g in genes:
        txts = db.children(g, featuretype = 'transcript')
        l = []
        for t in txts:
            length = db.children_bp(t, child_featuretype='exon')
            l.append((length, t.id))
        longest = sorted(l)[-1][1]
        longest_txts.add(longest)
    return longest_txts

def write_longest_txts(longest_txts, txt_file, outname):
    '''write fasta file with the longest transcript for each gene'''
    
    found = 0
    notfound = 0
    notfoundset = set()
    txt_seqs = Fasta(txt_file)
    with open(f'{outname}.fa', 'w') as g:
        for record in longest_txts:
            if record in txt_seqs:
                g.write('>%s\n%s\n' % (record, str(txt_seqs[record])))
                found += 1
            else:
                notfound += 1
                notfoundset.add(record)
    print(f'found: {found}')
    print(f'notfound: {notfound}')
    with open(f'{outname}_notfound.txt', 'w') as k:
        for t in notfoundset:
            k.write(f'{t}\n')

def main():

    db_file, txt_file, outname = sys.argv[1:]
    longest_txts = find_longest_txts(db_file)
    write_longest_txts(longest_txts, txt_file, outname)

if __name__ == '__main__':
    main()
