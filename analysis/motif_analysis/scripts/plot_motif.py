import argparse
parser = argparse.ArgumentParser(prog='plot_motif.py', description='''
    Plot motifs appear in all files
''')
parser.add_argument('--inputs', nargs = '+', help='''
    a list of input files containing motif name
''')
parser.add_argument('--dbnames', help='''
    name of motif db
''')
parser.add_argument('--dbpaths', help='''
    path of motif db
''')
parser.add_argument('--out_dir', help='''
    output file name
''')
args = parser.parse_args()

from ntpath import basename
import pandas
import pandas as pd

db_dic = dict(zip(args.dbnames, args.dbpaths))

for f in args.inputs:
    (from_info, name_info) = basename(f).split('.')[1 : 3]
    if from_info == 'from_enrichment':
        db = name_info
    elif from_info == 'from_predefined_list':
        db = name_info.split('__')[1]
    try:
        df = pd.read_table(args.f, sep = '\t')
    except pd.errors.EmptyDataError:
        continue
    motifs = df['']
