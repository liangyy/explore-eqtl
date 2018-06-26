import argparse
parser = argparse.ArgumentParser(prog='compute_ld.py', description='''
    Compute LD for the given list of variant pairs (r-squared)
''')
parser.add_argument('--pair_list', help='''
    variant pair list containing motif_snp and strong_eqtl with corresponding meta data
''')
parser.add_argument('--genotype', help='''
    individual level genotype data
''')
parser.add_argument('--output', help='''
    output file name
''')
args = parser.parse_args()

import os
import sys
import pandas as pd
import numpy as np
import subprocess
import io

def ld(snp, geno):
    snp1 = snp['motif_snp_id']
    snp2 = snp['strong_eqtl_id']
    geno1 = geno[geno.ix[:, 0] == snp1].ix[:, 1:].values[0]
    geno2 = geno[geno.ix[:, 0] == snp2].ix[:, 1:].values[0]
    r = np.corrcoef(geno1, geno2)[0, 1]
    return np.power(r, 2)

pattern = '$1=="{snp_id}"'
try:
    df = pd.read_table(args.pair_list, sep = '\t')
except pd.errors.EmptyDataError:
    cmd = 'touch {filename}'.format(filename = args.output)
    os.system(cmd)
    sys.exit()

snp_list = df['strong_eqtl_id'].tolist() + df['motif_snp_id'].tolist()

conds = []
for snp in set(snp_list):
    conds.append(pattern.format(snp_id = snp))
if_statement = ' || '.join(conds)
awk_cmd = '''awk '{{if({if_statement}) print $0}}' '''.format(if_statement = if_statement)
cmd = 'zcat < {geno} | {awk}'.format(geno = args.genotype, awk = awk_cmd)

subset_str = subprocess.check_output(cmd, shell = True)
handle = io.StringIO(subset_str.decode())
df_geno = pd.read_csv(handle, sep = '\t', header = None)

df['r2'] = df.apply(ld, axis = 1, geno = df_geno)
df.to_csv(args.output, sep = '\t', index = False)
