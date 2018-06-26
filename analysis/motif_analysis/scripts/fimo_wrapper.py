import argparse
parser = argparse.ArgumentParser(prog='fimo_wrapper.py', description='''
    Wrapper script to run fimo on a list of motif
''')
parser.add_argument('--motif_list', help='''
    motif list in TXT
''')
parser.add_argument('--motif_database', help='''
    motif database in MEME
''')
parser.add_argument('--fa', help='''
    sequence in FASTA
''')
parser.add_argument('--out_dir', help='''
    output directory name
''')
args = parser.parse_args()

import os

def get_motif_cmd(txt):
    cmd = ''
    with open(txt, 'r') as f:
        for i in f:
            i = i.strip()
            if i == '':
                continue
            cmd += '--motif {i} '.format(i = i)
    return cmd

motif_cmd = get_motif_cmd(args.motif_list)
if cmd.strip() != '':
    cmd = 'fimo --oc {out_dir} --verbosity 1 {motif_cmd} {motif_database} {fa}'.format(out_dir = args.out_dir, motif_cmd = motif_cmd, motif_database = args.motif_database, fa = args.fa)
else:
    cmd = 'touch {out_dir}/fimo.tsv'.format(out_dir = args.out_dir)
    sys.exit()
os.system(cmd)
