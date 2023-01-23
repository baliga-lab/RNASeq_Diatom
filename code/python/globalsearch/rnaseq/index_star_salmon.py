#!/usr/bin/env python3

import argparse
import os
from .run_star_salmon import create_genome_index


DESCRIPTION = """index_star_salmon.py - Create genome index for SLURM"""

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=DESCRIPTION)
    parser.add_argument('genomedir', help='genome directory')
    parser.add_argument('--genome_fasta', help='genome FASTA file')
    args = parser.parse_args()
    if args.genome_fasta is not None and os.path.exists(args.genome_fasta):
        genome_fasta = args.genome_fasta
    else:
        genome_fasta = glob.glob('%s/*.fasta' % (args.genomedir))[0]
    create_genome_index(args.genomedir, genome_fasta)
