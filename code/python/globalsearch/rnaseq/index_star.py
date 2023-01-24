#!/usr/bin/env python3

"""
Module to encapsulate STAR indexing.

It takes the genome directory and one or more FASTA files
and passes them to STAR to generate a genome index within
the genome directory.
If the index exists, it will skip the generation, to avoid
wasting time as this is a very costly step.
"""
import argparse
import os


DESCRIPTION = """index_star_salmon.py - Create genome index using STAR"""


####################### Create STAR index ###############################
### This should be specific for the organism
### Use the equation file maybe another script to create references
def create_genome_index(genome_dir, genome_fasta):
    index_command = ['STAR', '--runMode', 'genomeGenerate',
                     '--runThreadN', '32',
                     '--genomeDir', genome_dir,
                     '--genomeFastaFiles', genome_fasta,
                     '--genomeChrBinNbits', '16',
                     '--genomeSAindexNbases', '12']
    index_cmd = ' '.join(index_command)
    print(index_cmd)

    print ("\033[34m %s Indexing genome... \033[0m")
    if os.path.exists('%s/SAindex' % (genome_dir)):
        print ('Genome indexes exist. Not creating!')
    else:
        print ('Creating genome indexes')
        compl_proc = subprocess.run(index_command, check=True, capture_output=False, cwd=genome_dir)


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
