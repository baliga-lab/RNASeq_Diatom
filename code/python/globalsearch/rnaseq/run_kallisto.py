#!/usr/bin/env python3

"""
RNASeq Analysis pipeline using Kallisto
"""
import glob, sys, os, re, subprocess
from trim_galore import trim_galore, collect_trimmed_data, create_result_dirs
import argparse

# data and results directories
# These were the original examples
#RUN_DIR = "/proj/omics4tb2/Global_Search"
#DATA_DIR = "%s/Pilot_Pass/X204SC21081158-Z01-F003/raw_data" % RUN_DIR
#GENOME_DIR = "%s/reference_genomes/acerv_smic-reefGenomics" % RUN_DIR
#TRANSCRIPTOME_FILE = "%s/Acerv_Smic-ReefGenomics_merged.fasta" % GENOME_DIR

############# Functions ##############
####################### Run Kalisto ###############################
def run_kallisto(index_path, results_dir, pair_files):
    print('\033[33mRunning kallisto! \033[0m')
    # flatten the pair_files list into an input file list
    input_files = []
    for pair_file in pair_files:
        input_files.extend(pair_file)

    command = ['kallisto',
              'quant', '-i', index_path]
    command.extend(input_files)
    command.extend(['-o', results_dir,
                    '-b', '100', '--bias', '-t', '4', '--rf-stranded'])
    kallisto_cmd = ' '.join(command)
    print('Kallisto run command: "%s"' % kallisto_cmd)
    compl_proc = subprocess.run(command, check=True, capture_output=False)

 ####################### Create Kallisto index ###############################
def kallisto_index(index_path, transcriptome_path):
    print('\033[33mRunning kallisto index! \033[0m')
    if os.path.exists(index_path):
        print("index path '%s' exists, skipping the indexing" % index_path)
        return

    command = ['kallisto',
               'index',
               '-i', index_path,
               transcriptome_path]
    kallistoindex_cmd = ' '.join(command)
    print('kallisto index command: %s' % kallistoindex_cmd)
    compl_proc = subprocess.run(command, check=True, capture_output=False)

####################### Running the Pipeline ###############################
def run_pipeline(data_folder, results_folder, genome_dir, transcriptome_file, args):
    folder_count = 1

    # Loop through each data folder
    folder_name = data_folder.split('/')[-1]
    print('\033[33mProcessing Folder: %s\033[0m' %(folder_name))

    # Get the list of first file names in paired end sequences
    first_pair_files = glob.glob('%s/*_1.fq*' %(data_folder))
    print(first_pair_files)

    # Program specific results directories
    organism = os.path.basename(genome_dir)
    data_trimmed_dir = os.path.join(results_folder, folder_name, "trimmed")
    fastqc_dir = os.path.join(results_folder, folder_name, "fastqc_results")
    results_dir = os.path.join(results_folder, folder_name, organism)
    htseq_dir = os.path.join(results_folder, "htseqcounts")

    # Run create directories function to create directory structure
    create_result_dirs(data_trimmed_dir,fastqc_dir,results_dir, htseq_dir)

    index_path = os.path.join(results_folder, "%s_kallistoindex" % organism)
    kallisto_index(index_path, transcriptome_file)

    # Loop through each file and create filenames
    file_count = 1
    for first_pair_file in first_pair_files:
        first_file_name_full = first_pair_file.split('/')[-1]

        second_pair_file = first_pair_file.replace('_1.fq', '_2.fq')
        second_file_name_full = second_pair_file.split('/')[-1]
        file_ext = first_pair_file.split('.')[-1]

        print ('\033[32m Processing File: %s of %s (%s)\033[0m' %(file_count, len(first_pair_files), first_file_name_full ))

        first_file_name = re.split('.fq|.fq.gz',first_file_name_full)[0]
        second_file_name = re.split('.fq|.fq.gz',second_file_name_full)[0]
        print('first_file_name:%s, second_file_name:%s' %(first_file_name,second_file_name))

        # Collect Sample attributes
        exp_name = folder_name
        print("exp_name: %s" %(exp_name))
        lane = first_file_name.split("_")[-1]
        print("Lane: %s" %(lane))
        sample_id = re.split('.fq|.fq.gz', first_file_name)[0]
        print("sample_id: %s" %(sample_id))

        # 01. Run TrimGalore
        # WW: Activate/deactivate for debugging
        trim_galore(first_pair_file, second_pair_file, folder_name, sample_id,
                    file_ext, data_trimmed_dir, fastqc_dir)
        file_count += 1

        # Run folder level salmon analysis
        first_pair_group, second_pair_group, pair_files = collect_trimmed_data(data_trimmed_dir,file_ext)
        run_kallisto(index_path, results_dir, pair_files)

        folder_count += 1

    return data_trimmed_dir,fastqc_dir,results_dir


DESCRIPTION = """run_kallisto.py - run Kallisto pipeline"""

if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=DESCRIPTION)
    parser.add_argument('genomedir', help='genome directory')
    parser.add_argument('dataroot', help="parent of input directory")
    parser.add_argument('indir', help="input directory without the input dir part")
    parser.add_argument('transcriptome_file', help="path to transcriptome_file")
    parser.add_argument('outdir', help='output directory')
    args = parser.parse_args()

    print("Processing directory %s" % args.indir)
    data_folder = os.path.join(args.dataroot, args.indir)
    data_trimmed_dir, fastqc_dir, results_dir = run_pipeline(data_folder, args.outdir, args.genomedir,
                                                             args.transcriptome_file, args)
