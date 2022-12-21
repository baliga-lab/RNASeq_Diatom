#!/usr/bin/env python3
import glob
import sys
import os
import re
from trim_galore import trim_galore, collect_trimmed_data, create_result_dirs

# data and results directories
RUN_DIR = "/proj/omics4tb2/Global_Search"
DATA_DIR = "%s/Pilot_Pass/X204SC21081158-Z01-F003/raw_data" % RUN_DIR
GENOME_DIR = "%s/reference_genomes/acerv_smic-reefGenomics" % RUN_DIR
TRANSCRIPTOME_FILE = "%s/Acerv_Smic-ReefGenomics_merged.fasta" % GENOME_DIR
#RESULTS_FOLDER = "%s/Pilot_Pass/X204SC21081158-Z01-F003/RNASeq_Analysis" % RUN_DIR
RESULTS_FOLDER = "/proj/omics4tb2/wwu/Global_Search/kallisto-out"

############# Functions ##############
####################### Run Kalisto ###############################
def run_kallisto(first_pair_group,second_pair_group,results_dir,kallisto_input_files):
    print
    print('\033[33mRunning kallisto! \033[0m')
    kallisto_cmd = '/users/sturkars/kallisto/kallisto quant -i %s/Acerv_Smic-reefGenomics_kallistoindex %s -o %s  -b 100 --bias -t 4 --rf-stranded' % (GENOME_DIR, kallisto_input_files, results_dir)
    print('Kallisto run command: "%s"' % kallisto_cmd)
    os.system(kallisto_cmd)

 ####################### Create Kallisto index ###############################
def kallisto_index():
    print
    print('\033[33mRunning kallisto index! \033[0m')
    kallistoindex_cmd = '/users/sturkars/kallisto/kallisto index -i %s/Acerv_Smic-reefGenomics_kallistoindex %s' %(GENOME_DIR, TRANSCRIPTOME_FILE)
    print('kallisto index command:%s' %(kallistoindex_cmd))
    #os.system(kallistoindex_cmd)

####################### Running the Pipeline ###############################
def run_pipeline(data_folder):
    folder_count = 1

    # Loop through each data folder
    folder_name = data_folder.split('/')[-1]
    print
    print
    print('\033[33mProcessing Folder: %s\033[0m' %(folder_name))

    # Get the list of first file names in paired end sequences
    first_pair_files = glob.glob('%s/*_1.fq*' %(data_folder))
    print(first_pair_files)

    # Program specific results directories
    data_trimmed_dir = "%s/%s/trimmed" % (RESULTS_FOLDER, folder_name)
    fastqc_dir = "%s/%s/fastqc_results" % (RESULTS_FOLDER, folder_name)
    results_dir = "%s/%s/results_Kallisto_Acerv_Smic-reefGenomics" %(RESULTS_FOLDER, folder_name)
    htseq_dir = "%s/htseqcounts" % RESULTS_FOLDER

    # Run create directories function to create directory structure
    create_result_dirs(data_trimmed_dir,fastqc_dir,results_dir, htseq_dir)

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
        trim_galore(first_pair_file,second_pair_file,folder_name,sample_id,file_ext,data_trimmed_dir,fastqc_dir)

        file_count += 1

        # Run folder level salmon analysis
        print("FILE_EXT: ", file_ext)
        first_pair_group,second_pair_group,kallisto_input_files = collect_trimmed_data(data_trimmed_dir,file_ext)
        print("FOUND FILES: ", kallisto_input_files)
        run_kallisto(first_pair_group,second_pair_group,results_dir,kallisto_input_files)

        folder_count += 1

    return data_trimmed_dir,fastqc_dir,results_dir

if __name__ == '__main__':
    folder_name = str(sys.argv[1])
    print(folder_name)
    data_folder = "%s/%s" % (DATA_DIR, folder_name)
    print(data_folder)
    data_trimmed_dir,fastqc_dir,results_dir = run_pipeline(data_folder)
