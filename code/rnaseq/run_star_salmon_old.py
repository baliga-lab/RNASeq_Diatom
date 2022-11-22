#!/usr/bin/env python3

#############################################################
##### RNASeq Analysis Pipeline with STAR                #####
##### Last update: 01/13/2022 Serdar Turkarslan         #####
##### Institute for Systems Biology                     #####
############################################################
import glob, sys, os, string, datetime, re
import argparse

DESCRIPTION = """run_STAR_SALMON.py - run STAR and Salmon"""

####################### Create results directories ###############################
def create_dirs(data_trimmed_dir, fastqc_dir, results_dir, htseq_dir):
    dirs = [data_trimmed_dir, fastqc_dir, results_dir, htseq_dir]
    for dir in dirs:
        # create results folder
        #print(dir)
        if not os.path.exists('%s' %(dir)):
            os.makedirs('%s' %(dir))
        else:
            print('\033[31m %s directory exists. Not creating. \033[0m' %(dir))


####################### Trimgalore for quality and trimming ###############################
def trim_galore(first_pair_file, second_pair_file, folder_name, sample_id, file_ext, data_trimmed_dir,
                fastqc_dir):
    #print("1stpair:%s, 2ndpair:%s, folder_name:%s, sample_name:%s")%(first_pair_file,second_pair_file,folder_name,sample_name)
    print
    print ("\033[34m Running TrimGalore \033[0m")
    # create sample spepcific trimmed directory
    if not os.path.exists('%s' %(data_trimmed_dir)):
        os.makedirs('%s' %(data_trimmed_dir))
    # create sample spepcific fastqcdirectory
    if not os.path.exists('%s' %(fastqc_dir)):
        os.makedirs('%s' %(fastqc_dir))
    # run Command
    cmd = 'trim_galore --fastqc_args "--outdir %s/" --paired --output_dir %s/ %s %s' %(fastqc_dir,data_trimmed_dir,first_pair_file, second_pair_file)
    print
    print( '++++++ Trimgalore Command:', cmd)
    print
    os.system(cmd)


####################### Collect trimmed data files ###############################
def collect_trimmed_data(data_trimmed_dir, file_ext):
    # define result files
    if file_ext == "gz":
        first_pair_trimmed = glob.glob('%s/*_val_1.fq.gz'%(data_trimmed_dir))
        second_pair_trimmed = glob.glob('%s/*_val_2.fq.gz'%(data_trimmed_dir))
    else:
        first_pair_trimmed = glob.glob('%s/*_val_1.fq'%(data_trimmed_dir))
        second_pair_trimmed = glob.glob('%s/*_val_2.fq'%(data_trimmed_dir))
    print('Trimmed Files:\n 1st:%s \n 2nd:%s' %(first_pair_trimmed,second_pair_trimmed))
    print
    first_pair_group = ' '.join(first_pair_trimmed)
    second_pair_group = ' '.join(second_pair_trimmed)
    pair_files = []
    for file in first_pair_trimmed:
        mate_file = file.replace('_1_val_1.fq','2_val_2.fq')
        paired_mates = file + ' ' + mate_file
        pair_files.append(paired_mates)

    star_input_files = ' '.join(pair_files)

    return first_pair_group,second_pair_group, star_input_files


####################### Run STAR #####################################
def run_star(first_pair_group, second_pair_group, results_dir, star_input_files,
             folder_name, genome_dir):
    print
    print('\033[33mRunning STAR! \033[0m')

    outfile_prefix = '%s/%s_star_' %(results_dir, folder_name)
    star_options ="--runThreadN 16 --outSAMattributes All --genomeLoad LoadAndRemove --outFilterType Normal  --outFilterMismatchNoverLmax 0.03 --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2"

    cmd = 'STAR --genomeDir %s %s --readFilesIn %s %s --outFileNamePrefix %s' % (genome_dir, star_options,
                                                                               first_pair_group, second_pair_group, outfile_prefix)
    print('STAR run command:%s' %cmd)
    os.system(cmd)

####################### Run Salmon Count ###############################
def run_salmon_quant(results_dir, folder_name, genome_fasta):
    print
    print('\033[33mRunning salmon-quant! \033[0m')
    salmon_input = '%s/%s_star_Aligned.out.bam' % (results_dir,folder_name)
    cmd = 'salmon quant -t %s -l A -a %s -o %s/salmon_quant' % (genome_fasta, salmon_input, results_dir)
    print('salmon-count run command:%s' %cmd)
    os.system(cmd)


####################### Run HTSEq Count ###############################
def run_htseq(htseq_dir, results_dir, folder_name, genome_gff):
    print
    print('\033[33mRunning htseq-count! \033[0m')
    htseq_input = '%s/%s_star_Aligned.sortedByCoord.out.bam' %(results_dir, folder_name)
    cmd = 'htseq-count -s "reverse" -t "exon" -i "Parent" -r pos --max-reads-in-buffer 60000000 -f bam %s %s > %s/%s_htseqcounts.txt' %(htseq_input,
                                                                                                                                        genome_gff,htseq_dir,folder_name)
    print('htseq-count run command:%s' %cmd)
    #os.system(cmd)

####################### Create STAR index ###############################
def create_genome_index(genome_dir, genome_fasta):
    index_cmd = 'STAR --runMode genomeGenerate --runThreadN 4 --genomeDir %s --genomeFastaFiles %s --genomeChrBinNbits 16' % (genome_dir,
                                                                                                                            genome_fasta)
    print(index_cmd)

    print ("\033[34m %s Indexing genome... \033[0m")
    if os.path.exists('%s/SAindex' % (genome_dir)):
        print ('Genome indexes exist. Not creating!')
    else:
        print ('Creating genome indexes')
        os.system(index_cmd)


####################### Running the Pipeline ###############################

def run_pipeline(data_folder, results_folder, genome_dir, genome_fasta, genome_gff):
    folder_count = 1

    # Loop through each data folder
    #for data_folder in data_folders:
    folder_name = data_folder.split('/')[-1]
    print
    print
    print('\033[33mProcessing Folder: %s\033[0m' % (folder_name))

    # Get the list of first file names in paired end sequences
    DATA_SEARCH1 = '%s/*_1.fq*' % data_folder
    print("SEARCHING FIRST PAIRS IN: ", DATA_SEARCH1)
    first_pair_files = glob.glob('%s/*_1.fq*' % (data_folder))
    #second_pair_files = glob.glob('%s/_R2*.fastq*' %(data_folder))

    # Program specific results directories
    data_trimmed_dir = "%s/%s/trimmed" % (results_folder,folder_name)
    fastqc_dir = "%s/%s/fastqc_results" % (results_folder,folder_name)

    results_dir = "%s/%s/results_STAR_Salmon" %(results_folder, folder_name)
    htseq_dir = "%s/htseqcounts" % (results_dir)

    # Run create directories function to create directory structure
    create_dirs(data_trimmed_dir, fastqc_dir, results_dir, htseq_dir)

    print("FIRST_PAIR_FILES: ", first_pair_files)

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
        #sys.exit()

        # Run TrimGalore
        trim_galore(first_pair_file,second_pair_file,folder_name,sample_id,file_ext,data_trimmed_dir,fastqc_dir)
        file_count = file_count + 1

        # Collect Trimmed data for input into STAR
        first_pair_group,second_pair_group,star_input_files = collect_trimmed_data(data_trimmed_dir,file_ext)

        # Run STAR
        run_star(first_pair_group,second_pair_group,results_dir,star_input_files, folder_name, genome_dir)

        # Run Salmon Quant
        run_salmon_quant(results_dir, folder_name, genome_fasta)

        # Run HTSeq count
        run_htseq(htseq_dir, results_dir, folder_name, genome_gff)

        folder_count += 1

    return data_trimmed_dir, fastqc_dir, results_dir


if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=DESCRIPTION)
    parser.add_argument('genomedir', help='genome directory')
    parser.add_argument('dataroot', help="parent of input directory")
    parser.add_argument('indir', help="input directory (R<somenumber>)")
    parser.add_argument('outdir', help='output directory')
    args = parser.parse_args()

    now = datetime.datetime.now()
    timeprint = now.strftime("%Y-%m-%d %H:%M")
    data_folder = "%s/%s" % (args.dataroot, args.indir)

    genome_fasta = "%s/Past_Smic_merged_CDS_suffixed.fasta" % args.genomedir
    genome_gff = "%s/past_smic.genome.annotation.gff3" % args.genomedir

    create_genome_index(args.genomedir, genome_fasta)
    data_trimmed_dir,fastqc_dir,results_dir = run_pipeline(data_folder, args.outdir, args.genomedir, genome_fasta, genome_gff)
