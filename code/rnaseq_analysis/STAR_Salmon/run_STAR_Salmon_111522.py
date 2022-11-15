#!/usr/bin/env python3

#############################################################
##### RNASeq Analysis Pipeline with STAR                #####
##### Last update: 11/15/2022 Serdar Turkarslan         #####
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


####################### Reorder Fastq files By read Name ###############################
### Remove this function
def order_fq(first_pair_file, second_pair_file, data_folder, sample_id):
    #print("1stpair:%s, 2ndpair:%s, folder_name:%s, sample_name:%s")%(first_pair_file,second_pair_file,folder_name,sample_name)
    print
    print ("\033[34m Running Order Fq \033[0m")
    new_sample_id = sample_id.split("_")[0]
    
    # Convert 1st file to SAM by using picard FastqToSam tool
    cmd1 = 'picard FastqToSam F1=%s O=/proj/omics4tb2/Global_Search/Pilot_Fail_New/raw_data/%s/%s_unaligned_reads_1.sam SM=%s' %(first_pair_file, new_sample_id,new_sample_id,new_sample_id)
    print
    print( '++++++  FastqToSam Command for 1st File:', cmd1)
    print
    #os.system(cmd1)

    # Convert 2nd file to SAM by using picard FastqToSam tool
    cmd2 = 'picard FastqToSam F1=%s O=/proj/omics4tb2/Global_Search/Pilot_Fail_New/raw_data/%s/%s_unaligned_reads_2.sam SM=%s' %(second_pair_file, new_sample_id,new_sample_id,new_sample_id)
    print
    print( '++++++  FastqToSam Command for 2nd File:', cmd2)
    print
    #os.system(cmd2)

# Merge Sam files
    cmd3 = 'picard MergeSamFiles I=/proj/omics4tb2/Global_Search/Pilot_Fail_New/raw_data/%s/%s_unaligned_reads_1.sam I=/proj/omics4tb2/Global_Search/Pilot_Fail_New/raw_data/%s/%s_unaligned_reads_2.sam O=/proj/omics4tb2/Global_Search/Pilot_Fail_New/raw_data/%s/%s_merged_reads.sam SORT_ORDER=queryname' %(new_sample_id,new_sample_id,new_sample_id,new_sample_id,new_sample_id,new_sample_id)
    print
    print( '++++++  MergeSamFiles Command:', cmd3)
    print
    #os.system(cmd3)

    # UnMerge Sam files to fastq files
    cmd4 = 'picard SamToFastq I=/proj/omics4tb2/Global_Search/Pilot_Fail_New/raw_data/%s/%s_merged_reads.sam FASTQ=/proj/omics4tb2/Global_Search/Pilot_Fail_New/raw_data/%s/%s_reordered_1.fastq SECOND_END_FASTQ=/proj/omics4tb2/Global_Search/Pilot_Fail_New/raw_data/%s/%s_reordered_2.fastq UNPAIRED_FASTQ=/proj/omics4tb2/Global_Search/Pilot_Fail_New/raw_data/%s/%s_reordered_unpaired.fastq' %(new_sample_id,new_sample_id,new_sample_id,new_sample_id,new_sample_id,new_sample_id,new_sample_id,new_sample_id)
    print
    print( '++++++  UnMergeSamFiles Command:', cmd4)
    print
    #os.system(cmd4)


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
    #os.system(cmd)


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
### We need to add Read GRoup info
### --outSAMattrRGline ID:${i%_TF_R1_val_1.fq.gz}
### https://github.com/BarshisLab/danslabnotebook/blob/main/CBASSAS_GenotypeScreening.md

def run_star(first_pair_group, second_pair_group, results_dir, star_input_files,
             folder_name, genome_dir):
    print
    print('\033[33mRunning STAR! \033[0m')

    outfile_prefix = '%s/%s_%s_' %(results_dir, folder_name, args.starPrefix)
    
    star_options ="--runThreadN 32 --outSAMattributes All --genomeLoad LoadAndKeep --outFilterType Normal  --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM Unsorted --limitBAMsortRAM 5784458574 --readFilesCommand zcat --outReadsUnmapped Fastx --outFilterMismatchNmax %s --outFilterMismatchNoverLmax %s --outFilterScoreMinOverLread %s --outFilterMatchNmin %s" % (args.outFilterMismatchNmax, args.outFilterMismatchNoverLmax,args.outFilterScoreMinOverLread, args.outFilterMatchNmin)

    cmd = 'STAR --genomeDir %s %s --readFilesIn %s %s --outFileNamePrefix %s' % (genome_dir, star_options,first_pair_group, second_pair_group, outfile_prefix)

    print('STAR run command:%s' %cmd)
    #os.system(cmd)

####################### Deduplication ###############################
def dedup(results_dir,folder_name):
    print
    print('\033[33mRunning Deduplication! \033[0m')
    outfile_prefix = '%s/%s_%s_' %(results_dir, folder_name, args.starPrefix)

    aligned_bam = '%sAligned.out.bam' % (outfile_prefix)
    fixmate_bam = '%sFixmate.out.bam' % (outfile_prefix)
    ordered_bam = '%sOrdered.out.bam' % (outfile_prefix)
    markdup_bam = '%sMarkedDup.out.bam' % (outfile_prefix)
    markdupSTAR_bam = '%sProcessed.out.bam' % (outfile_prefix)
    nosingleton_bam = '%sNoSingleton.out.bam' % (outfile_prefix)
    nosingletonCollated_bam = '%sNoSingletonCollated.out.bam' % (outfile_prefix)
    
    
    # Add ms and MC tags for markdup to use later:
    fixmate_cmd = 'samtools fixmate -m %s %s' %(aligned_bam,fixmate_bam)
    # position order
    sort_cmd = 'samtools sort -o %s %s' %(ordered_bam,fixmate_bam)
    # mark duplicates
    markdup_cmd = 'samtools markdup -r -s %s %s' %(ordered_bam,markdup_bam)
    # removesingletons
    rmsingletons_cmd = 'samtools view -@ 8 -F 0x08 -b %s > %s' %(markdup_bam,nosingleton_bam)
    # STAR mark duplicates
    star_markdup_cmd = 'STAR --runThreadN 32 --runMode inputAlignmentsFromBAM --bamRemoveDuplicatesType UniqueIdenticalNotMulti --inputBAMfile %s --outFileNamePrefix %s' % (aligned_bam,outfile_prefix)
    # removesingletons fron STAR
    rmsingletonsSTAR_cmd = 'samtools view -@ 8 -b -F 0x400 %s > %s' %(markdupSTAR_bam,nosingleton_bam)
    # Collate reads by name
    collatereadsSTAR_cmd = 'samtools sort -o %s -n -@ 8 %s' %(nosingletonCollated_bam, nosingleton_bam)

    ## Samtools based BAM duplicate removal
    print()
    print('samtools fixmate run command:%s' %fixmate_cmd)
    #os.system(fixmate_cmd)

    print()
    print('samtools sort run command:%s' %sort_cmd)
    #os.system(sort_cmd)

    print()
    print('samtools mark diuplicates run command:%s' %markdup_cmd)
    #os.system(markdup_cmd)

    print()
    print('samtools rm singletons run command:%s' %rmsingletons_cmd)
    #os.system(rmsingletons_cmd)

    ## STAR based BAM duplicate removal
    # Mark duplicates with STAR
    print()
    print('STAR mark duplicates run command:%s' %star_markdup_cmd)
    os.system(star_markdup_cmd)

    # Remove marked duplicates withh samtools
    print()
    print('Samtools  STAR Dedup Remove run command:%s' %rmsingletonsSTAR_cmd)
    os.system(rmsingletonsSTAR_cmd)
    
    # Remove marked duplicates withh samtools
    print()
    print('Samtools  Collate reads by read name run command:%s' %collatereadsSTAR_cmd)
    os.system(collatereadsSTAR_cmd)


####################### Run Salmon Count ###############################
def run_salmon_quant(results_dir, folder_name, genome_fasta):
    outfile_prefix = '%s/%s_%s_' %(results_dir, folder_name, args.starPrefix)
    print(outfile_prefix)
 
    print
    print('\033[33mRunning salmon-quant! \033[0m')
    # check if we are performing deduplication
    if args.dedup == True:
        salmon_input = '%sNoSingletonCollated.out.bam' % (outfile_prefix)
    else:
        salmon_input = '%sAligned.out.bam' % (outfile_prefix)

    cmd = 'salmon quant -t %s -l A -a %s -o %s/%s_salmon_quant' % (genome_fasta, salmon_input, results_dir,args.salmonPrefix)
    print('salmon-count run command:%s' %cmd)
    os.system(cmd)


####################### Run HTSEq Count ###############################
#### We can remove this since we are using salmon quant
def run_htseq(htseq_dir, results_dir, folder_name, genome_gff):
    print
    print('\033[33mRunning htseq-count! \033[0m')
    htseq_input = '%s/%s_star_Aligned.sortedByCoord.out.bam' %(results_dir, folder_name)
    cmd = 'htseq-count -s "reverse" -t "exon" -i "Parent" -r pos --max-reads-in-buffer 60000000 -f bam %s %s > %s/%s_htseqcounts.txt' %(htseq_input,
                                                                                                                                        genome_gff,htseq_dir,folder_name)
    print('htseq-count run command:%s' %cmd)
    #os.system(cmd)

####################### Create STAR index ###############################
### This should be specific for the organism
### Use the equation file maybe another script to create references
def create_genome_index(genome_dir, genome_fasta):
    index_cmd = 'STAR --runMode genomeGenerate --runThreadN 32 --genomeDir %s --genomeFastaFiles %s --genomeChrBinNbits 16 --genomeSAindexNbases 12'% (genome_dir,genome_fasta)
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
    ## We need to make sure we capture fastq data files
    
    DATA_SEARCH1 = '%s/*_concat_fixed_1.fq*' % data_folder # test
    #DATA_SEARCH1 = '%s/*_1.fq*' % data_folder
    print("SEARCHING FIRST PAIRS IN: ", DATA_SEARCH1)
    first_pair_files = glob.glob('%s/*_concat_fixed_1.fq*' % (data_folder)) # testing concat
    #first_pair_files = glob.glob('%s/*_1.fq*' % (data_folder))
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

        order_fq(first_pair_file, second_pair_file, data_folder, sample_id)
        #sys.exit()

        # Run TrimGalore
        trim_galore(first_pair_file,second_pair_file,folder_name,sample_id,file_ext,data_trimmed_dir,fastqc_dir)
        file_count = file_count + 1
        #sys.exit()

        # Collect Trimmed data for input into STAR
        first_pair_group,second_pair_group,star_input_files = collect_trimmed_data(data_trimmed_dir,file_ext)

        # Run STAR
        run_star(first_pair_group,second_pair_group,results_dir,star_input_files, folder_name, genome_dir)
        #sys.exit()
        # Run Deduplication
        if args.dedup == True:
            print
            print('\033[33mRunning Deduplication: \033[0m')
            dedup(results_dir,folder_name)

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
    parser.add_argument('--dedup', action='store_true', help='should we deduplicate bam files (True or False)')
    parser.add_argument('--starPrefix', help="STAR output file name prefix")
    parser.add_argument('--salmonPrefix', help="Salmon output folder name prefix")
    parser.add_argument('--outFilterMismatchNmax', nargs='?', const=10, type=int)
    parser.add_argument('--outFilterMismatchNoverLmax', nargs='?', const=0.3, type=float)
    parser.add_argument('--outFilterScoreMinOverLread', nargs='?', const=0.66, type=float)
    parser.add_argument('--outFilterMatchNmin', nargs='?', const=0, type=int)

    #### Add argument for running star in two pass mode
    ### Kate to contribute relevant code
    

    args = parser.parse_args()

    now = datetime.datetime.now()
    timeprint = now.strftime("%Y-%m-%d %H:%M")
    data_folder = "%s/%s" % (args.dataroot, args.indir)

    genome_fasta = glob.glob('%s/*.fasta' % (args.genomedir))[0]
    genome_gff = "%s/past_smic.genome.annotation.gff3" % args.genomedir

    create_genome_index(args.genomedir, genome_fasta)
    data_trimmed_dir,fastqc_dir,results_dir = run_pipeline(data_folder, args.outdir, args.genomedir, genome_fasta, genome_gff)
