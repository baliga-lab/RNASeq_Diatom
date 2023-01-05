#!/usr/bin/env python3

#############################################################
##### RNASeq Analysis Pipeline with STAR                #####
##### Last update: 11/15/2022 Serdar Turkarslan         #####
##### Institute for Systems Biology                     #####
############################################################
import glob, sys, os, string, datetime, re
import argparse
import subprocess

from .find_files import find_fastq_files
from .trim_galore import trim_galore, collect_trimmed_data, create_result_dirs

DESCRIPTION = """run_STAR_SALMON.py - run STAR and Salmon"""

####################### Create results directories ###############################


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
    os.system(cmd1)

    # Convert 2nd file to SAM by using picard FastqToSam tool
    cmd2 = 'picard FastqToSam F1=%s O=/proj/omics4tb2/Global_Search/Pilot_Fail_New/raw_data/%s/%s_unaligned_reads_2.sam SM=%s' %(second_pair_file, new_sample_id,new_sample_id,new_sample_id)
    print
    print( '++++++  FastqToSam Command for 2nd File:', cmd2)
    print
    os.system(cmd2)

# Merge Sam files
    cmd3 = 'picard MergeSamFiles I=/proj/omics4tb2/Global_Search/Pilot_Fail_New/raw_data/%s/%s_unaligned_reads_1.sam I=/proj/omics4tb2/Global_Search/Pilot_Fail_New/raw_data/%s/%s_unaligned_reads_2.sam O=/proj/omics4tb2/Global_Search/Pilot_Fail_New/raw_data/%s/%s_merged_reads.sam SORT_ORDER=queryname' %(new_sample_id,new_sample_id,new_sample_id,new_sample_id,new_sample_id,new_sample_id)
    print
    print( '++++++  MergeSamFiles Command:', cmd3)
    print
    os.system(cmd3)

    # UnMerge Sam files to fastq files
    cmd4 = 'picard SamToFastq I=/proj/omics4tb2/Global_Search/Pilot_Fail_New/raw_data/%s/%s_merged_reads.sam FASTQ=/proj/omics4tb2/Global_Search/Pilot_Fail_New/raw_data/%s/%s_reordered_1.fastq SECOND_END_FASTQ=/proj/omics4tb2/Global_Search/Pilot_Fail_New/raw_data/%s/%s_reordered_2.fastq UNPAIRED_FASTQ=/proj/omics4tb2/Global_Search/Pilot_Fail_New/raw_data/%s/%s_reordered_unpaired.fastq' %(new_sample_id,new_sample_id,new_sample_id,new_sample_id,new_sample_id,new_sample_id,new_sample_id,new_sample_id)
    print
    print( '++++++  UnMergeSamFiles Command:', cmd4)
    print
    os.system(cmd4)


####################### Trimgalore for quality and trimming ###############################




####################### Run STAR #####################################
### We need to add Read GRoup info
### --outSAMattrRGline ID:${i%_TF_R1_val_1.fq.gz}
### https://github.com/BarshisLab/danslabnotebook/blob/main/CBASSAS_GenotypeScreening.md

def run_star(first_pair_group, second_pair_group, results_dir, folder_name, genome_dir):
    print('\033[33mRunning STAR! \033[0m')
    outfile_prefix = '%s/%s_%s_' %(results_dir, folder_name, args.starPrefix)
    star_options = ["--runThreadN", "32",
                    "--outFilterType", "Normal",
                    "--outSAMstrandField", "intronMotif",
                    "--outFilterIntronMotifs", "RemoveNoncanonical",
                    "--outSAMtype", "BAM", "Unsorted",
                    "--limitBAMsortRAM", "5784458574",
                    "--readFilesCommand", "zcat",
                    "--outReadsUnmapped", "Fastx",
                    "--outFilterMismatchNmax", str(args.outFilterMismatchNmax),
                    "--outFilterMismatchNoverLmax", str(args.outFilterMismatchNoverLmax),
                    "--outFilterScoreMinOverLread", str(args.outFilterScoreMinOverLread),
                    "--outFilterMatchNmin", str(args.outFilterMatchNmin)]
    if args.twopassMode:
        star_options.extend(["--twopassMode", "Basic"])
        star_options.extend(["--genomeLoad", "NoSharedMemory"])
    else:
        star_options.extend(["--genomeLoad", "LoadAndKeep"])

    command = ["STAR", "--genomeDir", genome_dir]
    command += star_options
    command += [ "--readFilesIn", first_pair_group,
                 second_pair_group,
                 "--outFileNamePrefix", outfile_prefix]
    if args.outSAMattributes != "Standard":
        out_sam_attrs = args.outSAMattributes.split()
        command.append('--outSAMattributes')
        command += out_sam_attrs

    cmd = ' '.join(command)
    print('STAR run command:%s' % cmd)
    compl_proc = subprocess.run(command, check=True, capture_output=False)

####################### Deduplication (not in _old) ###############################
def dedup(results_dir,folder_name):
    print('\033[33mRunning Deduplication! \033[0m')
    outfile_prefix = '%s/%s_%s_' %(results_dir, folder_name, args.starPrefix)

    aligned_bam = '%sAligned.out.bam' % (outfile_prefix)
    fixmate_bam = '%sFixmate.out.bam' % (outfile_prefix)
    ordered_bam = '%sOrdered.out.bam' % (outfile_prefix)
    markdup_bam = '%sMarkedDup.out.bam' % (outfile_prefix)
    markdupSTAR_bam = '%sProcessed.out.bam' % (outfile_prefix)
    nosingleton_bam = '%sNoSingleton.out.bam' % (outfile_prefix)
    nosingletonCollated_bam = '%sNoSingletonCollated.out.bam' % (outfile_prefix)

    # STAR mark duplicates
    star_markdup_command = ['STAR', '--runThreadN', '32',
                            '--runMode',
                            'inputAlignmentsFromBAM',
                            '--bamRemoveDuplicatesType', 'UniqueIdenticalNotMulti',
                            '--inputBAMfile', aligned_bam,
                            '--outFileNamePrefix', outfile_prefix]
    star_markdup_cmd = ' '.join(star_markdup_command)

    # removesingletons from STAR
    rmsingletonsSTAR_command = ['samtools', 'view', '-@', '8',
                                '-b', '-F', '0x400', markdupSTAR_bam,
                                '>', nosingleton_bam]
    rmsingletonsSTAR_cmd = ' '.join(rmsingletonsSTAR_command)

    # Collate reads by name
    collatereadsSTAR_command = ['samtools', 'sort', '-o',
                                nosingletonCollated_bam,
                                '-n', '-@', '8', nosingleton_bam]
    collatereadsSTAR_cmd = ' '.join(collatereadsSTAR_command)

    ## STAR based BAM duplicate removal
    # Mark duplicates with STAR
    print('STAR mark duplicates run command:%s' % star_markdup_cmd)
    compl_proc = subprocess.run(star_markdup_command, check=True, capture_output=False)

    # Remove marked duplicates withh samtools
    print('Samtools  STAR Dedup Remove run command:%s' % rmsingletonsSTAR_cmd)
    os.system(rmsingletonsSTAR_cmd)  # WW: use subprocess.run(), test shell interaction

    # Remove marked duplicates withh samtools
    print('Samtools  Collate reads by read name run command:%s' % collatereadsSTAR_cmd)
    compl_proc = subprocess.run(collatereadsSTAR_command, check=True, capture_output=False)


####################### Run Salmon Count ###############################
# WW: Check the names of the input files they will be different from _out
def run_salmon_quant(results_dir, folder_name, genome_fasta):
    outfile_prefix = '%s/%s_%s_' %(results_dir, folder_name, args.starPrefix)
    print(outfile_prefix)
    print
    print('\033[33mRunning salmon-quant! \033[0m')
    # check if we are performing deduplication
    if args.dedup:
        salmon_input = '%sNoSingletonCollated.out.bam' % (outfile_prefix)
    else:
        salmon_input = '%sAligned.out.bam' % (outfile_prefix)

    cmd = 'salmon quant -t %s -l A -a %s -o %s/%s_salmon_quant' % (genome_fasta, salmon_input, results_dir, args.salmonPrefix)
    print('salmon-count run command: %s' % cmd)
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
    os.system(cmd)

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
        compl_proc = subprocess.run(index_command, check=True, capture_output=False)


####################### Running the Pipeline ###############################

def run_pipeline(data_folder, results_folder, genome_dir, genome_fasta, genome_gff, args):
    folder_count = 1

    # Loop through each data folder
    folder_name = data_folder.split('/')[-1]
    print
    print
    print('\033[33mProcessing Folder: %s\033[0m' % (folder_name))

    # Get the list of first file names in paired end sequences
    ## We need to make sure we capture fastq data files
    first_pair_files = find_fastq_files(data_folder, args.fastq_patterns.split(','))

    # Program specific results directories
    data_trimmed_dir = "%s/%s/trimmed" % (results_folder,folder_name)
    fastqc_dir = "%s/%s/fastqc_results" % (results_folder,folder_name)

    results_dir = "%s/%s/results_STAR_Salmon" %(results_folder, folder_name)
    htseq_dir = "%s/htseqcounts" % (results_dir)

    # Run create directories function to create directory structure
    create_result_dirs(data_trimmed_dir, fastqc_dir, results_dir, htseq_dir)

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

        #order_fq(first_pair_file, second_pair_file, data_folder, sample_id)

        # Run TrimGalore
        trim_galore(first_pair_file,second_pair_file,folder_name,sample_id,file_ext,data_trimmed_dir,fastqc_dir)
        file_count = file_count + 1

        # Collect Trimmed data for input into STAR
        first_pair_group, second_pair_group, pair_files = collect_trimmed_data(data_trimmed_dir, file_ext)

        # Run STAR
        run_star(first_pair_group, second_pair_group, results_dir, folder_name, genome_dir)

        # Run Deduplication
        if args.dedup:
            print
            print('\033[33mRunning Deduplication: \033[0m')
            dedup(results_dir,folder_name)

        # Run Salmon Quant
        run_salmon_quant(results_dir, folder_name, genome_fasta)

        # Run HTSeq count
        if not genome_gff is None and os.path.exists(genome_gff):
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
    parser.add_argument('--fastq_patterns', help="FASTQ file patterns", default="*_{{pairnum}}.fq.*")
    parser.add_argument('--genome_gff', help='genome GFF file')
    parser.add_argument('--genome_fasta', help='genome FASTA file')
    parser.add_argument('--dedup', action='store_true', help='should we deduplicate bam files (True or False)')
    parser.add_argument('--twopassMode', action='store_true', help='run STAR in two-pass mode')
    parser.add_argument('--starPrefix', help="STAR output file name prefix")
    parser.add_argument('--salmonPrefix', help="Salmon output folder name prefix")
    parser.add_argument('--outFilterMismatchNmax', nargs='?', const=10, type=int)
    parser.add_argument('--outFilterMismatchNoverLmax', nargs='?', const=0.3, type=float)
    parser.add_argument('--outFilterScoreMinOverLread', nargs='?', const=0.66, type=float)
    parser.add_argument('--outFilterMatchNmin', nargs='?', const=0, type=int)
    parser.add_argument('--outSAMattributes', nargs='?', type=str, default="Standard")

    #### Add argument for running star in two pass mode
    ### Kate to contribute relevant code
    args = parser.parse_args()

    now = datetime.datetime.now()
    timeprint = now.strftime("%Y-%m-%d %H:%M")
    data_folder = "%s/%s" % (args.dataroot, args.indir)
    if args.genome_fasta is not None and os.path.exists(args.genome_fasta):
        genome_fasta = args.genome_fasta
    else:
        genome_fasta = glob.glob('%s/*.fasta' % (args.genomedir))[0]

    create_genome_index(args.genomedir, genome_fasta)
    data_trimmed_dir,fastqc_dir,results_dir = run_pipeline(data_folder, args.outdir, args.genomedir, genome_fasta,
        args.genome_gff, args)
