#!/usr/bin/env python3

#############################################################
##### RNASeq Analysis Pipeline - SplAdder               #####
##### Last update: 2023/03/01 by Yaqiao Li              #####
##### Institute for Systems Biology                     #####
#############################################################
import glob, sys, os, string, datetime, re, errno, time
import argparse
from multiprocessing import Pool, cpu_count, Manager
from subprocess import *
#from pathos.multiprocessing import ProcessingPoll as Pool
from functools import partial
#from find_files import find_fastq_files

DESCRIPTION = """run_spladder.py - run spladder"""

####################### Create results directories ###############################
def create_dirs(spladder_work_dir, spladder_out_dir, parsed_event_dir):
    dirs = [spladder_work_dir, spladder_out_dir, parsed_event_dir]
    for dir in dirs:
        if not os.path.exists('%s' %(dir)):
            os.makedirs('%s' %(dir))
        else:
            print('\033[31m %s directory exists. Not creating. \033[0m' %(dir))



####################### SplAdder Step1 Single graphs #############################
def SplAdder_step1_single_graphs(input_bam, genome_annotation, spladder_out_dir):
    command = ['spladder', 'build',
                '-o', spladder_out_dir,
                '-a', genome_annotation,
                '-b', input_bam,
                '--merge-strat', 'single',
                '--no-extract-ase',
                '--parallel', '2',
                '-v']
    cmd = ' '.join(command)
    print('SplAdder_Step1_single_graphs_run_command:\n%s' % cmd)
    outOut = open('stdout.out','a')
    errOut = open('ErrorReport.txt','a')
    errOut.write('\nSplAdder_Step1_single_graphs_run_command:\n%s\n' % cmd)
    spladderStep1Proc = Popen(cmd, shell=True, stdout=outOut, stderr=errOut)
    stdoutput = spladderStep1Proc.communicate()[0]
    outOut.close()
    errOut.close()
    
####################### SplAdder Step2 Merged graphs #############################
def SplAdder_step2_merged_graphs(genome_annotation, spladder_out_dir, input_bam_list):
    command = ['spladder', 'build',
                '-o', spladder_out_dir,
                '-a', genome_annotation,
                '-b', input_bam_list,
                '--merge-strat', 'merge_graphs',
                '--no-extract-ase',
                '--parallel', '40',
                '-v']
    cmd = ' '.join(command)
    print('SplAdder_Step2_merging_splice_graphs_run_command:\n%s' % cmd)
    outOut = open('stdout.out','a')
    errOut = open('ErrorReport.txt','a')
    errOut.write('\nSplAdder_Step2_merging_splice_graphs_run_command:\n%s\n' % cmd)
    spladderStep2Proc = Popen(cmd, shell=True, stdout=outOut, stderr=errOut)
    stdoutput = spladderStep2Proc.communicate()[0]
    outOut.close()
    errOut.close()
    #os.system(cmd)
    #compl_proc = subprocess.run(command, check=True, capture_output=False)

####################### SplAdder Step3 Quantification #############################
def SplAdder_step3_quantification(input_bam, genome_annotation, spladder_out_dir):
    command = ['spladder', 'build',
                '-o', spladder_out_dir,
                '-a', genome_annotation,
                '-b', input_bam,
                '--merge-strat', 'merge_graphs',
                '--no-extract-ase',
                '--quantify-graph',
                '--qmode', 'single',
                '--parallel', '2',
                '-v']
    cmd = ' '.join(command)
    print('SplAdder_Step3_splice_graphs_run_command:\n%s' % cmd)
    outOut = open('stdout.out','a')
    errOut = open('ErrorReport.txt','a')
    errOut.write('\nSplAdder_Step3_splice_graphs_run_command:\n%s\n' % cmd)
    spladderStep3Proc = Popen(cmd, shell=True, stdout=outOut, stderr=errOut)
    stdoutput = spladderStep3Proc.communicate()[0]
    outOut.close()
    errOut.close()

####################### SplAdder Step4 Aggregate Quantification####################
def SplAdder_step4_aggregate_quantification(genome_annotation, spladder_out_dir, input_bam_list):
    command = ['spladder', 'build',
                '-o', spladder_out_dir,
                '-a', genome_annotation,
                '-b', input_bam_list,
                '--merge-strat', 'merge_graphs',
                '--no-extract-ase',
                '--quantify-graph',
                '--qmode', 'collect',
                '--parallel', '40',
                '-v']
    cmd = ' '.join(command)
    print('SplAdder_Step4_aggregate_quantification_run_command:\n%s' % cmd)
    outOut = open('stdout.out','a')
    errOut = open('ErrorReport.txt','a')
    errOut.write('\nSplAdder_Step4_aggregate_quantification_run_command:\n%s\n' % cmd)
    spladderStep4Proc = Popen(cmd, shell=True, stdout=outOut, stderr=errOut)
    stdoutput = spladderStep4Proc.communicate()[0]
    outOut.close()
    errOut.close()
    #os.system(cmd)

####################### SplAdder Step5 Call Events ###############################
def SplAdder_step5_call_events(genome_annotation, spladder_out_dir, input_bam_list):
    command = ['spladder', 'build',
                '-o', spladder_out_dir,
                '-a', genome_annotation,
                '-b', input_bam_list,
                '--parallel', '40',
                '-v']
    cmd = ' '.join(command)
    print('SplAdder_Step5_call_events_run_command:\n%s' % cmd)
    outOut = open('stdout.out','a')
    errOut = open('ErrorReport.txt','a')
    errOut.write('\nSplAdder_Step5_call_events_run_command:\n%s\n' % cmd)
    spladderStep5Proc = Popen(cmd, shell=True, stdout=outOut, stderr=errOut)
    stdoutput = spladderStep5Proc.communicate()[0]
    outOut.close()
    errOut.close()
    #os.system(cmd)

####################### SplAdder Step6 Contrast Tests ############################
def SplAdder_step6_contrast_test(contrs, spladder_out_dir, contrast_dir):
    [var1, var2, const] = contrs.split()
    conditionA = contrast_dir + '/host_' + var1 + '_' + const + '.txt'
    conditionB = contrast_dir + '/host_' + var2 + '_' + const + '.txt'
    labelA = var1 + const
    labelB = var2 + const

    command = ['spladder', 'test',
                '-o', spladder_out_dir,
                '--out-tag', 'Rem_Prob',
                '--conditionA', conditionA,
                '--conditionB', conditionB,
                '--labelA', labelA,
                '--labelB', labelB,
                '--diagnose-plots',
                '--parallel', '5',
                '-v']
    cmd = ' '.join(command)
    print('SplAdder_Step6_contrast_test_run_command:\n%s' % cmd)
    outOut = open('stdout.out','a')
    errOut = open('ErrorReport.txt','a')
    errOut.write('\nSplAdder_Step6_contrast_test_run_command:\n%s\n' % cmd)
    spladderStep6Proc = Popen(cmd, shell=True, stdout=outOut, stderr=errOut)
    stdoutput = spladderStep6Proc.communicate()[0]
    outOut.close()
    errOut.close()

####################### SplAdder Step7 Parsing Statistic #########################
def SplAdder_step7_parsing_statistic(contrs, parseRscript):
    [var1, var2, const] = contrs.split()

    command = ['Rscript',
                parseRscript,
                var1,
                var2,
                const]
    cmd = ' '.join(command)
    print('SplAdder_Step7_parsing_statistic_run_command:\n%s' % cmd)
    outOut = open('stdout.out','a')
    errOut = open('ErrorReport.txt','a')
    errOut.write('\nSplAdder_Step7_parsing_statistic_run_command:\n%s\n' % cmd)
    spladderStep7Proc = Popen(cmd, shell=True, stdout=outOut, stderr=errOut)
    stdoutput = spladderStep7Proc.communicate()[0]
    outOut.close()
    errOut.close()


####################### Check Error Report  #####################################
def Check_Error_Report():
    fileopen = open('ErrorReport.txt', 'r')
    fileread = fileopen.read().replace('FileExistsError', '')
    fileopen.close()
    errOut = open('ErrorReport.txt','a')
    errorKeyWord = ['Execution halted', 'Error', 'error']
    cleanUp = True
    cf = open('cleanUp.sh','w')
    cf.write('#Remove temporary files\n\n')
    for e in errorKeyWord:
        if e in fileread:
            errOut.write('\nFind [ ' + e + ' ] error in report, please check it. For clean up manually please run cleanUp.sh\n')
            cleanUp = False
            for root, dirs, files in os.walk(spladder_out_dir):
                for f in files:
                    if f.endswith('pickle'):
                        cf.write('rm ' + os.path.join(root, f) + '\n')
    errOut.close()
    cf.close()
    return cleanUp

####################### Clean Up temporary files  ###############################
def Clean_Up(spladder_out_dir):
    for root, dirs, files in os.walk(spladder_out_dir):
        for f in files:
            if f.endswith('pickle'):
                os.remove(os.path.join(root, f))

####################### Run SplAdder ############################################
def run_spladder(spladder_work_dir, spladder_out_dir, parsed_event_dir, input_bam_list, genome_annotation, all_contrasts, contrast_dir, sampleType, args):
    # Run create directories function to create directory structure
    create_dirs(spladder_work_dir, spladder_out_dir, parsed_event_dir)
    
    #---run step 1
    bamlist = open(input_bam_list, 'r')
    bams = []
    for line in bamlist:
        single_bam = line.strip()
        bams.append(single_bam)
    print(bams)
    bamlist.close()
   
    partial_work_step1 = partial(SplAdder_step1_single_graphs, genome_annotation=genome_annotation, spladder_out_dir=spladder_out_dir)
    cpus = 16
    pool = Pool(processes=cpus)
    pool.map(partial_work_step1, bams)
    pool.close()
    pool.join()

    #---run step 2
    SplAdder_step2_merged_graphs(genome_annotation, spladder_out_dir, input_bam_list)

    #---run step 3
    partial_work_step3 = partial(SplAdder_step3_quantification, genome_annotation=genome_annotation, spladder_out_dir=spladder_out_dir)
    cpus = 16
    pool = Pool(processes=cpus)
    pool.map(partial_work_step3, bams)
    pool.close()
    pool.join()

    #---run step 4
    SplAdder_step4_aggregate_quantification(genome_annotation, spladder_out_dir, input_bam_list)

    #---run step 5
    SplAdder_step5_call_events(genome_annotation, spladder_out_dir, input_bam_list)

    #---run step 6
    contrast_list = open(all_contrasts, 'r')
    contrs = []
    for line in contrast_list:
        contrs.append(line.strip())
    contrast_list.close()

    partial_work_step6 = partial(SplAdder_step6_contrast_test, spladder_out_dir=spladder_out_dir, contrast_dir=contrast_dir)
    cpus = 16
    pool = Pool(processes=cpus)
    pool.map(partial_work_step6, contrs)
    pool.close()
    pool.join()



    #---run step 7
    parseRscript = spladder_out_dir + '/run_parse_spladder_output.R'

    strR1 = r'''args<-commandArgs(TRUE)
var1<-args[1]
var2 <-args[2]
const <-args[3]

# read in the event files for each contrast
setwd(paste0("'''

    strR2 = spladder_out_dir

    strR3 = r'''/testing_",var1,const,"_vs_",var2,const,"_Rem_Prob"))

ldf <- list()
listcsv <- dir(pattern = "*gene_unique.tsv")
for (k in 1:length(listcsv)){
    ldf[[k]] <- read.table(file = listcsv[k], header = T)
}

# combine into one DF
one_ldf <- do.call(rbind, ldf)

#add a column with the constant
one_ldf$const <- paste0(const)

#add a column with the contrast variables
one_ldf$vrs <- paste0(var1,"vs",var2)

# add a column for species name
one_ldf$spec <- c("'''

    strR3_1 = r'''")

#write out the results as one table
write.table(one_ldf, paste0("'''

    if(sampleType == 'host'):
        strR4 = parsed_event_dir + '/host_all_events'
    elif(sampleType == 'sym'):
        strR4 = parsed_event_dir + '/sym_all_events'
    else:
        strR4 = parsed_event_dir + '/other_all_events'
    strR5 = r'",var1,const,"_vs_",var2,const,".tsv"), quote = F, row.names = F)'
    strR = strR1 + strR2 + strR3 + sampleType + strR3_1 + strR4 + strR5 + '\n'
    parseR = open(parseRscript, 'w')
    parseR.write(strR)
    parseR.close()

    partial_work_step7 = partial(SplAdder_step7_parsing_statistic, parseRscript=parseRscript)
    cpus = 16
    pool = Pool(processes=cpus)
    pool.map(partial_work_step7, contrs)
    pool.close()
    pool.join()    

    #---run spladder program finished
    




####################### Main ####################################################
if __name__ == '__main__':
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=DESCRIPTION)
    parser.add_argument('-w', '--WorkDir', help='work directory')
    parser.add_argument('-i', '--inputBamFileList', help='input .bam file name list')
    parser.add_argument('-a', '--genomeAnnotation', help='genome GFF/GTF file')
    parser.add_argument('-c', '--allContrasts', help='sample contrasts information: all_contrasts.txt')
    parser.add_argument('-d', '--contrastDir', help='sample contrasts information: contrast_file')
    parser.add_argument('-s', '--sampleType', help='sample is host or symbiont or Other')
    
    args = parser.parse_args()
    
    work_dir = args.WorkDir
    input_bam_list = args.inputBamFileList
    genome_annotation = args.genomeAnnotation
    all_contrasts = args.allContrasts
    contrast_dir = args.contrastDir
    sampleType = args.sampleType
    
    spladder_work_dir = work_dir + '/' + sampleType + '_spladder_jobs'
    spladder_out_dir = work_dir + '/' + sampleType + '_spladder_jobs' + '/array_spladder_out'
    if sampleType == 'host' or sampleType == 'sym':
        parsed_event_dir = work_dir + '/host_sym_parsed_event_files'
    else:
        parsed_event_dir = work_dir + '/' + sampleType + '_parsed_event_files'

    errOut = open('ErrorReport.txt','w')
    starttime = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    errOut.write('Start: ' + starttime + '\n\n')
    errOut.close()
    run_spladder(spladder_work_dir, spladder_out_dir, parsed_event_dir, input_bam_list, genome_annotation, all_contrasts, contrast_dir, sampleType, args)
    errOut = open('ErrorReport.txt','a')
    endtime = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
    errOut.write('\nEnd: ' + endtime + '\n\n')
    errOut.close()
    
    # check errors
    cleanUp = Check_Error_Report()
    
    # clean up
    errOut = open('ErrorReport.txt','a')
    if cleanUp:
        Clean_Up(spladder_out_dir)
        cleantime = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        errOut.write('\nCleanedUp: ' + cleantime + '\n\n')
    else:
        cleantime = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
        errOut.write('\nCan not cleanedUp: ' + cleantime + '\n\n')     
    errOut.close()




