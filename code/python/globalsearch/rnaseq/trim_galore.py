import os
import glob

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
    #cmd = 'trim_galore --fastqc_args "--outdir %s/" --paired --output_dir %s/ %s %s' %(fastqc_dir,data_trimmed_dir,first_pair_file, second_pair_file)
    command = ['trim_galore',
                '--fastqc_args "--outdir %s/"' % fastqc_dir,
                '--paired',
                '--output_dir', '%s/' % data_trimmed_dir,
                first_pair_file, second_pair_file]
    cmd = ' '.join(command)
    print
    print( '++++++ Trimgalore Command:', cmd)
    print
    # TODO: check by subprocess.run() does not work here !!
    #compl_proc = subprocess.run(command, check=True, capture_output=False)
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
    first_pair_group = ' '.join(first_pair_trimmed)
    second_pair_group = ' '.join(second_pair_trimmed)
    pair_files = []
    for file in first_pair_trimmed:
        mate_file = file.replace('_1_val_1.fq','2_val_2.fq')
        paired_mates = file + ' ' + mate_file
        pair_files.append(paired_mates)

    star_input_files = ' '.join(pair_files)

    return first_pair_group,second_pair_group, star_input_files
