import os
from fs.osfs import OSFS
import fs

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

GZ_PATTERN = '/*_val_%d.fq.gz'
FQ_PATTERN = '/*_val_%d.fq'

def collect_trimmed_data(data_trimmed_dir, file_ext, rootfs=OSFS("/")):
    filesys = rootfs.opendir(data_trimmed_dir)
    # define result files
    if file_ext == "gz":
        first_pair_trimmed, second_pair_trimmed = [filesys.glob(GZ_PATTERN % i) for i in [1, 2]]
    else:
        first_pair_trimmed, second_pair_trimmed = [filesys.glob(FQ_PATTERN % i) for i in [1, 2]]

    first_pair_trimmed = [fs.path.combine(data_trimmed_dir, match.path)
                          for match in first_pair_trimmed]
    second_pair_trimmed = [fs.path.combine(data_trimmed_dir, match.path)
                           for match in second_pair_trimmed]

    first_pair_group = ' '.join(first_pair_trimmed)
    second_pair_group = ' '.join(second_pair_trimmed)
    pair_files = []

    for file in first_pair_trimmed:
        mate_file = file.replace('_1_val_1.fq','_2_val_2.fq')
        #paired_mates = "%s %s" % (file, mate_file)
        pair_files.append((file, mate_file))

    return first_pair_group, second_pair_group, pair_files

def create_result_dirs(data_trimmed_dir, fastqc_dir, results_dir, htseq_dir):
    dirs = [data_trimmed_dir, fastqc_dir, results_dir, htseq_dir]
    for dir in dirs:
        # create results folder
        #print(dir)
        if not os.path.exists('%s' %(dir)):
            os.makedirs('%s' %(dir))
        else:
            print('\033[31m %s directory exists. Not creating. \033[0m' %(dir))
