#!/usr/bin/env python3

"""
find_files_test.py - Unit tests for the globalsearch.rnaseq.find_files module
"""

import unittest
import xmlrunner
import os, sys
import fs
import globalsearch.rnaseq.find_files as find_files


class FindFilesTest(unittest.TestCase):

    def __make_input_folder(self, num_samples, pattern="R%d_%d.fq.gz"):
        for i in range(1, num_samples + 1):
            dir = "inputdata/R%d" % i
            self.mem_fs.makedir(dir)
            for j in range(2):
                fname = pattern % (i, j + 1)
                self.mem_fs.touch(fs.path.combine(dir, fname))

    def __make_include_file(self, path, includes):
        self.mem_fs.writetext(path, '\n'.join(includes))


    def setUp(self):
        self.mem_fs = fs.open_fs("mem://")
        self.mem_fs.makedir("inputdata")

    def tearDown(self):
        self.mem_fs.close()

    def test_rnaseq_data_folder_list_scandir(self):
        self.__make_input_folder(2)
        config = {'input_dir': '/inputdata'}
        result = find_files.rnaseq_data_folder_list(config, filesys=self.mem_fs)
        self.assertEqual(result, ['R1', 'R2'])

    def test_rnaseq_data_folder_list_includes(self):
        self.__make_input_folder(2)
        config = {'input_dir': '/inputdata', 'includes': ['R11', 'R12']}
        result = find_files.rnaseq_data_folder_list(config, filesys=self.mem_fs)
        self.assertEqual(result, ['R11', 'R12'])

    def test_rnaseq_data_folder_list_include_file(self):
        self.__make_input_folder(2)
        self.__make_include_file("/my_include.txt", ['R22', 'R23'])
        #config = {'input_dir': '/inputdata', 'includes': ['R11', 'R12']}
        config = {'input_dir': '/inputdata', 'include_file': '/my_include.txt'}
        result = find_files.rnaseq_data_folder_list(config, filesys=self.mem_fs)
        self.assertEqual(result, ['R22', 'R23'])

    def test_rnaseq_data_folder_list_includes_and_include_file(self):
        self.__make_input_folder(2)
        self.__make_include_file("/my_include.txt", ['R22', 'R23'])
        config = {
            'input_dir': '/inputdata',
            'includes': ['R11', 'R12'],
            'include_file': '/my_include.txt'
        }
        result = find_files.rnaseq_data_folder_list(config, filesys=self.mem_fs)
        self.assertEqual(sorted(result), ['R11', 'R12', 'R22', 'R23'])

    def test_find_fastq_files(self):
        self.__make_input_folder(2)
        result = find_files.find_fastq_files("/inputdata/R1", ["*_{{readnum}}.fq.*"], filesys=self.mem_fs)
        self.assertEqual(result, [('/inputdata/R1/R1_1.fq.gz', '/inputdata/R1/R1_2.fq.gz')])

    def test_find_fastq_files_multi_pattern(self):
        self.__make_input_folder(2)
        result = find_files.find_fastq_files("/inputdata/R1",
                                             ["*_{{readnum}}.fq.*", "*_{{readnum}}.fastq.*"],
                                             filesys=self.mem_fs)
        self.assertEqual(result, [('/inputdata/R1/R1_1.fq.gz', '/inputdata/R1/R1_2.fq.gz')])


if __name__ == '__main__':
    SUITE = []
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(FindFilesTest))
    if len(sys.argv) > 1 and sys.argv[1] == 'xml':
        xmlrunner.XMLTestRunner(output='test-reports').run(unittest.TestSuite(SUITE))
    else:
        unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(SUITE))
