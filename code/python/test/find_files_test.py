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

    def setUp(self):
        self.mem_fs = fs.open_fs("mem://")
        self.mem_fs.makedir("inputdata")

    def tearDown(self):
        self.mem_fs.close()

    def test_rnaseq_data_folder_list(self):
        self.assertEqual(self.mem_fs.listdir("/")[0], "inputdata")

    def test_find_fastq_files(self):
        pass

if __name__ == '__main__':
    SUITE = []
    SUITE.append(unittest.TestLoader().loadTestsFromTestCase(FindFilesTest))
    if len(sys.argv) > 1 and sys.argv[1] == 'xml':
        xmlrunner.XMLTestRunner(output='test-reports').run(unittest.TestSuite(SUITE))
    else:
        unittest.TextTestRunner(verbosity=2).run(unittest.TestSuite(SUITE))
