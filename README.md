# Global Search Code Repository

![Python Tests](https://github.com/baliga-lab/Global_Search/actions/workflows/python-package.yml/badge.svg)

## Description

This project manages all the software components related to the
Coral Reef Global Search project.

## Current pipeline tools

```gs_submit <config-file>```

Starting point point for the pipeline

## System requirements

  * Slurm cluster
  * STAR (https://github.com/alexdobin/STAR)
  * salmon 0.13.1 (https://combine-lab.github.io/salmon/)
  * [Kallisto (https://pachterlab.github.io/kallisto/)]
  * TrimGalore (https://github.com/FelixKrueger/TrimGalore)
  * samtools (http://www.htslib.org/)
  * htseq (https://github.com/htseq/htseq)
  * [Picard tools (https://broadinstitute.github.io/picard/)]
  * Python 3 (==3.10), or Anaconda 3

## Installation

Python:

```
$ pip install globalsearch
```

Within R:

```
$ library('devtools')
$ devtools::install_github('https://github.com/baliga-lab/Global_Search.git', ref="main", subdir="code/rpackage")
```

## Configuration format description

Configuration files are in JSON format of the following form

```
{
  "organisms": [<organism 1>, ...],
  "input_dir": <input directory>,
  "genome_dir": <genome file directory>,
  "output_dir": <output directory>,
  "postrun_output_dir": "<post-run output directory>",
  "log_dir": <log directory>,
  "genome_gff": <GFF file>,
  "genome_fasta": <FASTA file path>,
  "fastq_patterns": ["*_{{readnum}}.fq.*", "*_{{readnum}}.fastq.*"],
  "includes": [<directory name],
  "include_file": <path to file containing included directories>,
  "deduplicate_bam_files": false,
  "rnaseq_algorithm": "star_salmon",
  "star_options": {
     "outFilterMismatchNmax": 10,
     "outFilterMismatchNoverLmax": 0.3,
     "outFilterScoreMinOverLread": 0.66,
     "outFilterMatchNmin": 0,
     "twopassMode": false
  },
  "star_index_options": {
     "runThreadN": 32,
     "genomeChrBinNbits": 16,
     "genomeSAindexNbases": 12
  },
  "sbatch_options": {
     <pipeline_step_name>: {
       "options": [
          <sbatch options>
       ],
       "extras": [
          <additional lines for slurm job script>
       ]
     }
   }
}
```
