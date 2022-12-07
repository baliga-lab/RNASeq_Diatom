# Global Search Code Repository

## Description

This project manages all the software components related to the
Coral Reef Global Search project.

## Current pipeine tools

```gs_submit.sh <config-file>```

Starting point point for the pipeline

## System requirements

  * Slurm cluster
  * STAR (https://github.com/alexdobin/STAR)
  * salmon (https://combine-lab.github.io/salmon/)
  * samtools (http://www.htslib.org/)
  * htseq (https://github.com/htseq/htseq)
  * [Picard tools (https://broadinstitute.github.io/picard/)]
  * Python 3 (> 3.10)

## Installation

The commands in this repository can be used right away. Assume the 
repository is in `<yourpath>/Global_Search`. Add this to the .bashrc

```
export GS_HOME=<yourpath/Global_Search>
export PATH=$PATH:Global_Search/code/bin
```

## Configuration format description

Configuration files are in JSON format of the following form

```
{
  "input_dir": <input directory>,
  "genome_dir": <input directory>,
  "genome_gff": <GFF file>,
  "genome_fasta": <FASTA file>,
  "includes": [<directory name],
  "deduplicate_bam_files": false,
  "star_options": {
     "outFilterMismatchNmax": 10,
     "outFilterMismatchNoverLmax": 0.3,
     "outFilterScoreMinOverLread": 0.66,
     "outFilterMatchNmin": 0,
     "twopassMode": false
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
