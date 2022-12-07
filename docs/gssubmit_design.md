# gs_submit Design

## Description

gs_submit is the top-level command for submitting compute jobs for Global Search.


## Flow

### gs_prepare

gs_prepare is the first stage. It analyzes the configuration parameters for integrity
and ensures that everything is installed

### make_star_salmon_job

This tool takes the configuration parameters and creates a temporary job file

### run_star_salmon

The control script is the main control script to run the compute job

