"""
Analysis pipeline for paired end reads

Author: Ida Kups (iku)
"""
import os
import pathlib
import pandas as pd
import metadata_df


configfile: 'config_evaxion.yaml'

localrules: all, ok

PROJECTS = ['DNMT3A']
wgbs_df = metadata_df.DNMT_wgbs

prj2ref = {
     'DNMT3A':'mouse',
     'TET2':'mouse',
     'NZM':'human'
}

rule all:
    input: 'results/OK'

rule genome2cytosine_raw:
    #input: 'results/bismark/methylation_call_raw/{project}/{patient_id}_{technique}_{sra_id}/{sra_id}_bismark_bt2.bismark.cov.gz'
    input: 'results/bismark/methylation_call_raw/{project}/{patient_id}_{technique}_{sra_id}/{sra_id}.bismark.cov.gz'
    output: 'results/bismark/methylation_call_raw/{project}/{patient_id}_{technique}_{sra_id}/cytosine_report.txt.CpG_report.txt'
    threads: 8
    params:
        index = lambda wc: config[prj2ref[wc.project]]['bismark'],
    singularity: config['singularity']
    shell:
        """
        coverage2cytosine --dir $(dirname {output}) --genome_folder {params.index} -o cytosine_report.txt {input}
        """ 


use rule genome2cytosine_raw as genome2cytosine_trimmed with:
    #input: 'results/bismark/methylation_call_trimmed/{project}/{patient_id}_{technique}_{sra_id}/{sra_id}_trimmed_bismark_bt2.bismark.cov.gz'
    input: 'results/bismark/methylation_call_trimmed/{project}/{patient_id}_{technique}_{sra_id}/{sra_id}.bismark.cov.gz'
    output: 'results/bismark/methylation_call_trimmed/{project}/{patient_id}_{technique}_{sra_id}/cytosine_report.txt.CpG_report.txt'
    threads: 8


rule ok: # You only have to expand on the last rule in the branch and all depedencies will also be run, i.e. only expand on the methylatio_call, rsem and fastqc 
    """Collects all output from above and writes the OK file for `rule all`."""
    input:
        c = expand(rules.genome2cytosine_raw.output, zip, project=PROJECTS * len(wgbs_df), technique=wgbs_df['Assay Type'].tolist(), patient_id=wgbs_df['patient_id'].tolist(), sra_id=wgbs_df['Run'].tolist()),
        d = expand(rules.genome2cytosine_trimmed.output, zip, project=PROJECTS * len(wgbs_df), technique=wgbs_df['Assay Type'].tolist(), patient_id=wgbs_df['patient_id'].tolist(), sra_id=wgbs_df['Run'].tolist()),
    output: 
        temp(touch('results/OK'))