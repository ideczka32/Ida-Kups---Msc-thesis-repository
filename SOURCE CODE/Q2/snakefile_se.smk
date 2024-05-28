"""
Analysis pipeline for single reads

Author: Ida Kups (iku)
"""
import os
import pathlib
import pandas as pd
import metadata_df

configfile: 'config_evaxion.yaml'

localrules: all, ok

rule all:
    input: 'results/OK'

TMPDIR="/scratch"


prj2ref = {
     'DNMT3A':'mouse',
     'TET2':'mouse',
     'NZM':'human'
}


#FILF THE DF WITH CORRESPONDING VALUES
#IF ONLY ONE TYPE OF SEQUENCING DATA DON'T FORGET TO COMMENT OUT THE CORRESPONDING LINES IN EXPAND IN RULE OK
wgbs_df = metadata_df.nzm_wgbs_df
rnaseq_df = metadata_df.nzm_rnaseq_df
PROJECTS = ['NZM'] #list with a projectname




rule fetch_sra: #works the same regardless whether paired or single end i think cant use multiple threads
    output:
        # {sra_id} wc used globaly instead of "sample" so that naming scheme
        # follows SRA-ID to ensure correct function of prefetch so that the
        # whole directory is downlaed with potential reference sequences
        # this also allows for correct execution of fasterq-dump using the
        # SRA-id directory as input e.g "SRRXXXX" with resulting in full
        # download of all files again (which is default for fasterq-dump
        # if it cannot find the correct ID)
        sra=temp(directory("results/00_sra/{project}/{patient_id}/{technique}/{sra_id}")),
        #lambda wc: [temp(directory("results/00_sra/{project}/{patient_id}/{technique}/{sra_id}".format(sra_id))) for sra_id in patients2sra[wc.patient_id]]
    params:
        prefix=lambda wc, output: str(pathlib.Path(output[0]).parent),  #computes the prefix (parent directory) for the prefetch command
        ngc= lambda wc: os.path.realpath(config["ngc_token"])
    shadow:
        "minimal"
        #reduce amount of detailed information in the log file
    threads: 4
    resources:
        singleton=15,
        ntasks=10,
        cpus=50,
        mem_gb=50,
    log:
        "results/00_sra/{project}/{technique}/{patient_id}__{sra_id}.sra.log",
    shell:
        """
        prefetch --progress --ngc {params.ngc} \
        -v {wildcards.sra_id} -X 150G -O {params.prefix} \
         --resume yes 2> >(tee '{log}' >&2)
        """
        #-X max cache size --resume yes for resuming the interrupted downloads, standard arror 2> is captured and logged 


rule sra2fastq: #outputs only one read 
    input:
        sra = rules.fetch_sra.output,
    output:
        read="results/01_fastq/{project}/{patient_id}/{technique}/{sra_id}.fastq", 
    params:
        prefix=lambda wc, output: str(pathlib.Path(output.read).parent),
    threads: 8  #6-8 cores recommend from documantion of SRAToolkit 
    priority: 2
    resources:
        ntasks=10,
        cpus=50,
        mem_gb=50,
    shell:
        """
        fasterq-dump -O {params.prefix} --progress \
        --threads {threads} {input.sra} 
        """


rule trim: #input and output are now single files
    input:
        rules.sra2fastq.output,
    output:
        trimmed = "results/trimmed/{project}/{patient_id}/{technique}/{sra_id}_trimmed.fastq",
    resources:
        mem_gb=10,
    threads: 8
    shell: #using default settings as they did in the paper
        """
        cutadapt -j {threads} -o {output.trimmed} {input}
        """  



# rule star: #untrimmed doesn't throw ann error, but just uses one thread
#     input:
#         fastq = rules.sra2fastq.output
#     output:
#         tbam = temp('results/02_star/raw/{project}/{patient_id}_{technique}_{sra_id}.Aligned.toTranscriptome.out.bam'),
#         gbam = temp('results/02_star/raw/{project}/{patient_id}_{technique}_{sra_id}.Aligned.out.bam'),
#     params:
#         index = lambda wc: config[prj2ref[wc.project]]['star'], 
#         prefix = lambda wc, output: output.gbam.replace("Aligned.out.bam",""),
#         tmp_dir = lambda wc, output: output.gbam.replace(".Aligned.out.bam","") + "_tmpdir",
#     singularity: config['singularity']
#     resources:
#         tmpdir=TMPDIR,
#         ntasks=10,
#         cpus=50,
#         mem_gb=64,
#     threads: 32
#     shell: #--readFilesCommand gzip -cd \
#         '''
#         STAR --readFilesIn {input.fastq} \
#         --runThreadN {threads} \
#         --genomeDir {params.index} \
#         --outFileNamePrefix {params.prefix} \
#         --outSAMtype BAM Unsorted \
#         --outSAMunmapped Within KeepPairs \
#         --quantMode TranscriptomeSAM \
#         --outSAMattributes NH HI AS nM NM \
#         --outSAMattrRGline ID:{wildcards.patient_id} SM:{wildcards.patient_id} LB:{wildcards.patient_id} PL:ILLUMINA \
#         --outTmpDir {params.tmp_dir}
#         '''

rule star: #untrimmed doesn't throw ann error, but just uses one thread
    input:
        fastq = rules.sra2fastq.output
    output:
        tbam = 'results/02_star/raw/{project}/{patient_id}_{technique}_{sra_id}.Aligned.toTranscriptome.out.bam',
        gbam = 'results/02_star/raw/{project}/{patient_id}_{technique}_{sra_id}.Aligned.out.bam',
    params:
        index = lambda wc: config[prj2ref[wc.project]]['star'], 
        prefix = lambda wc, output: output.gbam.replace("Aligned.out.bam",""),
        tmp_dir = lambda wc, output: output.gbam.replace(".Aligned.out.bam","") + "_tmpdir",
    singularity: config['singularity']
    resources:
        tmpdir=TMPDIR,
        ntasks=10,
        cpus=50,
        mem_gb=64,
    threads: 32
    shell: #--readFilesCommand gzip -cd \
        '''
        STAR --readFilesIn {input.fastq} \
        --runThreadN {threads} \
        --genomeDir {params.index} \
        --outFileNamePrefix {params.prefix} \
        --outSAMtype BAM Unsorted \
        --outSAMunmapped Within KeepPairs \
        --quantMode TranscriptomeSAM \
        --outSAMattributes NH HI AS nM NM \
        --outSAMattrRGline ID:{wildcards.patient_id} SM:{wildcards.patient_id} LB:{wildcards.patient_id} PL:ILLUMINA \
        --outTmpDir {params.tmp_dir}
        '''


use rule star as star_trimmed with:
    input:
        fastq=rules.trim.output
    output:
        tbam = temp('results/02_star/trimmed/{project}/{patient_id}_{technique}_{sra_id}.Aligned.toTranscriptome.out.bam'),
        gbam = temp('results/02_star/trimmed/{project}/{patient_id}_{technique}_{sra_id}.Aligned.out.bam'),    

rule rsem:
    input: rules.star.output.tbam
    output:
        texpr = 'results/03_rsem/raw/{project}/{patient_id}_{technique}_{sra_id}.isoforms.results', #expression on transcriptome counts
        gexpr = 'results/03_rsem/raw/{project}/{patient_id}_{technique}_{sra_id}.genes.results',
    params:
        index = lambda wc: config[prj2ref[wc.project]]['rsem'],
        prefix = lambda wc, output: output.texpr.replace(".isoforms.results",""),
        tmp_dir = lambda wc, output: output.texpr.replace(".isoforms.results","") + "_tmpdir",
        threads = 28
    resources:
        tmpdir=TMPDIR,
        ntasks=10,
        cpus=50,
        mem_gb=48,
    singularity: config['singularity']
    threads: 32
    shell:
        '''
        rsem-calculate-expression --alignments --estimate-rspd \
        --num-threads {threads} \
        --temporary-folder {params.tmp_dir} \
        {input} {params.index} {params.prefix} &&
        rm {params.prefix}.transcript.bam
        '''

use rule rsem as rsem_trimmed with:
    input:
        rules.star_trimmed.output.tbam
    output:
        texpr = 'results/03_rsem/trimmed/{project}/{patient_id}_{technique}_{sra_id}.isoforms.results', #expression on transcriptome counts
        gexpr = 'results/03_rsem/trimmed/{project}/{patient_id}_{technique}_{sra_id}.genes.results',


rule fastqc: #on bisulfite data
    input: 
        reads=rules.sra2fastq.output.read
    output: touch("results/fastqc/{project}/{patient_id}/{technique}/{sra_id}.OK")
        #html_R1="results/fastqc/{project}/{patient_id}/{technique}/{sra_id}_fastqc.zip"
    params:
        extra = "--quiet",
        #tmp_dir = "results/fastqc/{project}/{patient_id}/{technique}"
    log:
        "results/fastqc/{project}/{patient_id}/logs/{technique}/{sra_id}.log"
    threads: 2
    resources:
        ntasks=10,
        cpus=50,
        mem_gb=5,
    singularity: config['fastqc_sif']
    shell:
        """
        fastqc -t {threads} {params.extra} -o $(dirname {output}) {input}
        """
        # mkdir -p {params.tmp_dir} && \
        # source /work/iku/exp1019-cancer-epigenetics-and-ervs/virt/bin/activate && \
        # /work/iku/exp1019-cancer-epigenetics-and-ervs/virt/bin/FastQC/fastqc -t {threads} -o {params.tmp_dir} {input.reads}


# rule TEcount: 
#     input:
#         rules.star.output.gbam
#     output:
#         "results/03_TEcount/{project}/{patient_id}_{technique}/{sra_id}.cntTable"
#     params:
#         GTF= lambda wc: config[prj2ref[wc.project]]['GTF'],
#         TE_GTF= lambda wc: config[prj2ref[wc.project]]['TE_GTF']
#     singularity:
#         config['TEtranscripts_sif']
#     shell:
#         """
#         TEcount -b {input} --format BAM --mode multi --project {wildcards.sra_id}  --GTF {params.GTF} --TE {params.TE_GTF} --outdir $(dirname {output})
#         """



        
rule fastqc_trim:
    input:
        reads=rules.trim.output
    output: touch("results/fastqc_trimmed/{project}/{patient_id}/{technique}/{sra_id}.OK")
        #html_R1="results/fastqc/{project}/{patient_id}/{technique}_trimmed/{sra_id}_fastqc.zip"
    params:
        extra = "--quiet",
        tmp_dir = "results/fastqc_trimmed/{project}/{patient_id}/{technique}"
    log:
        "results/fastqc/{project}/{patient_id}/logs/{technique}_trimmed/{sra_id}.log"
    threads: 2
    singularity: config['fastqc_sif']
    resources:
        ntasks=10,
        cpus=50,
        mem_gb=5,
    shell:
        """
        fastqc -t {threads} {params.extra} -o $(dirname {output}) {input}
        """

rule bismark_raw: #run bismark on untrimmed reads #should i use multicores
    input:
        reads=rules.sra2fastq.output.read
    output: bam= "results/bismark/raw/{project}/{technique}/{patient_id}_{sra_id}/{sra_id}_bismark_bt2.bam"
    params:
        index = lambda wc: config[prj2ref[wc.project]]['bismark'],
    singularity: config['singularity']
    resources:
        ntasks=10,
        cpus=50,
        mem_gb=64,
    threads: 32
    shell: #i think id doesnt work work with threads wildcard for some reason
        """
        output_file="{output.bam}"
        if [ ! -f "$output_file" ]; then
            bismark -q -o $(dirname "{output.bam}") --temp_dir ../scratch \
            --non_directional -p 4 --genome_folder {params.index} \
            --multicore 4 {input.reads}
        else
            echo "Output file $output_file already exists, skipping bismark execution."
        fi
        """  

rule methylation_call_raw: #get cytosine report for untrimmed reads 
    input:
        rules.bismark_raw.output
    output: touch("results/bismark/methylation_call_raw/{project}/{patient_id}_{technique}_{sra_id}/{sra_id}.OK")
    params:
        #prefix = "results/bismark/methylation_call_raw/{project}/{patient_id}_{technique}_{sra_id}",
        genome_folder = lambda wc: config[prj2ref[wc.project]]['genome_of_interest_dir'],
    singularity: config['singularity']
    resources:
        ntasks=10,
        cpus=50,
        mem_gb=24,
    threads: 6
    shell:
    #only Cpg context in the report
        """
        bismark_methylation_extractor -s --bedGraph --gzip \
        --cytosine_report --parallel {threads}  --genome_folder {params.genome_folder} \
        --comprehensive -o $(dirname {output}) {input} 
        """   

#same but for trimmed 
use rule bismark_raw as bismark_trimmed with:
    input:
        reads = rules.trim.output.trimmed
    output: bam = "results/bismark/trimmed/{project}/{technique}/{patient_id}_{sra_id}/{sra_id}_trimmed_bismark_bt2.bam"


use rule methylation_call_raw as methylation_call_trimmed with:
    input: rules.bismark_trimmed.output
    output: touch("results/bismark/methylation_call_trimmed/{project}/{patient_id}_{technique}_{sra_id}/{sra_id}.OK")


rule ok:
    """Collects all output from above and writes the OK file for `rule all`."""
    input:
        a = expand(rules.rsem.output, zip, project=PROJECTS * len(rnaseq_df), technique=rnaseq_df['Assay Type'].tolist(), patient_id=rnaseq_df['patient_id'].tolist(), sra_id=rnaseq_df['Run'].tolist()),
        b = expand(rules.rsem_trimmed.output, zip, project=PROJECTS * len(rnaseq_df), technique=rnaseq_df['Assay Type'].tolist(), patient_id=rnaseq_df['patient_id'].tolist(), sra_id=rnaseq_df['Run'].tolist()),
        c = expand(rules.methylation_call_raw.output, zip, project=PROJECTS * len(wgbs_df), technique=wgbs_df['Assay Type'].tolist(), patient_id=wgbs_df['patient_id'].tolist(), sra_id=wgbs_df['Run'].tolist()),
        d = expand(rules.methylation_call_trimmed.output, zip, project=PROJECTS * len(wgbs_df), technique=wgbs_df['Assay Type'].tolist(), patient_id=wgbs_df['patient_id'].tolist(), sra_id=wgbs_df['Run'].tolist()),
        # Add your outputs here, preferably using the correct expand syntax
    # rules.my_rule.output
    output: 
        temp(touch('results/OK'))



#a = expand(rules.rsem.output, zip, project=PROJECTS * len(rnaseq_df), technique=rnaseq_df['Assay Type'].tolist(), patient_id=rnaseq_df['patient_id'].tolist(), sra_id=rnaseq_df['Run'].tolist()),
#b = expand(rules.rsem_trimmed.output, zip, project=PROJECTS * len(rnaseq_df), technique=rnaseq_df['Assay Type'].tolist(), patient_id=rnaseq_df['patient_id'].tolist(), sra_id=rnaseq_df['Run'].tolist()),
        #c = expand(rules.methylation_call_raw.output, zip, project=PROJECTS * len(wgbs_df), technique=wgbs_df['Assay Type'].tolist(), patient_id=wgbs_df['patient_id'].tolist(), sra_id=wgbs_df['Run'].tolist()),
        #d = expand(rules.methylation_call_trimmed.output, zip, project=PROJECTS * len(wgbs_df), technique=wgbs_df['Assay Type'].tolist(), patient_id=wgbs_df['patient_id'].tolist(), sra_id=wgbs_df['Run'].tolist()),
        #e = expand(rules.fastqc.output, zip , project=PROJECTS* len(wgbs_df), technique=wgbs_df['Assay Type'].tolist(), patient_id=wgbs_df['patient_id'].tolist(), sra_id=wgbs_df['Run'].tolist()),
        #f = expand(rules.fastqc_trim.output, zip , project=PROJECTS* len(wgbs_df), technique=wgbs_df['Assay Type'].tolist(), patient_id=wgbs_df['patient_id'].tolist(), sra_id=wgbs_df['Run'].tolist()),