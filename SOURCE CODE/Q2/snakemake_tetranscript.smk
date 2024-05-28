"""
Analysis pipeline for single reads

Author: Ida Kups (iku)
"""
import os
import pathlib
import pandas as pd
import metadata_df

configfile: 'config_evaxion2.yaml'

localrules: all, ok

rule all:
    input: 'results/OK'

TMPDIR="/scratch"


prj2ref = {
     'DNMT3A':'mouse',
     'TET2':'mouse',
     'NZM':'human',
     'DNMT1':'human'
}


#FILF THE DF WITH CORRESPONDING VALUES
#IF ONLY ONE TYPE OF SEQUENCING DATA DON'T FORGET TO COMMENT OUT THE CORRESPONDING LINES IN EXPAND IN RULE OK
#wgbs_df = metadata_df.nzm_wgbs_df
# rnaseq_df = metadata_df.tet2_rnaseq_df
# PROJECTS = ['TET2'] #list with a projectname

# rnaseq_df = metadata_df.nzm_rnaseq_df
# PROJECTS = ['NZM'] #list with a projectname


rnaseq_df = metadata_df.DNMT_rnaseq
PROJECTS = ['TET2'] #list with a projectname

# rule TEtranscript: 
#     input:
#         m1 = '/work/iku/exp1019-cancer-epigenetics-and-ervs/01_BS_and_ERV/results/02_star/raw/TET2/1_RNA-Seq_SRR21643414.Aligned.out.bam',
#         m2 = '/work/iku/exp1019-cancer-epigenetics-and-ervs/01_BS_and_ERV/results/02_star/raw/TET2/2_RNA-Seq_SRR21643410.Aligned.out.bam',
#         m3 = '/work/iku/exp1019-cancer-epigenetics-and-ervs/01_BS_and_ERV/results/02_star/raw/TET2/3_RNA-Seq_SRR21643406.Aligned.out.bam',
#         m4 = '/work/iku/exp1019-cancer-epigenetics-and-ervs/01_BS_and_ERV/results/02_star/raw/TET2/4_RNA-Seq_SRR21643402.Aligned.out.bam',
#         m5 = '/work/iku/exp1019-cancer-epigenetics-and-ervs/01_BS_and_ERV/results/02_star/raw/TET2/5_RNA-Seq_SRR21643398.Aligned.out.bam',

#     output:
#        touch("results/04_TEtranscripts/TET2/tmp.OK") 
#     params:
#         GTF= lambda wc: config[prj2ref['TET2']]['GTF'],
#         TE_GTF= lambda wc: config[prj2ref['TET2']]['TE_GTF']
#     singularity:
#         config['TEtranscripts_sif']
#     shell:
#         """
#         TEtranscripts -t {input.m3} {input.m4} {input.m5} -c {input.m1} {input.m2} --mode multi --project TET2 --GTF {params.GTF} --TE {params.TE_GTF} --outdir $(dirname {output})
#         """
# nr 1 DONE
# nr 2 DONE
#nr 3 running
#nr 4 


# rule star: #untrimmed
#     input:
#         R1='/work/iku/exp1019-cancer-epigenetics-and-ervs/01_BS_and_ERV/results/01_fastq/DNMT1/4/RNA-Seq/SRR12849892_1.fastq',
#         R2='/work/iku/exp1019-cancer-epigenetics-and-ervs/01_BS_and_ERV/results/01_fastq/DNMT1/4/RNA-Seq/SRR12849892_2.fastq'
#     output:
#         tbam = 'results/022_star/raw/DNMT1/4.Aligned.toTranscriptome.out.bam',
#         gbam = 'results/022_star/raw/DNMT1/4.Aligned.out.bam',
#     params:
#         index = lambda wc: config[prj2ref['DNMT1']]['star'], 
#         prefix = lambda wc, output: output.gbam.replace("Aligned.out.bam",""),
#         temp_dir = TMPDIR 
#     singularity: config['singularity']
#     resources:
#         mem_gb=64,  
#         mem_mb=64000,
#         tmpdir=TMPDIR 
#     threads: 32
#     shell:
#         '''
#         STAR_TMPDIR=$(mktemp -u -p {params.temp_dir} star_tmp_XXXXXXXX)

#         STAR --readFilesIn {input.R1} {input.R2} \
#         --readFilesCommand cat \
#         --runThreadN {threads} \
#         --genomeDir {params.index} \
#         --outFileNamePrefix {params.prefix} \
#         --outSAMtype BAM Unsorted \
#         --outSAMunmapped Within KeepPairs \
#         --quantMode TranscriptomeSAM \
#         --outSAMattributes NH HI AS nM NM \
#         --outSAMattrRGline ID:SRR12849892 SM:SRR12849892 LB:SRR12849892 PL:ILLUMINA \
#         --outTmpDir $STAR_TMPDIR 
        
#         '''
#ZMIENIAC SRTA ID 

#ready to run 
rule TEtranscript: #DNMT1
    input:
        m1 = '/work/iku/exp1019-cancer-epigenetics-and-ervs/01_BS_and_ERV/results/022_star/raw/DNMT1/1.Aligned.out.bam',
        m2 = '/work/iku/exp1019-cancer-epigenetics-and-ervs/01_BS_and_ERV/results/022_star/raw/DNMT1/2.Aligned.out.bam',
        m3 = '/work/iku/exp1019-cancer-epigenetics-and-ervs/01_BS_and_ERV/results/022_star/raw/DNMT1/3.Aligned.out.bam',
        m4 = '/work/iku/exp1019-cancer-epigenetics-and-ervs/01_BS_and_ERV/results/022_star/raw/DNMT1/4.Aligned.out.bam',
        

    output:
       touch("results/04_TEtranscripts/DNMT1/tmp.OK") 
    params:
        GTF= lambda wc: config[prj2ref['DNMT1']]['GTF'],
        TE_GTF= lambda wc: config[prj2ref['DNMT1']]['TE_GTF']
    singularity:
        config['TEtranscripts_sif']
    shell:
        """
        TEtranscripts -t {input.m3} {input.m4}  -c {input.m1} {input.m2} --mode multi --project DNMT1 --GTF {params.GTF} --TE {params.TE_GTF} --outdir $(dirname {output})
        """


# rule TEtranscript: #nzm
#     input:
#         m1 = '/work/iku/exp1019-cancer-epigenetics-and-ervs/01_BS_and_ERV/results/02_star/raw/NZM/1_RNA-Seq_SRR12121871.Aligned.out.bam',
#         m2 = '/work/iku/exp1019-cancer-epigenetics-and-ervs/01_BS_and_ERV/results/02_star/raw/NZM/2_RNA-Seq_SRR12121872.Aligned.out.bam',
#         m3 = '/work/iku/exp1019-cancer-epigenetics-and-ervs/01_BS_and_ERV/results/02_star/raw/NZM/3_RNA-Seq_SRR12121873.Aligned.out.bam',
#         m4 = '/work/iku/exp1019-cancer-epigenetics-and-ervs/01_BS_and_ERV/results/02_star/raw/NZM/4_RNA-Seq_SRR12121874.Aligned.out.bam',
#         m5 = '/work/iku/exp1019-cancer-epigenetics-and-ervs/01_BS_and_ERV/results/02_star/raw/NZM/5_RNA-Seq_SRR12121875.Aligned.out.bam',
#         m6 = '/work/iku/exp1019-cancer-epigenetics-and-ervs/01_BS_and_ERV/results/02_star/raw/NZM/6_RNA-Seq_SRR12121876.Aligned.out.bam',
#         m7 = '/work/iku/exp1019-cancer-epigenetics-and-ervs/01_BS_and_ERV/results/02_star/raw/NZM/7_RNA-Seq_SRR12121877.Aligned.out.bam',
#         m8 = '/work/iku/exp1019-cancer-epigenetics-and-ervs/01_BS_and_ERV/results/02_star/raw/NZM/8_RNA-Seq_SRR12121878.Aligned.out.bam',
#         m9 = '/work/iku/exp1019-cancer-epigenetics-and-ervs/01_BS_and_ERV/results/02_star/raw/NZM/9_RNA-Seq_SRR12121879.Aligned.out.bam',
#     output:
#        touch("results/04_TEtranscripts/NZM/tmp.OK") 
#     params:
#         GTF= lambda wc: config[prj2ref['NZM']]['GTF'],
#         TE_GTF= lambda wc: config[prj2ref['NZM']]['TE_GTF']
#     singularity:
#         config['TEtranscripts_sif']
#     shell:
#         """
#         TEtranscripts -t {input.m1} {input.m2} {input.m3} {input.m6} {input.m7} -c {input.m4} {input.m5} {input.m8} {input.m9} --mode multi --project NZM --GTF {params.GTF} --TE {params.TE_GTF} --outdir $(dirname {output})
# #         """
# rule TEtranscript: #dnmt3a +/- vs +/+ 
#     input: 
#         m1 = '/work/iku/exp1019-cancer-epigenetics-and-ervs/01_BS_and_ERV/results/02_star/raw/DNMT3A/1_RNA-Seq_SRR2297229.Aligned.out.bam',
#         m2 = '/work/iku/exp1019-cancer-epigenetics-and-ervs/01_BS_and_ERV/results/02_star/raw/DNMT3A/2_RNA-Seq_SRR2297230.Aligned.out.bam',
#         m3 = '/work/iku/exp1019-cancer-epigenetics-and-ervs/01_BS_and_ERV/results/02_star/raw/DNMT3A/3_RNA-Seq_SRR2297231.Aligned.out.bam',
#         m4 = '/work/iku/exp1019-cancer-epigenetics-and-ervs/01_BS_and_ERV/results/02_star/raw/DNMT3A/4_RNA-Seq_SRR2297232.Aligned.out.bam',
#         m6 = '/work/iku/exp1019-cancer-epigenetics-and-ervs/01_BS_and_ERV/results/02_star/raw/DNMT3A/6_RNA-Seq_SRR2297234.Aligned.out.bam',
#     output:
#        touch("results/04_TEtranscripts/DNMT3A/tmp.OK") 
#     params:
#         GTF= lambda wc: config[prj2ref['DNMT3A']]['GTF'],
#         TE_GTF= lambda wc: config[prj2ref['DNMT3A']]['TE_GTF']
#     singularity:
#         config['TEtranscripts_sif']
#     shell:
#         """
#         TEtranscripts -t {input.m1} {input.m2}  -c {input.m3} {input.m4} {input.m6} --mode multi --project DNMT3A --GTF {params.GTF} --TE {params.TE_GTF} --outdir $(dirname {output})
#         """
# rule TEtranscript: #dnmt3a +/m vs +/+
#     input:
#         m3 = '/work/iku/exp1019-cancer-epigenetics-and-ervs/01_BS_and_ERV/results/02_star/raw/DNMT3A/3_RNA-Seq_SRR2297231.Aligned.out.bam',
#         m4 = '/work/iku/exp1019-cancer-epigenetics-and-ervs/01_BS_and_ERV/results/02_star/raw/DNMT3A/4_RNA-Seq_SRR2297232.Aligned.out.bam',
#         m5 = '/work/iku/exp1019-cancer-epigenetics-and-ervs/01_BS_and_ERV/results/02_star/raw/DNMT3A/5_RNA-Seq_SRR2297233.Aligned.out.bam',
#         m6 = '/work/iku/exp1019-cancer-epigenetics-and-ervs/01_BS_and_ERV/results/02_star/raw/DNMT3A/6_RNA-Seq_SRR2297234.Aligned.out.bam',
#         m7 = '/work/iku/exp1019-cancer-epigenetics-and-ervs/01_BS_and_ERV/results/02_star/raw/DNMT3A/7_RNA-Seq_SRR2297235.Aligned.out.bam',
#         m8 = '/work/iku/exp1019-cancer-epigenetics-and-ervs/01_BS_and_ERV/results/02_star/raw/DNMT3A/8_RNA-Seq_SRR2297236.Aligned.out.bam',
#         m9 = '/work/iku/exp1019-cancer-epigenetics-and-ervs/01_BS_and_ERV/results/02_star/raw/DNMT3A/9_RNA-Seq_SRR2297237.Aligned.out.bam',
#     output:
#        touch("results/04_TEtranscripts/{project}/{patient_id}/tmp.OK") 
#     params:
#         GTF= lambda wc: config[prj2ref[wc.project]]['GTF'],
#         TE_GTF= lambda wc: config[prj2ref[wc.project]]['TE_GTF']
#     singularity:
#         config['TEtranscripts_sif']
#     shell:
#         """
#         TEtranscripts -t {input.m5} {input.m7} {input.m8} {input.m9} -c {input.m3} {input.m4} {input.m6} --mode multi --project {wildcards.project} --GTF {params.GTF} --TE {params.TE_GTF} --outdir $(dirname {output})
#         """
rule ok:
    """Collects all output from above and writes the OK file for `rule all`."""
    input: 
        #a = expand(rules.TEtranscript.output, zip, project=PROJECTS * len(rnaseq_df), patient_id=rnaseq_df['patient_id'].tolist()),  
        rules.TEtranscript.output
        #rules.star.output
        
    # rules.my_rule.output
    output: 
        temp(touch('results/OK'))