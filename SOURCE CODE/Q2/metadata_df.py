import pandas as pd

def extract_patient_id(sample_name):
    return sample_name[-1]

def pop_rows(df, values_in_column, column_name):
    indexes = df[df[column_name].isin(values_in_column)].index
    dropped = df.drop(indexes, axis = 0)
    return dropped 



#NZM CELL LINE
#rna-seq PAIRED layout Bisulfite-seq SINGLE
#N.. , Skmel -> melanoma cell lines
#Hem -> normal melanocytes
NZM_df = pd.read_csv('./datasets_metadata/NZM_cell_line.txt', sep=',')
unpaired_experiments_NZM = ['Metastatic melanoma cell line NZM15_I', 'Metastatic melanoma cell line NZM9_II', 'Melanocyte cell line MelSt', 'Melanocyte cell line HemALP']
#remove unpaired experiments rows
NZM_df = pop_rows(NZM_df, unpaired_experiments_NZM, 'cell_line')

#rename incorrect cell
NZM_df.loc[9,'cell_line'] = 'Melanocyte cell line HemNLP'

nzm_wgbs_df = NZM_df[NZM_df['Assay Type'] == 'Bisulfite-Seq']
nzm_rnaseq_df=NZM_df[NZM_df['Assay Type'] == 'RNA-Seq']

# #add patient ids 
start_id = 1
end_id = start_id + len(nzm_wgbs_df)

nzm_wgbs_df.insert(0, 'patient_id', range(start_id, end_id))
nzm_rnaseq_df.insert(0, 'patient_id', range(start_id, end_id))

# #ADD INFO ABOUT WHAT GENOME WE'RE MAPPING TO 
# #TET2 DATASET
TET2_df = pd.read_csv('./datasets_metadata/TET2_metadata.txt')
TET2_df['patient_id'] = TET2_df['Sample Name'].apply(extract_patient_id)
unpaired_experiments_TET2 = ['GSM6594850']
TET2_df = pop_rows(TET2_df, unpaired_experiments_TET2, 'Sample Name')

# REMOVE MULTIPLE RUNS - can be skipped
TET2_df.drop_duplicates(subset =['BioSample'], keep='first', inplace=True, ignore_index=True)

tet2_wgbs_df=TET2_df[TET2_df['Assay Type'] == 'Bisulfite-Seq']
tet2_rnaseq_df=TET2_df[TET2_df['Assay Type'] == 'RNA-Seq']


# #DNMT3A DATASET
DNMT_wgbs = pd.read_csv('./datasets_metadata/dnmt_WGBS.txt')
DNMT_rnaseq = pd.read_csv('./datasets_metadata/dnmt_RNASEQ.txt')

# #add patient ids 
start_id = 1
end_id = start_id + len(DNMT_wgbs)


DNMT_wgbs.insert(0, 'patient_id', range(start_id, end_id))
DNMT_rnaseq.insert(0, 'patient_id', range(start_id, end_id))


#DNMT1 & DNMT3B 
DNMT1_rnaseq = pd.read_csv('./datasets_metadata/DNMT1_RnaSeq.txt')
DNMT1_wgbs = pd.read_csv('./datasets_metadata/DNMT1_WGBS.txt')

#pop patients without wgbs data
DNMT1_rnaseq.drop([2,5], inplace = True, )
DNMT1_rnaseq.reset_index(drop=True, inplace=True)

start_id = 1
end_id = start_id + len(DNMT1_wgbs)


DNMT1_wgbs.insert(0, 'patient_id', range(start_id, end_id))
DNMT1_rnaseq.insert(0, 'patient_id', range(start_id, end_id))

