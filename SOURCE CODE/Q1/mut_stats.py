import os
import pandas as pd
import gzip
import re

def get_vcf(project_dir): 
    #get the vcf file of interest
    vcf_files = [f for f in os.listdir(project_dir) if f.endswith('.vcf.gz')]
    if(len(vcf_files) == 0):
        print("NO VCF IN THE FOLDER")
        return None
    return vcf_files[0]


def get_df_names(vcf_path):
    #gets the column names that will be used in the dataframe
    #get the information for the df header
    with gzip.open(vcf_path, "rt") as ifile:
          for line in ifile:
            if line.startswith("#CHROM"):
                  vcf_names = [x for x in line.split('\t')]
                  break
    ifile.close()
    return vcf_names

def get_vcf_df(project_dir,filename): 
    #creates the df out of the vcf file
    file_path = os.path.join(project_dir, filename)
    header_names = get_df_names(file_path)
    vcf_df = pd.read_csv(file_path, compression = 'gzip', comment = '#', chunksize=10000, delim_whitespace=True, header=None, names=header_names) 
    return vcf_df

def extract_effect(info): 
    #extracts the effect of the mutation from the info column 
    pattern = r'CSQ=.*?\|.*?\|(\w*)\|.*?'
    match = re.search(pattern, info)
    if match:
        effect = match.group(1)
        return effect
    else: 
        print('no effect info')
        return None



def get_mut_stats(project_dir, filename):
    #creates a dataframe stroing the information about the number of mutations in each of the gene on interest in a patient
    vcf_df_reader = get_vcf_df(project_dir, filename)

    mut_counts = {
        'Dnmt3a': 0,'Dnmt3a_HIGH': 0,'Dnmt3a_MODERATE': 0,'Dnmt3a_LOW': 0,'Dnmt3a_MODIFIER': 0,
        'Dnmt1': 0,'Dnmt1_HIGH': 0,'Dnmt1_MODERATE': 0,'Dnmt1_LOW': 0,'Dnmt1_MODIFIER': 0,
        'Dnmt3b': 0,'Dnmt3b_HIGH': 0,'Dnmt3b_MODERATE': 0,'Dnmt3b_LOW': 0,'Dnmt3b_MODIFIER': 0,
        'Dnmt3L': 0,'Dnmt3L_HIGH': 0,'Dnmt3L_MODERATE': 0,'Dnmt3L_LOW': 0,'Dnmt3L_MODIFIER': 0,
        'Tet1': 0,'Tet1_HIGH': 0,'Tet1_MODERATE': 0,'Tet1_LOW': 0,'Tet1_MODIFIER': 0,
        'Tet2': 0,'Tet2_HIGH': 0,'Tet2_MODERATE': 0,'Tet2_LOW': 0,'Tet2_MODIFIER': 0,
        'Tet3': 0,'Tet3_HIGH': 0,'Tet3_MODERATE': 0,'Tet3_LOW': 0,'Tet3_MODIFIER': 0,
        'total_mut' : 0
    }

    # Iterate through chunks
    iter = 0
    for vcf_df in vcf_df_reader:
        print(vcf_df.head())
        print(iter)
        iter += 1
        vcf_df['POS'] = vcf_df['POS'].astype(int)
        #add column with information about the effect
        vcf_df['mut_effect'] = vcf_df['INFO'].apply(extract_effect)
        
        #get the total count of mutations to see if mutations in our genes are above the average
        Dnmt3a = vcf_df[(vcf_df['#CHROM'] == 'chr2') & (vcf_df['POS'].between(25227874, 25341925))]
        Dnmt1 = vcf_df[(vcf_df['#CHROM'] == 'chr19') & (vcf_df['POS'].between(10133346, 10194953))]
        Dnmt3b = vcf_df[(vcf_df['#CHROM'] == 'chr20') & (vcf_df['POS'].between(32762385, 32809356))]
        Dnmt3L = vcf_df[(vcf_df['#CHROM'] == 'chr21') & (vcf_df['POS'].between(44246339, 44261897))]
        Tet1 = vcf_df[(vcf_df['#CHROM'] == 'chr10') & (vcf_df['POS'].between(68560337, 68694487))]
        Tet2 = vcf_df[(vcf_df['#CHROM'] == 'chr4') & (vcf_df['POS'].between(105146876, 105279803))]
        Tet3 = vcf_df[(vcf_df['#CHROM'] == 'chr2') & (vcf_df['POS'].between(73984910, 74108176))]

        mut_counts['total_mut'] += vcf_df['POS'].nunique()
        #Dnmt3a
        mut_counts['Dnmt3a'] += len(Dnmt3a)
        mut_counts['Dnmt3a_HIGH'] = len(Dnmt3a[Dnmt3a['mut_effect'] == 'HIGH'])
        mut_counts['Dnmt3a_MODERATE'] = len(Dnmt3a[Dnmt3a['mut_effect'] == 'MODERATE'])
        mut_counts['Dnmt3a_LOW'] = len(Dnmt3a[Dnmt3a['mut_effect'] == 'LOW'])
        mut_counts['Dnmt3a_MODIFIER'] = len(Dnmt3a[Dnmt3a['mut_effect'] == 'MODIFIER'])
        #Dnmt1
        mut_counts['Dnmt1'] += len(Dnmt1)
        mut_counts['Dnmt1_HIGH'] = len(Dnmt1[Dnmt1['mut_effect'] == 'HIGH'])
        mut_counts['Dnmt1_MODERATE'] = len(Dnmt1[Dnmt1['mut_effect'] == 'MODERATE'])
        mut_counts['Dnmt1_LOW'] = len(Dnmt1[Dnmt1['mut_effect'] == 'LOW'])
        mut_counts['Dnmt1_MODIFIER'] = len(Dnmt1[Dnmt1['mut_effect'] == 'MODIFIER'])
        #Dnmt3b
        mut_counts['Dnmt3b'] += len(Dnmt3b)
        mut_counts['Dnmt3b_HIGH'] = len(Dnmt3b[Dnmt3b['mut_effect'] == 'HIGH'])
        mut_counts['Dnmt3b_MODERATE'] = len(Dnmt3b[Dnmt3b['mut_effect'] == 'MODERATE'])
        mut_counts['Dnmt3b_LOW'] = len(Dnmt3b[Dnmt3b['mut_effect'] == 'LOW'])
        mut_counts['Dnmt3b_MODIFIER'] = len(Dnmt3b[Dnmt3b['mut_effect'] == 'MODIFIER'])
        #Dnmt3L
        mut_counts['Dnmt3L'] += len(Dnmt3L)
        mut_counts['Dnmt3L_HIGH'] = len(Dnmt3L[Dnmt3L['mut_effect'] == 'HIGH'])
        mut_counts['Dnmt3L_MODERATE'] = len(Dnmt3L[Dnmt3L['mut_effect'] == 'MODERATE'])
        mut_counts['Dnmt3L_LOW'] = len(Dnmt3L[Dnmt3L['mut_effect'] == 'LOW'])
        mut_counts['Dnmt3L_MODIFIER'] = len(Dnmt3L[Dnmt3L['mut_effect'] == 'MODIFIER'])
        #Tet1
        mut_counts['Tet1'] += len(Tet1)
        mut_counts['Tet1_HIGH'] = len(Tet1[Tet1['mut_effect'] == 'HIGH'])
        mut_counts['Tet1_MODERATE'] = len(Tet1[Tet1['mut_effect'] == 'MODERATE'])
        mut_counts['Tet1_LOW'] = len(Tet1[Tet1['mut_effect'] == 'LOW'])
        mut_counts['Tet1_MODIFIER'] = len(Tet1[Tet1['mut_effect'] == 'MODIFIER'])
        #Tet2
        mut_counts['Tet2'] += len(Tet2)
        mut_counts['Tet2_HIGH'] = len(Tet2[Tet2['mut_effect'] == 'HIGH'])
        mut_counts['Tet2_MODERATE'] = len(Tet2[Tet2['mut_effect'] == 'MODERATE'])
        mut_counts['Tet2_LOW'] = len(Tet2[Tet2['mut_effect'] == 'LOW'])
        mut_counts['Tet2_MODIFIER'] = len(Tet2[Tet2['mut_effect'] == 'MODIFIER'])
        #Tet3
        mut_counts['Tet3'] += len(Tet3)
        mut_counts['Tet3_HIGH'] = len(Tet3[Tet3['mut_effect'] == 'HIGH'])
        mut_counts['Tet3_MODERATE'] = len(Tet3[Tet3['mut_effect'] == 'MODERATE'])
        mut_counts['Tet3_LOW'] = len(Tet3[Tet3['mut_effect'] == 'LOW'])
        mut_counts['Tet3_MODIFIER'] = len(Tet3[Tet3['mut_effect'] == 'MODIFIER'])        
        
        

    return mut_counts


def add_annotations(mut_stats, annotation_df, itr):
    #adds info about the patient id, tumor id and sampling date
    mut_stats['patient_ID'] = annotation_df.loc[itr, 'id']
    mut_stats['tumor_ID'] = annotation_df.loc[itr, 'entity_id']
    mut_stats['sampling_date'] = annotation_df.loc[itr, 'created_datetime']
    return mut_stats





