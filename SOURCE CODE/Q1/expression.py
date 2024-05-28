import pandas as pd
import os
import string
import re 
import json 

expression_data_folder = '/work/iku/exp1019-cancer-epigenetics-and-ervs/00_mut_and_ERV/expression_data'
subfolder_list = os.listdir(expression_data_folder)
metadata_path = '/work/iku/exp1019-cancer-epigenetics-and-ervs/00_mut_and_ERV/metadata.cart.2024-02-12.json'
with open(metadata_path, 'r') as file:
    metadata = json.load(file)
    

def read_expression_data(file_path):
    #reads expression data into dataframe
    expression_df = pd.read_csv(file_path, sep = '\t')
    return expression_df




def get_patient_ID(filename, metadata):
    #extracts patient id from the metadata 
    patient_metadata = [patient_data for patient_data in metadata if patient_data.get('file_name') == filename]
    pattern = r"^([^\-]+\-[^\-]+\-[^\-]+)"

    if patient_metadata:
        if 'annotations' in patient_metadata[0]:
            patient_dict = patient_metadata[0]
            annotations = patient_dict['annotations']
            annotations = annotations[0]
            match = re.match(pattern, annotations['entity_submitter_id'])
            return match.group(1)
        elif 'associated_entities' in patient_metadata[0]:
            patient_dict = patient_metadata[0]
            associated_entities = patient_dict['associated_entities']
            associated_entities = associated_entities[0]
            entity_submitter_id = associated_entities['entity_submitter_id']
            
            match = re.match(pattern, entity_submitter_id)
            if match:
                return match.group(1)
            else:
                return None
        else:
            print('no annotations for patient')
            print(patient_metadata)
            return None
    else:
        print('no metadata for the patient')
        return None

def get_expression_folder(filename, expression_data_folder): 
    #finds an appropriate expression folder
    pattern = r'^([^.]+)'
    match = re.match(pattern, filename)

    if match:
        result = match.group(1)
        result = result.replace('_', '-')
        exp_folder_path = os.path.join(expression_data_folder, result)
        if os.path.exists(exp_folder_path) and os.path.isdir(exp_folder_path):
            return exp_folder_path
        else:
            print('folder doesnt exist')
            return None
    else: 
        print('folder not found')
        return None

def find_expression_file(patient_ID, expression_folder):
    #finds expression file within the expression folder
    for filename in os.listdir(expression_folder):
        if filename.startswith(patient_ID):
            return os.path.join(expression_folder, filename)
    print('no expression file: ', patient_ID)
    return None

def count_breadth(expression_df):
    #calcultes patient's ERV breadth
    breadth = len(expression_df[expression_df['gene_id'].str.startswith('Hsap38') & (expression_df['TPM'] > 1) ] )    
    return breadth

def count_abundance(expression_df):
    #calculates patient's ERV abundance
    expression_df = expression_df[expression_df['gene_id'].str.startswith('Hsap38')]
    abundance = expression_df['TPM'].sum()
    return abundance


