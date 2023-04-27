import os
import pandas as pd
import numpy as np
from scipy.stats import norm, chi2, percentileofscore
import json
import sys
import math
import matplotlib.pyplot as plt

# Load the DB file
# df_db : Data frame of accumulated Experimental result information - Abundance
path_db = os.path.abspath('') + "/input/db_abundance.csv"
df_db = pd.read_csv(path_db)

# Load the Experiment result file
# df_exp : Data frame of Experimental result information - Abundance
path_exp = os.path.abspath('') + "/input/experiment_result_abundance.csv"
df_exp = pd.read_csv(path_exp)

# Load the merged Valencia output file
# df_valencia : Data frame of merged Valencia output
path_valencia = os.path.abspath('') + "/input/VALENCIA_output_merged.csv"
df_valencia = pd.read_csv(path_valencia)

# Insert data into DB - Merge the data frame df_db & df_exp

try: 
    df_db = pd.merge(df_db, df_exp, how='outer',on='taxa', suffixes=['', '_right']) 
    df_db = df_db.fillna(0)
    df_db = df_db.filter(regex='^(?!.*_right).*') # Eliminate duplicate columns

    df_db_rev = df_db.set_index(keys=['taxa'], inplace=False, drop=True)    
    df_db_rev.to_csv(path_db)
    
except:
    print("Check the Experiment result file")
    sys.exit()

    
# Delete the diversity, observed rows
if (list(df_exp['taxa'][0:2]) == ['diversity', 'observed']) & (list(df_db['taxa'][0:2]) == ['diversity', 'observed']):
    df_exp = df_exp.iloc[2:,:]
    df_db = df_db.iloc[2:,:]
else:
    print("Check the diversity & observed rows in the exp file or db file")
    sys.exit()


# Load the Phenotype-Microbiome file
# df_beta : Data frame of of Phenotype-Microbiome information
path_beta = os.path.abspath('') + "/input/phenotype_microbiome.xlsx"
df_beta = pd.read_excel(path_beta)
df_beta.rename(columns = {"Disease": "phenotype", "NCBI name": "ncbi_name", "MIrROR name": "microbiome", "Health sign": "beta", "subtract": "microbiome_subtract"}, inplace=True)
df_beta = df_beta[["phenotype", "ncbi_name", "microbiome", "beta","microbiome_subtract"]]
df_beta['beta'] = df_beta['beta'].replace({'증가': 1, '감소': -1})

li_new_sample_name = list(df_exp.columns)[1:]  
li_phenotype = list(dict.fromkeys(df_beta['phenotype']))

## Top 5 NCBI name 
li_phenotype_ncbi_name = []

for idx, row in df_beta.iterrows(): 
    if [row['phenotype'], row['ncbi_name']] not in li_phenotype_ncbi_name:
        li_phenotype_ncbi_name.append([row['phenotype'], row['ncbi_name']])

json_abundance = []

for i in range(len(li_new_sample_name)):
    for j in range(len(li_phenotype_ncbi_name)):
        
        condition_phen = (df_beta.phenotype == li_phenotype_ncbi_name[j][0]) & (df_beta.ncbi_name == li_phenotype_ncbi_name[j][1]) & (df_beta.beta == 1) 

        abundance = 0 
        for idx_beta, row_beta in df_beta[condition_phen].iterrows(): 
            if row_beta['microbiome'][:3] in ['s__', 'g__']:
                condition = (df_exp.taxa == row_beta['microbiome'])
                if len(df_exp[condition]) > 0:
                    abundance += df_exp[condition][li_new_sample_name[i]].values[0]

                    if (pd.isna(row_beta['microbiome_subtract']) is False):
                        li_micro_sub = row_beta['microbiome_subtract'].split('\n')
                        for micro_sub in li_micro_sub:
                            condition_sub = (df_exp.taxa == micro_sub)
                            if len(df_exp[condition_sub]) > 0:
                                 abundance -= df_exp[condition_sub][li_new_sample_name[i]].values[0]
                            
                json_abundance.append({"sample_name" : li_new_sample_name[i], "phenotype" : li_phenotype_ncbi_name[j][0], "ncbi_name" : li_phenotype_ncbi_name[j][1], "abundance" : abundance})
                
df_abundance = pd.DataFrame.from_dict(json_abundance)   

df_abundance = df_abundance.drop_duplicates(['sample_name', 'phenotype', 'ncbi_name'], keep='last')

df_top_five = pd.DataFrame(columns = ["sample_name", "phenotype", "ncbi_name","abundance"])

for i in range(len(li_new_sample_name)):
    for j in range(len(li_phenotype)):
    
        condition = (df_abundance.sample_name == li_new_sample_name[i]) & (df_abundance.phenotype == li_phenotype[j])
        df_new = df_abundance[condition].sort_values(by=['abundance'], ascending=False).head(5)
        df_top_five = pd.concat([df_top_five,df_new])

df_top_five = df_top_five.set_index(keys=['sample_name'], inplace=False, drop=True)           
df_top_five.to_excel('/home/kbkim/vaginal_microbiome/output/top5.xlsx')


# Subtract the abundance - df_exp

for idx_beta, row_beta in df_beta.iterrows(): 
    li_micro_sub = []

    if pd.isna(row_beta['microbiome_subtract']) is False:
        li_micro_sub = row_beta['microbiome_subtract'].split('\n')
        
        for micro_sub in li_micro_sub:
            condition = (df_exp.taxa == row_beta['microbiome'])
            condition_sub = (df_exp.taxa == micro_sub)
            
            if len(df_exp[condition_sub]) > 0:
                
                for sample_name in li_new_sample_name:
                    df_exp.loc[condition, sample_name] -= df_exp[condition_sub][sample_name].values[0]

# Calculate the MRS (Microbiome Risk Score)
# li_phenotype : Phenotype list 
# df_mrs : Data frame of MRS corresponding to specific phenotype and sample
df_mrs = pd.DataFrame(index = li_new_sample_name, columns = li_phenotype)
df_mrs = df_mrs.fillna(0) 

for i in range(len(li_new_sample_name)):
    for j in range(len(li_phenotype)):
        condition_phen = (df_beta.phenotype == li_phenotype[j])   
        mrs = 0
        
        for idx_beta, row_beta in df_beta[condition_phen].iterrows():
            condition_micro = (df_exp.taxa == row_beta['microbiome'])
            
            if (len(df_exp[condition_micro]) > 0):      
                x_i = df_exp[condition_micro][li_new_sample_name[i]].values[0]
                ln_x_i = math.log(x_i + 1e-15)  
                mrs += ln_x_i * row_beta['beta']
            
            elif (len(df_exp[condition_micro]) == 0):      
                x_i = 0
                ln_x_i = math.log(x_i + 1e-15)  
                mrs += ln_x_i * row_beta['beta']   
            
        mrs /= len(df_beta[condition_phen])       
        df_mrs.loc[li_new_sample_name[i], li_phenotype[j]] = mrs

# Load the mrs_db dataframe
# df_mrs_db : Data frame of mrs in the DB
path_mrs_db = os.path.abspath('') + "/input/vaginal_mrs.xlsx"
df_mrs_db = pd.read_excel(path_mrs_db, index_col=0)        

# Calculation - Percentile Rank
# Create an empty data frame with the same index and columns as the df_mrs data frame
df_percentile_rank = pd.DataFrame(index = li_new_sample_name, columns = li_phenotype)

# Loop through all samples and phenotypes and calculate the percentile rank
for i in range(len(li_new_sample_name)):
    for j in range(len(li_phenotype)):
        df_percentile_rank.loc[li_new_sample_name[i], li_phenotype[j]] = (percentileofscore(list(df_mrs_db[li_phenotype[j]]), df_mrs.loc[li_new_sample_name[i], li_phenotype[j]], kind='mean')).round(1)

# Outliers
# Replace percentile ranks that are less than or equal to 5 with 5, and those that are greater than or equal to 95 with 95
for i in range(len(li_phenotype)):
    df_percentile_rank.loc[df_percentile_rank[li_phenotype[i]]<=5, li_phenotype[i]] = 5.0
    df_percentile_rank.loc[df_percentile_rank[li_phenotype[i]]>=95, li_phenotype[i]] = 95.0

# Replace missing values with the string 'None'    
df_percentile_rank = df_percentile_rank.fillna('None')
        
# Merge the valencia output file & df_percentile_rank

df_percentile_rank = pd.merge(df_percentile_rank, df_valencia[['sampleID', 'subCST', 'score', 'CST']], how='left',left_index=True, right_on = 'sampleID') 
df_percentile_rank = df_percentile_rank.set_index(keys=['sampleID'], inplace=False, drop=True)    

# Save the output file - Percentile Rank of the samples

path_vaginal_percentile_rank_output = os.path.abspath('') + "/output/vaginal_percentile_rank.xlsx"
df_percentile_rank.to_excel(path_vaginal_percentile_rank_output)

print('Analysis Complete')       
        