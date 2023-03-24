import os
import pandas as pd
import numpy as np
from scipy.stats import norm
from scipy.stats import chi2
import json
import sys
import math
import matplotlib.pyplot as plt

# Load the DB file
# df_db : Data frame of accumulated Experimental result information - Abundance
path_db = os.path.dirname(os.path.abspath(__file__)) + "/input/db_abundance.csv"
df_db = pd.read_csv(path_db)

# Load the Experiment result file
# df_exp : Data frame of Experimental result information - Abundance
path_exp = os.path.dirname(os.path.abspath(__file__)) + "/input/experiment_result_abundance.csv"
df_exp = pd.read_csv(path_exp)

# Load the merged Valencia output file
# df_valencia : Data frame of merged Valencia output
path_valencia = os.path.dirname(os.path.abspath(__file__)) + "/input/VALENCIA_output_merged.csv"
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

# Update the statistic - <<Uncomment only when updating json_grs_statistics.json>>     
df_exp = df_db     
 
# Delete the diversity, observed rows
if (list(df_exp['taxa'][0:2]) == ['diversity', 'observed']) & (list(df_db['taxa'][0:2]) == ['diversity', 'observed']):
    df_exp.drop(df_exp.index[0:2], inplace=True)
    df_db.drop(df_db.index[0:2], inplace=True)
else:
    print("Check the diversity & observed rows in the exp file or db file")
    sys.exit()


# Load the Phenotype-Microbiome file
# df_beta : Data frame of of Phenotype-Microbiome information
path_beta = os.path.dirname(os.path.abspath(__file__)) + "/input/phenotype_microbiome.xlsx"
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
        
        condition_phen = (df_beta.phenotype == li_phenotype_ncbi_name[j][0]) & (df_beta.ncbi_name == li_phenotype_ncbi_name[j][1])

        abundance = 0 
        for idx_beta, row_beta in df_beta[condition_phen].iterrows():
             
            if (row_beta['beta'] == 1) & (row_beta['microbiome'][:3] in ['s__', 'g__']):
                condition = (df_exp.taxa == row_beta['microbiome'])
                if len(df_exp[condition]) > 0:
                    abundance += df_exp[condition][li_new_sample_name[i]].values[0]
                
            if (row_beta['beta'] == 1) & (pd.isna(row_beta['microbiome_subtract']) is False):
                li_micro_sub = row_beta['microbiome_subtract'].split('\n')
                    
                for micro_sub in li_micro_sub:
                    condition_sub = (df_exp.taxa == micro_sub)
                    if len(df_exp[condition_sub]) > 0:
                        abundance -= df_exp[condition_sub][li_new_sample_name[i]].values[0]
     
            
        json_abundance.append({"sample_name" : li_new_sample_name[i], "phenotype" : li_phenotype_ncbi_name[j][0], "ncbi_name" : li_phenotype_ncbi_name[j][1], "abundance" : abundance})
df_abundance = pd.DataFrame.from_dict(json_abundance)   

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

# Calculate the GRS 
# li_phenotype : Phenotype list 
# df_grs : Data frame of grs corresponding to specific phenotype and sample

df_grs = pd.DataFrame(index = li_phenotype, columns = li_new_sample_name)
df_grs = df_grs.fillna(0) 

for i in range(len(li_phenotype)):
    for j in range(len(li_new_sample_name)):
        condition_phen = (df_beta.phenotype == li_phenotype[i])   
        grs = 0
        
        for idx_beta, row_beta in df_beta[condition_phen].iterrows():
            condition_micro = (df_exp.taxa == row_beta['microbiome'])
            
            if (len(df_exp[condition_micro]) > 0):      
                x_i = df_exp[condition_micro][li_new_sample_name[j]].values[0]
                clr_x_i = math.log(x_i + 1e-15)  
                grs += clr_x_i * row_beta['beta']
                
        grs /= len(df_beta[condition_phen])       
        df_grs.loc[li_phenotype[i], li_new_sample_name[j]] = grs

# Histogram Plot - GRS - <<Uncomment only when updating json_grs_statistics.json>>
def save_histograms_to_file(df, filename):
    num_rows = df.shape[0]
    fig, axs = plt.subplots(num_rows, 1, figsize=(8, 6*num_rows))
    
    for i in range(num_rows):
        axs[i].hist(df.iloc[i,:], bins=10)
        axs[i].set_title(df.index.to_list()[i])
    
    plt.tight_layout()
    plt.savefig(filename)
    
save_histograms_to_file(df_grs, '/home/kbkim/vaginal_microbiome/output/grs_hist.png')    

# Sample Estimation - Population Standard deviation & Mean - <<Uncomment only when updating json_grs_statistics.json>>
# json_grs_statistics : Json of GRS statistics - phenotype, num_sample, mean_grs, std_grs

json_grs_statistics= []

for i in range(len(li_phenotype)):
  
    li_grs = list(df_grs.iloc[i])
    li_grs = [x for x in li_grs if np.isnan(x) == False]
  
    num_sample = len(li_grs)
    if num_sample > 0:
        std_grs_sample = np.std(li_grs)
        mean_grs_sample = np.mean(li_grs)
    
        chi_1 = (chi2.ppf(1-.025, df=num_sample-1))**0.5
        chi_2 = (chi2.ppf(1-.975, df=num_sample-1))**0.5
    
        std_grs_population = std_grs_sample * ((num_sample)**0.5) * (1/chi_1 + 1/chi_2) * 0.5
        json_grs_statistics.append({"phenotype" : li_phenotype[i], "num_sample" : num_sample, "mean_grs" : mean_grs_sample, "std_grs" : std_grs_population})
  
    else:
        json_grs_statistics.append({"phenotype" : li_phenotype[i], "num_sample" : 0, "mean_grs" : 'none', "std_grs" : 'none'})

# Save the list of the grs statistics in the DB - <<Uncomment only when updating json_grs_statistics.json>>

path_statistics = os.path.dirname(os.path.abspath(__file__)) + "/input/json_grs_statistics.json"

with open(path_statistics, 'w') as f:
    json.dump(json_grs_statistics, f, indent=4, ensure_ascii=False)

print('Update Complete')       
     
        
        
        
        
        
        
        