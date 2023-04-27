##<Usage: python Script.py {path_exp}>
### ex) python vaginal_percentile_rank_new.py "/home/kbkim/vaginal_microbiome/input/experiment_result_abundance.csv"

import os
import pandas as pd
from scipy.stats import percentileofscore, pearsonr
import sys
import math
import matplotlib.pyplot as plt
import numpy as np
from skbio.stats.composition import multiplicative_replacement, clr

# path_exp : Path of Merged Proportion file to analyze
path_exp = sys.argv[1] 

###################################
# MainClass
###################################
class VaginalDisease:
    def __init__(self, path_exp):
        """
        Initializes a VaginalDisease object.

        Parameters:
        path_exp (str): Path of Merged Proportion file to analyze.
        """        
        self.path_exp = path_exp
        
        curdir = os.path.abspath('')
        self.path_beta = f"{curdir}/input/phenotype_microbiome.xlsx"
        self.path_db = f"{curdir}/input/db_abundance.csv"
        self.path_mrs_db = f"{curdir}/input/vaginal_mrs.xlsx"  
        self.path_vaginal_percentile_rank_output = f"{curdir}/output/vaginal_percentile_rank.xlsx"
        self.path_topfive = f"{curdir}/output/top5.xlsx"
        self.path_valencia = f"{curdir}//input/VALENCIA_output_merged.csv"
        
        self.df_beta = None
        self.df_db = None
        self.df_exp = None
        self.df_mrs = None
        self.df_db_rev = None       
        self.df_mrs_db = None        
        self.df_percentile_rank = None
        self.df_abundance = None
        self.df_top_five = None
        self.df_new = None
        self.df_valencia = None

        self.li_new_sample_name = None
        self.li_phenotype = None
        self.li_phenotype_ncbi_name = None

    
    # Load the DB file
    # df_beta : Data frame of of Phenotype-Microbiome information
    # df_db : Data frame of accumulated Experimental result information - Abundance
    # df_exp : Data frame of Experimental result information - Abundance    
    def ReadDB(self):        
        rv = True
        rvmsg = "Success"
        
        try:       
            self.df_beta = pd.read_excel(self.path_beta)
            self.df_db = pd.read_csv(self.path_db)
            self.df_exp = pd.read_csv(self.path_exp)
            self.df_mrs_db = pd.read_excel(self.path_mrs_db, index_col=0) 
            self.df_valencia = pd.read_csv(self.path_valencia)

            self.df_beta.rename(columns = {"Disease": "phenotype", "NCBI name": "ncbi_name", "MIrROR name": "microbiome", "Health sign": "beta", "subtract": "microbiome_subtract"}, inplace=True)
            self.df_beta = self.df_beta[["phenotype", "ncbi_name", "microbiome", "beta", "microbiome_subtract"]]
            self.df_beta['beta'] = self.df_beta['beta'].replace({'증가': 1, '감소': -1})    
            
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print("Error has occurred in the ReadDB process")
            sys.exit()            
        return rv, rvmsg

    def InsertDataDB(self): 
        """
        Inserts data into the database by merging the data frames df_db and df_exp.

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """        
        rv = True
        rvmsg = "Success"
        
        try: 
            self.df_db = pd.merge(self.df_db, self.df_exp, how='outer',on='taxa', suffixes=['', '_right']) 
            self.df_db = self.df_db.fillna(0)
            self.df_db = self.df_db.filter(regex='^(?!.*_right).*') # Eliminate duplicate columns
                    
            self.df_db_rev = self.df_db.set_index(keys=['taxa'], inplace=False, drop=True)    
            self.df_db_rev.to_csv(self.path_db)

        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print("Error has occurred in the InsertDataDB process")
            sys.exit()
        return rv, rvmsg
    
    def TopFive(self): 
        try:     
            # li_new_sample_name : Sample name list 
            # li_phenotype : Phenotype list 
            self.li_new_sample_name = list(self.df_exp.columns)[1:]  
            self.li_phenotype = list(dict.fromkeys(self.df_beta['phenotype']))
            
            self.li_phenotype_ncbi_name = []

            for idx, row in self.df_beta.iterrows(): 
                if [row['phenotype'], row['ncbi_name']] not in self.li_phenotype_ncbi_name:
                    self.li_phenotype_ncbi_name.append([row['phenotype'], row['ncbi_name']])

            json_abundance = []

            for i in range(len(self.li_new_sample_name)):
                for j in range(len(self.li_phenotype_ncbi_name)):

                    condition_phen = (self.df_beta.phenotype == self.li_phenotype_ncbi_name[j][0]) & (self.df_beta.ncbi_name == self.li_phenotype_ncbi_name[j][1]) & (self.df_beta.beta == 1) 

                    abundance = 0 
                    for idx_beta, row_beta in self.df_beta[condition_phen].iterrows(): 
                        if row_beta['microbiome'][:3] in ['s__', 'g__']:
                            condition = (self.df_exp.taxa == row_beta['microbiome'])
                            if len(self.df_exp[condition]) > 0:
                                abundance += self.df_exp[condition][self.li_new_sample_name[i]].values[0]

                                if (pd.isna(row_beta['microbiome_subtract']) is False):
                                    li_micro_sub = row_beta['microbiome_subtract'].split('\n')
                                    for micro_sub in li_micro_sub:
                                        condition_sub = (self.df_exp.taxa == micro_sub)
                                        if len(self.df_exp[condition_sub]) > 0:
                                             abundance -= self.df_exp[condition_sub][self.li_new_sample_name[i]].values[0]

                            json_abundance.append({"sample_name" : self.li_new_sample_name[i], "phenotype" : self.li_phenotype_ncbi_name[j][0], "ncbi_name" : self.li_phenotype_ncbi_name[j][1], "abundance" : abundance})

            self.df_abundance = pd.DataFrame.from_dict(json_abundance)   

            self.df_abundance = self.df_abundance.drop_duplicates(['sample_name', 'phenotype', 'ncbi_name'], keep='last')

            self.df_top_five = pd.DataFrame(columns = ["sample_name", "phenotype", "ncbi_name","abundance"])

            for i in range(len(self.li_new_sample_name)):
                for j in range(len(self.li_phenotype)):

                    condition = (self.df_abundance.sample_name == li_new_sample_name[i]) & (self.df_abundance.phenotype == li_phenotype[j])
                    self.df_new = self.df_abundance[condition].sort_values(by=['abundance'], ascending=False).head(5)
                    self.df_top_five = pd.concat([self.df_top_five,self.df_new])

            self.df_top_five = self.df_top_five.set_index(keys=['sample_name'], inplace=False, drop=True)           
            self.df_top_five.to_excel(self.path_topfive)    
    
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print("Error has occurred in the TopFive process")
            sys.exit()
    
        return rv, rvmsg    
    
    
    def SubtractAbundance(self): 
        rv = True
        rvmsg = "Success"
        
        try: 
            # Delete the diversity, observed rows
            if (list(self.df_exp['taxa'][0:2]) == ['diversity', 'observed']) & (list(self.df_db['taxa'][0:2]) == ['diversity', 'observed']):
                self.df_exp = self.df_exp.iloc[2:,:]
                self.df_db = self.df_db.iloc[2:,:]
            
            # li_new_sample_name : Sample name list 
            # li_phenotype : Phenotype list 
            self.li_new_sample_name = list(self.df_exp.columns)[1:]  
            self.li_phenotype = list(dict.fromkeys(self.df_beta['phenotype']))
            
            # Subtract the abundance - df_exp
            for idx_beta, row_beta in self.df_beta.iterrows(): 
                li_micro_sub = []

                if pd.isna(row_beta['microbiome_subtract']) is False:
                    li_micro_sub = row_beta['microbiome_subtract'].split('\n')

                    for micro_sub in li_micro_sub:
                        condition = (self.df_exp.taxa == row_beta['microbiome'])
                        condition_sub = (self.df_exp.taxa == micro_sub)

                        if len(self.df_exp[condition_sub]) > 0:

                            for sample_name in self.li_new_sample_name:
                                self.df_exp.loc[condition, sample_name] -= self.df_exp[condition_sub][sample_name].values[0]                
           
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print("Check the diversity & observed rows in the exp file or db file")
            sys.exit()
    
        return rv, rvmsg

    def CalculateMRS(self): 
        rv = True
        rvmsg = "Success"
        
        try:                
            # df_mrs : Data frame of MRS corresponding to specific phenotype and sample
            self.df_mrs = pd.DataFrame(index = self.li_new_sample_name, columns = self.li_phenotype)
            self.df_mrs = self.df_mrs.fillna(0) 

            for i in range(len(self.li_new_sample_name)):
                for j in range(len(self.li_phenotype)):
                    condition_phen = (self.df_beta.phenotype == self.li_phenotype[j])   
                    mrs = 0

                    for idx_beta, row_beta in self.df_beta[condition_phen].iterrows():
                        condition_micro = (self.df_exp.taxa == row_beta['microbiome'])

                        if (len(self.df_exp[condition_micro]) > 0):      
                            abundance = self.df_exp[condition_micro][self.li_new_sample_name[i]].values[0] 
                            mrs += row_beta['beta'] * math.log10(100*abundance + 1)  

                    mrs /= len(self.df_beta[condition_phen])       
                    self.df_mrs.loc[self.li_new_sample_name[i], self.li_phenotype[j]] = mrs

        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print("Error has occurred in the CalculateMRS process")
            sys.exit()
    
        return rv, rvmsg                        

    
    def CalculatePercentileRank(self): 
        rv = True
        rvmsg = "Success"
        
        try:      
            # Create an empty data frame with the same index and columns as the df_mrs data frame
            self.df_percentile_rank = pd.DataFrame(index = self.li_new_sample_name, columns = self.li_phenotype)
            # Loop through all samples and phenotypes and calculate the percentile rank
            for i in range(len(self.li_new_sample_name)):
                for j in range(len(self.li_phenotype)):
                    self.df_percentile_rank.loc[self.li_new_sample_name[i], self.li_phenotype[j]] = (percentileofscore(list(self.df_mrs_db[self.li_phenotype[j]]), self.df_mrs.loc[self.li_new_sample_name[i], self.li_phenotype[j]], kind='mean')).round(1)
                 
            # Outliers
            # Replace percentile ranks that are less than or equal to 5 with 5, and those that are greater than or equal to 95 with 95
            for i in range(len(self.li_phenotype)):
                self.df_percentile_rank.loc[self.df_percentile_rank[self.li_phenotype[i]]<=5, self.li_phenotype[i]] = 5.0
                self.df_percentile_rank.loc[self.df_percentile_rank[self.li_phenotype[i]]>=95, self.li_phenotype[i]] = 95.0        
        
            # Replace missing values with the string 'None'    
            self.df_percentile_rank = self.df_percentile_rank.fillna('None')
            
            # Merge the valencia output file & df_percentile_rank

            self.df_percentile_rank = pd.merge(self.df_percentile_rank, self.df_valencia[['sampleID', 'subCST', 'score', 'CST']], how='left',left_index=True, right_on = 'sampleID') 
            self.df_percentile_rank = self.df_percentile_rank.set_index(keys=['sampleID'], inplace=False, drop=True)  
            
            # Save the output file - Percentile Rank of the samples
            self.df_percentile_rank.to_excel(self.path_vaginal_percentile_rank_output)

            print('Analysis Complete')         
            
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print("Error has occurred in the CalculatePercentileRank process")
            sys.exit()
    
        return rv, rvmsg                  


    
####################################
# main
####################################
if __name__ == '__main__':
    
    vaginal = VaginalDisease(path_exp)
    vaginal.ReadDB()
    vaginal.InsertDataDB()
    vaginal.SubtractAbundance()
    vaginal.CalculateMRS()    
    vaginal.CalculatePercentileRank()
    