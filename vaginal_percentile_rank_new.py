##<Usage: python Script.py {path_exp}>
### ex) python vaginal_percentile_rank_new.py "/home/kbkim/vaginal_microbiome/input/experiment_result_abundance.csv"

import os
import pandas as pd
from scipy.stats import percentileofscore, pearsonr
import sys
import math
import matplotlib.pyplot as plt
import numpy as np

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
        self.path_mrs_db = f"{curdir}/input/vaginal_mrs.xlsx"  
        self.path_vaginal_percentile_rank_output = f"{curdir}/output/vaginal_percentile_rank.xlsx"
        self.path_valencia = f"{curdir}/input/VALENCIA_output_merged.csv"
        self.path_harmful = f"{curdir}/output/harmful.xlsx"
        self.path_beneficial = f"{curdir}/output/beneficial.xlsx"    
        
        self.df_beta = None
        self.df_exp = None
        self.df_mrs = None
        self.df_mrs_db = None        
        self.df_percentile_rank = None
        self.df_abundance = None
        self.df_beneficial = None
        self.df_harmful = None        
        self.df_new = None
        self.df_valencia = None

        self.li_new_sample_name = None
        self.li_phenotype = None
        self.li_phenotype_ncbi_name = None

    
    # Load the DB file
    # df_beta : Data frame of of Phenotype-Microbiome information
    # df_exp : Data frame of Experimental result information - Abundance    
    def ReadDB(self):   
        """
        Read the data.

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """            
        rv = True
        rvmsg = "Success"
        
        try:       
            self.df_beta = pd.read_excel(self.path_beta)
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
    
    def SubtractAbundance(self): 
        """
        Subtract the abundance for each microbiome in the df_exp.

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """          
        rv = True
        rvmsg = "Success"
        
        try: 
            # Delete the diversity, observed rows
            if (list(self.df_exp['taxa'][0:2]) == ['diversity', 'observed']):
                self.df_exp = self.df_exp.iloc[2:,:]
            
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

    def BeneficialMicrobiome(self):     
        """
        Save the list of Beneficial Microbiome as an Excel file.

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """          
        rv = True
        rvmsg = "Success"
        
        try:              
            
            self.li_phenotype_ncbi_name = []
            for idx, row in self.df_beta.iterrows(): 
                if [row['phenotype'], row['ncbi_name']] not in self.li_phenotype_ncbi_name:
                    self.li_phenotype_ncbi_name.append([row['phenotype'], row['ncbi_name']]) 
                    
            json_abundance = []

            for i in range(len(self.li_new_sample_name)):
                for j in range(len(self.li_phenotype_ncbi_name)):

                    condition_phen = (self.df_beta.phenotype == self.li_phenotype_ncbi_name[j][0]) & (self.df_beta.ncbi_name == self.li_phenotype_ncbi_name[j][1]) & (self.df_beta.beta == -1) 

                    abundance = 0 
                    for idx_beta, row_beta in self.df_beta[condition_phen].iterrows(): 
                        condition = (self.df_exp.taxa == row_beta['microbiome'])
                        if len(self.df_exp[condition]) > 0:
                            abundance += self.df_exp[condition][self.li_new_sample_name[i]].values[0]

                        json_abundance.append({"sample_name" : self.li_new_sample_name[i], "phenotype" : self.li_phenotype_ncbi_name[j][0], "ncbi_name" : self.li_phenotype_ncbi_name[j][1], "abundance" : abundance})

            self.df_abundance = pd.DataFrame.from_dict(json_abundance)   

            self.df_abundance = self.df_abundance.drop_duplicates(['sample_name', 'phenotype', 'ncbi_name'], keep='last')

            self.df_beneficial = pd.DataFrame(columns = ["sample_name", "phenotype", "ncbi_name","abundance"])

            for i in range(len(self.li_new_sample_name)):
                for j in range(len(self.li_phenotype)):

                    condition = (self.df_abundance.sample_name == self.li_new_sample_name[i]) & (self.df_abundance.phenotype == self.li_phenotype[j])
                    self.df_new = self.df_abundance[condition].sort_values(by=['abundance'], ascending=False)
                    self.df_beneficial = pd.concat([self.df_beneficial,self.df_new])

            self.df_beneficial = self.df_beneficial.set_index(keys=['sample_name'], inplace=False, drop=True)           
            self.df_beneficial.to_excel(self.path_beneficial)    
    
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print("Error has occurred in the BeneficialMicrobiome process")
            sys.exit()
    
        return rv, rvmsg 
    
    def HarmfulMicrobiome(self):
        """
        Save the list of Harmful Microbiome as an Excel file.

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """  
        rv = True
        rvmsg = "Success"
        
        try:     
            json_abundance = []

            for i in range(len(self.li_new_sample_name)):
                for j in range(len(self.li_phenotype_ncbi_name)):

                    condition_phen = (self.df_beta.phenotype == self.li_phenotype_ncbi_name[j][0]) & (self.df_beta.ncbi_name == self.li_phenotype_ncbi_name[j][1]) & (self.df_beta.beta == 1) 

                    abundance = 0 
                    for idx_beta, row_beta in self.df_beta[condition_phen].iterrows(): 
                        condition = (self.df_exp.taxa == row_beta['microbiome'])
                        if len(self.df_exp[condition]) > 0:
                            abundance += self.df_exp[condition][self.li_new_sample_name[i]].values[0]

                        json_abundance.append({"sample_name" : self.li_new_sample_name[i], "phenotype" : self.li_phenotype_ncbi_name[j][0], "ncbi_name" : self.li_phenotype_ncbi_name[j][1], "abundance" : abundance})

            self.df_abundance = pd.DataFrame.from_dict(json_abundance)   

            self.df_abundance = self.df_abundance.drop_duplicates(['sample_name', 'phenotype', 'ncbi_name'], keep='last')

            self.df_harmful = pd.DataFrame(columns = ["sample_name", "phenotype", "ncbi_name","abundance"])

            for i in range(len(self.li_new_sample_name)):
                for j in range(len(self.li_phenotype)):

                    condition = (self.df_abundance.sample_name == self.li_new_sample_name[i]) & (self.df_abundance.phenotype == self.li_phenotype[j])
                    self.df_new = self.df_abundance[condition].sort_values(by=['abundance'], ascending=False)
                    self.df_harmful = pd.concat([self.df_harmful,self.df_new])

            self.df_harmful = self.df_harmful.set_index(keys=['sample_name'], inplace=False, drop=True)           
            self.df_harmful.to_excel(self.path_harmful)    
    
        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print("Error has occurred in the HarmfulMicrobiome process")
            sys.exit()
    
        return rv, rvmsg         
    def CalculateMRS(self): 
        """
        Calculate the MRS (Microbiome Risk Score).

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """ 
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
        """
        Calculate the Percentile Rank and Save the Percentile Rank data as an Excel file.

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """         
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
    vaginal.SubtractAbundance()
    vaginal.BeneficialMicrobiome()    
    vaginal.HarmfulMicrobiome() 
    vaginal.CalculateMRS()    
    vaginal.CalculatePercentileRank()
    