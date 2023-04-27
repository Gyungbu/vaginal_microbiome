##<Usage: python Script.py {path_exp}>
### ex) python vaginal_update_mrs_new.py "/home/kbkim/vaginal_microbiome/input/experiment_result_abundance.csv"

import os
import pandas as pd
import numpy as np
from scipy.stats import norm
from scipy.stats import chi2
import json
import sys
import math
import matplotlib.pyplot as plt

# path_exp : Path of Merged Proportion file to analyze
path_exp = sys.argv[1]

# Common Functions #
# Histogram Plot - mrs 
def save_histograms_to_file(df, filename):
    num_rows = df.shape[1]
    fig, axs = plt.subplots(num_rows, 1, figsize=(8, 6*num_rows))
    
    for i in range(num_rows):
        axs[i].hist(df.iloc[:,i], bins=10)
        axs[i].set_title(df.columns.to_list()[i])
    
    plt.tight_layout()
    plt.savefig(filename)
    
###################################
# MainClass
###################################
class VaginalDiseaseUpdateMRS:
    def __init__(self, path_exp):
        """
        Initializes a VaginalDiseaseUpdateMRS object.

        Parameters:
        path_exp (str): Path of Merged Proportion file to analyze.  
        """   
        
        self.path_exp = path_exp
        
        curdir = os.path.abspath('')
        self.path_beta = f"{curdir}/input/phenotype_microbiome.xlsx"
        self.path_db = f"{curdir}/input/db_abundance.csv"
        self.path_valencia = f"{curdir}/input/VALENCIA_output_merged.csv"             
        self.path_mrs_db = os.path.abspath('') + "/input/vaginal_mrs.xlsx"  
        self.path_hist = f"{curdir}/output/mrs_hist.png"

        self.df_beta = None
        self.df_db = None
        self.df_exp = None
        self.df_mrs = None
        self.df_db_rev = None

        self.li_new_sample_name = None
        self.li_phenotype = None
        self.li_microbiome = None            

            
    # Load the DB file
    # df_beta : Data frame of of Phenotype-Microbiome information
    # df_db : Data frame of accumulated Experimental result information - Abundance
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
            self.df_db = pd.read_csv(self.path_db)
            self.df_exp = pd.read_csv(self.path_exp)
    
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


    # Insert data into DB - Merge the data frame df_db & df_exp
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

            # Update the data - Convert df_exp to df_db
            self.df_exp = self.df_db            
            
            self.df_db_rev = self.df_db.set_index(keys=['taxa'], inplace=False, drop=True)    
            self.df_db_rev.to_csv(self.path_db)

        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print("Error has occurred in the InsertDataDB process")
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

    
    def UpdateMRS(self): 
        """
        Save the MRS data as an Excel file & Save the histogram as an png file.

        Returns:
        A tuple (success, message), where success is a boolean indicating whether the operation was successful,
        and message is a string containing a success or error message.
        """ 
        rv = True
        rvmsg = "Success"
        
        try:                 
            # Save the df_mrs
            self.df_mrs.to_excel(self.path_mrs_db)
   
            # Histogram Plot - mrs 
            save_histograms_to_file(self.df_mrs, self.path_hist)    
            print('Update Complete') 

        except Exception as e:
            print(str(e))
            rv = False
            rvmsg = str(e)
            print("Error has occurred in the UpdateMRS process")
            sys.exit()
    
        return rv, rvmsg       
    
    
####################################
# main
####################################
if __name__ == '__main__':
    
    vaginal = VaginalDiseaseUpdateMRS(path_exp)
    vaginal.ReadDB()
    vaginal.InsertDataDB()
    vaginal.SubtractAbundance()
    vaginal.CalculateMRS()       
    vaginal.UpdateMRS() 