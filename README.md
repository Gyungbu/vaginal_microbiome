# vaginal_microbiome : Calculate the MRS(Microbiome Risk Score) and percentile ranks of vaginal samples.

![Python](https://img.shields.io/badge/Python-v3.9.0-blue.svg?style=flat&logo=python)&nbsp;
![Pandas](https://img.shields.io/badge/pandas-v2.0.3-blue.svg?style=flat&logo=pandas)&nbsp;
![Numpy](https://img.shields.io/badge/NumPy-v1.25.0-blue.svg?style=flat&logo=numpy)&nbsp;
![Scipy](https://img.shields.io/badge/SciPy-v1.11.1-blue.svg?style=flat&logo=scipy)&nbsp;

## Installation

You can install the vaginal_microbiome with following command.
	
	git clone https://github.com/Gyungbu/vaginal_microbiome.git
 
The list of required packages for `script` is shown in the `requirements.txt` file. When you run the command below, the required packages will be downloaded. (version : `python 3.9.0`)
	
	conda create -n env_vaginal_microbiome
	conda activate env_vaginal_microbiome
	conda install pip  
	conda install python=3.9.0
	pip install -r ./vaginal_microbiome/requirements.txt 



# vaginal_update_mrs : (Optional) Update Reference Data
## How to use

### 1. Prepare Merged Proportion File
Place the csv file of your Merged Proportion File in the `./vaginal_microbiome/input/` folder.

Caveats: 

1. The first column of the Merged Proportion File is in the order of taxa, diversity, observed, and microbiome list.
2. From the second column of the Merged Proportion File, the sample name, the diversity value of the sample, the number of species in the sample, and the relative abundance value of each microbiome should be listed in the order.

### 2. Run vaginal_update_mrs
To run vaginal_update_mrs,
 
Run the command below:
  
    python ./vaginal_microbiome/vaginal_update_mrs.py {path_exp}
    ### ex) python vaginal_update_mrs.py "/home/kbkim/vaginal_microbiome/input/experiment_result_abundance.csv"

When vaginal_update_mrs is executed as above, the file `db_abundance.csv`, `vaginal_mrs.xlsx` will be created or modified in the `./vaginal_microbiome/input/` folder.
And the file `mrs_hist.png` will be created or modified in the `./vaginal_microbiome/output/` folder.

# vaginal_percentile_rank : Calculate the MRS(Microbiome Risk Score) and percentile ranks of vaginal samples.
## How to use

### 1. Prepare Input data
Place the csv file of your proportion file in the `./vaginal_microbiome/input/` folder.
Caveats: 

1. The first column of the Merged Proportion File is in the order of taxa, diversity, observed, and microbiome list.
2. From the second column of the Merged Proportion File, the sample name, the diversity value of the sample, the number of species in the sample, and the relative abundance value of each microbiome should be listed in the order.

### 2. Run vaginal_percentile_rank
To run vaginal_percentile_rank,
 
Run the command below:

    python ./vaginal_microbiome/vaginal_percentile_rank.py {path_exp}
    ### ex) python vaginal_percentile_rank.py "/home/kbkim/vaginal_microbiome/input/experiment_result_abundance.csv"

When vaginal_percentile_rank is executed as above, the file `beneficial.xlsx`, `harmful.xlsx` and `vaginal_percentile_rank.xlsx` will be created or modified in the `./vaginal_microbiome/output/` folder.



