# Microbiome
UConn Otolaryngology Microbiome Studies

This folder is the repository to contain all the files necessary for data analysis for cleaned 
bacterial specimen for studies pertaining to the microbiome in the shared labs of Dr. Daniel
Roberts and Dr. Yanjiao Zhou. 

`1-sample-data.csv` lists unique patient IDs as they are recorded in clinic during sample intake. 
It also holds the unique IDs as they are processed during bacterial sequencing. It then lists 
the status of the sequenced samples. Some are still raw, having not been preprocessed for further 
analysis, while some are clean, and ready to be analyzed. 

`2-processed-data.csv` contains the set of all such clean files, along with the presence of each 
bacterial group. 

`3-plots.Rmd` is the final file containing the code used to process the files. 
