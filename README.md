Microbiome Lab
---
### Dr. Roberts and Dr. Zhou
***Also Mustafa and Avery***

This repository has all the data and scripts for the figures and analysis for Dr. Robert's microbiome research from the summer of 2025. 

The main home directory, the one with this README file, has a bunch of short *python scripts* and *csv files* labelled from a to f that were involved in taking microbiome data from two separate time periods
and combining them, then cleaning the data. 

It was labelled with letters for personal organization, you can ignore that. 

---


**f1-cleaned-ready-general.csv** is the current best file if you want to analyze microbiome data, its preprocessed, turned into relative values, low abundance bacterium have been dropped or aggregated. 

**f2-non-aggregate.csv** is the same thing, without the low abundance aggregation. Some metrics in microbiome research, like richness, are ruined by aggregation so you want to have that too. 

---

Folder **Middle Ear** contains the matlab script **b1-middle-ear.m** which uses files in that directory to generate middle ear analysis.

Folder **Sinus** contains the script **b1-sinus-analysis.m** doing this for sinus. 

And Folder **All** contains **a2-analysis.m** doing this for some more broad figures. 

---

**b2-demographics.csv** in folder **Middle Ear** contains most recent anonomized demographics info which may help with future data stratification. 

and finally,

**Analysis.md** has a preview of all figures with current data analyzed in Summer 2025. 
