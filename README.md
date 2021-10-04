# LCD Data Repository and Analysis Routines
Data repository and analysis routines used to study prion-like low-complexity domains as described in https://www.biorxiv.org/content/10.1101/2021.01.01.425046v1

This directory contains multiple data files in the Data-Files directory:

SAXS_Data.xlsx - small-angle X-ray scattering data for variants described in the study.

CD_Data.xlsx - circular dichroism data for variants described in the study.

Phase_Separation_Data.xlsx - dilute and dense phase concentrations of phase separated systems described in the study.

pH_Dependence_Measurements.xlsx - pH-dependent csat measurements of the +7K+12D variant.

vant_Hoff_Data.xlsx - thermodynamic values inferred from the csat measurements of the variants used in this study.

hnRNPA1_Homologs.fasta - hnRNPA1 homolog sequences used in the study.

hnRNPA1_Homologs_IDRs.fasta - hnRNPA1 homolog IDR sequences used in the study.

FUS_Fam_Homologs.fasta - FUS family homolog sequences used in the study.

FUS_Fam_Homologs_IDRs.fasta - FUS family homolog IDR sequences used in the study.

hnRNPA1_Homologs_IDRs_csat.xlsx - Estimated csat values of hnRNPA1 homolog IDRs.

DIC_Images_1.zip and DIC_Images_2.zip - zip files containing DIC microscopy images of condensates formed by many of the variants used in this study.

In addition, this directory contains four Python modules used for data analysis:

Rescaled_Csat.py - a module with functions to rescale csat values given a mean-field model and to plot the rescaled values. This module can also be used to parameterize a mean-field model.

vant_Hoff_Csat.py - a module with functions to perform a thermodynamic van't Hoff analysis given csat values of different proteins and to plot the linear regressions from this analysis. 

PLCD_Analyzer.py - a master module that takes advantage of the other modules to output a custom plot of rescaled csat values and a thermodynamic analysis of csat values.

data_set_example.py - a module that defines some settings and data required for PLCD_Analyzer.py to function properly. All data are taken from the preprint listed above.

In order to use these analysis methods, one should first verify that running PLCD_Analyzer.py shows a 2-panel plot of rescaled csat values and a thermodynamic analysis.

Next, one should reformat data_set_example.py or create a new data file to contain their personal csat data. Specifically, one can follow the template set in data_set_example.py or use an alternative method. For an alternative method, the modules expect a dictionary, seq_dic, with the form seq_dic[_Variant Name_]["csat"][_Temperature_] = X, where anything in italics should be replaced with an appropriate string. The modules also expect seq_dic to have other keywords, which can be discerned from the example module. In addition, the data_set module will require a list of variant names, a list of temperatures from experiments, and a list of marker types to be used with matplotlib. The latter can be taken directly from the example module.

After formatting the data_set module, one should format PLCD_Analyzer.py to contain the correct variant names to be analyzed as well as the correct temperature list to use in the mean-field model.

After accounting for the above changes, the analysis routines should be fully functional. Further adjustments can be made to customize, for example, the mean-field model, whether variants with fewer than 5 csat measurements are included in the thermodynamic analysis, and the format of the plots/data by modifying the respective modules.

Please direct any questions to Mina Farag, MD/PhD candidate in the Pappu Lab (minafarag@wustl.edu). 
