# LCD_Analyzer
Analysis routines used to study prion-like low-complexity domains as described in https://www.biorxiv.org/content/10.1101/2021.01.01.425046v1

This directory contains four files:
Rescaled_Csat.py - a module with functions to rescale csat values given a mean-field model and to plot the rescaled values. This module can also be used to parameterize a mean-field model.
vant_Hoff_Csat.py - a module with functions to perform a thermodynamic van't Hoff analysis given csat values of different proteins and to plot the linear regressions from this analysis. 
PLCD_Analyzer.py - a master module that takes advantage of the other modules to output a custom plot of rescaled csat values and a thermodynamic analysis of csat values.
data_set_example.py - a module that defines some settings and data required for PLCD_Analyzer.py to function properly. All data are taken from the preprint listed above.

In order to use these analysis methods, one should first verify that running PLCD_Analyzer.py correctly shows a 2-panel plot of rescaled csat values and a thermodynamic analysis.

Next, one should reformat data_set_example.py or create a new data file to contain their personal csat data. Specifically, one can follow the template set in data_set_example.py or use an alternative method. For an alternative method, the modules expect a dictionary, seq_dic, with the form seq_dic[<Variant Name>]["csat"][<Temperature>] = X, where anything surrounded by <> should be replaced with an appropriate string. The modules also expect seq_dic to have other keywords, which can be discerned from the example module. In addition, the data_set module will require a list of variant names, a list of temperatures from experiments, and a list of marker types to be used with matplotlib. The latter can be take directly from the example module.

After formatting the data_set module, one should format PLCD_Analyzer.py to contain the correct variant names to be analyzed as well as the correct temperature list to use in the mean-field model.

After accounting for the above changes, the analysis routines should be fully functional. Further adjustments can be made to customize, for example, the mean-field model, whether variants with fewer than 5 csat measurements are included in the thermodynamic analysis, and the format of the plots/data by modifying the respective modules.

Please direct any questions to Mina Farag, MD/PhD candidate in the Pappu Lab (minafarag@wustl.edu). 
