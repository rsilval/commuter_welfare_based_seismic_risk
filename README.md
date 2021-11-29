# commuter_welfare_based_seismic_risk

This repository contains code that allows partial reproduction of the results shown in the article titled: "Commuter welfare-based probabilistic seismic risk assessment of regional road networks"

The are three main codes to generate the results of the article:
  1- ModuleT_Single_Scenario.py: Code that runs the traffic model and network disruption for a seismic scenario, using LODES data. The output of this code is a nested dictionary of tuples, where the keys are origin and destination and the tuple indicates number of users and travel time.
  
  2- Postprocessing_Outputs_Regular.py: Code that uses the results of ModuleT_Single_Scenario.py to compute welfare metrics such as expected annual welfare loss disaggregated per income.
  
  3- generate_plots_paper.py: Code that generate te figures in the paper using outputs of the previous code for several seismic scenarions


Important Notes:

- To run the code it is necessary to create an output folder.
- Trrafic data information, OD matrix, has to be downloaded from: https://lehd.ces.census.gov/data/#:~:text=ca_od_main_JT05_2019.csv.gz
