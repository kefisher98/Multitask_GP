## Contents of Directory

* Training Data (files are described as they will appear after calling Multitask_GP/retrieve_data.jl)

  - 3b_data.xyz - ASE readable extended xyz file with all water trimer geometries

  - trimer_dft.csv - DFT and CCSD(T) predictions for ~3000 water trimers. Entries in Index column correspond to 3b_data.xyz

  - trimer_ccsdt.csv - CCSD(T) predictions for ~6000 water trimers. Entries in Index column correspond to 3b_data.xyz

  - organic_molecules.xyz - ASE readable extended xyz file with all small organic molecules geometries

  - organic_ip_data.csv - DFT and CCSD(T) predictions for ~3000 small organic molecules. Entries in Index column correspond to organic_molecules.xyz

* Special training data indices

  - trimer_inds.csv - rows of trimer_dft.csv set aside to be used for hyperparameter optimization

  - water_cost_inds.csv - rows of trimer_dft.csv used for the water trimer cost model

  - organic_inds.csv - rows of organic_ip_data.csv set aside to be used for hyperparameter optimization

  - organic_cost_inds.csv - rows of organic_ip_data.csv used for the organic molecules cost model

* Optimized hyperparameters

  - water_all_params.csv - all optimized kernel hyperparameters for water trimer example obtained by initializing optimization with different random seeds

  - water_params.csv - representative optimized kernel hyperparameters for water trimer example

  - organic_params.csv - optimized kernel hyperparameters for organic example

* Inference model predictions

  - water_target<number of test molecules>_seed<random seed>.csv - ouput of all water inference models
  
  - organic_target<number of test molecules>_seed<random seed>.csv - output of all organic molecules inference models

* Representative output files from psi4 calculations  (in directory, psi4_out)
  
  - <molecule index>_<quantum chemistry method>_water_cost.out - full output from psi4 for water three body energy calculations; the molecule index matches trimer_dft.csv and 3b_data.xyz
  
