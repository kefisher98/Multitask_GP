## Contents

This directory contains the data and code necessary for reproducing the results in *TBD*.

The paper considers two examples:

1. Prediction of the 3 body energy of water trimers

2. Prediction of the ionization potential of small organic molecules

For each example, we fit inference models to make more efficient predictions of our quantities of interest. In this section, we will summarize the contents of the directory and the role they play in this task. The contents of the data subdirectory are described within their own README file. In the following sections, we will provide more detail on how to run each example. Within this directory, we have

* Scripts to train and test inference models

  - run_water.jl

  - run_organic.jl

* Scripts to produce plots

  - make_plots.jl

* Scripts to optimize kernel hyperparameters

  - optimize_water.jl

  - optimize_organic.jl

* Code to support inference, optimization, and visualization

  - common.jl

  - plotSupport.jl

* Scripts to generate training data as well as cost models of training data generation

  - calc_water_dft.py

  - calc_water_cost.py

  - calc_organic_cost.py

* Code to support quantum chemistry calculations

  - q_chem.py

**NOTE:** The results of each of the following sections are saved in this repository, so it is not necessary to run the scripts in any particular order. Furthermore, if the Julia scripts are run as they are, their existing results will be replaced with new results of the same name.

## Training Inference Models and Reproducing Figures

The code in this section requires an installation of
[Julia 1.9.4](https://julialang.org/downloads/#current_stable_release),
of [DScribe](https://singroup.github.io/dscribe/1.1.x/install.html)
and of [ASE](https://wiki.fysik.dtu.dk/ase/).

With this setup, we train inference models and reproduce the plots of the paper by running: 

```bash
julia --project=@. -e "import Pkg; Pkg.instantiate()"  # Install dependencies
julia run_water.jl # generate data from water trimer 3 body energy example
julia run_organic.jl      # generate data from organic molecules ionization potential example
julia make_plots.jl  # Generate plots
```

**NOTE:** Tests are designed to demonstrate performance of inference methods for a range of datasets. As a consequence, this implementation is not the most efficient version of these methods. In particular, training with secondary data for target molecular systems trains a new GP for each target to isolate its impact. These redundant training steps would not be necessary in practice.


## Optimizing Hyperparameters

The code in this section requires an installation of
[Julia 1.9.4](https://julialang.org/downloads/#current_stable_release),
of [DScribe](https://singroup.github.io/dscribe/1.1.x/install.html)
and of [ASE](https://wiki.fysik.dtu.dk/ase/).

Optimized hyperparameters corresponding to each quantum chemisty method considered in this paper are saved in the data subdirectory and are ready for use in inference models. To regenerate these results, run

```bash
julia --project=@. -e "import Pkg; Pkg.instantiate()"  # Install dependencies
julia optimize_water.jl # create .CSV files with optimized kernel hyperparameters for the water trimer 3 body energy example
julia optimize_organic.jl      # create .CSV files with optimized kernel hyperparameters for the organic molecules ionization potential example
```

**NOTE:**  The optimization procedures are more extensive than strictly necessary. In particular, the optimization procedure for the water example is repeated for many random initial conditions to produce an idea of the distribution of estimated hyperparameters.



## Reproducing Training Data

The code in this section requires an installation of
[Python 3.9.16](https://www.python.org/downloads/) and of [ASE](https://wiki.fysik.dtu.dk/ase/).

The following code will reproduce density functional theory (DFT) training data for the water example. The final two scripts run 10 DFT and CCSD(T) computations for each example to produce a cost model.

```bash
pip3 install -r requirements.txt  # Install dependencies
python calc_water_dft.py # compute 3 body energy of water trimers in the training set using DFT
python calc_water_cost.py # compute 3 body energy for 10 representative trimers using DFT and CCSD(T) to estimate computation cost
python calc_organic_cost.py # compute energy of 10 representative organic molecules using DFT and CCSD(T) to estimate computation cost
```

The water example also uses CCSD(T) training data from [https://github.com/jmbowma/q-AQUA](https://github.com/jmbowma/q-AQUA). 

Training data for the organic molecule example was generated as described in 

C. Duan, S. Chen, M. G. Taylor, F. Liu, and H. J. Kulik, “Machine learning to tame divergent density functional approximations: a new path to consensus materials design principles,” Chem. Sci. 12, 13021–13036 (2021). [https://doi.org/10.1039/D1SC03701C ](https://doi.org/10.1039/D1SC03701C). 




## Parallelization

Many procedures in these scripts are repeated with multpile random initializations or are otherwise independent. For improved efficiency, the code can be run with multiple CPUs:

```bash
julia run_water.jl <Task ID> <Number of Tasks>
```

**All** of the above scripts can support these additional inputs in order to divide their tasks among multiple CPUs. 
