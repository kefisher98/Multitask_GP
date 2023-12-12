# import Julia libraries
using LinearAlgebra
using Random
using StatsBase
using Stheno
using CSV
using DataFrames
using KernelFunctions
using Distributions
using Plots
using Plots.Measures
using Combinatorics

# import Python libraries
using PyCall
d     = pyimport("dscribe")
dsd   = pyimport("dscribe.descriptors")
ase   = pyimport("ase.io")
atoms = pyimport("ase")


# import local files
path = @__DIR__
include( path*"/common.jl")
include( path*"/plotSupport.jl")

# ========================================================
# Define data sets to iterate over
# ========================================================

# List of random seeds where we initialize tests
seeds        = [21,22,23,24,25,26]

# Number of target data points
target       = 320

# List of sizes of core data set to consider
core         = [5,10,20,40,80,160,320]

# To perform tests of the single task method with larger training sets,
#   add training set sizes to this array
single_extra = [640,1280]



# List tasks to include
levels = ( primary   =  "CCSDT",
           secondary = ["PBE0", "PBE","PBE0_DH","BLYP"] )

# Subsets of Secondary levels (this is how we ask for scrambled delta tests)
multitask_subsets = [ [1], [2],  [1,2], [1,3,4,2]  ]
delta_subsets     = [ [1], [2],  [1,2], [2,1], [1,2,3], [1,3,2], [3,2,1], [1,2,3,4], [1,3,4,2] ]


# ========================================================
# Location of Input Files
# ========================================================

# data, xyz, params, augment data, optimization inds

files = ( data         = path*"/data/organic_ip_data.csv",    # prediction data for all molecules and levels of theory
          xyz          = path*"/data/organic_molecules.xyz",  # extended xyz file containing all molecular systems
          params       = path*"/data/organic_params.csv",     # GP hyperparameters (obtained via optimize_water.jl)
          opt_inds     = path*"/data/organic_inds.csv",       # indices in data that are set aside for optimization
          output_path  = path*"/data/organic_",               # initial path and string for output file names
          figure_path  = path*"/figures/organic_",             # initial path and string for figure file names
         )

# ========================================================
# Inference
# ========================================================

# local seed vector
if length(ARGS) > 0
    my_task_id   = parse(Int64, ARGS[1] ) + 1
    num_tasks    = parse(Int64, ARGS[2] )
    my_seeds     = seeds[ my_task_id:num_tasks:length(seeds) ]
else
    my_seeds     = seeds
end


Norganic = maximum( vcat( 7*core, single_extra ) ) + target


# set up 
for seed in my_seeds

    
    # set up
    rm( files.output_path*"target"*string(target)*"_seed"*string(seed)*".csv", force=true )
    GP_settings, Δ_settings, guide, all_params = set_up( ; files, levels, seed, target, 
                                                           rcut=4, sigma=0.4, nmax=8, lmax=8, species=["C","N","O","H"] )

    # read in features and data
    Random.seed!(seed);
    X, Y  = get_XY( GP_settings, Norganic )



    # inference
    for set in delta_subsets
        Random.seed!(seed);
        test_levels = ( primary   =  levels.primary,
                        secondary =  levels.secondary[set] )
        training_set_iteration(  Δ_settings, 
                                 set_delta_params( all_params, test_levels ),
                                 Y, X,
                                 merge( guide,  (s=test_levels.secondary,)   ),
                                 core, target,
                                [delta];
                                C_sets=[1.],
                                task_fractions=[0.5,1],
                                A_sets=[0.,1.,2.,3.,4.,5.,6.] )
    end

    for set in multitask_subsets
        Random.seed!(seed);
        test_levels = ( primary   =  levels.primary,
                        secondary =  levels.secondary[set] )
        training_set_iteration(  GP_settings, 
                                 set_multitask_params( all_params, test_levels ), 
                                 Y, X,
                                 merge( guide,  (s=test_levels.secondary,)   ),
                                 core, target,
                                [multitask, multitask_with_target];
                                 task_fractions=[0,0.5,1],
                                 A_sets=[0.,1.,2.,3.,4.,5.,6.])
    end

    Random.seed!(seed);
    test_levels = ( primary   =  levels.primary,
                    secondary =  [] )

    training_set_iteration(  GP_settings,
                             set_multitask_params( all_params, test_levels ),
                             Y, X,
                             guide,
                             core, target,
                            [single];
                             C_sets=[0.],
                             A_sets=[0.] )

    training_set_iteration(  GP_settings,
                             set_multitask_params( all_params, test_levels ),
                             Y, X,
                             guide,
                             single_extra, target,
                            [single];
                             C_sets=[0.],
                             A_sets=[0.] )

end

