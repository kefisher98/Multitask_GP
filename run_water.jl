
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
single_extra = [640,1280,2560]

# List tasks to include
levels = ( primary   =  "CCSDT",
           secondary = ["PBE","SCAN"] )


# ========================================================
# Location of Input Files
# ========================================================

# data, xyz, params, augment data, optimization inds

files = ( data         = path*"/data/trimer_dft.csv",       # prediction data for all molecules and levels of theory
          augment_data = path*"/data/trimer_ccsdt.csv",     # additional primary task prediction data
          xyz          = path*"/data/3b_data.xyz",          # extended xyz file containing all molecular systems
          params       = path*"/data/water_params.csv",     # GP hyperparameters (obtained via optimize_water.jl)
          opt_inds     = path*"/data/trimer_inds.csv",      # indices in data that are set aside for optimization
          output_path  = path*"/data/water_",               # initial path and string for output file names
          figure_path  = path*"/figures/water_",             # initial path and string for figure file names
         )

# ========================================================
# Inference 
# ========================================================

# number of data points to read from CSV file
Ntrimer = 7*maximum(core) + target
Nccsd   = maximum(single_extra) + target

# local seed vector
if length(ARGS) > 0
    my_task_id   = parse(Int64, ARGS[1] ) + 1
    num_tasks    = parse(Int64, ARGS[2] )
    my_seeds     = seeds[ my_task_id:num_tasks:length(seeds) ]
else 
    my_seeds     = seeds
end




for seed in my_seeds
    
    # set up
    rm( files.output_path*"target"*string(target)*"_seed"*string(seed)*".csv", force=true )
    GP_settings, Δ_settings, guide, all_params = set_up( ; files, levels, seed, target, species=["O","H"] )

    # read in features and data
    Random.seed!(seed);
    X, Y  = get_XY( GP_settings, Ntrimer )

    # Multitask Inference: One Secondary Task =====================================================

    for i in 1:length(levels.secondary)
        
        Random.seed!(seed);
        test_levels = ( primary   =  levels.primary,
                        secondary = [levels.secondary[i]] )
        training_set_iteration( GP_settings,
                                set_multitask_params( all_params, test_levels ),
                                Y, X,
                                merge( guide,  (s=test_levels.secondary,)   ),
                                core, target,
                               [multitask, multitask_with_target] )
    end
    
    # Multitask Inference: Two secondary Tasks ======================================================
    Random.seed!(seed);
    test_levels = levels
    training_set_iteration( GP_settings,
                            set_multitask_params( all_params, test_levels ),
                            Y, X,
                            guide, core, target,
                           [multitask, multitask_with_target] )

    test_levels = ( primary   = levels.primary,
                    secondary = reverse(levels.secondary) )
    training_set_iteration( GP_settings, 
                            set_multitask_params( all_params, test_levels ),
                            Y, X,
                            merge( guide, (s=reverse(levels.secondary),) ), 
                            core, target,
                           [multitask_with_target] )


    # Delta Method: 1 level =====================================================
    Random.seed!(seed);
    test_levels = ( primary   = "CCSDT_SCAN",
                    secondary = [] )
    training_set_iteration(  Δ_settings, 
                             set_delta_params(all_params, test_levels),
                             Y, X,
                             merge( guide, ( p="CCSDT", s=["SCAN"] )),
                             core, target,
                            [delta];
                             C_sets=[1.], A_sets=[0.] )


    Random.seed!(seed);
    test_levels = ( primary   = "CCSDT_PBE",
                    secondary = [] )
    training_set_iteration(  Δ_settings, 
                             set_delta_params(all_params, test_levels),
                             Y, X,
                             merge( guide, ( p="CCSDT", s=["PBE"] )),
                             core, target,
                            [delta];
                             C_sets=[1.], A_sets=[0.] )


    # Delta Method: 2 level =====================================================

    Random.seed!(seed);
    test_levels = ( primary   =  "CCSDT_SCAN",
                    secondary = ["SCAN_PBE"] )
    training_set_iteration(  Δ_settings, 
                             set_delta_params(all_params, test_levels),
                             Y, X,
                             merge( guide, ( p="CCSDT", s=["SCAN","PBE"] )),
                             core, target,
                            [delta];
                            C_sets=[1.] )


    Random.seed!(seed);
    test_levels = ( primary   =  "CCSDT_PBE",
                    secondary = ["PBE_SCAN"] )
    training_set_iteration(  Δ_settings, 
                             set_delta_params(all_params, test_levels),
                             Y, X,
                             merge( guide, ( p="CCSDT", s=["PBE", "SCAN"] )),
                             core, target,
                            [delta];
                             C_sets=[1.] )

    # Multitask of Deltas =======================================================

    Random.seed!(seed);
    test_levels = ( primary   =  "CCSDT_SCAN",
                    secondary = ["PBE_SCAN_CCSDT"] )
    training_set_iteration(  GP_settings, 
                             set_multitask_params(all_params, test_levels),
                             Y, X,
                             merge( guide, ( p="CCSDT_SCAN", s=["PBE_SCAN"] )),
                             core, target,
                            [multitask, multitask_with_target] )


    Random.seed!(seed);
    test_levels = ( primary   =  "CCSDT_PBE",
                    secondary = ["SCAN_PBE_CCSDT"] )
    training_set_iteration(  GP_settings, 
                             set_multitask_params(all_params, test_levels),
                             Y, X,
                             merge( guide, ( p="CCSDT_PBE", s=["SCAN_PBE"] )),
                             core, target,
                            [multitask, multitask_with_target] )

    Random.seed!(seed);
    test_levels = ( primary   =  "CCSDT_SCAN",
                    secondary = ["CCSDT_PBE_SCAN"] )
    training_set_iteration(  GP_settings, 
                             set_multitask_params(all_params, test_levels),
                             Y, X,
                             merge( guide, ( p="CCSDT_SCAN", s=["CCSDT_PBE"] )),
                             core, target,
                            [multitask] )


    Random.seed!(seed);
    test_levels = ( primary   =  "CCSDT_PBE",
                    secondary = ["CCSDT_SCAN_PBE"] )
    training_set_iteration(  GP_settings, 
                             set_multitask_params(all_params, test_levels),
                             Y, X,
                             merge( guide, ( p="CCSDT_PBE", s=["CCSDT_SCAN"] )),
                             core, target,
                            [multitask] )

    # CCSD(T) only =============================================================
    # single
    Random.seed!(seed);
    test_levels = ( primary   = levels.primary,
                    secondary = [] )
    training_set_iteration(  GP_settings,
                             set_multitask_params( all_params, test_levels ),
                             Y, X, guide,
                             core, target,
                            [single] ;
                             C_sets=[0.],
                             A_sets=[0.]  )



    # get data
    Random.seed!(seed);
    X, Y  = augment_XY( GP_settings, [Ntrimer, Nccsd] )

    # single, larger tests
    Random.seed!(seed);
    test_levels = ( primary   = levels.primary,
                    secondary = [] )
    training_set_iteration(  GP_settings,
                             set_multitask_params( all_params, test_levels ),
                             Y, X, guide,
                             single_extra, target,
                            [single] ;
                             C_sets=[0.],
                             A_sets=[0.]  )



end

