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

# size of each data set
core         = 500
additional   = 500
target       = 1500

# List tasks to include
levels = ( primary   =  "CCSDT",
           secondary = ["PBE"] )

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
          output_path  = path*"/data/calib_",                 # initial path and string for output file names
          figure_path  = path*"/figures/calib_",              # initial path and string for figure file names
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




# set up 
for seed in my_seeds

    
    # set up
    rm( files.output_path*"target"*string(target)*"_seed"*string(seed)*".csv", force=true )
    GP_settings, Δ_settings, guide, all_params = set_up( ; files, levels, seed, target, 
                                                           rcut=4, sigma=0.4, nmax=8, lmax=8, species=["C","N","O","H"],
                                                           custom_heading=[:index,:tasks,:mu,:sigma,:truth_ccsdt, :truth_pbe] )

    # read in features and data
    Random.seed!(seed);
    data    = CSV.read( GP_settings.files.data, DataFrame  )
    remove  = CSV.read( GP_settings.files.opt_inds, DataFrame)[!,1]
    rows    = 1:core+additional+target
    X,Y     = get_XY( GP_settings, data, rows, GP_settings.tasks; differences=false  )


    # Assign indices
    Random.seed!(seed);
    T = sample( rows, Int64(ceil(target)), replace=false  )
    C = sample( setdiff(rows,T), Int64(ceil(core)), replace=false  )
    A = setdiff(rows,vcat(C,T))

    
    # Single task ###########################################

    # set parameters
    test_levels = ( primary   =  levels.primary,
                    secondary =  [] )
    params = set_multitask_params( ( p=   GP_settings.prior_mean( Y[C,guide.p]), s= [] ), 
                                     set_multitask_params( all_params, test_levels ) )

    # infer
    Random.seed!(seed);
    μ,σ,dt = infer( GP_settings.model, params, GPPPInput( :P, ColVecs( X[:,C] )), Y[C,guide.p], GPPPInput( :P, ColVecs( X[:,T] ) ))
    
    # record
    addToFile( DataFrame( hcat( data[rows[T],:index], ones(length(T)), μ, σ, Y[T,guide.p], Y[T,guide.s[1]]   ), :auto  ), guide.filename)


    # Multitask #############################################

    # set parameters
    Random.seed!(seed);
    test_levels = levels
    params = set_multitask_params( all_params, test_levels )
    params = set_multitask_params( ( p=   GP_settings.prior_mean( Y[C,guide.p]), 
                                     s= [ GP_settings.prior_mean( Y[vcat(C,A),guide.s[1]]) -
                                         params.ρ[1]*GP_settings.prior_mean(data[C,guide.p]) ] ),
                                     params )


    # infer
    Random.seed!(seed);
    μ,σ,dt = infer( GP_settings.model, 
                    params, 
                    BlockData( GPPPInput( :P,    ColVecs( X[:,C] ) ),
                               GPPPInput( :S1,   ColVecs( X[:,vcat(C,A)])  ) ),
                    vcat( Y[C,guide.p], Y[vcat(C,A),guide.s[1]]  ),   
                    GPPPInput( :P, ColVecs( X[:,T] )  ) )

    # record
    addToFile( DataFrame( hcat( data[rows[T],:index], 2*ones(length(T)), μ, σ, Y[T,guide.p], Y[T,guide.s[1]]   ), :auto  ), guide.filename)


end

