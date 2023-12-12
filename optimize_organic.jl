using LinearAlgebra
using Statistics
using Distributions
using Random
using CSV
using DataFrames
using Combinatorics
using Stheno

using Zygote: gradient
using Optimization
using OptimizationOptimJL
using OptimizationEvolutionary

# Python
using PyCall
d     = pyimport("dscribe")
dsd   = pyimport("dscribe.descriptors")
ase   = pyimport("ase.io")
atoms = pyimport("ase")

# local files
path  = @__DIR__
include( path*"/common.jl")

# =========================================================
# Optimization Options
# =========================================================

optAlg  = CMAES(μ=40,λ=80)
grad    = Optimization.AutoZygote()
noise   = 0.015


# List tasks to include
levels = ( primary   =  "CCSDT",
           secondary = ["PBE0", "PBE","BLYP", "PBE0_DH"] )


# ========================================================
# Location of Input Files
# ========================================================

# data, xyz, params, augment data, optimization inds

files = ( data         = path*"/data/organic_ip_data.csv",    # prediction data for all molecules and levels of theory
          xyz          = path*"/data/organic_molecules.xyz",  # extended xyz file containing all molecular systems
          params       = path*"/data/organic_params.csv",     # final hyperparameters for inference
          opt_inds     = path*"/data/organic_inds.csv",       # indices in data that are set aside for optimization
         )


# =========================================================
# Read Data
# =========================================================

# optimization indices
inds     = CSV.read( files.opt_inds, DataFrame)[1:50,1]

# retrieve data
data   = CSV.read( files.data,   DataFrame )[inds, :]

# construct features
options  = ( globalize=average_feature, positions=getNonHydrogenPositions )
desc     = dsd.SOAP( species=["C","N","O","H"], periodic=false, rcut=4, sigma=0.4, nmax=8, lmax=8 )
ats      = [ ase.read( files.xyz, index=string(i) ) for i in data[:,:index] ]
X        = get_features( desc, ats ; settings=options )

println(size(X))

# make parameter files
rm( files.params, force=true )
makeFile( ["task","variance","lengthscale","rho"], files.params )



# =========================================================
# Record Optimization Results for MLE
# =========================================================

# local seed vector
c = collect( combinations(1:length(levels.secondary)+1,2) )
if length(ARGS) > 0
    my_task_id      = parse(Int64, ARGS[1] ) + 1
    num_tasks       = parse(Int64, ARGS[2] )
    my_secondary    = levels.secondary[ my_task_id:num_tasks:length(levels.secondary) ]
    my_combinations = c[ my_task_id:num_tasks:length(c) ]
else
    my_secondary    = levels.secondary
    my_combinations = c
end


# primary only
if levels.secondary[1] in my_secondary
    results = max_likelihood( X,
                              data[:,levels.primary],
                              noise;
                              optAlg,
                              grad )
    if results.converged
        addToFile( DataFrame( task=levels.primary,
                              variance=results.v,
                              lengthscale=results.ℓ,
                              rho=1. ), files.params  )
    end
end

# Secondary Tasks
for s in my_secondary
    ρ       = cor(  data[:,levels.primary], data[:,s]  )
    results = max_likelihood(   X,
                                data[:,s] - ρ*data[:,levels.primary],
                                noise;
                                optAlg,
                                grad )
    if results.converged
        addToFile( DataFrame( task=s,
                              variance=results.v,
                              lengthscale=results.ℓ,
                              rho=ρ ), files.params  )
    end
end
    

# Deltas
tasks = vcat( levels.primary, levels.secondary  )
for (a,b) in my_combinations
    results = max_likelihood(   X,
                                data[:,tasks[a]] - data[:,tasks[b]],
                                noise;
                                optAlg,
                                grad )
    if results.converged
        for n in [ tasks[a]*"_"*tasks[b], tasks[b]*"_"*tasks[a]  ]
            addToFile( DataFrame( task=n,
                                  variance=results.v,
                                  lengthscale=results.ℓ,
                                  rho=1. ), files.params  )
        end
    end
end




