using LinearAlgebra
using Statistics
using Distributions
using Random
using CSV
using DataFrames
using Stheno
using Combinatorics
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

optAlg  = LBFGS()
grad    = Optimization.AutoZygote()
noise   = 0.015
seeds   = collect(1:100)

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
          output       = path*"/data/test.csv",#path*"/data/water_all_params.csv", # full set of all hyperparameters obtained from optimization
          params       = path*"/data/water_params.csv",     # final parameters for inference
          opt_inds     = path*"/data/trimer_inds.csv",      # indices in data that are set aside for optimization
         )

# =========================================================
# Read Data
# =========================================================

# optimization indices
inds     = CSV.read( files.opt_inds, DataFrame)[:,1]

# retrieve data
data   = CSV.read( files.data,   DataFrame )[inds, :]

# construct features
options  = ( globalize=average_feature, positions=getNonHydrogenPositions )
desc     = dsd.SOAP( species=["O","H"], periodic=false, rcut=10, sigma=0.4, nmax=8, lmax=8 )
ats      = [ ase.read( files.xyz, index=string(i) ) for i in data[:,:index] ]
X        = get_features( desc, ats ; settings=options )


# make parameter files
rm( files.output, force=true )
makeFile( ["task","variance","variance_mle","lengthscale","rho","seed","objective_mle","objective_mae"], files.output )


# =========================================================
# Record Optimization Results for different shufflings of the data
# =========================================================

# local seed vector
if length(ARGS) > 0
    my_task_id   = parse(Int64, ARGS[1] ) + 1
    num_tasks    = parse(Int64, ARGS[2] )
    my_seeds     = seeds[ my_task_id:num_tasks:length(seeds) ]
else
    my_seeds     = seeds
end


tasks = []
for seed in my_seeds

    Random.seed!(seed)
    shuffled = sample(1:length(inds), length(inds), replace=false)
    #shuffled = sample(1:length(inds), 24, replace=false)

    # primary only
    results = mle_mae( X[:,shuffled], 
                       data[shuffled,levels.primary], 
                       noise, 
                       ( task=levels.primary, rho=1., seed=seed, name=files.output );
                       optAlg, 
                       grad )
    push!( tasks, levels.primary )

    # Secondary Tasks
    for s in levels.secondary
        ρ       = cor(  data[:,levels.primary], data[:,s]  )
        results = mle_mae(   X[:,shuffled],
                             data[shuffled,s] - ρ*data[shuffled,levels.primary],
                             noise,
                           ( task=s, rho=ρ, seed=seed, name=files.output );
                             optAlg,
                             grad )
        push!( tasks, s )

    end

    # Deltas: between primary and secondary
    for s in levels.secondary
        results = mle_mae(   X[:,shuffled],
                             data[shuffled,levels.primary] - data[shuffled,s],
                             noise,
                           ( task=levels.primary*"_"*s, rho=1., seed=seed, name=files.output );
                             optAlg,
                             grad )
        push!( tasks, levels.primary*"_"*s )

        if results.converged
            for z in filter( x->x!=s, levels.secondary )
                ρ       = cor(data[:,levels.primary] - data[:,s], data[:,:CCSDT] - data[:,z]  )
                addToFile( DataFrame( task=[ levels.primary*"_"*s*"_"*z ],
                                      variance=results.v,
                                      variance_mle=results.v_mle,
                                      lengthscale=results.ℓ,
                                      rho=ρ,
                                      seed=seed,
                                      objective_mle=results.mle,
                                      objective_mae=results.mae   ), files.output  )
                push!( tasks, levels.primary*"_"*s*"_"*z )
            end
        end
    end



    # Deltas: between secondary tasks
    for (a,b) in collect( combinations(1:length(levels.secondary),2) )

        s = levels.secondary[a]
        z = levels.secondary[b]

        results = mle_mae(   X[:,shuffled],
                             data[shuffled,s] - data[shuffled,z],
                             noise,
                           ( task=s*"_"*z, rho=1., seed=seed, name=files.output );
                             optAlg,
                             grad )
        push!( tasks, s*"_"*z )

        if results.converged
            addToFile( DataFrame( task=[z*"_"*s],
                                  variance=results.v,
                                  variance_mle=results.v_mle,
                                  lengthscale=results.ℓ,
                                  rho=1.,
                                  seed=seed,
                                  objective_mle=results.mle,
                                  objective_mae=results.mae   ), files.output  )
            push!( tasks, z*"_"*s )

            ρ       = cor( data[:,s] - data[:,z], data[:,levels.primary] - data[:,z]  )
            addToFile( DataFrame( task=[s*"_"*z*"_"*levels.primary],
                                  variance=results.v,
                                  variance_mle=results.v_mle,
                                  lengthscale=results.ℓ,
                                  rho=ρ,
                                  seed=seed,
                                  objective_mle=results.mle,
                                  objective_mae=results.mae   ), files.output  )
            push!( tasks, s*"_"*z*"_"*levels.primary )

            ρ       = cor( data[:,z] - data[:,s], data[:,levels.primary] - data[:,s]  )
            addToFile( DataFrame( task=[z*"_"*s*"_"*levels.primary],
                                  variance=results.v,
                                  variance_mle=results.v_mle,
                                  lengthscale=results.ℓ,
                                  rho=ρ,
                                  seed=seed,
                                  objective_mle=results.mle,
                                  objective_mae=results.mae   ), files.output  )
            push!( tasks, z*"_"*s*"_"*levels.primary )
        end
    end


end


# =========================================================
# Select median of variance parameters
# =========================================================

results = CSV.read( files.output, DataFrame )

rm( files.params, force=true )
makeFile(["task","variance","lengthscale","rho"], files.params )

for t in tasks
    t_params = multifilter( results, [:task], [(t,)] )
    v        = median(t_params[!,:variance])
    ℓ        = median(t_params[!,:lengthscale])
    ρ        = median(t_params[!,:rho])
    addToFile( DataFrame( task=t, variance=v, lengthscale=ℓ, rho=ρ ), files.params )
end
