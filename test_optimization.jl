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

# List tasks to include
levels = ( primary   =  "CCSDT",
           secondary = ["SCAN"] )


# ========================================================
# Location of Input Files
# ========================================================

# data, xyz, params, augment data, optimization inds

files = ( data         = path*"/data/trimer_dft.csv",       # prediction data for all molecules and levels of theory
          augment_data = path*"/data/trimer_ccsdt.csv",     # additional primary task prediction data
          xyz          = path*"/data/3b_data.xyz",          # extended xyz file containing all molecular systems
          output       = path*"/data/test_optimization.csv",# full set of all hyperparameters obtained from optimization
          params       = path*"/data/water_params.csv",     # final parameters for inference
          opt_inds     = path*"/data/trimer_inds.csv",      # indices in data that are set aside for optimization
         )

# =========================================================
# Read Data
# =========================================================

# optimization indices
inds     = CSV.read( files.opt_inds, DataFrame)[:,1]
test_n   = [10,45,80,115,150]

# retrieve optimization data
data   = CSV.read( files.data,   DataFrame )[inds, :]

# construct features
options  = ( globalize=average_feature, positions=getNonHydrogenPositions )
desc     = dsd.SOAP( species=["O","H"], periodic=false, rcut=10, sigma=0.4, nmax=8, lmax=8 )
ats      = [ ase.read( files.xyz, index=string(i) ) for i in data[:,:index] ]
X        = get_features( desc, ats ; settings=options )


# make parameter files
rm( files.output, force=true )
makeFile( ["task","variance","variance_mle","lengthscale","rho","seed","objective_mle","objective_mae","time_mle","time_mae","N"], files.output )


# =========================================================
# Testing Functions
# =========================================================


# time hyperparameter optimization
function time_opt( X, Y, test_n ; limit=20, ρ=1, task  )

    time_mle = 0 
    time_mae = 0

    for n in test_n
        found    = 0
        tried    = 0
        while found < limit 

            Random.seed!( tried+1 )
            shuffled = sample(1:length(inds), n, replace=false)

            # lengthscale and variance estimation: maximize likelihood
            t_mle     = @elapsed begin
                mle   = max_likelihood( X[:,shuffled], Y[shuffled], noise; θ_init=[1.], lower=[1e-10], upper=[100.], optAlg, grad )
            end
            time_mle += t_mle

            # proceed if MLE converged
            if mle.converged
                # variance estimation: minimize mean absolute error
                t_mae     = @elapsed begin
                    mae   = min_mae_var( X[:,shuffled], Y[shuffled], mle.ℓ, noise; optAlg, grad, maxtime=30  )
                end
                time_mae += t_mae

                # proceed if MAE converged
                if mae.converged && mae.v[1] < 1e3
                    addToFile( DataFrame( task=task,
                                          variance=mae.v,
                                          variance_mle=mle.v,
                                          lengthscale=mle.ℓ,
                                          rho=ρ,
                                          seed=tried+1,
                                          objective_mle=mle.objective,
                                          objective_mae=mae.objective,
                                          time_mle=time_mle,
                                          time_mae=time_mae,
                                          N=n ), files.output  )
            
                    time_mae = 0
                    found   += 1
                end
                time_mle = 0
            end
            tried += 1
        end
    end
end

# find MAE achieved by hyperparamter estimates
function opt_eval( X, Y, inds, results )

    n   = size(results,1)
    mae = zeros(n)
    mle = zeros(n)

    for r in inds

        Random.seed!( r+1 )
        train   = sample( 1:640, 320, replace=false)
        test    = setdiff(1:640, train)
        params  = (X=X, Y=Y, μ=mean(Y[train]), ℓ=results[r,:lengthscale], noise=noise, train=train, test=test)
        mle[r]  = mae_var( results[r,:variance_mle], params)
        mae[r]  = mae_var( results[r,:variance], params)

    end
    return mle, mae
end


# =========================================================
# Time optimization (run procedure 20 times)
# =========================================================

limit = 20
ρ     = cor(  data[:,levels.primary], data[:,levels.secondary[1]]  )

time_opt( X, data[:,levels.primary] - ρ*data[:,levels.secondary[1]],  test_n ; limit, ρ, task=levels.secondary[1]  )
time_opt( X, data[:,levels.primary],                                  test_n ; limit, ρ=1, task=levels.primary )


# =========================================================
# Test optimization results
# =========================================================


# retrieve data outside of optimization set
rows     = sample(setdiff(1:2992,inds), 640)
data     = CSV.read( files.data,   DataFrame )[rows, :]

# construct features
options  = ( globalize=average_feature, positions=getNonHydrogenPositions )
desc     = dsd.SOAP( species=["O","H"], periodic=false, rcut=10, sigma=0.4, nmax=8, lmax=8 )
ats      = [ ase.read( files.xyz, index=string(i) ) for i in data[:,:index] ]
X        = get_features( desc, ats ; settings=options )


# read in results
results  = CSV.read( files.output, DataFrame  )

# obtain optimization mae
mle_s, mae_s = opt_eval( X, data[:,levels.primary] - ρ*data[:,levels.secondary[1]], collect(161:180), results    )
mle_p, mae_p = opt_eval( X, data[:,levels.primary], collect(101:160), results     )

# augment dataframe
results[:,:mle_mae] = vcat(mle_s,mle_p)
results[:,:mae_mae] = vcat(mae_s,mae_p)
CSV.write( files.output, results )


