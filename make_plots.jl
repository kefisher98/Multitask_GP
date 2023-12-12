
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


# List tasks to include
levels = ( primary   =  "CCSDT",
           secondary = ["PBE","SCAN"] )


# ========================================================
# Location of Input Files
# ========================================================

# data, xyz, params, augment data, optimization inds

files = ( water_output   = path*"/data/water_",                    # location of water output data
          water_figures  = path*"/figures/water_",                 # initial path and string for water figures
          organic_output   = path*"/data/organic_",                # location of organic output data
          organic_figures  = path*"/figures/organic_",             # initial path and string for organic figures
         )

# ========================================================
# Seeds
# ========================================================

# local seed vector
if length(ARGS) > 0
    my_task_id   = parse(Int64, ARGS[1] ) + 1
    num_tasks    = parse(Int64, ARGS[2] )
    my_seeds     = seeds[ my_task_id:num_tasks:length(seeds) ]
else
    my_seeds     = seeds
end


# ========================================================
# Water Figures
# ========================================================


# create file list
results = [  files.water_output*"target"*string(target)*"_seed"*string(s)*".csv" for s in my_seeds]

# plot comparison of delta method, multitask method, and combination
scalefontsizes(2)

for set in ["C","A","CA"]
    p2  = plot_two_tasks( results; datasets=set )
    savefig( p2,  files.water_figures*"2level_"*set*".pdf" )

    p3  = plot_three_tasks( results; datasets=set )
    savefig( p3,  files.water_figures*"3level_"*set*".pdf" )

    p3t = plot_three_tasks( results; datasets=set*"T" )
    savefig( p3t,  files.water_figures*"3level_"*set*"T.pdf" )
end



# compare cases with secondary C data to equivalent cases with out
cs  = compare_S_power( results; datasets="CA", lowest="PBE", get_cost=water_cost, include_single=false, range=[0.,0.03] )
cst = compare_S_power( results; datasets="CAT", lowest="PBE", get_cost=water_cost, include_single=false, range=[0.,0.03] )
savefig( cs,  files.water_figures*"ca_versus_a.pdf" )
savefig( cst, files.water_figures*"cat_versus_at.pdf" )


# plot comparison of delta method, multitask method, and combination
p = plot_delta(results)
savefig( p, files.water_figures*"delta_CCSDT-PBE.pdf" )

p = plot_delta(results; quantity="CCSD(T)-SCAN" )
savefig( p, files.water_figures*"delta_CCSDT-SCAN.pdf" )

scalefontsizes(0.5)


# ========================================================
# Organic Figures
# ========================================================

# create file list
results = [  files.organic_output*"target"*string(target)*"_seed"*string(s)*".csv" for s in my_seeds]

# plot comparison of delta method, multitask method, and combination
scalefontsizes(2)


# layer plots
for a in [0.,0.5,1.]
    cs_layers = plot_layers( results; title="", datasets="CA", markersize=8,
                                      get_cost=organic_cost, lowest=["PBE","PBE","PBE"],
                                      levels=[5,3,2], c_size=[5,10,20,40,80,320],
                                      t_fraction=a, c_fraction=1.)
    cst_layers = plot_layers( results; title="", datasets="CAT", markersize=8,
                                      get_cost=organic_cost, lowest=["PBE","PBE","PBE"],
                                      levels=[5,3,2], c_size=[5,10,20,40,80,320],
                                      t_fraction=a, c_fraction=1.)
    savefig( cs_layers,  files.organic_figures*"ca_"*string(a)*"_layers.pdf" )
    savefig( cst_layers, files.organic_figures*"cat_"*string(a)*"_layers.pdf" )
end

scalefontsizes(0.5)

