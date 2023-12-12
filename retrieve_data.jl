using Downloads
using CSV
using DataFrames

# import Python libraries
using PyCall
ase   = pyimport("ase.io")
atoms = pyimport("ase")

# ====================================
# Settings
# ====================================

# file paths
path     = @__DIR__
remote   = ( trimer_xyz   = "https://github.com/jmbowma/q-AQUA/raw/main/3b_data.xyz",
             organic_xyz  = "https://github.com/kefisher98/IP_EA_deltaSCF/raw/main/organic_molecules.xyz",
             organic_data = "https://github.com/kefisher98/IP_EA_deltaSCF/raw/main/organic_ip_data.csv" )
home     = ( trimer_xyz   = path*"/data/3b_data.xyz",
             organic_xyz  = path*"/data/organic_molecules.xyz",
             organic_data = path*"/data/organic_ip_data.csv",
             ccsdt_dft    = path*"/data/trimer_dft.csv",
             ccsdt        = path*"/data/trimer_ccsdt.csv" )

# Hartree to eV
h2eV = 27.2114079527

# =====================================
# Download
# =====================================

Downloads.download( remote.trimer_xyz,   home.trimer_xyz   )
Downloads.download( remote.organic_xyz,  home.organic_xyz  )
Downloads.download( remote.organic_data, home.organic_data )

# =====================================
# Reformat Trimer data
# =====================================

# read in data
ccsdt     = CSV.read( home.ccsdt,     DataFrame )
ccsdt_dft = CSV.read( home.ccsdt_dft, DataFrame )
data      = [ parse(Float64, collect( ase.read( home.trimer_xyz, index=string(i) )["info"])[1] )  for i in ccsdt[!,:index] ]

# add to DataFrame
ccsdt[!,:CCSDT_Hartree] = data 
ccsdt[!,:CCSDT]         = data.*h2eV # conversion to eV
ccsdt_dft[:,:CCSDT]     = ccsdt[1:size(ccsdt_dft,1),:CCSDT]

# save
CSV.write( home.ccsdt,     ccsdt     )
CSV.write( home.ccsdt_dft, ccsdt_dft )

