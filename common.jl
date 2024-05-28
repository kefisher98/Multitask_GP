

##############################################################################
#
# Problem settings structure
#
##############################################################################


struct settings
    model           # function to define GP model; one of [model_delta, model_SE]
    tasks           # Array of strings identifying each task
    files           # namedTuple with file paths to data
    prior_mean      # function to calculate prior mean from training data; one of [mean, 0]
    globalize       # function to compute global features from local features; one of []
    positions       # function to locate atoms in molecule for feature construction; one of []
    species         # species encountered by SOAP feature
    rcut            # SOAP parameter: radius cut off
    sigma           # SOAP parameter: sigma of Gaussian neighborhood
    nmax            # SOAP parameter: n max of spherical harmonic expansion
    lmax            # SOAP parameter: l max of spherical harmonic expansion
end



##############################################################################
#
# DataFrame Filtering
#
##############################################################################


# Description: retrieve rows indices of dataframe where a given value appears in a given column

# Input: df  - dataframe to be filtered
#        col - column to screen
#        val - value to look for

function filter_inds( df, col, val  )

    return collect(1:size(df,1))[df[!,col] .=== val]

end

###########


# Description: retrieve rows of dataframe where columns equal some provided value

# Inputs: df        - dataframe to be filtered
#         filCols   - columns to screen ( ie [:color, :shape]   )
#         selectors - if these values appear in columns keep the dataframe row ( ie [ ("green", "circle" ), ("orange", "triangle"  )  ]  )
#         leftover  - if true, return all rows filtered out of the dataframe as a second result

# Note: calling the function with the above examples will return all rows of df where the color is green AND shape is circle OR the color is orange AND the shape is triangle

function multifilter(df, filCols, selectors; leftover=false)

    function filterRules(cols...)::Bool
        (cols) in selectors
    end

    return leftover ? (filter( filCols => filterRules , df ), filter( filCols => !filterRules , df )) : filter( filCols => filterRules , df )

end

###########

# Description: retrieve rows of a dataframe where values of given columns fall in a given range

# Inputs: df        - dataframe to be filtered
#         filCols   - columns to screen ( ie [:MAE, :RMSE]   )
#         lower     - lower threshold for values in filCols that correspond to rows we should keep ( ie 0.0 )
#         upper     - upper threshold for values in filCols that correspond to rows we should keep ( ie 0.1 )

function multifilter_range(  df, filCols, lower, upper )

    function filterRules(cols...)::Bool
        prod( [(lower <= c) for c in cols] )*prod( [(c<upper) for c in cols ] )
    end

    return filter( filCols => filterRules , df )

end


##############################################################################
#
# Array Manipulation
#
#############################################################################

### concatenation

# Description: apply vcat and hcat with a filter to remove any length 0 DataFrames

rhcat(A) =  hcat( [a for a in Iterators.filter( x->length(x)>0, A)]...  )

rvcat(A) =  vcat( [a for a in Iterators.filter( x->length(x)>0, A)]...  )



# Description: Divide elements of an array into fractional arrays (where the original ordering is preserved) and return a given fraction

# Inputs: arr    - array to divide
#         split  - number of fractional arrays to divide arr into
#         select - which fractional arrays to return
#
# Example: part( collect(1:18, 6, [2,3]  ) returns [4,5,6,7,8,9]

function part( arr, split, select  )
    section = div( length(arr), split )
    result  = split === select[end] ? vcat( [arr[ (s-1)*section + 1: s*section ] for s in select[1:end-1]]..., arr[ (select[end]-1)*section + 1:end  ]  ) :
                        vcat( [arr[ (s-1)*section + 1: s*section ] for s in select]...  )
    return result
end

#######

# Description: returns first part of an array

# Inputs: arr   - array to split
#         split - number of parts to split the array into

# Examples: fp( collect(1:18); split=3  ) returns [1,2,3,4,5,6]
#           fp( collect(1:18); split=2  ) returns [1,2,...,9]

function fp( arr; split=2  )
    return arr[ 1:div( length(arr), split ) ]
end

#######

# Description: returns last part of an array

# Inputs: arr   - array to split
#         split - number of parts to split the array into
#
# Examples: lp( collect(1:18); split=3  ) returns [13,14,15,16,17,18]
#           lp( collect(1:18); split=2  ) returns [10,11,...18]

function lp( arr; split=2  )
    return arr[ div( length(arr), split )+1:end ]
end

#######

# Description: replaces the last fraction of one array with the first fraction of a different array

# Inputs: C     - array for which we will remove the last elements in order to swap in a fraction of another array
#         S     - array from which we will grab the first fraction of elements and swap into C
#         split - the number of pieces we will split S into

# Note: The result should always be an array of length equal to length(S) with its first elements taken from C and final elements taken from S
#
# Examples: swap( [1,2,3,4], [11,12,13,14]; split=2  ) returns [1,2,11,12]
#           swap( [1,2,3,4,5,6,7,8] [11,12,13,14,15,16]; split=3 ) returns [1,2,3,4,11,12]
#
# If C is too short to create an array with the same length as S by swapping in the fraction of S, we keep all components of C, and add the first length(S) - length(C) components of S
#
# Example: swap( [1,2,3], [11,12,13,14,15,16]; split=3 )  ) returns [1,2,3,11,12,13] 

function swap( C, S; split=2  )

    if length(S)-div(length(S),split) > length(C)
        return vcat(C, S[1:end-length(C)])
    else
        return vcat(  C[ 1:length(S)-div(length(S),split) ],  fp(S;split)        )
    end
end

########

# Description: replace final element in array with new value
# Inputs: S   - array to modify
#         loc - new value to put in last element of S
#
# Example: pop( [5,6,7], 12) returns [5,6,12]  
function pop( S, loc)
    return vcat( S[1:end-1], loc  )
end

########

# Description: This function first calls the function swap, then calls the function pop, as defined above
function popswap( C, S, loc; split=2)
    return pop( swap( C, S; split ), loc  )
end


##############################################################################
#
# Tuple Manipulation
#
#############################################################################

############
# Dict to Tuple

# retrieve keys of dictionary
dictkeys(d::Dict)   = (collect(keys(d))...,)

# retrieve values of dictionary
dictvalues(d::Dict) = (collect(values(d))...,)

# take dictionary d and return tuple with the same keys and values
namedtuple(d::Dict{Symbol,T}) where {T} = NamedTuple{dictkeys(d)}(dictvalues(d))


#############################################################################
#
# Feature Construction
#
##############################################################################

# Description: returns the positions of all atoms that are not hydrogen for a given molecule
# Input: at - ase molecule 

function getNonHydrogenPositions(at)

    pos = at.get_positions()[ at.get_chemical_symbols() .!= "H"  ,:]
    return [pos[r,:] for r in 1:size(pos,1)]

end

##########

# Description: returns all positions of atoms in a given molecule
# Input:       at - ase molecule
# Note:        the point of this function is to provide a function in the same form as getNonHydrogenPositions()

getAllPositions(at) = at.get_positions()

##########

# Description: Features are originally generated as a vector of vectors. 
#              This function reformats that into a matrix where each column is the feature of a different molecule
# Input:       features - vector of feature vectors

molecule_columns( features  ) = reshape(collect(Iterators.flatten(features)), (length(features[1]),length(features)))

##########

# Note: the following two functions are designed to be stored in a settings struct,
#       and used as strategies to globalize the matrix of local features generated by desc.create(at)

# Description: Function to average local features
# Input: feature - matrix with a local feature at each row  

average_feature(  feature )   = mean( feature; dims=1 )

# Description: Function to concatenate all local features
# Input: feature - matrix with a local feature at each row

cat_feature( feature ) = reshape( feature', size(feature,1)*size(feature,2), 1)


# Description: Function to construct an average feature vector for each element in a given list of species

# Inputs: feature - matrix with a local feature at each row 
#         at      - ase molecule
#         species - vector of chemical symbols that local features correspond to 

function average_by_species( feature, at; species=["H","O"]  )
    
    # get chemical symbols
    symbols = at.get_chemical_symbols()

    # average over species the flatten
    return hcat( [ mean( feature[symbols .=== s,:]; dims=1 )  for s in species ]...  )

end

##########

# Description: This function takes a vector of features and formats them to be compatible with a Stheno GP model. There is a safeguard to exclude empty features.

# Inputs: label - Vector of GPPPInput labels corresponding to each feature. 
#                 See model_SE and model_delta for examples of labels (ie: :P, :S1, :Δ12, etc)
#         V     - Vector of features that will be used to implement GP model

# Output: GPPPInput features compatible with Stheno model

function featurize( label, V  )
    return BlockData( [ GPPPInput( label[i], ColVecs( V[i] ) )
                        for i in Iterators.filter( x->length(V[x])>0, 1:length(V) )  ]... )
end



##########

# Description: Function to construct a matrix where each column corresponds to a feature. 

# Input: desc     - dscribe descriptor 
#        ats      - vector of ase molecule to create features for 
#        settings - settings structure. the relevant fields are
#                       settings.positions in { getAllPositions, getNonHydrogenPositions   }
#                       settings.globalize in { average_feature, cat_feature }

function get_features(desc,ats::Vector; settings)

    # get average of non hydrogen local features; each molecule is represented by a column
    return molecule_columns( [settings.globalize( desc.create(at,  positions=settings.positions(at)) )  for at in ats[:,1]] )

end

# Description: Function to construct GPPPInput features

# Input: desc     - dscribe descriptor 
#        ats      - vector of ase molecule to create features for 
#        GP_ID    - GPPPInput label corresponding to the part of the GPPP model to which the feature corresponds  
#                   See model_SE and model_delta for examples of labels (ie: :P, :S1, :Δ12, etc)
#        settings - settings structure. the relevant fields are
#                       settings.positions in { getAllPositions, getNonHydrogenPositions   }
#                       settings.globalize in { average_feature, cat_feature }

function get_features(desc,ats::Vector, GP_ID; settings)

    # format as input to GPPP
    return GPPPInput( GP_ID,  ColVecs( get_features(desc,ats; settings)  ) )

end


###########

# Description: Function to construct a matrix where each column corresponds to a feature. Specifically, this function is used when the feature we are interested in is a difference between SOAP features.

# Input: desc     - dscribe descriptor
#        ats      - two column matrix ase molecule to create features for
#                   (a feature will be created for each row as the difference between the SOAP features for the molecules in the two columns)
#        settings - settings structure. the relevant fields are
#                       settings.positions in { getAllPositions, getNonHydrogenPositions   }
#                       settings.globalize in { average_feature, cat_feature }

function get_features(desc, ats; settings)
    
    return molecule_columns( [ settings.globalize( desc.create(ats[r,1],  positions=settings.positions(ats[r,1])) -
                                                   desc.create(ats[r,2],  positions=settings.positions(ats[r,2]))   )
                                                   for r in 1:size(ats,1) ] )

end


# Description: Function to construct GPPPInput features. Specifically, this function is used when the feature we are interested in is a difference between SOAP features.

# Input: desc     - dscribe descriptor 
#        ats      - two column matrix ase molecule to create features for
#                   (a feature will be created for each row as the difference between the SOAP features for the molecules in the two columns)
#        GP_ID    - GPPPInput label corresponding to the part of the GPPP model to which the feature corresponds  
#                   See model_SE and model_delta for examples of labels (ie: :P, :S1, :Δ12, etc)
#        settings - settings structure. the relevant fields are
#                       settings.positions in { getAllPositions, getNonHydrogenPositions   }
#                       settings.globalize in { average_feature, cat_feature }


function get_features(desc, ats, GP_ID; settings)
    
    # format as input to GPPP
    return GPPPInput( GP_ID,  ColVecs( get_features(desc, ats; settings)  ) )

end




##############################################################################
#
# GP Models
#
##############################################################################

# Description: Function to construct a squared exponential kernel. This function is meant to be used to construct Stheno GPPP.

# Inputs: v - variance hyperparameter
#         ℓ - lengthscale hyperparameter

SE( v, ℓ ) =  v * with_lengthscale(SqExponentialKernel(), ℓ )


############

# Description: Function that creates a GPPP model of a single or asymmetric multitask method which can incorporate up to six levels of theory.

# Inputs: θ - NamedTuple with parameters of the model. Fields include
#               v = Named Tuple with fields
#                   p = variance hyperparameter for primary task
#                   d = list of variance hyperparameters for secondary tasks 
#               l = Named Tuple with fields
#                   p = lengthscale hyperparameter for primary task
#                   d = list of lengthscale hyperparameters for secondary tasks 
#               ρ = list of correlation hyperparameters for secondary tasks 
#

# Output: A model which can be instantiated with training data.


function model_SE(θ::NamedTuple)
    return @gppp let

        # primary GP
        P  = GP(  θ.μ.p, SE( θ.v.p, θ.l.p )  )
        
        # Secondary GPs
        D1 = GP(  θ.μ.d[1], SE( θ.v.d[1], θ.l.d[1]  )   )
        S1 =   (  θ.ρ[1] * P)  +  D1

        D2 = GP(  θ.μ.d[2], SE( θ.v.d[2], θ.l.d[2]  )   )
        S2 =   (  θ.ρ[2] * P)  +  D2

        D3 = GP(  θ.μ.d[3], SE( θ.v.d[3], θ.l.d[3]  )   )
        S3 =   (  θ.ρ[3] * P)  +  D3

        D4 = GP(  θ.μ.d[4], SE( θ.v.d[4], θ.l.d[4]  )   )
        S4 =   (  θ.ρ[4] * P)  +  D4

        D5 = GP(  θ.μ.d[5], SE( θ.v.d[5], θ.l.d[5]  )   )
        S5 =   (  θ.ρ[5] * P)  +  D5

    end
end

##########

# Description: Function to format parameters for multitask GP model. This function will format all parameters provided without any reordering.


# Input: params - DataFrame with rows corresponding to each task. Columns should include [variance, lengthscale, rho]
# Output: NamedTuple of parameters formatted to be compatible as an input to model_SE

# Note: The dataframe should only include as many tasks as will be used in the model. While the model_SE function supports six tasks, extraneous tasks will be zeroed out.


function set_multitask_params( params::DataFrame )

    remainder = 6 - size( params  )[1]
        
    return ( v= ( p=params[1,:variance],    d=vcat( params[2:end,:variance],    [params[end,:variance]     for i in 1:remainder]...)),
             l= ( p=params[1,:lengthscale], d=vcat( params[2:end,:lengthscale], [params[end,:lengthscale]  for i in 1:remainder]...)  ),
             ρ= vcat( params[2:end,:rho], [params[end,:rho]  for i in 1:remainder]... )  )
end


##############

# Description: Function to format parameters for multitask GP model. This function will downselect parameters according to a provided subset and put parameters in the order of that subset.

# Input: all_params - DataFrame with rows corresponding to each task. Columns should include [variance, lengthscale, rho]
#        levels     - named tuple containing the tasks we want parameters for. Fields include
#                         primary   = the primary task
#                         secondary = the secondary tasks, in the order they are considered 

# Output: NamedTuple of parameters formatted to be compatible as an input to model_SE

# Note: The dataframe should only include as many tasks as will be used in the model. While the model_SE function supports six tasks, extraneous tasks will be zeroed out.


function set_multitask_params( all_params::DataFrame, levels::NamedTuple )

    params = multifilter( all_params, [:task], [(levels.primary,)] )
    for s in levels.secondary
        params = vcat( params, multifilter( params, [:task], [(s,)] ) )
    end

    remainder = 6 - size( params  )[1]
        
    return ( v= ( p=params[1,:variance],    d=vcat( params[2:end,:variance],    [params[end,:variance]     for i in 1:remainder]...)),
             l= ( p=params[1,:lengthscale], d=vcat( params[2:end,:lengthscale], [params[end,:lengthscale]  for i in 1:remainder]...)  ),
             ρ= vcat( params[2:end,:rho], [params[end,:rho]  for i in 1:remainder]... )  )
end

#########

# Description: Function to format parameters for multitask GP model. Specifically, this function adds prior mean parameters to a previously created parameter NamedTuple.

# Input: means  - a vector of prior means for each task
#        params - NamedTuple formatted as the output of set_multitask_params( params::DataFrame )
#
# Output: NamedTuple of parameters formatted to be compatible as an input to model_SE

# Note: The dataframe should only include as many tasks as will be used in the model. While the model_SE function supports six tasks, extraneous tasks will be zeroed out.


function set_multitask_params( means::NamedTuple, params::NamedTuple  )
    
    if haskey( means, :s  )
        return merge( params, (μ=( p=means.p, d=vcat(means.s, zeros(5-length(means.s) )  ) ),) )
    else
        return merge( params, (μ=( p=means.p, d=zeros(5)  ),) )
    end

end

############

# Description: Function that creates a GPPP model of the delta method which can incorporate up to six levels of theory.

# Inputs: θ - NamedTuple with parameters of the model. Fields include
#               μ = a vector of mean parameters for the delta between each pair of levels
#               v = a vector of variance parameters for the delta between each pair of levels
#               l = a vector of lengthscale parameters for the delta between each pair of levels

# Output: A model which can be instantiated with training data.

function model_delta(θ::NamedTuple)
    return @gppp let

        # paired levels
        Δ12  = GP(  θ.μ[1], SE( θ.v[1], θ.l[1] )  )
        Δ23  = GP(  θ.μ[2], SE( θ.v[2], θ.l[2] )  )
        Δ34  = GP(  θ.μ[3], SE( θ.v[3], θ.l[3] )  )
        Δ45  = GP(  θ.μ[4], SE( θ.v[4], θ.l[4] )  )
        Δ56  = GP(  θ.μ[5], SE( θ.v[5], θ.l[5] )  )


        # combination levels 
        Δ13  = Δ12 + Δ23
        Δ14  = Δ12 + Δ23 + Δ34
        Δ15  = Δ12 + Δ23 + Δ34 + Δ45
        Δ16  = Δ12 + Δ23 + Δ34 + Δ45 + Δ56


    end
end

##########

# Description: Function to format parameters for delta GP model. This function will format all parameters provided without any reordering.

# Input: params - DataFrame with rows corresponding to each task. Columns should include [variance, lengthscale]
# Output: NamedTuple of parameters formatted to be compatible as an input to model_delta

# Note: The dataframe should only include as many tasks as will be used in the model. While the model_delta function supports six tasks, extraneous tasks will be zeroed out.


function set_delta_params( params::DataFrame )

    remainder = 5 - size(params)[1]
    return ( v= vcat(  params[:,:variance],    [ params[end,:variance]     for i in 1:remainder   ]... ),
             l= vcat(  params[:,:lengthscale], [ params[end,:lengthscale]  for i in 1:remainder   ]... ) )
end

#########


# Description: Function to format parameters for delta GP model. This function will downselect parameters according to a provided subset and put parameters in the order of that subset.


# Input: params - DataFrame with rows corresponding to each task. Columns should include [variance, lengthscale]
#        levels     - named tuple containing the tasks we want parameters for. Fields include
#                         primary   = the primary task
#                         secondary = the secondary tasks, in the order they are considered
                         
# Output: NamedTuple of parameters formatted to be compatible as an input to model_delta

# Note: The dataframe should only include as many tasks as will be used in the model. While the model_delta function supports six tasks, extraneous tasks will be zeroed out.

function set_delta_params( all_params::DataFrame, levels::NamedTuple )

    params = multifilter( all_params, [:task], [(levels.primary,)] )
    for s in levels.secondary
        params = vcat( params, multifilter( params, [:task], [(s,)] ) )
    end

    remainder = 5 - size(params)[1]
    return ( v= vcat(  params[:,:variance],    [ params[end,:variance]     for i in 1:remainder   ]... ),
             l= vcat(  params[:,:lengthscale], [ params[end,:lengthscale]  for i in 1:remainder   ]... ) )
end



#########

# Description: Function to format parameters for delta GP model. Specifically, this function adds prior mean parameters to a previously created parameter NamedTuple.

# Input: means  - a vector of prior means for each task
#        params - NamedTuple formatted as the output of set_delta_params( params::DataFrame )
#
# Output: NamedTuple of parameters formatted to be compatible as an input to model_delta

# Note: The dataframe should only include as many tasks as will be used in the model. While the model_delta function supports six tasks, extraneous tasks will be zeroed out.

function set_delta_params( means::AbstractArray, params::NamedTuple  )

    return merge( params, ( μ=vcat(means, zeros(5-length(means))),)  )

end





##############################################################################
#
# Parameter Optimization
#
##############################################################################

# Description: Function which performs both mle estimation of lengthscale and minimum mae estimation of variance for a squared exponential kernel.

# Input: X         - matrix where each column is a feature
#        Y         - vector of outputs where each entry corresponds to a column of X
#        noise     - variance of data noise
#        context   - namedTuple with fields
#                       rho  = correlation parameter
#                       seed = seed used to initialize random number generator
#                       name = name of file to which to save results
#        θ_init    - initial value of variable that will be optimized
#        lower     - lower bound on variable that will be optimized
#        udata[shuffled,levels.primary]pper     - upper bound on variable that will be optimized
#        optAlg    - optimization algorithm
#        grad      - method for computing gradient

# Output: NamedTuple with variance and lengthscale hyperparameters, conervence status, and final objective value



function mle_mae( X, Y, noise, context; 
                  θ_init=[1.], lower=[1e-10], upper=[100.], optAlg=LBFGS(), grad=Optimization.AutoZygote(), maxtime=Inf )

    # lengthscale estimate (maximum likelihood)
    mle   = max_likelihood( X, Y, noise; θ_init, lower, upper, optAlg, grad )
    
    if mle.converged

        # variance estimation (minimize mean absolute error)
        mae   = min_mae_var( X, Y, mle.ℓ, noise; optAlg, grad, maxtime  )
        
        if mae.converged && mae.v[1] < 1e6
            addToFile( DataFrame( task=context.task, 
                                  variance=mae.v, 
                                  variance_mle=mle.v,
                                  lengthscale=mle.ℓ, 
                                  rho=context.rho, 
                                  seed=context.seed, 
                                  objective_mle=mle.objective,
                                  objective_mae=mae.objective ), context.name  )
            return ( converged=true, v=mae.v, v_mle=mle.v,  ℓ=mle.ℓ, mle=mle.objective, mae=mae.objective )
        end
    end
    return ( converged=false,)

end

# Description: Function to find MLE estimates of lengthscale, mean, and variance for the squared exponential kernel.

# Input: X         - matrix where each column is a feature
#        Y         - vector of outputs where each entry corresponds to a column of X
#        noise     - variance of data noise
#        θ_init    - initial value of variable that will be optimized
#        lower     - lower bound on variable that will be optimized
#        upper     - upper bound on variable that will be optimized
#        optAlg    - optimization algorithm
#        grad      - method for computing gradient
#        report    - boolean parameter which determines whether the optimization report is printed

# Output: NamedTuple with mean, variance, and lengthscale hyperparameters, conervence status, and final objective value

function max_likelihood( X, Y, noise; θ_init=[1.], lower=[1e-10], upper=[100.], optAlg=LBFGS(), grad=Optimization.AutoZygote(), report=true  )
    

    if typeof(optAlg) <: LBFGS || typeof(optAlg) <: BFGS
        convergence = Optim.converged
        min         = Optim.minimum
    elseif typeof(optAlg) <: CMAES
        convergence = Evolutionary.converged
        min         = Evolutionary.minimum
    else
        error("Requested optimization algorithm is unsupported. Use LBFGS, BFGS, or CMAES.")
    end


    # call optimization procedure
    ℓ = optParams(nlml_recurse_SE, (X=X, Y=Y, noise=noise), θ_init, lower, upper; optAlg, grad)

    # get μ and σ²
    μ, σ² = μ_and_σ2( X, Y ;ℓ )


    # report to audience
    if report
        @show ℓ.original
    end

    return ( μ=μ, v=σ², ℓ=ℓ.u[1], converged=convergence(ℓ.original), objective=min(ℓ.original) )
end


##########

# Descrition: Function to estimate variance for the squared exponential kernel by minimizing mean absolute error.

# Input: X         - matrix where each column is a feature
#        Y         - vector of outputs where each entry corresponds to a column of X
#        ℓ         - lengthscale hyperparameters
#        noise     - variance of data noise
#        θ_init    - initial value of variable that will be optimized
#        lower     - lower bound on variable that will be optimized
#        upper     - upper bound on variable that will be optimized
#        optAlg    - optimization algorithm
#        grad      - method for computing gradient
#        report    - boolean parameter which determines whether the optimization report is printed

# Output: NamedTuple with mean, variance, and lengthscale hyperparameters, conervence status, and final objective value

function min_mae_var( X, Y, ℓ, noise; 
                      θ_init=[1.], lower=[1e-10], upper=[100.], 
                      optAlg=LBFGS(), grad=Optimization.AutoZygote(), report=true, maxtime=Inf  )

    if typeof(optAlg) <: LBFGS || typeof(optAlg) <: BFGS
        convergence = Optim.converged
        min         = Optim.minimum
    elseif typeof(optAlg) <: CMAES
        convergence = Evolutionary.converged
        min         = Evolutionary.minimum
    else
        error("Requested optimization algorithm is unsupported. Use LBFGS, BFGS, or CMAES.")
    end

    # split data
    N     = 2*div( size(Y,1), 3 )
    train = sample(  1:size(Y,1), N, replace=false )
    test  = setdiff( 1:size(Y,1), train  )

    # optimize variance 
    params = (X=X, Y=Y, μ=mean(Y[train]), ℓ=ℓ, noise=noise, train=train, test=test)
    v      = optParams(mae_var, params, θ_init, lower, upper; optAlg, grad, maxtime)

    # report to audience
    if report
        @show v.original
    end

    return (μ=params.μ, v=exp.(v.u[1]), ℓ=ℓ, converged=convergence(v.original), objective=min(v.original) )


end


########## 

# Description: Function to set up and call optimization routine for generic objective, parameters, and algorithm.

# Inputs: objective - function to optimize
#         params    - NamedTuple of parameters necessary for objective
#         θ_init    - initial value of variable that will be optimized
#         lower     - lower bound on variable that will be optimized
#         upper     - upper bound on variable that will be optimized
#         optAlg    - optimization algorithm
#         grad      - method for computing gradient

# Output: solution of optimization problem

function optParams(objective, params, θ_init, lower, upper; optAlg=LBFGS(), grad=Optimization.AutoZygote(), maxtime=Inf)

    # Define optimization routing
    f    = OptimizationFunction(objective, grad)
    prob = OptimizationProblem(f, θ_init, params; lb=lower, ub=upper)

    # solve with optAlg algorithm
    sol = solve(prob, optAlg; maxtime)

    return sol


end

################################################################################
# Min MAE

# Description: Objective function for optimization. Specifically, for a given variance parameter, this function implements a Gaussian Process model and returns its mean absolute prediction error.

# Input: v      - variance hyperparameter
#        params - NamedTuple with fields
#                       X     - matrix where each column corresponds to a feature
#                       Y     - vector of output elements each corresponding to a column of X
#                       noise - data noise variance
#                       train - indices of columns of X/rows of Y which will be used to train the GP model
#                       test  - indices of columns of X/rows of Y which will be used to test the GP model
#                       μ     - mean of training data
#                       ℓ     - lengthscale hyperparameter

# Output: Mean absolute error of GP model for variance input


function mae_var( v, params)

    v = exp.(v)

    # construct K_ff
    dim  = length( params.train )
    K_ff = v.*recurse_Ψ( params.X[:,params.train]; col=1, dim, ℓ=params.ℓ  ) + diagm( fill( params.noise, dim )) 

    # construct K_*f
    dim  = length( params.test )
    K_sf = v.*recurse_Ψ( params.X[:,params.test], params.X[:,params.train]; col=1, dim, ℓ=params.ℓ  )
    
    # prediction
    star = K_sf*( K_ff\ ( params.Y[params.train] -  fill(params.μ, length(params.train))  )) + fill( params.μ, length(params.test)   )

    # return mae
    return mean( abs.( star - params.Y[params.test]) )

end

################################################################################
# Gaussian Maximum Likelihoods

# Description: Function to compute a row of a covariance matrix with unit variance using the squared exponential kernel.

# Input: Xcol - feature corresponding to the row of the covariance matrix which we are constructing
#        Xmat - matrix of all features included in the covariance matrix (each column is a feature)
#        ℓ    - lengthscale hyperparameter

# Note: features and ℓ should be column vectors

Ψrow( Xcol, Xmat;ℓ ) = exp.( sum( -(Xcol .- Xmat).^2 ./ 2(ℓ.^2), dims=1 ) )

##########

# Description: Recursively construct a covariance matrix with unit variance from features X using the squared exponential kernel.

# Inputs: X     - matrix where each column corresponds to a feature
#         col   - index of column of X that a particular call to this function compares with the entire matrix X
#                 ( Initial calls to this function should use 1, and recursive calls will increment col )
#         dim   - the total number of columns of X
#         ℓ     - lengthscale hyperparameter

# Output: matrix comparing features in X (starting with the feature located at col) to the entire matrix X using the squared exponential kernel


function recurse_Ψ( X; col, dim, ℓ   )
    if col < dim
        return vcat( Ψrow( X[:,col], X; ℓ ), recurse_Ψ( X; col=col+1, dim, ℓ   )   )
    else
        return Ψrow( X[:,col], X; ℓ )
    end
end

# Description: Recursively construct a matrix comparing features in X to features in Y using the squared exponential kernel.

# Inputs: X     - matrix where each column corresponds to a feature belonging to the first set of features that we are comparing
#         Y     - matrix where each column corresponds to a feature belonging to the second set of features that we are comparing
#         col   - index of column of X that a particular call to this function compares with Y. 
#                 ( Initial calls to this function should use 1, and recursive calls will increment col )
#         dim   - the total number of columns of X
#         ℓ     - lengthscale hyperparameter

# Output: matrix comparing features in X (starting with the feature located at col) to features in Y using the squared exponential kernel

function recurse_Ψ( X, Y; col, dim, ℓ   )
    if col < dim
        return vcat( Ψrow( X[:,col], Y; ℓ ), recurse_Ψ( X, Y; col=col+1, dim, ℓ   )   )
    else
        return Ψrow( X[:,col], Y; ℓ )
    end
end


########

# Description: Objective function for optimization. Specifically, this function implements the negative log likelihood with a recursive construction of the covariance matrix to be compatible with autozygote.

# Inputs: ℓ      - lengthscale hyperparameter 
#         params - NamedTuple with fields
#                       X     - matrix where each column corresponds to a feature
#                       Y     - vector of output elements each corresponding to a column of X
#                       noise - data noise variance

# Output: negative log likelihood

function nlml_recurse_SE(ℓ, params)

    # factorize Ψ
    dim  = size( params.X, 2  )
    Ψ    = recurse_Ψ( params.X; col=1, dim, ℓ  ) + diagm( params.noise*ones(dim))

    # variance
    μ    = sum( Ψ\params.Y )   / sum(  Ψ\ones(dim)  )
    σ²   = (params.Y .- μ)' *(  Ψ\(params.Y .- μ)  ) / dim

    # log likelihood
    lml     = -0.5*dim*log(σ²) - 0.5*log( abs( det(Ψ) )  )

    return -lml

end


########

# Description: Function to construct a covariance matrix with unit variance for set of training features X.

# Inputs: X     - matrix where each column corresponds to a feature
#         ℓ     - lengthscale hyperparameter of kernel
#         noise - data noise variance

# Output: Covariance matrix for features X with unit variance. (Outer functions will multiply this result by the variance parameter)

function Ψ_mat(X;ℓ,noise=0.015)

    # each column j of X is the feature of the jth observation
    mat = zeros( size(X,2), size(X,2) )

    println("construct...")
    for i in 1:size(X,2)

        # diagonal
        mat[i,i] = Ψ(X[:,i],X[:,i];ℓ) + noise

        # off diagonals
        for j in i+1:size(X,2)

            ij       = Ψ(X[:,i],X[:,j];ℓ)
            mat[i,j] = ij
            mat[j,i] = ij

        end
    end

    return mat

end



#########

# Description: Squared exponential kernel with unit variance.

# Inputs: X1, X2 - features to compare
#         ℓ      - lengthscale of kernel

# Output: kernel value comparing X1 and X2

Ψ( X1,X2;ℓ ) = exp( sum( -(X1-X2).^2 ./ 2(ℓ.^2) ) )

#########

# Description: Compute MLE mean and variance estimates given data, lengthscale, and noise.

# Inputs: X     - matrix where each column corresponds to a molecular feature
#         Y     - vector of predictions (each element corresponds to a column of X)
#         ℓ     - lengthscale of squared exponential kernel
#         noise - variance of data noise

# Outout: μ     - MLE mean
#         σ²    - MLE variance

function μ_and_σ2( X, Y; ℓ, noise=0.015)

    # factorize Ψ
    factors = lu!( Ψ_mat(X;ℓ,noise) )

    # mean and variance
    μ       = sum( factors\Y )   / sum(  factors\ones(size(X,2),1)  )
    σ²      = (Y .- μ)'*(  factors\(Y .- μ)  ) / size(X,2)

    return μ,σ²[1]

end



##############################################################################
#
# Inner loop Inference
#
##############################################################################

# Description: Function to train a Stheno GP model and make new predictions.

# Inputs: model          - GP model to use; one of [model_delta, model_SE]
#         θ              - NamedTuple with parameters of GP model
#         train_features - a GPPPInput containing the input data used to train the GP model
#         train_data     - vector of output data corresponeding to train_features which is used to train the GP model
#         test_features  - a GPPPInput containing the input data where the GP model will make predictions
#         noise          - variance of data noise

# Outputs: μ    - predicted mean for each data point
#          σ    - predicted standard deviation for each data point
#          dt   - a two element vector: [ training time, inference time  ]


function infer( model, θ, train_features, train_data, test_features; noise=0.015  )

    t1      = time()
    GP_fun  = model(θ)
    GP_vec  = GP_fun(train_features, noise)
    post    = posterior(GP_vec,  train_data  )
    t2      = time()
    star    = marginals( post(test_features) )
    dt      = [t2-t1, time()-t2]

    return mean.(star), std.(star), dt

end

###############

# Description: Function to train a Stheno GP model and record performance. Specifically, this model handles cases where we have a vector of test data points.

# Inputs: name           - file name (including path) where inference results will be saved
#         model          - GP model to use; one of [model_delta, model_SE]
#         θ              - NamedTuple with parameters of GP model
#         train_features - a GPPPInput containing the input data used to train the GP model
#         train_data     - vector of output data corresponeding to train_features which is used to train the GP model
#         test_features  - a GPPPInput containing the input data where the GP model will make predictions
#         test_data      - vector of output data corresponeding to test_features which is used to test the GP model
#         identifier     - identifier - a vector of specific details about the test case, including the number of training points from primary and secondary tasks
#                          ( this vector is constructed in the middle loop of inference, but may be augmented with user input parameters )
#         noise          - variance of data noise

# Outputs: Summary of performance of inference is saved to CSV file at name.

function infer( name, model, θ, train_features, train_data, test_features, test_data::AbstractArray, identifier; noise=0.015)

    # predictions
    μ, σ, dt = infer( model, θ, train_features, train_data, test_features; noise  )

    # error
    addToFile( identifier, dt, μ, σ, test_data, name  )

end

################


# Description: Function to train a Stheno GP model and record performance. Specifically, this model handles cases where we have a single test data point.

# Inputs: name           - file name (including path) where inference results will be saved
#         model          - GP model to use; one of [model_delta, model_SE]
#         θ              - NamedTuple with parameters of GP model
#         train_features - a GPPPInput containing the input data used to train the GP model
#         train_data     - vector of output data corresponeding to train_features which is used to train the GP model
#         test_features  - a GPPPInput containing the input data where the GP model will make a prediction
#         test_data      - the output data point corresponeding to the test_feature which is used to test the GP model
#         identifier     - identifier - a vector of specific details about the test case, including the number of training points from primary and secondary tasks
#                          ( this vector is constructed in the middle loop of inference, but may be augmented with user input parameters )
#         noise          - variance of data noise

# Outputs: Summary of performance of inference is saved to CSV file at name.

function infer( name, model, θ, train_features, train_data, test_features, test_data::Number, identifier; noise=0.015)

    # predictions
    μ, σ, dt = infer( model, θ, train_features, train_data, test_features; noise  )

    # error
    addToFile( identifier, dt, μ[1], σ[1], test_data, name  )

end


###############

# Description: Function to train a Stheno GP model and record performance. Specifically, this model handles cases where we both have a vector of test data points and are making a prediction of a difference/delta.
# Inputs: name           - file name (including path) where inference results will be saved
#         model          - GP model to use; one of [model_delta, model_SE]
#         θ              - NamedTuple with parameters of GP model
#         train_features - a GPPPInput containing the input data used to train the GP model
#         train_data     - vector of output data corresponeding to train_features which is used to train the GP model
#         test_features  - a GPPPInput containing the input data where the GP model will make predictions
#         test_data      - vector of output data corresponeding to test_features which is used to test the GP model
#         baseline       - value which we add to the mean prediction of a difference/delta to obtain final prediction to compare with test_data
#         factor         - factor by which we multiply prediction to obtain final prediction to compare with test_data
#         identifier     - identifier - a vector of specific details about the test case, including the number of training points from primary and secondary tasks
#                          ( this vector is constructed in the middle loop of inference, but may be augmented with user input parameters )
#         noise          - variance of data noise

# Outputs: Summary of performance of inference is saved to CSV file at name.



function infer( name, model, θ, train_features, train_data, test_features, test_data::AbstractArray, baseline::AbstractArray, factor, identifier; noise=0.015)

    # predictions
    μ, σ, dt = infer( model, θ, train_features, train_data, test_features; noise  )

    # error
    addToFile( identifier, dt, μ, σ, test_data, baseline, factor, name  )

end

################

# Description: Function to train a Stheno GP model and record performance. Specifically, this model handles cases where we both have a single test data point and are making a prediction of a difference/delta.

# Inputs: name           - file name (including path) where inference results will be saved
#         model          - GP model to use; one of [model_delta, model_SE]
#         θ              - NamedTuple with parameters of GP model
#         train_features - a GPPPInput containing the input data used to train the GP model
#         train_data     - vector of output data corresponeding to train_features which is used to train the GP model
#         test_features  - a GPPPInput containing the input data where the GP model will make a prediction
#         test_data      - the output data point corresponeding to the test_feature which is used to test the GP model
#         baseline       - value which we add to the mean prediction of a difference/delta to obtain final prediction to compare with test_data
#         factor         - factor by which we multiply prediction to obtain final prediction to compare with test_data
#                          ( this value will be 1 unless we are learning a scaled difference )
#         identifier     - identifier - a vector of specific details about the test case, including the number of training points from primary and secondary tasks
#                          ( this vector is constructed in the middle loop of inference, but may be augmented with user input parameters )
#         noise          - variance of data noise

# Outputs: Summary of performance of inference is saved to CSV file at name.



function infer( name, model, θ, train_features, train_data, test_features, test_data::Number, baseline::Number, factor, identifier; noise=0.015)

    # predictions
    μ, σ, dt = infer( model, θ, train_features, train_data, test_features; noise  )

    # error
    addToFile( identifier, dt, μ[1], σ[1], test_data, baseline, factor, name  )

end


#############################################################################
#
# Input
#
#############################################################################

# Description: Function that sets up inference settings structure and guide tuple to provide details on a particular test. This function also initializes the output file.

# Inputs: files  - NamedTuple containing paths to data files
#                  fields include:
#                      data        - path to prediction data used to train GP
#                      data_augmet - path to additional data for larger tests on primary task
#                      xyz         - path to exended xyz file with molecular configurations for training and testing
#                      params      - path to GP hyperparameter CSV file
#                      opt_inds    - path to CSV file with indices of data set aside for optimization
#                      output_path - initial path and string for output file names
#                      figure_path - initial path and string for output file names
#         levels - NamedTuple with names of tasks; fields include
#                      primary     - string with primary task name; ie "CCSDT"
#                      secondary   - vector of strings with secondary task names; ie ["PBE","SCAN"]
#         seed   - integer we use to intialize the random number generator for a given test
#         target - number of molecules with which we test an inference method
#         feature_method - function that will map local features to global feature; one of [average_feature, cat_feature]
#         atom_position  - function that will get the positions of atoms in a molecule that are used as centers in a feature;
#                          one of [getNonHydrogenPositions, getAllPositions]
#         species - all species that will be encountered in a SOAP neighborhood 
#         rcut    - SOAP parameter: radius cut off
#         sigma   - SOAP parameter: sigma of Gaussian neighborhood
#         nmax    - SOAP parameter: n max of spherical harmonic expansion
#         lmax    - SOAP parameter: l max of spherical harmonic expansion


# Output: GP_settings - settings struct with information about SOAP parameters and the GP model for single and multitask cases
#         Δ_settings  - settings struct with information about SOAP parameters and the GP model for delta cases
#         guide       - Named Tuple with informaion about the test case and the names of tasks; 
#                       details are at the last section of this script
#         params      - DataFrame with parameters in file located at files.params

function set_up( ; files::NamedTuple, levels::NamedTuple, seed::Int, target::Number,
                   feature_method=average_feature, atom_positions=getNonHydrogenPositions,
                   species, rcut=10, sigma=0.4, nmax=8, lmax=8, custom_heading=false  )

    # define settings features
    tasks       = vcat( levels.primary, levels.secondary )
    GP_settings = settings( model_SE,    tasks, files, mean, feature_method, atom_positions, species, rcut, sigma, nmax, lmax)
    Δ_settings  = settings( model_delta, tasks, files, mean, feature_method, atom_positions, species, rcut, sigma, nmax, lmax)

    # guide to test case and data
    guide = ( p          = levels.primary,
              s          = levels.secondary,
              case       =  "C",
              identifier = [],
              filename   = files.output_path*"target"*string(target)*"_seed"*string(seed)*".csv" )

    # initialize parameters
    params = CSV.read( files.params, DataFrame )

    # initialize output
    if custom_heading===false
        headers = Symbol.( vcat( [:levels, :lowest, :high_level, :low_level, :method, :datasets, :task_fractions,
                                  :train_time, :test_time],
                                 [[ "mae_", "median_ae_", "minimum_ae_", "maximum_ae_" ].*s
                                     for s in ["error", "relative", "std"]]...,
                                 [["pearson_", "spearman_","kendall_" ].*s
                                    for s in ["error", "std"]]... ) )
    else
        headers = Symbol.( custom_heading )
    end
    makeFile( headers, guide.filename )

    return GP_settings, Δ_settings, guide, params


end


##########

# Description: Function to retrieve X and Y data pairs for training and testing an inference model. This function randomly samples Ndata rows from the entire data frame, excluding rows that are earmarked for optimization.

# Inputs: settings    - settings struct with information about SOAP parameters and the GP model for single and multitask cases
#         Ndata       - the number of X and Y data pairs to retrieve
#         differences - boolean which is true if we want to return differences between all Y prediction methods
#                       ie. if two Y methods are "CCSDT" and "PBE", we return "CCSDT", "PBE", and "CCSDT-PBE"

# Outputs: X - matrix where each column corresponds to a molecular feature
#          Y - DataFrame with rows corresponding to different molecular systems and columns corresponding to tasks -- that is, each row corresponds to a column of X


function get_XY( settings::settings, Ndata::Number; differences=true )

    # trimers
    data    = CSV.read( settings.files.data, DataFrame  )

    # select rows
    remove     = CSV.read( settings.files.opt_inds, DataFrame)[!,1]
    rows       = sample( setdiff( 1:size(data)[1], remove ), Ndata, replace=false )
    
    return get_XY( settings, data, rows, settings.tasks; differences  )
end


###############

# Description: Function to retrieve X and Y data pairs for training and testing an inference model. This function is aimed at cases where we previously drew n data points for all tasks using get_XY and now want to draw m>n data points from a subset of the tasks. For instance, we may want more CCSD(T) data points than DFT data points to train larger single task GP models. 

# Inputs: settings    - settings struct with information about SOAP parameters and the GP model for single and multitask cases
#         Ndata       - a vector containing: [ n: the number of data points drawn previously, m: the total numbed of data points we want now  ]
#         differences - boolean which is true if we want to return differences between all Y prediction methods
#                       ie. if two Y methods are "CCSDT" and "PBE", we return "CCSDT", "PBE", and "CCSDT-PBE"
#         match       - this parameter aims at standardizing draws of the augmented data set with previous data sets obtained with get_XY();
#                       it is equal to the total data points available when we drew the previous data set. 
#                       (ie, if we drew n data points out of 3n data points, match should equal 3n)
#         tasks       - all tasks we want Y data for

# Outputs: X - matrix where each column corresponds to a molecular feature
#          Y - DataFrame with rows corresponding to different molecular systems and columns corresponding to tasks -- that is, each row corresponds to a column of X



function augment_XY( settings::settings, Ndata::AbstractArray; differences=true, match=2992, tasks=["CCSDT"] )

    # trimers
    data    = CSV.read( settings.files.augment_data, DataFrame  )

    # select rows
    remove     = CSV.read( settings.files.opt_inds, DataFrame)[!,1]
    rows       = sample( setdiff( 1:match, remove ), Ndata[1],          replace=false )
    if size(data)[1] > match && Ndata[2]-Ndata[1]>0 
        rows       = vcat( rows, sample( match+1:size(data)[1], Ndata[2]-Ndata[1], replace=false )  )
    end
    return get_XY( settings, data, rows, tasks; differences )
end


###############

# Description: Function to retrieve X and Y data pairs for training and testing an inference model. This function retrieves data from specified row indices and generates the corresponding features.

# Inputs: settings    - settings struct with information about SOAP parameters and the GP model for single and multitask cases
#         data        - full DataFrame with a different molecular system on each row and different task predictions on different columns
#         rows        - set of indices of molecular systems that we want X, Y pairs for
#         tasks       - all tasks we want Y data for (these tasks must be column names in data)
#         differences - boolean which is true if we want to return differences between all Y prediction methods
#                       ie. if two Y methods are "CCSDT" and "PBE", we return "CCSDT", "PBE", and "CCSDT-PBE"

# Outputs: X - matrix where each column corresponds to a molecular feature
#          Y - DataFrame with rows corresponding to different molecular systems and columns corresponding to tasks -- that is, each row corresponds to a column of X



function get_XY( settings::settings, data::DataFrame, rows::AbstractArray, tasks::AbstractArray; differences=true )

    Y = data[rows,tasks]
    if differences
        [ Y[!,tasks[a]*"_"*tasks[b]] = data[rows,tasks[a]]-data[rows,tasks[b]] 
          for (a,b) in collect(permutations(1:length(tasks),2)) ]
    end

    # features
    structures = [ ase.read( settings.files.xyz, index=string(i) ) for i in data[rows,:index] ]
    desc       = dsd.SOAP(species=settings.species,periodic=false,rcut=settings.rcut,nmax=settings.nmax,lmax=settings.lmax,sigma=settings.sigma)
    X          = get_features( desc, structures; settings )

    return X, Y
    
end



#############################################################################
#
# Output
#
############################################################################# 

# Description: Function to make a CSV file containing an empty DataFrame with column names headers.

# Inputs: headers - column names for DataFrame to be saved
#         name    - name (including path) of CSV file

# Output: The DataFrame is added to the CSV file given by name.


function makeFile(headers,name)

    df = DataFrame()
    [df[!,h] = Any[]  for h in headers ]
    addToFile(df,name)

end

###########

# Description: Function to add a DataFrame to a CSV file. If the CSV file exists, the function will append to it. If it does not exist already, the file will be created.

# Inputs: df   - DataFrame to add to CSV file
#         name - name (including path) of CSV file

# Output: The DataFrame is added to the CSV file given by name.

function addToFile(df,name)

    if isfile(name)
        CSV.write(name, df, append=true )
    else
        CSV.write(name, df )
    end

end

##############

# Description: Function to compute and record statistics of inference performance. 

# Inputs : identifier - a vector of specific details about the test case, including the number of training points from primary and secondary tasks
#                       ( this vector is constructed in the middle loop of inference, but may be augmented with user input parameters )
#          dt         - a two element vector: [ training time, inference time  ]
#          μ          - predicted mean for each data point
#          σ          - predicted standard deviation for each data point
#          test_data  - the true values to which we will compare predictions
#          name       - the name (including path) of the CSV file to which we will write the results 

# Output: This function calls a further addToFile function to write statistics to the CSV file given by name


function addToFile( identifier, dt, μ, σ, test_data::AbstractArray, name  )

    addToFile( DataFrame( reshape(  vcat( identifier,
                                          dt,
                                          abs_stats(  μ-test_data ),
                                          abs_stats( (μ-test_data) ./ test_data ),
                                          abs_stats(  σ  ),
                                          cors( μ, test_data ),
                                          cors( σ, abs.(μ-test_data) ) ), 1, length(identifier)+length(dt)+18 ),   :auto ), name )

end

##############

# Description: Function to compute and record statistics of inference performance. Specifically, this function is for cases where we have inferred a Delta/difference and want statistics with respect to the absolute value.

# Inputs : identifier - a vector of specific details about the test case, including the number of training points from primary and secondary tasks
#                       ( this vector is constructed in the middle loop of inference, but may be augmented with user input parameters )
#          dt         - a two element vector: [ training time, inference time  ]
#          μ          - predicted mean for each data point
#          σ          - predicted standard deviation for each data point
#          test_data  - the true values to which we will compare predictions
#          baseline   - vector of values which we add to μ to obtain final prediction to compare with test_data
#          factor     - factor by which we multiply μ+baseline to obtain final prediction to compare with test_data 
#                       (this will be 1 unless we have learned a scaled difference)
#          name       - the name (including path) of the CSV file to which we will write the results 

# Output: This function calls a further addToFile function to write statistics to the CSV file given by name

function addToFile( identifier, dt, μ, σ, test_data::AbstractArray, baseline::AbstractArray, factor, name  )

    addToFile( DataFrame( reshape(  vcat( identifier,
                                          dt,
                                          abs_stats(  μ-test_data ),
                                          abs_stats( (μ-test_data) ./ (factor*(baseline+test_data)) ),
                                          abs_stats(  σ  ),
                                          cors( factor*(baseline+μ), (factor*(baseline+test_data)) ),
                                          cors( σ, abs.(μ-test_data) ) ), 1, length(identifier)+length(dt)+18 ),   :auto ), name )


end

#############################################################################
#
# Error calulation
#
#############################################################################

# Description: this function calculates statistics on the absolute values of elements in a vector.

# Input: v - a vector of elements for which to compute summary statistics
# Output: a list: [ Mean Absolute Error (MAE), Median Absolute Error, Minimum Absolute Error, Maximum Absolute Error  ]

# Note: This function is used to evaluate the performance of inference methods.

function abs_stats( v  )

    return [mean(abs.(v)), median(abs.(v)), minimum(abs.(v)), maximum(abs.(v))]

end

#######

# Description: Function to calculate correlation coefficients between two vectors

# Inputs: v, w - the two vectors of data for which function computes the correlation coefficients 
# Output: a list: [ pearson's correlation coefficient, spearman's correlation coefficient, kendall's correlation coefficient ]

# Note: This function is used to evaluate the performance of inference methods. We evaluate the correlation of prediction to the truth and uncertainty predictions to the error.

function cors( v, w  )

    return  [ cor( v, w ), corspearman( v, w ), corkendall( v, w )  ]

end


#################################################################################
#
# Middle Loop Inference
# Implementation of single, multitask, and delta methods
#
################################################################################


# Description: function which performs single task inference.

# Inputs: settings     - settings structure containing inference model and method for computing prior mean
#         params       - NamedTuple with variance, lengthscale, and correlation parameters
#         data         - DataFrame with rows corresponding to different molecular systems and columns corresponding to tasks 
#         features     - matrix where each column corresponds to a molecular feature--that is, columns of this matrix correspond to rows of the data DataFrame 
#         rows         - NamedTuple with fields
#                               p = indices used to train the primary task
#                               s = indices used to train each secondary task
#                               T = indices of testing/target molecules
#         cols          - NamedTuple with names of tasks as they appear as column names in the data dataframe (this tuple is called guide in outer loop)

# Output: This function calls the infer function which writes inference results to the file saved under the key, 'filename', in cols.

function single( settings, params, data::DataFrame, features, rows, cols )

    # set parameters
    params = set_multitask_params( ( p=   settings.prior_mean(data[rows.p,cols.p]), s= [] ), params )

    # set test ID
    identifier = vcat( [ 1, cols.p, length(rows.p), 0, cols.p, cols.p, 0 ], cols.identifier )

    # infer. results are written to CSV within inference function
    infer(     cols.filename,
               settings.model,                                                  # inference model
               params,
               GPPPInput( :P,    ColVecs( features[:,rows.p] ) ),
               data[rows.p,cols.p],                                             # training data
               GPPPInput( :P,  ColVecs( features[:,rows.t] ) ),                 # testing features
               data[rows.t,cols.p],                                             # testing data
               identifier )

end



##################

# Description: function which performs multitask inference for an arbitrary number of tasks. The structure of the training set is determined in outer loops. 

# Inputs: settings     - settings structure containing inference model and method for computing prior mean
#         params       - NamedTuple with variance, lengthscale, and correlation parameters
#         data         - DataFrame with rows corresponding to different molecular systems and columns corresponding to tasks 
#         features     - matrix where each column corresponds to a molecular feature--that is, columns of this matrix correspond to rows of the data DataFrame 
#         rows         - NamedTuple with fields
#                               p = indices used to train the primary task
#                               s = indices used to train each secondary task
#                               T = indices of testing/target molecules
#         cols          - NamedTuple with names of tasks as they appear as column names in the data dataframe (this tuple is called guide in outer loop)

# Output: This function calls the infer function which writes inference results to the file saved under the key, 'filename', in cols.


function multitask( settings, params, data::DataFrame, features, rows, cols  )

    # set parameters
    params = set_multitask_params( ( p=   settings.prior_mean(data[rows.p,cols.p]),
                                     s= [ settings.prior_mean(data[rows.s[i],cols.s[i]]) - 
                                          params.ρ[i]*settings.prior_mean(data[rows.p,cols.p])     
                                          for i in 1:length(rows.s)  ]  ), params )

    # set test ID
    identifier = vcat( [ 1+length(cols.s), cols.s[end], length(rows.p), sum([length(s) for s in rows.s]), "Multitask", cols.case, cols.shared ], cols.identifier )
    
    # infer. results are written to CSV within inference function
    infer(     cols.filename,
               settings.model,                                                            # inference model
               params,
               BlockData( GPPPInput( :P,    ColVecs( features[:,rows.p] ) ),
                         [GPPPInput( Symbol( "S"*string(i) ), ColVecs( features[:, rows.s[i]]  )  ) for i in 1:length(rows.s)  ]... ),
               vcat( data[rows.p,cols.p], [data[rows.s[i],cols.s[i]] for i in 1:length(rows.s) ]... ),                              # training data
               GPPPInput( :P,  ColVecs( features[:,rows.t] ) ),                                  # testing features
               data[rows.t,cols.p],                                                           # testing data
               identifier )

end

#########


# Description: function which is designed to test the impact of including secondary task predictions of the target data in the training set. These secondary target training data points are added to the training set one by one. That is, we train a model for each target molecule and test its acccuracy predicting that target. In practice, if we want to use secondary target data, we would train one model on all points--a much more efficient approach. This function is designed to isolate the advantage of the method from the advantage of training on a large amount of data.

# Inputs: settings     - settings structure containing inference model and method for computing prior mean
#         params       - NamedTuple with variance, lengthscale, and correlation parameters
#         data         - DataFrame with rows corresponding to different molecular systems and columns corresponding to tasks
#         features     - matrix where each column corresponds to a molecular feature--that is, columns of this matrix correspond to rows of the data DataFrame
#         rows         - NamedTuple with fields
#                               p = indices used to train the primary task
#                               s = indices used to train each secondary task
#                               T = indices of testing/target molecules
#         cols          - NamedTuple with names of tasks as they appear as column names in the data dataframe (this tuple is called guide in outer loop)

# Output: This function calls the addToFile function which writes inference results to the file saved under the key, 'filename', in cols.



function multitask_with_target( settings, params, data::DataFrame, features, rows, cols; noise=0.015  )
    
    # set parameters
    params = set_multitask_params( ( p=   settings.prior_mean(data[rows.p,cols.p]),
                                     s= [ settings.prior_mean(data[rows.s[i],cols.s[i]]) - params.ρ[i]*settings.prior_mean(data[rows.p,cols.p])     
                                          for i in 1:length(rows.s)  ]  ), params )

    # set test ID
    identifier = vcat( [ 1+length(cols.s), cols.s[end], length(rows.p), sum([length(s) for s in rows.s])+1, "Multitask", cols.case*"T", cols.shared ], cols.identifier )


    # perform inference by training on one low level prediction for the target at a time
    output = [ infer( settings.model, 
                      params,
                      BlockData(  GPPPInput( :P,    ColVecs( features[:,rows.p] ) ),
                                 [GPPPInput( Symbol( "S"*string(i) ), ColVecs( features[:, rows.s[i]]  )  ) for i in 1:length(rows.s)-1  ]...,
                                 GPPPInput( Symbol( "S"*string(length(rows.s)) ), ColVecs( features[:, vcat( rows.s[end], rows.t[ind] )]  ) )   ),
                      vcat( data[rows.p,cols.p], [data[rows.s[i],cols.s[i]] for i in 1:length(rows.s) ]..., data[rows.t[ind],cols.s[end]] ),
                      GPPPInput( :P,  ColVecs( reshape(  features[:,rows.t[ind]],:,1 ) ) ); 
                      noise  )
               for ind in 1:length(rows.t)   ]

    μ      = vcat( [ result[1] for result in output  ]...)
    σ      = vcat( [ result[2] for result in output  ]...) 

    dt_sum = sum( [ result[3] for result in output  ]  ) 


    # error
    addToFile( identifier, dt_sum, μ, σ, data[rows.t,cols.p], cols.filename  )

end


#####################

# Description: function which implements the delta method for an arbitrary number of tasks. The structure of the training set is determined in outer loops.

# Inputs: settings     - settings structure containing inference model and method for computing prior mean
#         params       - NamedTuple with variance, lengthscale, and correlation parameters
#         data         - DataFrame with rows corresponding to different molecular systems and columns corresponding to tasks
#         features     - matrix where each column corresponds to a molecular feature--that is, columns of this matrix correspond to rows of the data DataFrame
#         rows         - NamedTuple with fields
#                               p = indices used to train the primary task
#                               s = indices used to train each secondary task
#                               T = indices of testing/target molecules
#         cols          - NamedTuple with names of tasks as they appear as column names in the data dataframe (this tuple is called guide in outer loop)

# Output: This function calls the infer function which writes inference results to the file saved under the key, 'filename', in cols.



function delta( settings, params, data::DataFrame, features, rows, cols  )


    # clear non-overlapping indices
    rows = ( t=rows.t, d = vcat( [intersect( rows.p, rows.s[1] )],   [ [intersect( rows.s[i-1], rows.s[i] )] for i in 2:length(rows.s)  ]...  ),  ) 



    # set parameters
    params = set_delta_params( vcat(   settings.prior_mean( data[rows.d[1], cols.p] - data[rows.d[1], cols.s[1]] ), 
                                     [ settings.prior_mean( data[rows.d[i], cols.s[i-1]] - data[rows.d[i], cols.s[i]] ) for i in 2:length(cols.s)  ]... ),
                               params )

    
    identifier = vcat( [ 1+length(cols.s), cols.s[end], length(rows.d[1]), sum([length(d) for d in rows.d])+1, "Delta", cols.case*"T", cols.shared ], cols.identifier )
    infer(      cols.filename, 
                settings.model,                                                 # inference model
                params,
                BlockData( GPPPInput( :Δ12,    ColVecs( features[:,rows.d[1]] ) ),
                          [GPPPInput( Symbol( "Δ"*string(i)*string(i+1) ), ColVecs( features[:, rows.d[i]]  )  ) for i in 2:length(rows.d)  ]... ),
                vcat(  data[rows.d[1],cols.p]           - data[rows.d[1],cols.s[1]], 
                      [data[rows.d[i],cols.s[i-1]] - data[rows.d[i],cols.s[i]]  for i in 2:length(cols.s)  ]... ),
                GPPPInput( Symbol( "Δ1"*string(length(cols.s)+1) ), ColVecs( features[:,rows.t]  )  ),
                data[rows.t, cols.p] - data[rows.t, cols.s[end]], 
                data[rows.t, cols.s[end]],
                1,
                identifier)                      # testing features

end




############################################################################################################
#
# Outer Loop Inference
#
############################################################################################################


# Description: step function centered at x=0
step(n) = Int64(n>0)

##################

# Description: Function to manage inference for training sets of a range of sizes and constructions according to the CAT framework. The function takes a fixed list of sizes for the C set of the primary task. The sizes of the C and A sets of the secondary tasks are determined by lists of ratios (keyword inputs, C_setsand A_sets, provided by the user) to the size of the primary C set. If there are multiple secondary tasks, it is also possible to control the fraction of the A set that is shared between the secondary tasks through another keyword input (A_fractions).

# NOTE: see bottom of page for more detail on the fields of named Tuples

# Inputs: settings     - settings structure containing inference model and method for computing prior mean
#         params       - NamedTuple with variance, lengthscale, and correlation parameters
#         data         - DataFrame with rows corresponding to different molecular systems and columns corresponding to tasks
#         features     - matrix where each column corresponds to a molecular feature--that is, columns of this matrix correspond to rows of the data DataFrame
#         guide        - NamedTuple with output filename, extra test identification data, and levels of theory to include in inference
#                        (the list will be used to select columns from the data DataFrame, so they should be given the same name)
#         core         - list of all sizes of C set to use for the primary task
#         target       - number of molecules to include in testing data set
#         tests        - list of inference methods to apply to each training data set
#                        may include: ["single", "multitask", "multitask_with_target", "delta"]
#         C_sets       - list of numbers > 0; for each entry r of this array and c of core, we will consider a secondary C set of size r*c
#         A_sets       - list of numbers > 0; for each entry r of this array and c of core, we will consider a secondary A set of size r*c
#         A_fractions  - list of numbers in [0,1]; for each entry f of this array, we consider a case where fraction f of the A set molecules 
#                        are used to train every secondary task, and the remaining molecules are split among the secondary tasks

# Outputs: Functions called within this function write inference results to a CSV with name given by the 'filename' field of guide.


function training_set_iteration( settings::settings, 
                                 params::NamedTuple, 
                                 data::DataFrame, 
                                 features::AbstractArray, 
                                 guide::NamedTuple, 
                                 core::AbstractArray,
                                 target::Number, 
                                 tests::AbstractArray; 
                                 C_sets=[0.,1.], 
                                 A_sets=[0.,0.5,1.,2.,3.,4.,5.,6.],
                                 task_fractions=[1.] )

   

    # select indices of targets
    T = sample( 1:size(data,1), Int64(ceil(target)), replace=false  )
    C = sample( setdiff(1:size(data,1),T), Int64(ceil(maximum(core))), replace=false  )
    A = sample( setdiff(1:size(data,1),vcat(C,T)),  Int64(ceil(maximum(A_sets)*maximum(core))), replace=false )

    # safeguard against redundant calculations
    task_fractions = length(guide.s) > 1 ? task_fractions : [1.]
    if C_sets == A_sets == [0.]
        cases = Iterators.product( core, C_sets, A_sets, task_fractions )
    else
        cases = Iterators.filter( x -> x[[2,3]] !== (0.,0.),
                                                           Iterators.product( core, C_sets, A_sets, task_fractions ))
    end

    # iterate over training sets
    for (cn, cRatio, aRatio, tCommon) in cases

        println((cn, cRatio, aRatio, tCommon))

        C_inds     = split_indices( cRatio*cn, length(guide.s), tCommon )
        A_inds     = split_indices( aRatio*cn, length(guide.s), tCommon ) 
        test_guide = merge( guide, ( case="C"[1:step(cRatio)]*"A"[1:step(aRatio)], shared=tCommon )  )
        rows       = ( t=T,
                       p=C[1:cn] ,
                       s=[ vcat( C[C_inds[i]], A[A_inds[i]] )  for i in 1:length(guide.s)  ]  )
        
        println(rows)
        println()

        [ test( settings, params, data, features, rows, test_guide  ) for test in tests  ]

    end
end

#################

# Description: This function creates a set of indices corresponding to each secondary tasks. Those indices indicate which molecules in the A set that particular secondary task has access to.

# Inputs: n_data  - the total number of molecules in the secondary task
#         n_tasks - the total numbed of secondary tasks
#         shared  - the fraction of the A set molecules that all secondary tasks will train on 
#                   (all other molecules will be split among the secondary tasks)

# Outputs: Vector of vectors. Each element of the outer vector is the set of A indices used for training a given secondary task.

function split_indices( n_data, n_tasks, shared )

    # get indices which all tasks will have access to
    common_ind = 1:Int64(ceil( shared*n_data ))

    inds = fill( collect(common_ind), n_tasks, 1 )
    if shared < 1
        start    = maximum(vcat(common_ind,[0])) + 1
        n_unique = Int64(round( (n_data - length(common_ind)) / n_tasks ))
        for i in 1:n_tasks-1
            inds[i] = vcat( inds[i], start:start+n_unique-1 )
            start  += n_unique 
        end
        inds[n_tasks] = vcat( inds[n_tasks], start:n_data )
    end
    return inds

end

############################################################################################################
#
# Tuple Descriptions
#
############################################################################################################

################################
# guide tuple

# fields : p          = name of primary task. This will be the column name corresponding to this task in the data DataFrame
#          s          = list of names of secondary tasks. These correspond to column names in the data data frame
#          filename   = the name (including directory path) of the CSV file where we print output
#          case       = initialized at "C". This field will be updated to include datasets used in training a particular inference model based on the CST framework
#          identifier = if extra columns are added to the output file, we can include values to print to these columns here

 
################################
# params tuple

# for multitask cases
# fields : v = Named Tuple with fields
#               p = variance hyperparameter for primary task
#               d = list of variance hyperparameters for secondary tasks (order should correspond to order of secondary parameters in guide)
#          l = Named Tuple with fields
#               p = lengthscale hyperparameter for primary task
#               d = list of lengthscale hyperparameters for secondary tasks (order should correspond to order of secondary parameters in guide)
#          ρ = list of correlation hyperparameters for secondary tasks (order should correspond to order of secondary parameters in guide)
#
# for delta cases
# fields:  μ = a vector of mean parameters for the delta between each pair of levels
#          v = a vector of variance parameters for the delta between each pair of levels
#          l = a vector of lengthscale parameters for the delta between each pair of levels


