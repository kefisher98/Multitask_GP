


# =================================================================
# Cost Model
# =================================================================

function water_cost( df )
    return df[!,:high_level]*(24550.3) + df[!,:low_level]*(665.0)
end

function organic_cost( df )
    return df[!,:high_level]*(13288.3) + df[!,:low_level]*(54.5)
end


# =================================================================
# Read Results for Visualization
# =================================================================


function read_and_avg( files, filtercols, filtervals, skiprows, avgcol  )

    results = multifilter( CSV.read(files[1],DataFrame), filtercols, filtervals )[skiprows,avgcol]
    for f in files[2:end]
        results += multifilter( CSV.read(f,DataFrame), filtercols, filtervals )[skiprows,avgcol]
    end

    return results ./ length(files)
end

#########

function read_and_std( files, filtercols, filtervals, skiprows, avgcol  )

    results = multifilter( CSV.read(files[1],DataFrame), filtercols, filtervals )[skiprows,avgcol]
    for f in files[2:end]
        results = hcat( results, multifilter( CSV.read(f,DataFrame), filtercols, filtervals )[skiprows,avgcol] )
    end

    return std( results, dims=2  )
end

#######

function get_low_level_size( c_size; dims, code, max_ratio=6 )

    num = 0
    if cmp( code[end:end], "T" )==0 
        num += 1 
    end
        
    if cmp( code[1:1], "C" ) == 0  
        # shared tasks
        num += Int64(ceil( c_size * dims.t_shared  )) * dims.tasks

        # remainder
        num += c_size  - Int64(ceil( c_size * dims.t_shared ))
    end
    
    if cmp( code[1:1], "A" )==0  || cmp( code[2:2],"A" )==0 

        # shared tasks
        num += Int64(ceil( c_size * dims.t_shared * max_ratio )) * dims.tasks

        # remainder
        num += Int64(ceil( c_size * max_ratio )) - Int64(ceil( c_size * dims.t_shared * max_ratio ))
    end

    return num
end


function report_differences( files, lowest, datasets )
    
    # read in data
    ref   = read_and_avg( files, [:method], [("CCSDT",)], 1:7, :mae_error  )
    mt_1  = read_and_avg( files, [:method,:datasets,:levels,:lowest], [("Multitask",datasets[1],2,lowest)]  , :, :mae_error  )
    mt_2  = read_and_avg( files, [:method,:datasets,:levels,:lowest], [("Multitask",datasets[2],2,lowest)]  , :, :mae_error  )

    # tile reference for comparison
    ref   = repeat( ref, 7,1)

    return ref-mt_1, ref-mt_2, mt_1-mt_2
end

#=

 Single to A: 0.0005  - 0.006 (median: 0.004)
 A to CA:     -3.9e-5 - 0.002 (median: 0.0007)

 Single to AT: 0.002 - 0.01   (median: 0.006)
 AT to CAT:    0.0006 - 0.007 (median: 0.002)

=#

# =================================================================
# Visualize: water or organic
# =================================================================

function plot_two_tasks( files; 
                         title="", datasets="CA", lowest="SCAN", arrow=false,
                         markersize=8, get_cost=water_cost, yticks=[10^-2.25,10^-2,10^-1.75] )

    default(;fontfamily="Computer Modern");
    p = plot(; ixrotation=50, size=(750,750), 
               xlabel="Cost [s]", ylabel="Mean Absolute Error [eV]", 
               legend=:topright, margins=15mm)

    # no target
    cost  = get_cost( multifilter( CSV.read(files[1], DataFrame), [:method,:datasets,:levels,:lowest], [("Multitask",datasets,2,lowest)]  ) )
    mt    = read_and_avg( files, [:method,:datasets,:levels,:lowest], [("Multitask",datasets,2,lowest)]  , :, :mae_error  ) 
    scatter!(  cost, mt,
               label=datasets, color=:slateblue, markerstrokecolor=:slateblue; markersize, markerstrokealpha=0, markerstrokewidth=0 )


    # target
    cost  = get_cost( multifilter( CSV.read(files[1], DataFrame), [:method,:datasets,:levels,:lowest], [("Multitask",datasets*"T",2,lowest)]  ) )
    mt    = read_and_avg( files, [:method,:datasets,:levels,:lowest], [("Multitask",datasets*"T",2,lowest)]  , :, :mae_error  ) 
    scatter!(  cost, mt,
               label=datasets*"T", color=:gold3, markerstrokecolor=:gold3; markersize, markerstrokealpha=0, markerstrokewidth=0 )

    # reference
    cost  = get_cost( multifilter( CSV.read(files[1], DataFrame), [:method], [("CCSDT",)]  ) )
    ref   = read_and_avg( files, [:method], [("CCSDT",)], :, :mae_error  )
    plot!(cost, ref, xaxis=:log, yaxis=:log, color=:black, label="CCSD(T)")

    # optional arrow
    if arrow
        c = :lightcoral
        x1 = 2.8e6 
        x2 = 1.6e7
        y  = 0.0072
        plot!( [x1,x2], [y,y]; linewidth=6,c,label=false, alpha=0.8 )
        scatter!( [x1], [y]; markershape=:ltriangle, markerstrokecolor=c,linewidth=4,c,label=false, markersize=13 )
        scatter!( [x2], [y]; markershape=:rtriangle, markerstrokecolor=c,linewidth=4,c,label=false, markersize=13 )
    end
    plot!( ;title)
    plot!(minorgrid=true,minorgridalpha=0.2)
    plot!(yticks=yticks)

    return p

end

#######


function plot_three_tasks( files; title="", datasets="CA", markersize=8, get_cost=water_cost, lowest=["SCAN","PBE"], yticks=[10^-2.25,10^-2,10^-1.75])

    default(;fontfamily="Computer Modern");
    p = plot(; ixrotation=50, size=(750,750), 
               xlabel="Cost [s]", ylabel="Mean Absolute Error [eV]", legend=:topright, margins=15mm)

    # plot multitask 3 levels
    colors = [:slateblue, :turquoise, :skyblue]
    if datasets in ["C","A","CA"]

        cost  = get_cost( multifilter( CSV.read(files[1], DataFrame), [:method,:datasets,:levels,:lowest], [("Multitask",datasets,3,lowest[1])]  ) )
        mt    = read_and_avg( files, [:method,:datasets,:levels,:lowest], [("Multitask",datasets,3,lowest[1])]  , :, :mae_error  )
        scatter!(  cost, mt,
                   label="3 levels", color=colors[1], markerstrokecolor=colors[1]; markersize, markerstrokealpha=0, markerstrokewidth=0 )

    else
        for (i,lw) in enumerate( lowest)
            cost  = get_cost( multifilter( CSV.read(files[1], DataFrame), [:method,:datasets,:levels,:lowest], [("Multitask",datasets,3,lw)]  ) )
            mt    = read_and_avg( files, [:method,:datasets,:levels,:lowest], [("Multitask",datasets,3,lw)]  , :, :mae_error  ) 
            scatter!(  cost, mt,
                       label="3 levels ("*lw*")", color=colors[i], markerstrokecolor=colors[i]; markersize, markerstrokealpha=0, markerstrokewidth=0 )
        end
    end



    # plot multitask 2 levels 
    colors = [:gold3, :darkorange4, :antiquewhite]
    for (i,lw) in enumerate(lowest)
        cost  = get_cost( multifilter( CSV.read(files[1], DataFrame), [:method,:datasets,:levels,:lowest], [("Multitask",datasets,2,lw)]  ) )
        mt    = read_and_avg( files, [:method,:datasets,:levels,:lowest], [("Multitask",datasets,2,lw)]  , :, :mae_error  ) 
        scatter!(  cost, mt,
                   label="2 levels ("*lw*")", color=colors[i], markerstrokecolor=colors[i]; markersize, markerstrokealpha=0, markerstrokewidth=0 )
    end
    
    # reference
    cost  = get_cost( multifilter( CSV.read(files[1], DataFrame), [:method], [("CCSDT",)]  ) )
    ref   = read_and_avg( files, [:method], [("CCSDT",)], :, :mae_error  )
    plot!(cost, ref, xaxis=:log, yaxis=:log, color=:black, label="CCSD(T)")


    plot!( ;title)
    plot!(minorgrid=true,minorgridalpha=0.2)
    plot!(yticks=yticks)

    return p

end


function plot_layers( files; title="", datasets="CA", markersize=8, 
                             get_cost=water_cost, lowest=["PBE","PBE","PBE"], 
                             levels=[5,3,2], c_size=[5,10,20,40,80,320],
                             yticks=[10^-0.6,10^-0.3,10^(0)], t_fraction=1., c_fraction=1. )

    default(;fontfamily="Computer Modern");
    p = plot(; ixrotation=50, size=(750,750), xlabel="Cost [s]", ylabel="Mean Absolute Error [eV]", legend=:topright, margins=15mm)

    markershape = :star4
    colors      = [:slateblue,:palegreen,:gold3,:teal]

    for (i,lv) in enumerate(levels)

        t_shared = lv > 2 ? t_fraction : 1.

        # scatter plot
        cost  = get_cost( multifilter( CSV.read(files[1], DataFrame), [:method,:datasets,:levels,:task_fractions,:lowest], [("Multitask",datasets,lv,t_shared,lowest[i])]  ) )
        mt    = read_and_avg( files, [:method,:datasets,:levels,:task_fractions,:lowest], [("Multitask",datasets,lv,t_shared,lowest[i])]  , :, :mae_error  )
        scatter!( cost, mt,
                        color=colors[i], markerstrokecolor=colors[i], label=string(lv)*" Tasks";
                        markersize, markerstrokewidth=0, markerstrokealpha=0, markeralpha=0.7, markershape )
        
        # ribbon plot
        a_size_low  = get_low_level_size.( c_size; dims=(tasks=lv-1, c_fraction=c_fraction, t_shared=t_shared), code=datasets, max_ratio=1 )
        a_size_high = get_low_level_size.( c_size; dims=(tasks=lv-1, c_fraction=c_fraction, t_shared=t_shared), code=datasets, max_ratio=6 )


        low_sets    = [ ("Multitask", datasets, lv, lowest[i], t_shared, c, a)   for (c,a) in Iterators.zip(c_size, a_size_low) ]
        high_sets   = [ ("Multitask", datasets, lv, lowest[i], t_shared, c, a)   for (c,a) in Iterators.zip(c_size, a_size_high) ]


        cost        = get_cost( multifilter( CSV.read(files[1], DataFrame), [:method,:datasets,:levels,:lowest,:task_fractions,:high_level,:low_level], low_sets ))
        mt_low      = read_and_avg( files, [:method,:datasets,:levels,:lowest,:task_fractions,:high_level,:low_level], low_sets, :, :mae_error )
        mt_high     = read_and_avg( files, [:method,:datasets,:levels,:lowest,:task_fractions,:high_level,:low_level], high_sets, :, :mae_error )



        plot!( cost, mt_low,
               fillrange=(mt_low,mt_high), fillalpha=0.3,
               fillcolor=colors[i], color=colors[i],
               alpha=0.4, label=false)

    end

    # reference
    cost  = get_cost( multifilter( CSV.read(files[1], DataFrame), [:method], [("CCSDT",)]  ) )
    ref   = read_and_avg( files, [:method], [("CCSDT",)], :, :mae_error  )
    plot!(cost, ref, xaxis=:log, yaxis=:log, color=:black, label="CCSD(T)")

    plot!( ;title)
    plot!(minorgrid=true,minorgridalpha=0.2)
    plot!(yticks=yticks, ylims=[10^-0.75,10^0])

    return p
end




function compare_S_power( files; datasets="CAT", lowest="PBE", get_cost=water_cost, include_single, range=[0.,0.1] )

    markershape = :circle

    # read and average test results
    s_only_cost   = get_cost( multifilter( CSV.read(files[1], DataFrame), [:method,:datasets,:levels,:lowest], [("Multitask",datasets[2:end],2,lowest)]  ) )
    
    s_only = read_and_avg( files, [:method,:datasets,:levels,:lowest], [("Multitask",datasets[2:end],2,lowest)]  , :, :mae_error  )
    with_c = read_and_avg( files, [:method,:datasets,:levels,:lowest], [("Multitask",datasets,2,lowest)]  , :, :mae_error  )

    default(;fontfamily="Computer Modern");
    ref_cost  = get_cost( multifilter( CSV.read(files[1], DataFrame), [:method], [("CCSDT",)]  ) )
    ref       = read_and_avg( files, [:method], [("CCSDT",)], :, :mae_error  )
    
    p = plot( range, range, label="y=x", color=:green, linewidth=2, xlabel=datasets[2:end]*" MAE [eV]", size=(750,750), margin=5mm, legend=:bottomright )

    if include_single
        scatter!( s_only, repeat( ref[1:7], 7, 1 ), label="CCSD(T)", color=:black, linewidth=2, markersize=11  )
        label  = datasets
        ylabel = "MAE [eV]"
    else
        label  = false
        ylabel = datasets*" MAE [eV]"
    end
    scatter!( s_only, with_c,
              color=:slateblue, markerstrokecolor=:slateblue;
              ylabel, label,
              markersize=11, markerstrokewidth=0, markerstrokealpha=0, markeralpha=0.7, markershape )

    return p
end



# =================================================================
# Visualize: water or organic
# =================================================================


function level_difference( files; lowest, markersize=8)


    colors=[:navyblue,
            :antiquewhite,
            :slateblue,
            :gold3,
            :deepskyblue,
            :darkorange4,
            :paleturquoise1]


    default(;fontfamily="Computer Modern");
    p = plot(; title=lowest, ixrotation=50, size=(750,750), 
               xlabel="Cost [s]",  ylabel="Percent Improvement in MAE", legend=:outertopright, margins=15mm)

    for (i,ds) in enumerate(["C","CT","A","AT","CA","CAT"])
        for Ct in [5,10,20,40,80,160,320]

            inds  = ds in ["C","CT"] ? collect(1:1) : collect(1:3)
            n_dft = read_and_avg( files, [:method, :lowest, :datasets, :levels, :high_level], [("Multitask",lowest, ds, 2, Ct)] , :, :low_level  )[inds]
            l2    = read_and_avg( files, [:method, :lowest, :datasets, :levels, :high_level], [("Multitask",lowest, ds, 2, Ct)] , :, :mae_error  )[inds]

            inds  = ds in ["C","CT"] ? collect(1:1)  : [3,4,6]
            low3  = ds in ["C","A","CA"] ? "SCAN" : lowest
            l3    = read_and_avg( files, [:method, :lowest, :datasets, :levels, :high_level], [("Multitask",low3, ds, 3, Ct)] , :, :mae_error  )[inds]

            lb    = Ct==5 ? ds : false
            scatter!(  n_dft, 100*(l2-l3) ./ l2,
                       label=lb, color=colors[i], markerstrokecolor=:black;
                       markersize, markerstrokewidth=0.5 )
        end
    end
    hline!([0.], color=:black, label=false,xaxis=:log)

    plot!(minorgrid=true,minorgridalpha=0.2)

    return p

end

########

function plot_delta( files; quantity="CCSD(T)-PBE", multitask=true, markersize=8, yticks=[10^-2.5,10^-2,10^-1.5])

    # assign filtering parameters so that plots are consistent
    # (consistent in the sense that they try to estimate the same quantity with the same DFT target data)
    if cmp(quantity, "CCSD(T)-PBE")==0
        cases = ["PBE","SCAN_PBE","PBE","PBE"]
    elseif cmp(quantity, "CCSD(T)-SCAN")==0
        cases = ["SCAN","PBE_SCAN","SCAN","SCAN"]
    else
        error("Unsupported input for keyword quantity. Use 'CCSD(T)-PBE' or 'CCSD(T)-SCAN'.")
    end

    default(;fontfamily="Computer Modern");
    p = plot(; ixrotation=50, size=(750,750), xlabel="Cost [s]", ylabel="Mean Absolute Error [eV]", legend=:topright, margins=15mm, ylims=(0.0009,0.035))

    # plot multitask
    if multitask
        cost  = water_cost( multifilter( CSV.read(files[1], DataFrame), [:method,:datasets,:levels,:lowest], [("Multitask","CAT",3,cases[1])]  ) )
        mt    = read_and_avg( files, [:method,:datasets,:levels,:lowest], [("Multitask","CAT",3,cases[1])]  , :, :mae_error  )
        scatter!(  cost, mt,
                   label="Multitask", color=:darkorange4, markerstrokecolor=:darkorange4; markersize, markerstrokealpha=0, markerstrokewidth=0 )
    end
    
    # delta 3 level
    cost  = water_cost( multifilter( CSV.read(files[1], DataFrame), [:method,:levels,:datasets,:lowest], [("Delta",3,"CAT",cases[4])]  ) )
    dt    = read_and_avg( files, [:method,:levels,:datasets,:lowest], [("Delta",3,"CAT",cases[4])], :, :mae_error  )
    scatter!(  cost, dt,
               label="Δ hML", color=:paleturquoise1, markerstrokecolor=:paleturquoise1; markersize, markerstrokealpha=0, markerstrokewidth=0 )

    # delta 2 level
    cost  = water_cost( multifilter( CSV.read(files[1], DataFrame), [:method,:levels,:lowest], [("Delta",2,cases[3])]  ) )
    dt    = read_and_avg( files, [:method,:levels,:lowest], [("Delta",2,cases[3])], :, :mae_error  )
    scatter!(  cost, dt,
               label="Δ-learning", color=:slateblue, markerstrokecolor=:slateblue; markersize, markerstrokealpha=0, markerstrokewidth=0 )



    # plot multitask Δ
    cost  = water_cost( multifilter( CSV.read(files[1], DataFrame), [:method,:datasets,:lowest], [("Multitask","CAT",cases[2])]  ) )
    mt    = read_and_avg( files, [:method,:datasets,:lowest], [("Multitask","CAT",cases[2])], :, :mae_error  )
    scatter!(  cost, mt,
               label="Multitask Δ", color=:gold3, markerstrokecolor=:gold3; markersize, markerstrokealpha=0, markerstrokewidth=0 )


    cost  = water_cost( multifilter( CSV.read(files[1], DataFrame), [:method], [("CCSDT",)]  ) )
    ref   = read_and_avg( files, [:method], [("CCSDT",)], :, :mae_error  )
    plot!(cost, ref, xaxis=:log, yaxis=:log, color=:black, label="CCSD(T)")


    plot!(minorgrid=true,minorgridalpha=0.2)
    plot!(yticks=yticks)

    return p
end

######

# ===================================================================================
# Calibration
# ===================================================================================

# keys has fields truth, mu, sigma
# rows has fields calibrate, target
function conform( data, keys, rows; α=0.05  )
    # create scores
    scores = sort( abs.( data[rows.calibrate,keys.mu] - data[rows.calibrate,keys.truth]  )./ data[rows.calibrate,keys.sigma] )
    
    # find quantile
    n      = length(rows.calibrate)
    q_hat  = scores[ Int64(ceil( (n+1)*(1-α) )) ]

    # return dataframe with rescaled distributions
    return DataFrame(  hcat( Matrix(data[rows.target,[keys.mu,keys.sigma,keys.truth]]), q_hat*data[rows.target,keys.sigma]   ), 
                      [keys.mu,keys.sigma,keys.truth,:conformal] )
end

# scatter plot
#   original error+ std and adjusted 
function error_std( data, keys; color=:cornflowerblue,markerstrokecolor=:slateblue, markersize=3  )

    ymax  = maximum( data[:,:conformal]  ) 
    ylims = (0,ymax+0.03)

    # original
    s1 = plot( ; size=(300,300), xlabel="absolute error", ylabel="standard deviation" )
    scatter!( abs.(data[:,keys.truth] - data[:,keys.mu]), data[:,keys.sigma]; 
              color, markerstrokecolor, markersize, label=false, ylims, title="Initial Prediction"  )

    # adjusted
    s2 = plot( ; size=(300,300), xlabel="absolute error", ylabel="standard deviation" )
    scatter!( abs.(data[:,keys.truth] - data[:,keys.mu]), data[:,:conformal]; 
              color, markerstrokecolor, markersize, label=false, ylims, title="Conformalized"   )

    # together
    s  = plot!( s1, s2, layout=(1,2), size=(600,300)  )

    return s
end



# original and adjusted intervals for a subset of cases
function intervals( data, keys; nplot=20, ylabel="QOI [eV]", adjust=0.1, ylims=(6,13)  )

    # calculate coverage
    before =  sum( ( data[:,keys.truth] .> data[:,keys.mu] - data[:,keys.sigma]  ) .* 
                   ( data[:,keys.truth] .< data[:,keys.mu] + data[:,keys.sigma]  )     ) / size(data,1)
    after  =  sum( ( data[:,keys.truth] .> data[:,keys.mu] - data[:,:conformal]  ) .* 
                   ( data[:,keys.truth] .< data[:,keys.mu] + data[:,:conformal]  )     ) / size(data,1)
    
    # initialize plot, coverage in title
    fig = plot( ; size=(400,300), 
                  title="Intial Coverage: "*string(before)*"\nConformalized Coverage: "*string(after), 
                  xlabel="molecule index", 
                  ylabel,
                  ylims )

    
    for n = 1:nplot
        
        # line
        vline!([n]; color=:gray, label=false, alpha=0.2)

        # final intervals
        label = n==1 ? "conformalized" : false
        plot!([n-adjust,n+adjust], data[n,keys.mu]*[1,1],ribbon=data[n,:conformal]*[1,1]; label, color=:cornflowerblue , alpha=1 )
    
        # original intervals
        label = n==1 ? "inital interval" : false
        plot!([n-adjust,n+adjust], data[n,keys.mu]*[1,1],ribbon=data[n,keys.sigma]; label, color=:darkorange3 , alpha=0.5 )

        # truth
        label = n==1 ? "truth" : false
        scatter!( [n], [data[n,keys.truth]]; label, color=:forestgreen, markerstrokecolor=:forestgreen, markershape=:hline  )
    end
    return fig
end

# =====================================================================================
# Hyperparameter Tests
# =====================================================================================

function plot_errors(file, N_vec;
                     tasks=["CCSDT","SCAN"],
                     titles=["CCSD(T)","SCAN-ρ[CCSD(T)]"],
                     y=[0.04,0.004],
                     colors=[:slateblue,
                             :gold3,
                             :deepskyblue,
                             :darkorange4,
                             :paleturquoise1] )

    data = CSV.read( file, DataFrame  )

    figs = []
    for (j,t) in enumerate(tasks)
        default(;fontfamily="Computer Modern");
        fig    = plot( [0.0, y[j]], [0.0, y[j]], seriestype = :straightline,
                       color=:gray, linewidth=2, label="y=x", size=(750,750),
                       xlabel="Test MAE for MLE objective", ylabel="Test MAE for MAE objective" )

        for (i,N) in enumerate(N_vec)
            select =  multifilter( data, [:task,:N], [(t,N)] )
            scatter!( select[:,:mle_mae], select[:,:mae_mae],
                      color=colors[i], markerstrokecolor=:black, alpha=0.7, label="M="*string(N), markersize=8  )
        end
        push!(figs,fig)
    end
    return figs
end



# =====================================================================================
# Heat Maps
# =====================================================================================

function SOAPheat(path; colname="MAE ip GP, Average Features", title="Error", 
                        nmax=[6,8,10,12],
                        lmax=[3,6,8],
                        sigma=[0.3,0.4,0.5],
                        r_cut=[3,4,5],
                        color=:blues, ctitle="", size=(1000,800), clims=nothing)

    # read data
    data = CSV.read(path,DataFrame)

    # dimensions per map
    n_num = length(nmax)
    l_num = length(lmax)

    # plot settings
    default(;fontfamily="Computer Modern");
    hmaps    = []
    subtitle = vcat(["σ="*string(s) for s in sigma ], fill("",6,1))
    ytitle   = collect( Iterators.flatten( [ ["r cutoff ="*string(r)*"\n\n", "", ""]  for r in r_cut] ) )
    clims    = isnothing(clims) ? (minimum(data[!,colname]), maximum(data[!,colname])) : clims

    # iterate through SOAP parameters, reframe, make heatmap
    loc = 1
    for r in r_cut
        for σ in sigma
            fil = multifilter(data,[:r_cut,:sigma],[(r,σ)])
            mat = reshape(fil[!,colname], l_num, n_num )
            push!(hmaps, heatmap(nmax,lmax,mat; color, colorbar=false, clims,
                                 title=subtitle[loc], xlabel="n max", ylabel=ytitle[loc]*"ℓ max"))
            loc+=1
        end
    end

    # add colorbar
    cb = scatter([0,0], [0,1], zcolor=[0,3], clims=clims,
                 xlims=(1,1.1), label="", c=color, colorbar_title=ctitle, framestyle=:none)

    # combined figure
    l    = @layout [grid(3,3) a{0.035w}]
    all  = plot(hmaps[1], hmaps[2], hmaps[3], hmaps[4], hmaps[5], hmaps[6], hmaps[7], hmaps[8], hmaps[9], cb,
                layout=l, link=:all, margins=2.5mm, size=size  )

    return all
end

# =====================================================================================
# Kernel Distances
# =====================================================================================

function getNonHydrogenPositions(at)

    pos = at.get_positions()[ at.get_chemical_symbols() .!= "H"  ,:]
    return [pos[r,:] for r in 1:size(pos,1)]

end


######


function getAvgFeatures(desc,ats; settings=nothing)

    # compute features for each atom
    features  =  [mean( desc.create(at,  positions=getNonHydrogenPositions(at)); dims=1 )  for at in ats[:,1]]

    # reformat as Matrix with each atom's global features corresponding to a column
    fMat      = reshape(collect(Iterators.flatten(features)), (length(features[1]),length(features)))

    # format as input to GPPP
    return fMat


end


########


function getAvgFeatures(desc,ats,species::NamedTuple; avgApproach=avgBySpecies, settings=nothing)

    # compute features for each atom
    features  =  [avgApproach(desc.create(at,  positions=getNonHydrogenPositions(at)), at, desc, species)  for at in ats]


    # reformat as Matrix with each atom's global features corresponding to a column
    fMat      = reshape(collect(Iterators.flatten(features)), (length(features[1]),length(features)))

    # format as input to GPPP
    return fMat


end

#########

getNormedFeatures(desc,ats) =  [normalize(desc.create(at,  positions=getNonHydrogenPositions(at)) )  for at in ats]

#########


function avgBySpecies( feature, at, desc, species )


    # average over species the flatten
    avg_feature = zeros(size(species.rollcall,1),size(feature,2))


    for (ind,specie) in enumerate(species.rollcall)

        # remove Hydrogen from the symbol list
        symbols = at.get_chemical_symbols()[at.get_chemical_symbols().!="H"]

        # average feature where (if anywhere) this specie is in the atom
        avg_feature[ind,:] = mean( feature[symbols .=== specie,:]; dims=1 )

        # input solitary atom if species was not present in main compound
        if sum(symbols .=== specie)===0
            avg_feature[ind,:] = desc.create( species.solitary[ind], positions=getNonHydrogenPositions(species.solitary[ind]) )
        end
    end

    return collect(Iterators.flatten(avg_feature'))

end

#######


function KernelDistanceTriple(at_x,at_y,desc,location; degree=2., species=["C","N","O","H"] )

    # information about species under consideration
    rc                                = species[species .!= "H"]
    avgBySpeciesInfo                  = (rollcall=rc, solitary=[ase.read(location*"isolated_molecules/"*at*".xyz") for at in rc])


    # average
    x      = getAvgFeatures(desc,[at_x])
    y      = getAvgFeatures(desc,[at_y])
    avg    = dsk.AverageKernel(metric="linear")
    k_avg  = (avg.create( [x',y'] )[1,2]).^degree

    # average by species
    x      = getAvgFeatures(desc,[at_x],avgBySpeciesInfo)
    y      = getAvgFeatures(desc,[at_y],avgBySpeciesInfo)
    k_spe  = (avg.create( [x',y'] )[1,2]).^degree

    # REMatch
    x      = getNormedFeatures(desc,[at_x])[1]
    y      = getNormedFeatures(desc,[at_y])[1]

    rem    = dsk.REMatchKernel(metric="linear", alpha=0.5, threshold=1e-6)
    k_rem  = (rem.create( [x,y] )[1,2]).^degree

    return (sqrt(2 - 2k_avg), sqrt(2 - 2k_spe), sqrt(2 - 2k_rem))
end


function distancePairs(ats, location; rcut=3, nmax=3, lmax=3, sigma=0.5, degree=2.0, species=["C","N","O","H"])

    # define SOAP
    desc = dsd.SOAP(species=species,periodic=false,rcut=rcut,nmax=nmax,lmax=lmax,sigma=sigma)


    # distance pairs
    dists = zeros(binomial(size(ats,1),2),3)
    ind   = 1

    for i in 1:size(ats,1)
        for j in i+1:size(ats,1)
            dists[ind,:]  = collect( KernelDistanceTriple(ats[i],ats[j],desc,location; degree, species) )
            ind          +=1
        end
    end

    return dists

end




function scatterDistance(dists; title="", size=(800,400),color=:red)

    # limits
    lims = (0, maximum(dists)+0.05*maximum(dists) )

    default(;fontfamily="Computer Modern");
    # first pair
    ρ_1 = cor(dists[:,1],dists[:,2])
    s_1 = scatter(dists[:,1],dists[:,2],
                  xlabel="Average",ylabel="Average By Species",
                  markerstrokecolor=color,
                  margin=3mm,
                  #title="\nρ="*string(round(ρ_1,digits=4)),
                  xlims=lims,ylims=lims,
                  markersize=0.5,markeralpha=0.5,markerstrokealpha=0,color=color,
                  aspect_ratio=:equal,legend=false)


    # Second pair
    ρ_2 = cor(dists[:,1],dists[:,3])
    s_2 = scatter(dists[:,1],dists[:,3],
                  xlabel="Average",ylabel="REMatch",
                  margin=3mm,
                  markerstrokecolor=color,
                  #title=title*"\nρ="*string(round(ρ_2,digits=4)),
                  xlims=lims,ylims=lims,
                  markersize=0.5,markeralpha=0.5,markerstrokealpha=0, color=color,
                  aspect_ratio=:equal, legend=false)



    # Third pair
    ρ_3 = cor(dists[:,2],dists[:,3])
    s_3 = scatter(dists[:,2],dists[:,3],
                  xlabel="Average By Species",ylabel="REMatch",
                  margin=3mm,
                  markerstrokecolor=color,
                  #title="\nρ="*string(round(ρ_3,digits=4)),
                  xlims=lims,ylims=lims,
                  markersize=0.5,markeralpha=0.5,markerstrokealpha=0, color=color,
                  aspect_ratio=:equal, legend=false)


    s_all = plot(s_1,s_2,s_3,layout=(1,3),size=size, margin=3mm)

    return s_all, s_1, s_2, s_3, ρ_1, ρ_2, ρ_3
end

