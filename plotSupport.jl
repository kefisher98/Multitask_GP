


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

function plot_two_tasks( files; title="", datasets="CA", lowest="SCAN", markersize=8, get_cost=water_cost, yticks=[10^-2.25,10^-2,10^-1.75] )

    p = plot(; ixrotation=50, size=(900,900), xlabel="Cost [s]", ylabel="Mean Absolute Error [eV]", legend=:topright, margins=15mm)

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
    plot!(cost, ref, xaxis=:log, yaxis=:log, color=:black, label="CCSD(T) Only")


    plot!( ;title)
    plot!(minorgrid=true,minorgridalpha=0.2)
    plot!(yticks=yticks)

    return p

end

#######


function plot_three_tasks( files; title="", datasets="CA", markersize=8, get_cost=water_cost, lowest=["SCAN","PBE"], yticks=[10^-2.25,10^-2,10^-1.75])

    p = plot(; ixrotation=50, size=(900,900), xlabel="Cost [s]", ylabel="Mean Absolute Error [eV]", legend=:topright, margins=15mm)

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
    plot!(cost, ref, xaxis=:log, yaxis=:log, color=:black, label="CCSD(T) Only")


    plot!( ;title)
    plot!(minorgrid=true,minorgridalpha=0.2)
    plot!(yticks=yticks)

    return p

end


function plot_layers( files; title="", datasets="CA", markersize=8, 
                             get_cost=water_cost, lowest=["PBE","PBE","PBE"], 
                             levels=[5,3,2], c_size=[5,10,20,40,80,320],
                             yticks=[10^-0.6,10^-0.3,10^0], t_fraction=1., c_fraction=1. )

    p = plot(; ixrotation=50, size=(900,900), xlabel="Cost [s]", ylabel="Mean Absolute Error [eV]", legend=:topright, margins=15mm)

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
    plot!(cost, ref, xaxis=:log, yaxis=:log, color=:black, label="CCSD(T) Only")

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

    ref_cost  = get_cost( multifilter( CSV.read(files[1], DataFrame), [:method], [("CCSDT",)]  ) )
    ref       = read_and_avg( files, [:method], [("CCSDT",)], :, :mae_error  )
    
    p = plot( range, range, label="y=x", color=:green, linewidth=2, xlabel=datasets[2:end]*" MAE [eV]", size=(600,600), margin=5mm, legend=:bottomright )

    if include_single
        scatter!( s_only, repeat( ref[1:7], 7, 1 ), label="CCSD(T) Only", color=:black, linewidth=2, markersize=11  )
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


    p = plot(; title=lowest, ixrotation=50, size=(1200,800), xlabel="Cost [s]",  ylabel="Percent Improvement in MAE", legend=:outertopright, margins=15mm)

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

    p = plot(; ixrotation=50, size=(900,900), xlabel="Cost [s]", ylabel="Mean Absolute Error [eV]", legend=:topright, margins=15mm)

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
               label="Δ (3 levels)", color=:paleturquoise1, markerstrokecolor=:paleturquoise1; markersize, markerstrokealpha=0, markerstrokewidth=0 )

    # delta 2 level
    cost  = water_cost( multifilter( CSV.read(files[1], DataFrame), [:method,:levels,:lowest], [("Delta",2,cases[3])]  ) )
    dt    = read_and_avg( files, [:method,:levels,:lowest], [("Delta",2,cases[3])], :, :mae_error  )
    scatter!(  cost, dt,
               label="Δ (2 levels)", color=:slateblue, markerstrokecolor=:slateblue; markersize, markerstrokealpha=0, markerstrokewidth=0 )



    # plot multitask Δ
    cost  = water_cost( multifilter( CSV.read(files[1], DataFrame), [:method,:datasets,:lowest], [("Multitask","CAT",cases[2])]  ) )
    mt    = read_and_avg( files, [:method,:datasets,:lowest], [("Multitask","CAT",cases[2])], :, :mae_error  )
    scatter!(  cost, mt,
               label="Multitask Δ", color=:gold3, markerstrokecolor=:gold3; markersize, markerstrokealpha=0, markerstrokewidth=0 )


    cost  = water_cost( multifilter( CSV.read(files[1], DataFrame), [:method], [("CCSDT",)]  ) )
    ref   = read_and_avg( files, [:method], [("CCSDT",)], :, :mae_error  )
    plot!(cost, ref, xaxis=:log, yaxis=:log, color=:black, label="CCSD(T) Only")


    plot!(minorgrid=true,minorgridalpha=0.2)
    plot!(yticks=yticks)

    return p
end

######
    
