import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Nanchen2006")

@time begin

    import Chemostat_Nanchen2006
    const ChN = Chemostat_Nanchen2006

    const Nd = ChN.NanchenData # experimental data

    import UtilsJL
    const UJL = UtilsJL

    using Plots
end

## -------------------------------------------------------------------
fileid = "1.1"
function mysavefig(p, pname; params...) 
    fname = UJL.mysavefig(p, string(fileid, "_", pname), Nd.plotsdir(); params...)
    @info "Plotting" fname
end

## -------------------------------------------------------------------
# D vs X
let
    p = plot(; title = "Kayser", 
        xlabel = string("X (", Nd.unit(:X), ")"), 
        ylabel = string("D (", Nd.unit(:D), ")")
    )
    plot!(p, Nd.val(:D), Nd.val(:X); 
        label = "", ls = :dash, lw = 3, color = :black, 
        alpha = 0.6
    )
    scatter!(p, Nd.val(:D), Nd.val(:X); 
        label = "", m = 8, color = :black
    )
    mysavefig(p, "X_vs_D") 
end

## -------------------------------------------------------------------
# BALANCE
let
    ps = Plots.Plot[]
    for met in [:GLC]
        p = plot(; title = string("Balance: ", met), 
            xlabel = "feed", 
            ylabel = "exch + drain" 
        )

        exps = 1:17

        feed = Nd.cval.(met, exps) .* Nd.val.(:D, exps)
        exch = Nd.uval.(met, exps) .* Nd.val.(:X, exps) .|> abs
        drain = zeros(length(exch)) # This experiments hase sg == 0

        scatter!(p, feed, exch .+ drain; 
            label = "", m = 8, color = :black
        )
        vals = [feed; exch .+ drain] |> sort
        plot!(p, vals, vals;
            label = "", ls = :dash, lw = 3, alpha = 0.6, color = :black
        )
        push!(ps, p)
    end
    mysavefig(ps, "Balances") 
end