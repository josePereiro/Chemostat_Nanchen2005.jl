import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Nanchen2006")

@time begin
    using Serialization

    import Chemostat_Nanchen2006
    const ChN = Chemostat_Nanchen2006
    const iJR = ChN.iJR904

    using Plots
    import GR
    GR.inline("png")

    import UtilsJL
    const UJL = UtilsJL
    using Serialization
    using Base.Threads

end

## ----------------------------------------------------------------------------
fileid = "2.1"
mysavefig(p, pname; params...) = 
    UJL.mysavefig(p, string(fileid, "_", pname), iJR.plotsdir(); params...)

## ----------------------------------------------------------------------------
let
    mon = UJL.OnDiskMonitor(iJR.cachedir(), "monitor.jld2")
    live_proves = Dict()

    UJL.watch(mon; wt = 15.0) do ddat

        isempty(ddat) && return
        
        for (exp, tdat) in ddat
            method = get(tdat, :method, "")
            
            # Check activity
            old_lprove = get!(live_proves, exp, -1.0)
            new_lprove = get(tdat, :live_prove, -1.0)
            old_lprove == new_lprove && continue
            live_proves[:exp] = new_lprove

            for datk in [:round, :gd]
                kdat = get!(tdat, datk, Dict())
                ps = Plots.Plot[]

                # means
                for (avk, limk) in [
                            (:vg_avPME, :cgD_X), 
                            (:biom_avPME, :exp_growth), 
                        ]
                    avdat = get(kdat, avk, [])
                    p = plot(avdat;  title = string(avk), 
                        xlabel = "iter", ylabel = string(avk), 
                        lw = 3, label = string(avk)
                    )
                    limdat = get(tdat, limk, 0.0)
                    hline!(p, [limdat]; label = string(limk), 
                        lw = 3, ls = :dash, color = :black, 
                    )
                    plot!(p; legend = :topleft)
                    push!(ps, p)
                end
                
                # betas
                for bk in [:vg_beta, :biom_beta]
                    dat = get(kdat, bk, [])
                    p = plot(dat;  title = string(bk), 
                        xlabel = "iter", ylabel = string(bk), 
                        lw = 3, label = string(bk)
                    )
                    plot!(p; legend = :topleft)
                    push!(ps, p)
                end

                mysavefig(ps, "monitor"; datk, exp, method)
                @info("Done", exp, datk)
            end
        end
    end
end