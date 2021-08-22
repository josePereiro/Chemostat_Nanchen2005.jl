import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Nanchen2006")

@time begin
    import SparseArrays
    import Base.Threads: @threads, threadid, SpinLock

    using Serialization

    import Chemostat_Nanchen2006
    const ChN = Chemostat_Nanchen2006

    const iJR = ChN.iJR904
    const Nd = ChN.NanchenData # experimental data
    const Bd = ChN.BegData    # cost data

    import Chemostat
    import Chemostat.LP: MathProgBase
    const Ch = ChN.Chemostat
    const ChU = Ch.Utils
    const ChSS = Ch.SteadyState
    const ChLP = Ch.LP
    const ChEP = Ch.MaxEntEP
    const ChSU = Ch.SimulationUtils

    using Statistics

    import UtilsJL
    const UJL = UtilsJL
    using Serialization
    using Base.Threads
    UJL.set_cache_dir(iJR.cachedir())

    using Plots
    import GR
    !isinteractive() && GR.inline("png")
end

## ----------------------------------------------------------------------------
function extract_flxs(model, epout)
    # inner flxs
    idermap = iJR.load_inners_idermap()
    iJR_iders, avs, vas = [], [], []
    for (iJR_ider, model_iders) in idermap
        # flxs
        ep_av = ChU.av(model, epout, model_iders[1])
        ep_std = sqrt(ChU.va(model, epout, model_iders[1]))
        if length(model_iders) == 2 # reversible
            # r = r+ - r-
            ep_av -= ChU.av(model, epout, model_iders[2])
            # ep_std += sqrt(ChU.va(model, epout, model_iders[2]))
            ep_std = NaN
        end
        push!(avs, ep_av)
        push!(vas, ep_std)
        push!(iJR_iders, iJR_ider)
    end
    return iJR_iders, avs, vas
end

## ----------------------------------------------------------------------------
let
    model = iJR.load_model("max_model")
    exglc_ider = iJR.EX_GLC_IDER
    Nd_rxns_map = iJR.load_rxns_map() # inner reacts
    iJR_ider_subs = iJR.load_inner_rxns_subs()

    fpath = iJR.procdir("MSE_study_params.jls")
    exp, method, exp_biom_beta, exp_vg_beta, biom_betas, vg_betas = deserialize(fpath)

    p = plot(;xlabel = "vg_beta", ylabel = "MSE")
    for biom_beta in biom_betas

        fname = UJL.mysavename("MSE_study_epouts", "jls"; biom_beta)
        fpath = iJR.procdir(fname)
        !isfile(fpath) && continue
        epouts = deserialize(fpath)

        MSEs = []
        for vg_beta in vg_betas
            
            epout = epouts[vg_beta]
            exp_exglc = Nd.uval("GLC", exp)

            # iJR_iders, avs, vas = extract_flxs(model, epout)

            sum = 0.0
            iter = zip(extract_flxs(model, epout)...)
            for (iJR_ider, iJR_av, iJR_va) in iter
                Nd_ider = Nd_rxns_map[iJR_ider]
                exp_av = Nd.val(Nd_ider, exp)
                sum += (abs(iJR_av/exp_exglc) - abs(exp_av/exp_exglc))^2
            end
            MSE = sum / length(iter)
            push!(MSEs, MSE)
        end # for vg_beta

        plot!(p, log10.(vg_betas), log10.(MSEs); label = "")
    end # for biom_beta
    p
end