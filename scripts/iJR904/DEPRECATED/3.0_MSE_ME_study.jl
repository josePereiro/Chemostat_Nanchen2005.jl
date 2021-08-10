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
end

## ----------------------------------------------------------------------------
# ME data
DAT_FILE_PREFFIX =  "maxent_ep_dat"
function dat_file(;kwargs...)
    fname = UJL.mysavename(DAT_FILE_PREFFIX, "jls"; kwargs...)
    iJR.procdir(fname)
end

## ----------------------------------------------------------------------------
let
    # ME dat
    exp = 1
    method = :ME_MAX_POL
    datfile = dat_file(;method, exp)
    # !isfile(datfile) && continue
    
    dat = deserialize(datfile)
    
    ## ----------------------------------------------------------------------------
    # betas
    base = 1.3
    offset = 2.0
    bins = 100
    exp_biom_beta, exp_vg_beta = dat[:exp_beta]

    biom_log_range = range(0.0, log(base, exp_biom_beta * offset); length = bins)
    biom_betas = [0.0; base .^ (biom_log_range); exp_biom_beta] |> unique! |> sort!
    println() 

    vg_log_range = range(0.0, log(base, exp_vg_beta * offset); length = bins)
    vg_betas = [0.0; base .^ (vg_log_range); exp_vg_beta] |> unique! |> sort!
    println() 

    # save ranges
    pfpath = iJR.procdir("MSE_study_params.jls")
    serialize(pfpath, (;exp, method, exp_biom_beta, 
        exp_vg_beta, biom_betas, vg_betas)
    )

    ## ----------------------------------------------------------------------------
    # globals
    wlk = ReentrantLock()
    model0 = dat[:model] |> ChU.compressed_model
    M, N = size(model0)
    biom_idx = ChU.rxnindex(model0, iJR.BIOMASS_IDER)
    vg_idx = ChU.rxnindex(model0, iJR.EX_GLC_IDER)
    
    last_epout = nothing

    # Feed jobs
    Ch = Channel(nthreads()) do ch
        for (biom_betai, biom_beta) in enumerate(biom_betas)
            put!(ch, (biom_betai, biom_beta))
        end
    end
    
    # ----------------------------------------------------------------------------
    @threads for _ in 1:nthreads()
        thid = threadid()
        thid == 1 && continue

        model = deepcopy(model0) |> ChU.uncompressed_model
        beta_vec = zeros(N)
        # epmodel = ChU.EPModel(model; beta_vec, alpha = Inf)
        
        for (biom_betai, biom_beta) in Ch
            
            # ----------------------------------------------------------------------------
            epouts = Dict()

            # ----------------------------------------------------------------------------
            # check cache
            fname = UJL.mysavename("MSE_study_epouts", "jls"; biom_beta)
            fpath = iJR.procdir(fname)
            if isfile(fpath)
                epouts = deserialize(fpath)
                last_epout = deepcopy(epouts[first(vg_betas)])
                continue
            end

            # ----------------------------------------------------------------------------
            epout = deepcopy(last_epout)
            
            # ----------------------------------------------------------------------------
            for (vg_betai, vg_beta) in enumerate(vg_betas)
                
                up = rem(vg_betai, 10) == 0.0 || vg_betai == 1 || vg_betai == bins
                up && lock(wlk) do
                    @info("Doing  $("-"^50)", 
                        bins,
                        (biom_betai, vg_betai),
                        (exp_biom_beta, exp_vg_beta),
                        (biom_beta, vg_beta), 
                        thid
                    ); println()
                end
                
                # Update
                # ChU.update_solution!(epmodel, epout)
                # epmodel.beta_vec[biom_idx] = biom_beta
                # epmodel.beta_vec[vg_idx] = vg_beta
                beta_vec[biom_idx] = biom_beta
                beta_vec[vg_idx] = vg_beta

                # EP
                # epout = ChEP.converge_ep!(epmodel; 
                #     epsconv = 1e-4, verbose = false
                # )
                epout = ChEP.maxent_ep(model; 
                    alpha = Inf, beta_vec,
                    epsconv = 1e-2, verbose = false
                )

                # TODO: check err and do recovery

                # save
                lock(wlk) do
                    epouts[vg_beta] = deepcopy(epout)
                end

            end # for vg_beta

            # save
            serialize(fpath, epouts)
            
            last_epout = epouts[first(vg_betas)]

        end # for biom_beta

    end # for thid

end