function do_max_pol_b0(method, model_key)
    

    objider = iJR.BIOMASS_IDER

    # handle cache
    datfile = dat_file(;method)
    check_cache(;method) && return

    # setup
    model = iJR.load_model(model_key)

    lock(WLOCK) do
        @info("Doing... ", 
            method, threadid()
        ); println()
    end

    # maxent
    epout = ChEP.maxent_ep(model; 
        alpha = Inf, damp = 0.9, epsconv = 1e-4, 
        verbose = true, maxiter = 5000
    )
        
    # storing
    lock(WLOCK) do
        # Storing
        dat = Dict()
        dat[:exp_beta] = 0.0
        dat[:epouts] = Dict(0.0 => epout)
        dat[:model] = model |> ChU.compressed_model

        # caching
        serialize(datfile, dat)

        @info("Finished ", threadid())
        println()
    end
end

## ------------------------------------------------------------------------
let
    for (method, model_key) in [
        # (ME_MAX_POL_B0, "max_model"),
        (ME_MAX_POL_B0_COSTLESS, "max_model_costless")
    ]
        maxent_max_pol(method, model_key)
    end
end