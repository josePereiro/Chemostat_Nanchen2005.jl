# initial approach
function do_z_expected_ug_bounded(method, model_key)

    # orig model
    iterator = Nd.val(:D) |> enumerate |> collect 
    @threads for (exp, D) in iterator
        thid = threadid()

        ## -------------------------------------------------------------------
        # ep
        me_params = lglob(iJR, :maxent, :params)
        @extract me_params: alpha epsconv maxiter damp maxvar minvar

        # approach
        convth = 0.01
        converr = nothing

        ## -------------------------------------------------------------------
        # handle cache
        datfile = dat_file(;method, exp)
        is_cached(;method, exp) && continue

        # prepare model
        model = iJR.load_model(model_key, exp)
        biomidx = ChU.rxnindex(model, iJR.BIOMASS_IDER)
        M, N = size(model)
        exp_growth = Nd.val(:D, exp)
        growth_ub = ChU.ub(model, iJR.BIOMASS_IDER)
        feasible = exp_growth <= growth_ub
        if !feasible
            lock(WLOCK) do
                @info("Not feasible (skipping)", 
                    exp, method, 
                    growth_ub ,exp_growth, 
                    thid
                ); println()
            end
            continue
        end
        ChU.ub!(model, iJR.BIOMASS_IDER, growth_ub * 1.1) # open a bit helps EP

        lock(WLOCK) do
            nzabs_range = ChU.nzabs_range(model.S)
            @info("Starting... ", 
                exp, method,
                size(model), nzabs_range, 
                feasible,
                thid
            ); println()
        end

        # simulation
        dat = ldat(() -> Dict(), iJR, datfile)
        epouts = get!(dat, :epouts, Dict())
        init_len = length(epouts)
        beta_vec = zeros(N)
        approach_status = get!(dat, :approach_status, :running)
        approach_status == :finished && continue
        
        # log approach
        epout_seed = isempty(epouts) ? nothing : epouts[maximum(keys(epouts))]
        betas = [0.0; 10.0.^(3:0.05:15)]
        nan_beta = first(betas)

        for approach in [:log_approach, :linear_approach]
            
            lock(WLOCK) do
                @info("Starting", approach, 
                    exp, method,  
                    length(epouts),
                    thid
                ); println()
            end

            for beta in betas

                nan_beta = beta
                haskey(epouts, beta) && continue

                beta_vec[biomidx] = beta
                epout = nothing
                try
                    epout = ChEP.maxent_ep(model; 
                        beta_vec, alpha, damp, maxiter,
                        epsconv, maxvar, minvar, 
                        verbose = false,
                        solution = epout_seed
                    )
                catch err; end

                # stop conditions
                biom_avPME = isnothing(epout) ? 0.0 : ChU.av(model, epout, biomidx)
                fail = isnothing(epout) || isnan(biom_avPME) || biom_avPME == 0.0 
                converr = abs(biom_avPME - exp_growth)/exp_growth

                # info
                lock(WLOCK) do
                    @info("Results", approach, 
                        exp, method, beta, 
                        exp_growth, growth_ub, biom_avPME, 
                        converr, convth,
                        length(epouts),
                        thid
                    ); println()
                end

                # error
                fail && break

                # updating
                epout_seed = epouts[beta] = epout

                # convergence
                conv = converr < convth || biom_avPME > exp_growth 
                conv && break

            end # for betas

            # Catching
            update = init_len != length(epouts)
            update && lock(WLOCK) do
                sdat(iJR, dat, datfile)
                @info("Catching", approach, 
                    exp, method,  
                    length(epouts),
                    basename(datfile),
                    thid
                ); println()
            end

            # lineal approach
            last_beta = maximum(keys(epouts))
            step = abs(last_beta - nan_beta) / 100.0
            iszero(step) && break
            betas = range(last_beta, 1.0e15; step)

        end # for approach

        # saving
        lock(WLOCK) do
            dat[:approach_status] = :finished
            dat[:model] = model |> ChU.compressed_model
            dat[:exp_beta] = maximum(keys(epouts))
            sdat(iJR, dat, datfile)
            @info("Finished", exp, method,  
                length(epouts),
                thid
            ); println()
        end
       
    end # for (exp, D)

    do_z_expected_ug_bounded_further_conv(method)

end


## ----------------------------------------------------------------------------
# Further convergence
function do_z_expected_ug_bounded_further_conv(method)

    # ep
    me_params = lglob(iJR, :maxent, :params)
    @extract me_params: alpha epsconv maxiter damp maxvar minvar

    iterator = Nd.val(:D) |> enumerate |> collect 
    @threads for (exp, D) in iterator

        ## -------------------------------------------------------------------
        # thread globals
        
        # handle cache
        datfile = dat_file(;method, exp)
        dat = ldat(iJR, datfile)
        model = ChU.uncompressed_model(dat[:model])
        epouts = dat[:epouts]

        exp_beta = dat[:exp_beta]
        exp_epout = epouts[exp_beta]

        lock(WLOCK) do
            @info("Converging...", exp, method,
                exp_beta, exp_epout.status, 
                threadid()
            ); println()
        end
        converg_status = get!(dat, :converg_status, :undone)
        converg_status == :done && continue

        new_epout = nothing
        if exp_epout.status == :unconverged
            try;
                biomidx = ChU.rxnindex(model, iJR.BIOMASS_IDER)
                beta_vec = zeros(size(model, 2)); 
                beta_vec[biomidx] = exp_beta
                new_epout = ChEP.maxent_ep(model; 
                    beta_vec, alpha, damp, epsconv, 
                    maxvar, minvar, verbose = false, 
                    solution = exp_epout, maxiter
                )
            catch err; @warn("ERROR", err); println() end

            biom_avPME = isnothing(new_epout) ? 0.0 : ChU.av(model, new_epout, iJR.BIOMASS_IDER)
            fail = isnan(biom_avPME) || biom_avPME == 0.0 
            epouts[exp_beta] = fail ? exp_epout : new_epout
        end
        
        # Saving
        lock(WLOCK) do
            @info("Saving...", exp, exp_beta, 
                exp_epout.status,
                threadid()
            ); println()
        end
        dat[:converg_status] = :done
        sdat(iJR, dat, datfile)
    end
end
