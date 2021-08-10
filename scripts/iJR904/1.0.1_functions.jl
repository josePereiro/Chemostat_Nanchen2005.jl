## -------------------------------------------------------------------
function load_raw_model()
    model = iJR.load_model("raw_model")
    model_size = size(model)
    nz_abs = ChU.nzabs_range(model.S)
    @info("MAT MODEL LOADED", model_size, nz_abs)
    ChN.test_fba(model, iJR.BIOMASS_IDER; summary = false)
    return model
end

## -------------------------------------------------------------------
# rescale bounds bounds
# The abs maximum bounds will be set to 100
function rescale_bounds!(model)
    @info("CLAMP BOUNDS", iJR.ABS_MAX_BOUND)
    foreach(model.rxns) do ider
            ChU.isfixxed(model, ider) && return # fixxed reaction are untouched

            old_ub = ChU.ub(model, ider)
            new_ub = old_ub == 0.0 ? 0.0 : iJR.ABS_MAX_BOUND
            ChU.ub!(model, ider, new_ub)

            old_lb = ChU.lb(model, ider)
            new_lb = old_lb == 0.0 ? 0.0 : -iJR.ABS_MAX_BOUND
            ChU.lb!(model, ider, new_lb)
    end
    return model
end

## -------------------------------------------------------------------
# CLOSING EXCHANGES
function close_exchanges!(model, exchs)
    
    @info("CLOSE EXCANGES", length(exchs))
    # Close, for now, all ChU.exchanges for avoiding it to be in revs
    # The reversible reactions will be splited for modeling cost
    # Exchanges have not associated cost, so, we do not split them
    foreach(exchs) do rxn
        ChU.ub!(model, rxn, 0.0) # Closing all outtakes
        ChU.lb!(model, rxn, 0.0) # Closing all intakes
    end
    return model
end

## -------------------------------------------------------------------
# ENZYMATIC COST
function add_cost(model; cost_ub = 1.0)
    ## -------------------------------------------------------------------
    # ENZYMATIC COST INFO
    # The cost will be introduced as a reaction, we follow the same cost models as 
    # Beg et al. (2007): https://doi.org/10.1073/pnas.0609845104.
    # A new balance equations is then added:
    #        Σ(rᵢ*costᵢ) + tot_cost = 0
    #    Because the cost coefficients (costᵢ) < 0 (it resamble a reactant), the system must allocate 
    #    the fluxes (rᵢ) so that Σ(rᵢ*costᵢ) = tot_cost, and tot_cost
    #    are usually bounded [0.0, 1.0]
    cost_info = Dict()
    fwd_ider(rxn) = string(rxn, ChU.FWD_SUFFIX);
    bkwd_ider(rxn) = string(rxn, ChU.BKWD_SUFFIX);
    for rxn in model.rxns
        # The ChU.exchanges, the atpm and the biomass are synthetic reactions, so, 
        # they have should not have an associated enzimatic cost 
        any(startswith.(rxn, ["EX_", "DM_"])) && continue
        rxn == iJR.BIOMASS_IDER && continue
        rxn == iJR.ATPM_IDER && continue
            
        # Only the internal, non reversible reactions have an associated cost
        # We will split the rev reactions, so we save the cost for both versions (fwd, bkwd)
        if ChU.isrev(model, rxn)
            cost_info[fwd_ider(rxn)] = -iJR.beg_enz_cost(rxn)
            cost_info[bkwd_ider(rxn)] = -iJR.beg_enz_cost(rxn)
        else
            cost_info[rxn] = -iJR.beg_enz_cost(rxn)
        end
    end

    ## -------------------------------------------------------------------
    # SPLITING REVS
    nfwd_suffix = ChU.FWD_SUFFIX
    bkwd_suffix = ChU.BKWD_SUFFIX
    @info("SPLITING REVS", nfwd_suffix, bkwd_suffix)
    model = ChU.split_revs(model;
        get_fwd_ider = fwd_ider,
        get_bkwd_ider = bkwd_ider,
    );

    ## -------------------------------------------------------------------
    # ADDING COST REACCION
    cost_met_id = "cost"
    cost_exch_id = iJR.COST_IDER

    # info
    to_add = cost_info |> length
    min_abs_cost = cost_info |> values .|> abs |> minimum
    max_abs_cost = cost_info |> values .|> abs |> maximum
    @info("ADDING COST", 
        to_add, min_abs_cost, max_abs_cost, cost_met_id, cost_exch_id
    )

    M, N = size(model)
    cost_met = ChU.Met(cost_met_id, S = collect(values(cost_info)), rxns = collect(keys(cost_info)), b = 0.0)
    model = ChU.expanded_model(model, M + 1, N + 1)
    ChU.set_met!(model, ChU.findempty(model, :mets), cost_met)
    cost_exch = ChU.Rxn(cost_exch_id, S = [1.0], mets = [cost_met_id], lb = -iJR.ABS_MAX_BOUND, ub = 0.0, c = 0.0)
    ChU.set_rxn!(model, ChU.findempty(model, :rxns), cost_exch);

    # tot_cost is the exchange that controls the bounds of the 
    # enzimatic cost contraint, we bound it to [0, 1.0]
    ChU.lb!(model, cost_exch_id, 0.0);
    ChU.ub!(model, cost_exch_id, cost_ub);

    return model
end

## -------------------------------------------------------------------
# SET BASE EXCHANGE
function reset_exchanges(model, exchs)

    model = ChU.fix_dims(model)

    @info("SETTING EXCHANGES") 
    # To control the intakes just the metabolites defined in the 
    # base_intake_info (The minimum medium) will be opened.
    # The base model will be constraint as in a cultivation with 
    # experimental minimum xi
    # see Cossios paper (see README)

    foreach(exchs) do rxn
        (rxn in [iJR.COST_IDER]) && return
        ChU.ub!(model, rxn, iJR.ABS_MAX_BOUND) # Opening all outakes
        ChU.lb!(model, rxn, 0.0) # Closing all intakes
    end

    # see Cossios paper (see README) for details in the Chemostat bound constraint
    xi = minimum(Nd.val(:xi))
    intake_info = iJR.load_base_intake_info()
    ChSS.apply_bound!(model, xi, intake_info; emptyfirst = true)

    ## -------------------------------------------------------------------
    ChN.test_fba(model, iJR.BIOMASS_IDER, iJR.COST_IDER)
    return model
end

## -------------------------------------------------------------------
function scale_model(model, scale_factor)

    base_nzabs_range = ChU.nzabs_range(model.S)
    base_size = size(model)
    
    # Scale model (reduce S ill-condition)
    scl_model = iszero(scale_factor) ? 
        deepcopy(model) :
        ChU.well_scaled_model(model, scale_factor)
    
    scl_size = size(scl_model)
    scl_nzabs_range = ChU.nzabs_range(scl_model.S)

    @info("Model", 
        base_size, base_nzabs_range, 
        scl_size, scl_nzabs_range
    ); println()

    return scl_model
end

## -------------------------------------------------------------------
function make_fva_model(model, exp, D; scale_factor = 1000.0)
    @info("DOING FVA", exp, D)

    ## -------------------------------------------------------------------
    # prepare model
    model0 = scale_model(deepcopy(model), scale_factor)

    M, N = size(model0)
    exp_xi = Nd.val(:xi, exp)
    intake_info = iJR.intake_info(exp)
    ChSS.apply_bound!(model0, exp_xi, intake_info; 
        emptyfirst = true)

    ChN.test_fba(exp, model0, iJR.BIOMASS_IDER, iJR.COST_IDER)
    fva_model = ChLP.fva_preprocess(model0, 
        check_obj = iJR.BIOMASS_IDER,
        verbose = true
    );
    ChN.test_fba(exp, fva_model, iJR.BIOMASS_IDER, iJR.COST_IDER)

    return fva_model
    
end

## -------------------------------------------------------------------
function make_max_model(model; scale_factor = 1000.0)
    # This model is bounded by the maximum rates found for EColi.
    # Data From:
    # Varma, (1993): 2465–73. https://doi.org/10.1128/AEM.59.8.2465-2473.1993.
    # Extract max exchages from FIG 3 to form the maximum polytope

    @info("DOING MAX MODEL")

    max_model = scale_model(deepcopy(model), scale_factor)
    
    # Biomass
    # 2.2 1/ h
    ChU.bounds!(max_model, iJR.BIOMASS_IDER, 0.0, 2.2)
    
    Fd_rxns_map = iJR.load_rxns_map() 
    # 40 mmol / gDW h
    # ChU.bounds!(max_model, Fd_rxns_map["GLC"], -40.0, 0.0)
    ChU.bounds!(max_model, Fd_rxns_map["GLC"], -8.0, 0.0)
    # 45 mmol/ gDW
    ChU.bounds!(max_model, Fd_rxns_map["AC"], 0.0, 40.0)
    # 55 mmol/ gDW h
    ChU.bounds!(max_model, Fd_rxns_map["FORM"], 0.0, 55.0)
    # 20 mmol/ gDW h
    ChU.bounds!(max_model, Fd_rxns_map["O2"], -20.0, 0.0)
    
    # fva
    max_model = ChLP.fva_preprocess(max_model, 
        check_obj = iJR.BIOMASS_IDER,
        verbose = true
    );

    ## -------------------------------------------------------------------
    test_model = deepcopy(max_model)
    for (exp, D) in Nd.val(:D) |> enumerate
        cgD_X = Nd.cval(:GLC, exp) * Nd.val(:D, exp) / Nd.val(:X, exp)
        ChU.lb!(test_model, iJR.EX_GLC_IDER, -cgD_X)
        fbaout = ChLP.fba(test_model, iJR.BIOMASS_IDER)
        biom = ChU.av(test_model, fbaout, iJR.BIOMASS_IDER)
        @info("Test", exp, cgD_X, D, biom); println()
    end

    return max_model
end