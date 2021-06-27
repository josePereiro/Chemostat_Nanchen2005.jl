import DrWatson: quickactivate
quickactivate(@__DIR__, "Chemostat_Nanchen2006")

@time begin

    import SparseArrays

    import Chemostat_Nanchen2006
    const ChN = Chemostat_Nanchen2006

    const iJR = ChN.iJR904
    const Nd = ChN.NanchenData # experimental data
    const Bd = ChN.BegData    # cost data

    import Chemostat
    const Ch = Chemostat
    const ChU = Ch.Utils
    const ChSS = Ch.SteadyState
    const ChLP = Ch.LP
end


## -------------------------------------------------------------------
# Helper function
include("1.0.1_functions.jl")

## -------------------------------------------------------------------
MODELS_FILE = iJR.procdir("base_models.bson")
MODELS = isfile(MODELS_FILE) ? 
    ChU.load_data(MODELS_FILE; verbose = false) : 
    Dict{Any, Any}();

## -------------------------------------------------------------------
# COST-CONSTRAINT MODELS
let
    compress(model) = ChU.compressed_model(model)

    ## --------------------------------------
    ChU.tagprintln_inmw("COST_CONSTRAINED MODELS\n")

    base_model = load_raw_model()
    exchs = ChU.exchanges(base_model)
    
    rescale_bounds!(base_model)
    close_exchanges!(base_model, exchs)
    base_model = add_cost(base_model)
    base_model = reset_exchanges(base_model, exchs)

    MODELS["base_model"] = compress(base_model)

    ## --------------------------------------
    scale_factor = 1000.0
    # fva_preprocessed models
    for (exp, D) in Nd.val(:D) |> enumerate

        DAT = get!(MODELS, "fva_models", Dict())
        haskey(DAT, exp) && continue # cached

        fva_model = make_fva_model(base_model, exp, D; 
            scale_factor
        )
        
        # storing
        DAT[exp] = compress(fva_model)

        # caching
        ChU.save_data(MODELS_FILE, MODELS);
        GC.gc()
    end

    ## --------------------------------------
    # MAX MODEL
    max_model = make_max_model(base_model; 
        scale_factor
    )

    # saving
    MODELS["max_model"] = compress(max_model)
    ChU.save_data(MODELS_FILE, MODELS)

    println()
end

## -------------------------------------------------------------------
# COST-LESS MODELS
let
    compress(model) = ChU.compressed_model(model)

    ## --------------------------------------
    ChU.tagprintln_inmw("COST-LESS MODELS\n")

    base_model = load_raw_model()
    exchis = ChU.exchanges(base_model)
    exchs = base_model.rxns[exchis]
    
    rescale_bounds!(base_model)
    close_exchanges!(base_model, exchs)
    base_model = add_cost(base_model; cost_ub = 1.0)
    base_model = reset_exchanges(base_model, exchs)

    MODELS["base_model_costless"] = compress(base_model)

    ## --------------------------------------
    scale_factor = 1000.0
    # fva_preprocessed models
    for (exp, D) in Nd.val(:D) |> enumerate

        DAT = get!(MODELS, "fva_models_costless", Dict())
        haskey(DAT, exp) && continue # cached

        fva_model = make_fva_model(base_model, exp, D; 
            scale_factor
        )
        
        # storing
        DAT[exp] = compress(fva_model)

        # caching
        ChU.save_data(MODELS_FILE, MODELS);
        GC.gc()

    end

    ## --------------------------------------
    # MAX MODEL
    max_model = make_max_model(base_model; 
        scale_factor
    )

    # saving
    MODELS["max_model_costless"] = compress(max_model)
    ChU.save_data(MODELS_FILE, MODELS)

end