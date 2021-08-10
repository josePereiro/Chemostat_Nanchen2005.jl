using ProjAssistant
@quickactivate 

## -------------------------------------------------------------------
@time begin

    import SparseArrays
    using ArgParse

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

## ----------------------------------------------------------------------------
# arg settings
ARGSET = ArgParseSettings()
@ArgParse.add_arg_table! ARGSET begin
    "--ignore-cached"
        help = "Ingnore on disk version of data"
        action = :store_true
end

ARGS_DICT = ArgParse.parse_args(ARGSET)
if isinteractive()
    # dev values
    ignore_cached = true
else
    ignore_cached = ARGS_DICT["ignore-cached"]
end
@info("ARGS", ignore_cached)

## -------------------------------------------------------------------
# Helper function
include("1.0.1_functions.jl")

## -------------------------------------------------------------------
MODELS_FILE = procdir(iJR, "base_models.bson")
MODELS = ldat(() -> Dict{Any, Any}(), iJR, MODELS_FILE);

## -------------------------------------------------------------------
# file globals
use_cached = false

## -------------------------------------------------------------------
let
    compress(model) = ChU.compressed_model(model)

    ## --------------------------------------
    @info "COST_CONSTRAINED MODELS";

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
        !ignore_cached && haskey(DAT, exp) && continue # cached

        fva_model = make_fva_model(base_model, exp, D; 
            scale_factor
        )
        
        # storing
        DAT[exp] = compress(fva_model)

        # caching
        sdat(iJR, MODELS, MODELS_FILE);
        GC.gc()
    end

    ## --------------------------------------
    # MAX MODEL
    max_model = make_max_model(base_model; scale_factor)

    # saving
    MODELS["max_model"] = compress(max_model)
    sdat(iJR, MODELS, MODELS_FILE)

    println()
end