using ProjAssistant
@quickactivate

@time begin
    import SparseArrays
    import Base.Threads: @threads, threadid, SpinLock

    import SimTools
    const SimT = SimTools

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

    using ArgParse
    using Statistics
    using Serialization
    using Base.Threads
    using ExtractMacro
end

## ----------------------------------------------------------------------------
# arg settings
ARGSET = ArgParse.ArgParseSettings()
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

## ----------------------------------------------------------------------------
# FILE GLOBALS
const WLOCK = ReentrantLock()
const IGNORE_EXPS = []

## -------------------------------------------------------------------
# UTILS
dat_file(;kwargs...) = procdir(iJR, "maxent_ep_dat", kwargs..., ".jls")

function is_cached(;kwargs...)
    ignore_cached && return false
    thid = threadid()
    datfile = dat_file(;kwargs...)
    isfile(datfile) && lock(WLOCK) do
        @info("Cache Found",
            datfile, thid
        ); println()
    end; isfile(datfile)
end

## -------------------------------------------------------------------
# MAXENT FUNS
include("2.0.1_ME_MAX_POL.jl")
include("2.0.2_ME_MAX_POL_B0.jl")
include("2.0.3_ME_Z_EXPECTED_G_BOUNDED.jl")

## -------------------------------------------------------------------
# ME GLOBALS
sglob(iJR, :maxent, :params) do
    (;
        alpha = Inf,
        epsconv = 9e-4,
        maxiter = 3000,
        damp = 0.9,
        maxvar = 1e50,
        minvar = 1e-50,
    )
end

## -------------------------------------------------------------------
# ME_MAX_POL
let
    method = iJR.ME_MAX_POL
    model_key = "max_model"
    
    maxent_max_pol(method, model_key)
end

## ----------------------------------------------------------------------------
# ME_MAX_POL_B0
let
    method = iJR.ME_MAX_POL_B0
    model_key = "max_model"
    
    do_max_pol_b0(method, model_key)
end

## ----------------------------------------------------------------------------
# ME_Z_EXPECTED_G_BOUNDED
let
    method = iJR.ME_Z_EXPECTED_G_BOUNDED
    model_key = "fva_models"
    
    do_z_expected_ug_bounded(method, model_key)
end
